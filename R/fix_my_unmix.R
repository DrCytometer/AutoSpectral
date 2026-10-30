# fix_my_unmix.R
#
# Simplifications: every fully stained file is unmixed with the one
# `af.spectra` library, and the spillover spread is not recomputed under the
# corrected spectra (re-run get.spectral.variants() with them for that).

#' @title Fix My Unmix
#'
#' @description
#' Corrects the component of a reference spectra error that lives inside the row
#' space of the spectra themselves, using fully stained samples and the fact
#' that a marker-negative population must read zero in every channel it is
#' negative for.
#'
#' This is the complement of `correct.unmixing.signatures()`, which corrects the
#' shape of individual rows from detector-space residuals and is structurally
#' blind to in-span error. Run this second, on the spectra that function
#' returns, because in unmixed space a row-shape error and a genuine spillover
#' error are indistinguishable.
#'
#' The correction runs in two phases. The first estimates the residual spillover
#' matrix pair by pair, from a robust line through the events that are negative
#' for the target channel: their apparent abundance there is spillover from the
#' source, so its slope against source abundance is the residual coefficient.
#' The slope is taken from the source's clearly positive events only: those
#' above both the spread-scaled unstained boundary and the upper bound of the
#' stained sample's own negative population in that channel (its mode plus
#' the mirrored lower tail, `negative.tail`), so the tilted tail of a negative
#' population collinear with other dyes never sets a slope. The
#' source-negative bulk fixes where the line passes (`estimator =
#' "anchored"`) but not its slope, and target negativity is a band around the
#' current line after a first one-sided pass, so the source's own
#' cell-to-cell variation is trimmed on both sides.
#'
#' Each coefficient earns a trust weight that damps the step it applies: from a
#' prior variance (spectral variants propagated through the unmixing matrix)
#' against the coefficient's standard error, from the agreement between the
#' binned envelope and median slopes (co-expression), from how many clearly
#' positive events the source has beyond the rate the unstained sample gives
#' (`source.events`), from its agreement with a joint lasso of the channel over
#' every other fluorophore (`glasso.check`), and from the share of the source's
#' own positive population it leaves below the channel's negative boundary. A
#' negative correlation between a source's abundance and a channel's corrected
#' value has no biological explanation the way co-expression can, so that last
#' check reaches a hard rejection at `max.hypernegative.frac`. Phase one stops
#' when every coefficient with at least `stop.trust` trust has settled
#' (`stop.rule = "settled"`), or when the change stalls.
#'
#' The second phase does not adopt the back-solved spectrum as it stands. It
#' uses the spillover matrix only to define clean populations, then re-measures
#' each changed fluorophore's signature directly from background-subtracted raw
#' data with `extract.raw.signature()`, fitted jointly with the dyes present on
#' those events, so the result is a physical spectrum rather than whichever
#' combination of the existing rows reproduces the required compensation. A
#' candidate is adopted only if it passes every gate. When the re-measured row
#' fails a gate that says nothing against the row itself and phase one's
#' evidence for it is strong, phase one's back-solved row is tried through the
#' same gates (`phase.two.rescue`), and a tandem whose parent is in the panel has
#' a last candidate, its starting row with only the parent fraction changed by
#' the tandem-to-parent coefficient phase one measured (`tandem.candidate`):
#' breakdown, or a reference with more breakdown than the sample.
#'
#' Every change is held to a plausible size: a detector of a row may change by
#' at most `plausibility.ratio` times its starting value or `plausibility.floor`
#' (peak = 1), whichever is greater. Real reference errors fit inside this box;
#' a dim marker co-expressed on a dye's cells, or autofluorescence taken up by
#' it, adds signal where the dye had none and does not. Phase one takes the
#' out-of-box part of each back-solved row back out of its coefficients; phase
#' two resets out-of-box detectors of an adopted row and rejects a row with
#' more than `plausibility.reject` of its change outside the box.
#'
#' With `smooth.updates` the correction runs in rounds, each on a fresh random
#' share of the stained events unmixed again with the current spectra, and each
#' moving the spectra only a step `1 / round` toward that round's result. A
#' single pass adopts its subsample's gate decisions, which are hard cuts on
#' noisy statistics; the rounds average them. With `smooth.mode =
#' "consensus"` every round instead starts from the input spectra, and a row is
#' adopted only when enough rounds corrected it and their results agree.
#'
#' Autofluorescence is by default removed per event (`bg.mode = "per.cell"`),
#' each event taking the single `af.spectra` row that minimises its apparent
#' fluorophore signal, as `unmix.autospectral()` does, so the correction is
#' measured on the unmixing it will be used with. `bg.mode = "af.deconv"` fits a
#' shared multi-component basis instead.
#'
#' @importFrom stats mad median qnorm quantile sd setNames lm.wfit
#'
#' @param spectra The spectral matrix, fluorophores x detectors, L-infinity
#'   normalised. Ideally the output of `correct.unmixing.signatures()`.
#' @param unstained.sample File path and name for a raw unstained sample,
#'   acquired the same day and matching the autofluorescence of the fully
#'   stained samples.
#' @param fully.stained.sample File path(s) and name(s) of one or more raw fully
#'   stained samples. Several are read and pooled, each contributing an equal
#'   share of `max.events`.
#' @param flow.control The flow.control list.
#' @param asp The AutoSpectral parameter list.
#' @param variants The variant list returned by `get.spectral.variants()`. Its
#'   `spillover.spread` sets abundance-dependent positivity boundaries, and its
#'   deltas the variant leakage prior (`leakage.prior`). If `variants` is
#'   `NULL`, positivity boundaries are flat: the hypernegative check then
#'   removes the trust of nearly every coefficient and large positive errors
#'   are hidden from the pair fit.
#' @param af.name Character or `NULL`, the name of an autofluorescence row in
#'   `spectra`. It is never treated as a panel fluorophore and never corrected.
#'   Default `"AF"`.
#' @param af.basis Optional matrix (components x detectors) from
#'   `get.af.basis()`. When `NULL` and `bg.mode = "af.deconv"`, it is built from
#'   the unstained sample. Default `NULL`.
#' @param af.n.pc Integer or `"auto"`, passed to `get.af.basis()`. Default
#'   `"auto"`.
#' @param af.spectra Matrix (n.af x detectors), a library of candidate
#'   autofluorescence spectra from `get.af.spectra()`, row 1 its population
#'   mean. Required, with at least two rows, when `bg.mode = "per.cell"`. Also
#'   the reference of the AF gate and of `af.likeness` in the logs.
#' @param bg.mode Character. `"per.cell"` (default) assigns each event the
#'   `af.spectra` row whose removal leaves the smallest total apparent
#'   fluorophore signal (`assign.af.fluorophores()`) and subtracts only that
#'   one; the library's mean row stands in for autofluorescence in the
#'   conditioning, leakage and prior calculations. `"af.deconv"` fits a
#'   multi-component autofluorescence basis jointly with the panel;
#'   `"af.row"` uses the single `af.name` row already in `spectra`;
#'   `"global.mean"` uses the mean unstained spectrum as a single background
#'   row; `"none"` fits no background. Under every mode but `"per.cell"`, a
#'   fluorophore inseparable from the background (`max.hotspot`) is frozen.
#' @param large.gate Logical, whether to use a large scatter gate. Default
#'   `TRUE`.
#' @param scatter.gate Logical, whether to gate on scatter at all before
#'   fitting. `FALSE` keeps every event, ignoring `large.gate`. Default `FALSE`.
#' @param landmark.quantile Numeric in `(0, 1)` or `NULL`: events above this
#'   quantile on either scatter parameter are kept in addition to whatever the
#'   gate polygon keeps. The threshold is computed once from
#'   `unstained.sample` and reused for the stained samples. Ignored when
#'   `scatter.gate = FALSE`. Default `NULL`.
#' @param max.events Integer, the most events read from the unstained sample,
#'   and the total read from the fully stained samples, split equally between
#'   them. Each file is subsampled with its own seed. `NULL` reads every event.
#'   With 200,000 events about half of the corrections on 40-50 colour panels
#'   varied between disjoint halves of the data. Default `400000`.
#' @param max.iter Integer, maximum spillover-matrix iterations. Default `20`.
#' @param downsample Logical or numeric. `FALSE` disables downsampling; a
#'   numeric gives the number of events to use, chosen by a stratified sample
#'   over each event's positive fluorophores so a dim or rare dye's own positive
#'   population is not thinned at the same rate as the background. Default
#'   `FALSE`.
#' @param downsample.background.frac Numeric in (0, 1), the share of
#'   `downsample` reserved for events positive for nothing. Default `0.3`.
#' @param downsample.min.stratum Integer, the floor below which a
#'   fluorophore's own positive population is kept whole by the stratified
#'   downsample. Default `2000`.
#' @param unstained.threshold Numeric in (0.5, 1), the quantile of the
#'   unstained events defining the upper (positivity) boundary. Default `0.95`.
#' @param lower.unstained.threshold Numeric in (0.5, 1); the lower
#'   (hypernegative) boundary is the `1 - lower.unstained.threshold` quantile
#'   of the unstained events. A correctly corrected population has about
#'   `1 - lower.unstained.threshold` of its events below it, so this must stay
#'   well under `max.hypernegative.frac` or the hypernegative check removes the
#'   trust of correct coefficients. Default `0.99`.
#' @param unstained.margin Numeric, multiplier applied to both boundaries.
#'   Default `1.3`.
#' @param spread.kappa Numeric, how many spillover-spread standard deviations
#'   beyond the flat boundary still count as negative. The Residual Model
#'   spread is narrower than the real cell-to-cell variation of stained cells,
#'   so the default is `3`.
#' @param negative.tail Numeric in (0, 0.5), the share of the stained sample's
#'   own negative population in a source channel allowed above its upper
#'   bound (the population's mode plus its mirrored lower tail). A source's
#'   clearly positive events are those above both that bound and the
#'   positivity boundary. Default `0.001`.
#' @param source.events Numeric pair. A source's clearly positive events,
#'   less the number the unstained sample gives at the same event count by the
#'   same rule, earn its coefficients no trust at or below the first value and
#'   full trust at or above the second, rising smoothly between. Partial trust
#'   accumulates over iterations, so a low ramp lets a few hundred events add
#'   error. Default `c( 200, 400 )`.
#' @param envelope.quantiles Numeric pair. The first is the binned envelope;
#'   the second a comparison quantile whose disagreement with the first reports
#'   co-expression. Default `c( 0.05, 0.5 )`.
#' @param disagreement.reference Character, the slope the median is compared
#'   with. `"envelope.corrected"` (default) moves each bin's envelope back to its
#'   centre by the bin's own lower-tail scale, measured at
#'   `envelope.scale.quantile`; the plain `"envelope"` falls with the spread as
#'   the source brightens, so it reports disagreement on clean populations;
#'   `"anchored"` compares with the fitted line, which moves with co-expression
#'   inside the band as the median does.
#' @param envelope.scale.quantile Numeric, the upper quantile of each bin's
#'   lower tail for `"envelope.corrected"`. Default `0.25`.
#' @param min.negative.frac Numeric, the smallest fraction of target-negative
#'   events the brightest source bins may contain before the pair is declared
#'   unidentifiable and left alone. Default `0.10`.
#' @param max.disagreement Numeric, envelope-versus-median slope disagreement
#'   above which the coefficient's trust weight is reduced. Default `0.5`.
#' @param min.hypernegative.frac Numeric, the fraction of a source's own
#'   positive population tolerated below the channel's negative boundary
#'   before trust begins to fall. Default `0.01`.
#' @param max.hypernegative.frac Numeric, the fraction at which trust reaches
#'   zero; between the two, trust falls linearly. The fraction checked is
#'   absolute, not a change against what already sat past the boundary.
#'   Default `0.05`.
#' @param min.hypernegative.events Integer, the fewest events a source's
#'   population may contain before its hypernegative fraction is estimated.
#'   Default `200`.
#' @param hypernegative.dominant Logical, whether only the source's clearly
#'   positive events that it also dominates are judged. An event another
#'   fluorophore dominates reads positive for the source only through that
#'   fluorophore's error, and removing the source's slope pushes it below
#'   zero. Requires `source.dominant`. Default `TRUE`.
#' @param leakage.prior Logical, whether to weight coefficients by the spillover
#'   a plausible spectral error could produce, from
#'   `get.variant.leakage.prior()`. When `FALSE`, or when no variants are
#'   supplied, the estimate supplies its own prior variance. Default `TRUE`.
#' @param span.fraction Numeric in (0, 1], passed to
#'   `get.variant.leakage.prior()`. Default `0.6`.
#' @param min.negative.events Integer, minimum events for a pair fit. Default
#'   `200`.
#' @param min.bin.negative Integer, minimum target-negative events an abundance
#'   bin must contain before its envelope quantile is used. Default `25`.
#' @param min.span Numeric, minimum source abundance span in units of that
#'   fluorophore's own threshold. Default `5`.
#' @param min.rise Numeric, the fitted rise across the source abundance span,
#'   in standard deviations of the target's negative population, below which a
#'   coefficient is treated as unresolvable. Default `1`.
#' @param n.levels Integer, the maximum number of abundance bins for the
#'   signature fit. Must exceed the panel size by at least three when
#'   `multivariate` is `TRUE`. Default `60`.
#' @param min.bin.events Integer, the fewest events an abundance bin may contain,
#'   passed to `extract.raw.signature()`. Default `50`.
#' @param n.levels.pair Integer, abundance bins for the pair estimator's
#'   envelope. Default `10`.
#' @param multivariate Logical, passed to `extract.raw.signature()`, whether
#'   the active abundances are fitted jointly rather than the nuisance rows
#'   being subtracted beforehand. Default `TRUE`.
#' @param ridge Numeric, ridge penalty for that joint fit. Default `1e-6`.
#' @param max.truncated.events Integer, cap on the events used for the robust
#'   pair fit and the lasso. Default `20000`.
#' @param max.mask.passes Integer, how many times the target-negative selection
#'   is recomputed with the current spillover estimate removed. Default `3`.
#' @param estimator Character, which pair estimator supplies the coefficients.
#'   `"anchored"` (default) is a Huber line through the clearly positive
#'   source events that are negative for the target, anchored at the robust
#'   centre of the source-negative bulk. `"truncated"` is a robust line through
#'   every target-negative event, bulk included. `"envelope"` is the binned
#'   lower envelope.
#' @param glasso.check Logical, whether each pair coefficient's trust is
#'   scaled by its agreement with a joint lasso of the channel over every
#'   other fluorophore, refitted without the penalty on the lasso's selection,
#'   `1 - |pair - lasso| / (|pair| + |lasso|)`, zero where the lasso selects
#'   nothing. Default `TRUE`.
#' @param glasso.warm Logical, whether the lasso is fitted after the pair
#'   slopes, with them removed from its first negativity mask. Without it, a
#'   source whose positive residual spillover lifts its positive events above
#'   the target's boundary is missing from the lasso's fit, which then returns
#'   zero or the wrong sign and removes the trust of a well-fitted slope.
#'   Default `TRUE`.
#' @param glasso.n.lambda,glasso.lambda.min.ratio,glasso.margin.frac Lasso
#'   path length, smallest penalty relative to the largest, and sparsity
#'   margin; see `correct.spectra.glasso()`. Defaults `40`, `1e-3`, `0.05`.
#' @param glasso.tolerance Numeric, lasso coefficients smaller than this are
#'   zeroed. Default `1e-4`.
#' @param af.gate Logical, whether a source's step is dropped when the change
#'   it implies for that row looks like autofluorescence: the largest squared
#'   cosine between the change and a row of `af.spectra` (or of the background
#'   basis) above `af.gate.max`. The same test rejects a phase-two change
#'   (`"af.like"`). Default `FALSE`.
#' @param af.gate.max Numeric in `[0, 1]`. Default `0.8`.
#' @param plausibility.ratio Numeric, how many times its starting value a
#'   detector of a row may change by. `0` turns the plausibility box off.
#'   Default `1`.
#' @param plausibility.floor Numeric, the change any detector may take
#'   regardless of its starting value (peak = 1). Default `0.003`.
#' @param plausibility.reject Numeric in (0, 1], the share of a phase-two
#'   change (Euclidean norm) outside the box above which the row is rejected
#'   (`"implausible"`). Default `0.5`.
#' @param plausibility.reset Logical, whether a detector of an adopted row
#'   outside the box keeps its starting value rather than taking the bound.
#'   Default `TRUE`.
#' @param phase.two.rows Character, `"flagged"` (default) re-measures only the
#'   rows phase one changed by at least `convergence.threshold`; the others
#'   are logged as `"not.flagged"`. `"all"` re-measures every row.
#' @param phase.two.active Character, `"present"` (default) fits a signature
#'   jointly with the dyes present in the row's own population only (see
#'   `phase.two.present`), so a population of a few hundred events still gets
#'   a joint fit; `"all"` fits every panel abundance.
#' @param phase.two.present Named numeric, the rules for a dye being present in
#'   a population: above its boundary in a share `frac` of its events, an
#'   abundance correlated with the target's at least `corr`, or a mean
#'   absolute abundance at least `mean` times its boundary. Default
#'   `c( frac = 0.02, corr = 0.3, mean = 0.5 )`.
#' @param phase.two.fallback Character, what a row whose population is too
#'   small for the joint fit is judged on. `"remeasured"` (default) is the
#'   fit with every other row held at its spectrum, which books their errors
#'   against the row but carries the out-of-span part of its own error;
#'   `"backsolved"` is phase one's back-solved row; `"reject"` rejects the row
#'   (`"not.joint"`).
#' @param phase.two.rescue Character, `"backsolved"` (default) tries phase
#'   one's back-solved row as a second candidate when the first was rejected
#'   for one of `rescue.reasons`, the row is flagged, its source weight is at
#'   least `rescue.min.weight` and the back-solved row is within `max.angle`;
#'   `"off"` does not.
#' @param rescue.reasons Character, the rejections that say nothing against the
#'   row itself. Default `c( "leakage", "conditioning", "no.fit", "not.joint" )`.
#' @param rescue.min.weight Numeric in `[0, 1]`. Default `0.5`.
#' @param leakage.channels Character, the off-target channels the leakage gate
#'   sums. `"negative"` (default) only those the row's own population is
#'   negative for (fewer than half of its events above the boundary), so a
#'   marker expressed on the same cells is not counted as leakage; `"all"`
#'   every channel.
#' @param tandem.candidate Logical, whether a tandem whose parent is in the
#'   panel (`tandem.parent` in `fluorophore_database.csv`) has its starting row
#'   with only the parent fraction changed, by phase one's tandem-to-parent
#'   coefficient, as a last candidate when the others fail. It passes the
#'   angle, conditioning and leakage gates but not the plausibility box, since
#'   breakdown adds parent signal where a tandem has little. Default `TRUE`.
#' @param tandem.min.coefficient Numeric, the smallest tandem-to-parent
#'   coefficient that gives a candidate. Default `0.01`.
#' @param source.dominant Logical, whether source-positive events dominated by
#'   another fluorophore are dropped from the pair fit. Default `TRUE`.
#' @param spread.addback Logical, whether the measured spillover spread is added
#'   back to the envelope before its slope is taken. Default `FALSE`.
#' @param anchor.weight Numeric, the weight ceiling for the zero-abundance
#'   anchor bin, in multiples of the best positive bin. Default `1`.
#' @param max.coefficient Numeric, the largest residual spillover coefficient
#'   accepted. Rows are L-infinity normalised, so a coefficient this large means
#'   one row is wrong by that fraction of another row's whole spectrum. Default
#'   `0.5`.
#' @param convergence.threshold Numeric, residual spillover coefficient below
#'   which phase one stops, and the change that flags a row for phase two.
#'   Default `0.003`.
#' @param convergence.quantile Numeric, the quantile of the judged off-diagonal
#'   errors the stop test uses. Default `1`, the largest.
#' @param stop.rule Character, the errors the stop test judges. `"settled"`
#'   (default) those of coefficients with at least `stop.trust` trust, below
#'   `convergence.threshold` for `stop.confirm` iterations in a row, so a slope
#'   too weakly trusted ever to finish does not hold phase one open; `"all"`
#'   every fitted error, once.
#' @param stop.trust Numeric in (0, 1]. Default `0.05`.
#' @param stop.confirm Integer. Default `2`.
#' @param coefficient.decay Numeric in `[0, 1)`, the maximum fraction by which
#'   an off-diagonal coefficient is pulled toward zero on one iteration,
#'   scaled by `1 - trust`. Default `0`.
#' @param update.spectra Logical, whether to run the raw-space signature phase.
#'   When `FALSE`, only the spillover and compensation matrices are returned,
#'   and smoothing does not apply. Default `TRUE`.
#' @param step Numeric, the fraction of each accepted signature change applied.
#'   Default `1`.
#' @param intercept Logical, whether the signature fit carries an intercept.
#'   Default `TRUE`.
#' @param min.explained,max.explained Numeric, the range of `explained.total`,
#'   the fraction of the brightest bin's signal the whole signature fit
#'   reproduces. Defaults `0.8`, `1.2`.
#' @param max.resid Numeric, maximum relative fit residual. Default `0.03`.
#' @param max.intercept Numeric, maximum relative intercept. Default `0.03`.
#' @param min.bg.align Numeric, minimum cosine between the fitted intercept and
#'   the candidate signature. Default `-0.9`.
#' @param max.angle Numeric, maximum angular change of a row in degrees. A large
#'   real correction (a tandem) can need about 20. Default `30`.
#' @param max.clamp.frac Numeric, maximum fraction of a candidate row's
#'   absolute mass that non-negativity clamping may remove. Default `0.15`.
#' @param max.anchor Numeric, maximum unremoved background relative to the
#'   brightest fitted signal. Default `0.10`.
#' @param max.vif Numeric, maximum variance inflation factor for the
#'   fluorophore's own abundance within its population. Default `500`.
#' @param max.condition.increase Numeric, the factor by which a single accepted
#'   row may increase the condition number of the unmixing design. Default
#'   `1.05`.
#' @param peak.shift.min.rel Numeric, a candidate row whose peak detector
#'   differs from the current row is only accepted if the current signature's
#'   value at the new peak, relative to its own peak, is at least this. Default
#'   `0.7`.
#' @param max.hotspot Numeric, hotspot scale above which a fluorophore is
#'   considered inseparable from the background basis and is frozen; not used
#'   under `bg.mode = "per.cell"`. Default `5`.
#' @param leakage.margin Numeric, the fraction by which held-out leakage must
#'   increase before a candidate is refused for it. Default `0.05`.
#' @param smooth.updates Logical, whether to run in rounds (see Description).
#'   Default `TRUE`.
#' @param smooth.rounds Integer, the number of rounds. Default `5`.
#' @param smooth.fraction Numeric in (0, 1], the share of the stained events
#'   each round draws. Default `0.5`.
#' @param smooth.mode Character, `"damped"` (default): each round starts from
#'   the current spectra and moves the rows it changed a step `1 / round`
#'   toward its result; `"consensus"`: each round starts from `spectra`, and a
#'   row is adopted as the mean of the rounds that accepted it when at least
#'   `consensus.min.rounds` of them did and their mean angle to that mean is at
#'   most `consensus.max.spread` of the move it makes.
#' @param consensus.min.rounds Numeric in (0, 1]. Default `0.6`.
#' @param consensus.max.spread Numeric. Default `0.5`.
#' @param n.threads Integer, OpenMP threads for the compiled pair estimator and
#'   the per-event library assignment. Default `1`.
#' @param figures Logical, whether to write the spillover heatmap. Default
#'   `TRUE`.
#' @param save Logical, whether to write the csv outputs. Default `TRUE`.
#' @param true.spectra Optional numeric matrix (fluorophores x detectors),
#'   independently known ground truth with row names matching `spectra`.
#'   Purely diagnostic: the returned `recovery` table reports the angular
#'   error against it before and after this run.
#' @param min.deg.start Numeric, degrees. Below this starting angular error,
#'   `recovered` is reported as `0`. Default `0.1`.
#' @param verbose Logical, controls messaging. Default `TRUE`.
#' @param keep.history Logical, whether to retain a per-iteration snapshot of
#'   phase one in `spillover.history`. Default `FALSE`.
#'
#' @return A named list:
#' \describe{
#'   \item{`spectra`}{The corrected spectra, every row of `spectra`. Rows that
#'     failed a gate are unchanged.}
#'   \item{`spectra.backsolved`}{The back-solved rows of phase one (of the last
#'     round). Diagnostic only: they reproduce the compensation exactly but are
#'     not constrained to be physical spectra.}
#'   \item{`spillover`, `compensation`, `trust`}{The residual spillover matrix
#'     (fluorophores x fluorophores), its inverse, and the trust weights of
#'     the last iteration.}
#'   \item{`coefficient.log`}{Per-pair fit diagnostics from the final
#'     iteration.}
#'   \item{`source.log`}{Per-source diagnostics from the final iteration: the
#'     negative-population bound, the clearly positive events, the number the
#'     unstained rate accounts for, the excess, the resulting weight, the AF
#'     likeness of the row change the step implied, whether the AF gate
#'     dropped it, and the share of its back-solved row's change outside the
#'     plausibility box (`plausible.project`).}
#'   \item{`signature.log`}{Per-fluorophore signature statistics and gate
#'     outcomes: `candidate` (`"remeasured"`, `"backsolved"` or `"parent"`
#'     for an adopted row), `rescued.from` (the rejection a later candidate
#'     replaced), `tandem.parent`, `tandem.coefficient`, `implausible.frac`.
#'     With smoothing, the last round's, with `rounds.accepted` and every row
#'     the rounds moved marked accepted.}
#'   \item{`convergence.log`, `stop.reason`}{Per-iteration delta history and
#'     why phase one stopped.}
#'   \item{`smooth.log`}{`NULL` without smoothing; else one row per round.}
#'   \item{`spillover.history`}{`NULL` unless `keep.history = TRUE`.}
#'   \item{`unmixed.final`, `residual.final`}{The fully stained events'
#'     abundance compensated with the final `compensation`, and their
#'     background-subtracted raw data; with smoothing, the last round's
#'     events.}
#'   \item{`thresholds.final`, `threshold.matrix.final`, `dominant.final`}{The
#'     final flat and spread-scaled positivity boundaries and each event's
#'     winning fluorophore (`0` for none).}
#'   \item{`neg.thresholds.final`, `neg.threshold.matrix.final`}{The negative
#'     boundaries the hypernegative check used, measured from the unstained
#'     population's own lower tail. Use these, not `-threshold.matrix.final`,
#'     to audit the result against a negative boundary.}
#'   \item{`af.basis`, `af.hotspot`, `af.frozen`}{The background rows of the
#'     design, their coupling to the panel, and the fluorophores frozen because
#'     of it.}
#'   \item{`af.assignment`}{`NULL` unless `bg.mode = "per.cell"`, in which case
#'     how often each `af.spectra` row was assigned and its mean abundance.}
#'   \item{`leakage.prior`}{The variant leakage prior, or `NULL`.}
#'   \item{`recovery`}{Angular errors against `true.spectra` at the start, of
#'     phase one's back-solved rows (single pass only) and after; `recovered`
#'     is `(deg.start - deg.after) / deg.start`. `NULL` without
#'     `true.spectra`.}
#' }
#'
#' @export

fix.my.unmix <- function(
    spectra,
    unstained.sample,
    fully.stained.sample,
    flow.control,
    asp,
    variants                  = NULL,
    af.name                   = "AF",
    af.basis                  = NULL,
    af.n.pc                   = "auto",
    af.spectra                = NULL,
    bg.mode                   = c( "per.cell", "af.deconv", "af.row", "global.mean", "none" ),
    large.gate                = TRUE,
    scatter.gate              = FALSE,
    landmark.quantile         = NULL,
    max.events                = 400000L,
    max.iter                  = 20L,
    downsample                = FALSE,
    downsample.background.frac = 0.3,
    downsample.min.stratum    = 2000L,
    unstained.threshold       = 0.95,
    lower.unstained.threshold = 0.99,
    unstained.margin          = 1.3,
    spread.kappa              = 3,
    negative.tail             = 0.001,
    source.events             = c( 200, 400 ),
    envelope.quantiles        = c( 0.05, 0.5 ),
    disagreement.reference    = c( "envelope.corrected", "envelope", "anchored" ),
    envelope.scale.quantile   = 0.25,
    min.negative.frac         = 0.10,
    max.disagreement          = 0.5,
    min.hypernegative.frac    = 0.01,
    max.hypernegative.frac    = 0.05,
    min.hypernegative.events  = 200L,
    hypernegative.dominant    = TRUE,
    leakage.prior             = TRUE,
    span.fraction             = 0.6,
    min.negative.events       = 200L,
    min.bin.negative          = 25L,
    min.span                  = 5,
    min.rise                  = 1,
    n.levels                  = 60L,
    min.bin.events            = 50L,
    n.levels.pair             = 10L,
    multivariate              = TRUE,
    ridge                     = 1e-6,
    max.truncated.events      = 20000L,
    max.mask.passes           = 3L,
    estimator                 = c( "anchored", "truncated", "envelope" ),
    glasso.check              = TRUE,
    glasso.warm               = TRUE,
    glasso.n.lambda           = 40L,
    glasso.lambda.min.ratio   = 1e-3,
    glasso.margin.frac        = 0.05,
    glasso.tolerance          = 1e-4,
    af.gate                   = FALSE,
    af.gate.max               = 0.8,
    plausibility.ratio        = 1,
    plausibility.floor        = 0.003,
    plausibility.reject       = 0.5,
    plausibility.reset        = TRUE,
    phase.two.rows            = c( "flagged", "all" ),
    phase.two.active          = c( "present", "all" ),
    phase.two.present         = c( frac = 0.02, corr = 0.3, mean = 0.5 ),
    phase.two.fallback        = c( "remeasured", "backsolved", "reject" ),
    phase.two.rescue          = c( "backsolved", "off" ),
    rescue.reasons            = c( "leakage", "conditioning", "no.fit", "not.joint" ),
    rescue.min.weight         = 0.5,
    leakage.channels          = c( "negative", "all" ),
    tandem.candidate          = TRUE,
    tandem.min.coefficient    = 0.01,
    source.dominant           = TRUE,
    spread.addback            = FALSE,
    anchor.weight             = 1,
    max.coefficient           = 0.5,
    convergence.threshold     = 0.003,
    convergence.quantile      = 1,
    stop.rule                 = c( "settled", "all" ),
    stop.trust                = 0.05,
    stop.confirm              = 2L,
    coefficient.decay         = 0,
    update.spectra            = TRUE,
    step                      = 1,
    intercept                 = TRUE,
    min.explained             = 0.8,
    max.explained             = 1.2,
    max.resid                 = 0.03,
    max.intercept             = 0.03,
    min.bg.align              = -0.9,
    max.angle                 = 30,
    max.clamp.frac            = 0.15,
    max.anchor                = 0.10,
    max.vif                   = 500,
    max.condition.increase    = 1.05,
    peak.shift.min.rel        = 0.7,
    max.hotspot               = 5,
    leakage.margin            = 0.05,
    smooth.updates            = TRUE,
    smooth.rounds             = 5L,
    smooth.fraction           = 0.5,
    smooth.mode               = c( "damped", "consensus" ),
    consensus.min.rounds      = 0.6,
    consensus.max.spread      = 0.5,
    n.threads                 = 1L,
    figures                   = TRUE,
    save                      = TRUE,
    true.spectra              = NULL,
    min.deg.start             = 0.1,
    verbose                   = TRUE,
    keep.history              = FALSE
) {

  bg.mode                <- match.arg( bg.mode )
  estimator              <- match.arg( estimator )
  disagreement.reference <- match.arg( disagreement.reference )
  phase.two.rows         <- match.arg( phase.two.rows )
  phase.two.active       <- match.arg( phase.two.active )
  phase.two.fallback     <- match.arg( phase.two.fallback )
  phase.two.rescue       <- match.arg( phase.two.rescue )
  leakage.channels       <- match.arg( leakage.channels )
  stop.rule              <- match.arg( stop.rule )
  smooth.mode            <- match.arg( smooth.mode )

  if ( length( source.events ) != 2 || !( source.events[ 1 ] >= 0 ) ||
       !( source.events[ 2 ] > source.events[ 1 ] ) )
    stop( "`source.events` must be a pair with 0 <= low < high.", call. = FALSE )

  if ( !( negative.tail > 0 && negative.tail < 0.5 ) )
    stop( "`negative.tail` must lie between 0 and 0.5.", call. = FALSE )

  if ( length( envelope.quantiles ) != 2 )
    stop( "`envelope.quantiles` must be a pair.", call. = FALSE )

  if ( !( unstained.threshold > 0.5 && unstained.threshold < 1 &&
          lower.unstained.threshold > 0.5 && lower.unstained.threshold < 1 ) )
    stop( "The unstained quantiles must lie between 0.5 and 1.", call. = FALSE )

  if ( plausibility.ratio < 0 || plausibility.floor < 0 ||
       !( plausibility.reject > 0 && plausibility.reject <= 1 ) )
    stop( paste0( "The plausibility box needs `plausibility.ratio` >= 0, ",
                  "`plausibility.floor` >= 0 and 0 < `plausibility.reject` <= 1." ),
          call. = FALSE )

  if ( !all( c( "frac", "corr", "mean" ) %in% names( phase.two.present ) ) )
    stop( "`phase.two.present` must name `frac`, `corr` and `mean`.", call. = FALSE )

  if ( !( rescue.min.weight >= 0 && rescue.min.weight <= 1 ) )
    stop( "`rescue.min.weight` must lie between 0 and 1.", call. = FALSE )

  if ( stop.rule == "settled" &&
       ( !( stop.trust > 0 && stop.trust <= 1 ) || stop.confirm < 1 ) )
    stop( "`stop.rule = \"settled\"` needs 0 < `stop.trust` <= 1 and `stop.confirm` >= 1.",
          call. = FALSE )

  if ( smooth.rounds < 1 || !( smooth.fraction > 0 && smooth.fraction <= 1 ) )
    stop( "Smoothing needs `smooth.rounds` >= 1 and 0 < `smooth.fraction` <= 1.",
          call. = FALSE )

  if ( !( consensus.min.rounds > 0 && consensus.min.rounds <= 1 ) ||
       consensus.max.spread < 0 )
    stop( "Consensus needs 0 < `consensus.min.rounds` <= 1 and `consensus.max.spread` >= 0.",
          call. = FALSE )

  if ( hypernegative.dominant && !source.dominant )
    stop( "`hypernegative.dominant` requires `source.dominant`.", call. = FALSE )

  seed <- if ( is.null( asp$bird.seed ) ) 0L else as.integer( asp$bird.seed )

  spectra <- as.matrix( spectra )

  if ( is.null( rownames( spectra ) ) )
    stop( "`spectra` must have fluorophore row names.", call. = FALSE )

  fluorophores  <- setdiff( rownames( spectra ), af.name )
  fluorophore.n <- length( fluorophores )

  if ( fluorophore.n < 2 )
    stop( "At least two panel fluorophores are required.", call. = FALSE )

  if ( multivariate && update.spectra && n.levels < fluorophore.n + 3L )
    stop( sprintf( paste0( "The joint signature fit needs `n.levels` of at least ",
                           "%d for %d fluorophores." ),
                   fluorophore.n + 3L, fluorophore.n ), call. = FALSE )

  if ( bg.mode == "per.cell" ) {

    if ( is.null( af.spectra ) )
      stop( paste0( "`bg.mode = \"per.cell\"` requires `af.spectra`; ",
                    "prepare using `get.af.spectra()`." ), call. = FALSE )

    if ( nrow( as.matrix( af.spectra ) ) < 2 )
      stop( "`af.spectra` must have at least two rows.", call. = FALSE )
  }

  if ( !is.null( af.spectra ) ) {

    af.spectra <- as.matrix( af.spectra )

    if ( is.null( rownames( af.spectra ) ) )
      rownames( af.spectra ) <- paste0( "AF", seq_len( nrow( af.spectra ) ) )

    if ( !identical( colnames( af.spectra ), colnames( spectra ) ) )
      stop( "`af.spectra` detectors do not match `spectra`.", call. = FALSE )
  }

  if ( save && ! dir.exists( asp$fix.unmixing.dir ) )
    dir.create( asp$fix.unmixing.dir, recursive = TRUE )

  spillover.spread <- variants$spillover.spread

  if ( is.null( spillover.spread ) )
    warning( paste0( "No `spillover.spread` in `variants`; positivity ",
                     "boundaries will be flat, so the hypernegative check ",
                     "removes the trust of most coefficients and large positive ",
                     "errors are hidden from the pair fit." ),
             call. = FALSE )

  if ( leakage.prior && is.null( variants ) )
    warning( paste0( "`leakage.prior` requires `variants`; falling back to ",
                     "the estimate's own variance." ), call. = FALSE )

  # ---------------------------------------------------------------------------
  # Read and gate
  # ---------------------------------------------------------------------------

  read.gated <- function( file.name, gate.polygon = NULL,
                          landmark.threshold = NULL, label,
                          cap = Inf, cap.seed = seed ) {

    if ( verbose )
      message( sprintf( "\033[34mReading %s.\033[0m", label ) )

    expr.data <- readFCS( file.name, columns = flow.control$scatter.and.channel.spectral )
    gate.data <- expr.data[ , flow.control$scatter.parameter ]

    if ( !scatter.gate ) {

      keep <- seq_len( nrow( expr.data ) )
      gate.polygon       <- NULL
      landmark.threshold <- NULL

    } else {

      if ( is.null( gate.polygon ) )
        gate.polygon <- do.gate(
          gate.data, viability.gate = FALSE, large.gate = large.gate,
          samp = label,
          scatter.and.channel.label = flow.control$scatter.and.channel.label,
          control.type = "cells", asp )

      if ( is.null( landmark.threshold ) && !is.null( landmark.quantile ) )
        landmark.threshold <- apply(
          gate.data, 2, stats::quantile, probs = landmark.quantile, names = FALSE )

      in.polygon <- .point.in.polygon(
        gate.data[ , 1 ], gate.data[ , 2 ], gate.polygon$x, gate.polygon$y ) != 0

      # A scatter-only rescue, applied identically to the unstained file that
      # defines `landmark.threshold` and the stained files that reuse it, so the
      # background and the correction it feeds come from the same population.
      in.landmark <- if ( is.null( landmark.threshold ) ) FALSE else
        gate.data[ , 1 ] > landmark.threshold[ 1 ] |
        gate.data[ , 2 ] > landmark.threshold[ 2 ]

      keep <- which( in.polygon | in.landmark )
    }

    if ( length( keep ) > cap ) {
      set.seed( cap.seed )
      keep <- sort( keep[ sample.int( length( keep ), cap ) ] )
    }

    list( data = expr.data[ keep, flow.control$spectral.channel, drop = FALSE ],
          gate = gate.polygon,
          landmark = landmark.threshold )
  }

  event.cap <- if ( is.null( max.events ) ) Inf else max.events

  stained.files <- as.character( fully.stained.sample )
  stained.share <- ceiling( event.cap / length( stained.files ) )

  unstained.in <- read.gated( unstained.sample, NULL, NULL, "unstained raw",
                              cap = event.cap, cap.seed = seed )

  stained.raw <- do.call( rbind, lapply( seq_along( stained.files ), function( i )
    read.gated( stained.files[ i ], unstained.in$gate, unstained.in$landmark,
                sprintf( "fully stained raw %d of %d", i, length( stained.files ) ),
                cap = stained.share, cap.seed = seed + i )$data ) )

  unstained.raw <- unstained.in$data

  if ( ncol( stained.raw ) != ncol( spectra ) )
    stop( "`spectra` and the spectral channels of the data do not match.",
          call. = FALSE )

  # ---------------------------------------------------------------------------
  # Projection and the run
  # ---------------------------------------------------------------------------

  projector <- .fix.projector(
    unstained.raw = unstained.raw,
    stained.raw   = stained.raw,
    spectra       = spectra,
    bg.mode       = bg.mode,
    af.name       = af.name,
    af.basis      = af.basis,
    af.n.pc       = af.n.pc,
    af.spectra    = af.spectra,
    n.threads     = n.threads,
    verbose       = verbose )

  spectra.fluor <- spectra[ fluorophores, , drop = FALSE ]

  tandems <- if ( tandem.candidate ) .fix.tandem.pairs( fluorophores ) else
    stats::setNames( character( 0 ), character( 0 ) )

  if ( verbose && length( tandems ) > 0 )
    message( sprintf( "\033[34mTandems with their parent in the panel: %s.\033[0m",
                      paste( names( tandems ), tandems, sep = " < ", collapse = ", " ) ) )

  core.data <- c( "spectra.fluor", "unstained.fit", "stained.fit",
                  "background.basis", "af.gate.spectra", "freeze.coupled",
                  "spillover.spread", "variants", "reference.spectra",
                  "tandems", "seed" )
  core.options <- mget( setdiff( names( formals( .fix.my.unmix.core ) ), core.data ) )

  run.core <- function( start, projection, round.seed ) {
    do.call( .fix.my.unmix.core, c( list(
      spectra.fluor     = start,
      unstained.fit     = projection$unstained,
      stained.fit       = projection$stained,
      background.basis  = projector$background.basis,
      af.gate.spectra   = projector$af.gate.spectra,
      freeze.coupled    = bg.mode != "per.cell",
      spillover.spread  = spillover.spread,
      variants          = variants,
      reference.spectra = spectra.fluor,
      tandems           = tandems,
      seed              = round.seed ), core.options ) )
  }

  smooth <- smooth.updates && smooth.rounds > 1 && update.spectra

  if ( smooth ) {

    smoothed <- .fix.smooth(
      spectra.fluor        = spectra.fluor,
      project              = projector$project,
      n.stained            = nrow( stained.raw ),
      run.core             = run.core,
      rounds               = as.integer( smooth.rounds ),
      fraction             = smooth.fraction,
      mode                 = smooth.mode,
      consensus.min.rounds = consensus.min.rounds,
      consensus.max.spread = consensus.max.spread,
      plausibility.ratio   = plausibility.ratio,
      plausibility.floor   = plausibility.floor,
      plausibility.reset   = plausibility.reset,
      seed                 = seed,
      verbose              = verbose )

    result     <- smoothed$result
    projection <- smoothed$projection

  } else {

    set.seed( seed )
    projection <- projector$project( spectra.fluor )
    result     <- run.core( spectra.fluor, projection, seed )
  }

  # ---------------------------------------------------------------------------
  # Recovery against ground truth, outputs
  # ---------------------------------------------------------------------------

  recovery <- NULL

  if ( !is.null( true.spectra ) ) {

    common <- intersect( fluorophores, rownames( true.spectra ) )

    if ( length( common ) > 0 ) {

      truth <- true.spectra[ common, colnames( spectra ), drop = FALSE ]

      deg.start  <- .signature.row.angle( spectra.fluor[ common, , drop = FALSE ], truth )
      deg.phase1 <- if ( smooth ) rep( NA_real_, length( common ) ) else
        .signature.row.angle( result$spectra.backsolved[ common, , drop = FALSE ], truth )
      deg.after  <- .signature.row.angle( result$spectra[ common, , drop = FALSE ], truth )

      recovery <- data.frame(
        fluorophore = common,
        deg.start   = deg.start,
        deg.phase1  = deg.phase1,
        deg.after   = deg.after,
        recovered   = ifelse( deg.start < min.deg.start, 0,
                              ( deg.start - deg.after ) / deg.start ),
        accepted    = result$accepted[ common ],
        row.names   = NULL )
    }
  }

  af.assignment <- NULL

  if ( bg.mode == "per.cell" ) {

    af.index     <- projection$stained$af.index
    af.abundance <- projection$stained$af.abundance

    tab <- table( af.index )
    idx <- as.integer( names( tab ) )

    af.assignment <- data.frame(
      af.name        = rownames( af.spectra )[ idx ],
      n.events       = as.integer( tab ),
      mean.abundance = vapply( idx, function( j )
        mean( af.abundance[ af.index == j ] ), numeric( 1 ) ),
      row.names        = NULL,
      stringsAsFactors = FALSE )

    af.assignment <- af.assignment[ order( -af.assignment$n.events ), ]

    if ( save )
      utils::write.csv( af.assignment, row.names = FALSE,
                        file = file.path( asp$fix.unmixing.dir,
                                          "FixMyUnmix_af_assignment.csv" ) )
  }

  spectra.new <- spectra
  spectra.new[ fluorophores, ] <- result$spectra[ fluorophores, , drop = FALSE ]

  if ( figures )
    spectral.heatmap(
      spectra = result$spillover,
      title = "FixMyUnmix_spillover_heatmap",
      plot.dir = asp$fix.unmixing.dir,
      legend.label = "Spillover",
      save = TRUE
    )

  if ( save ) {

    spillover.filename <- if ( is.null( asp$fix.spillover.filename ) )
      "FixMyUnmix_spillover.csv" else asp$fix.spillover.filename

    compensation.filename <- if ( is.null( asp$fix.compensation.filename ) )
      "FixMyUnmix_compensation.csv" else asp$fix.compensation.filename

    spectra.filename <- if ( is.null( asp$fix.spectra.filename ) )
      "FixMyUnmix_spectra.csv" else asp$fix.spectra.filename

    utils::write.csv( result$spillover,
                      file = file.path( asp$fix.unmixing.dir, spillover.filename ) )
    utils::write.csv( result$compensation,
                      file = file.path( asp$fix.unmixing.dir, compensation.filename ) )
    utils::write.csv( spectra.new,
                      file = file.path( asp$fix.unmixing.dir, spectra.filename ) )
  }

  list(
    spectra                    = spectra.new,
    spectra.backsolved         = result$spectra.backsolved,
    spillover                  = result$spillover,
    compensation               = result$compensation,
    trust                      = result$trust,
    coefficient.log            = result$coefficient.log,
    source.log                 = result$source.log,
    signature.log              = result$signature.log,
    convergence.log            = result$convergence.log,
    stop.reason                = result$stop.reason,
    smooth.log                 = result$smooth.log,
    spillover.history          = result$spillover.history,
    unmixed.final              = result$unmixed.final,
    residual.final             = result$residual.final,
    thresholds.final           = result$thresholds.final,
    threshold.matrix.final     = result$threshold.matrix.final,
    neg.thresholds.final       = result$neg.thresholds.final,
    neg.threshold.matrix.final = result$neg.threshold.matrix.final,
    dominant.final             = result$dominant.final,
    af.basis                   = result$af.basis,
    af.hotspot                 = result$af.hotspot,
    af.frozen                  = result$af.frozen,
    af.assignment              = af.assignment,
    leakage.prior              = result$leakage.prior,
    recovery                   = recovery
  )
}


#' Projection of the unstained and fully stained events onto the panel, with
#' the background removed.
#'
#' The background rows that stay fixed for a run are built once here. Under
#' `bg.mode = "per.cell"` each event takes the `af.spectra` row whose removal
#' leaves the smallest total apparent fluorophore signal (`.fix.assign.af()`),
#' and only that row is subtracted; the library's mean row (row 1) is the
#' background row of the design. Every other mode fits its background rows
#' jointly with the panel by least squares and subtracts their non-negative
#' part.
#'
#' @return A list: `project( spectra.fluor, keep = NULL )`, which unmixes the
#'   unstained events and the stained events `keep` (all when `NULL`) with
#'   `spectra.fluor` and returns `list( unstained, stained )`, each a list of
#'   `abundance` (events x fluorophores, uncompensated), `residual` (events x
#'   detectors, background removed) and, under per.cell, `af.index` and
#'   `af.abundance`; `background.basis`, the background rows of the design
#'   (`NULL` for none); and `af.gate.spectra`, the spectra the AF likeness of a
#'   row change is measured against.
#'
#' @noRd
.fix.projector <- function( unstained.raw, stained.raw, spectra, bg.mode,
                            af.name, af.basis, af.n.pc, af.spectra,
                            n.threads = 1L, verbose = TRUE ) {

  background.basis <- switch(
    bg.mode,
    per.cell = af.spectra[ 1, , drop = FALSE ],
    af.deconv = {
      if ( is.null( af.basis ) )
        af.basis <- get.af.basis( unstained.raw, n.pc = af.n.pc,
                                  verbose = verbose )
      as.matrix( af.basis )
    },
    af.row = {
      if ( is.null( af.name ) || ! af.name %in% rownames( spectra ) )
        stop( "`bg.mode = \"af.row\"` requires an `af.name` row in `spectra`.",
              call. = FALSE )
      spectra[ af.name, , drop = FALSE ]
    },
    global.mean = {
      bg <- colMeans( unstained.raw )
      matrix( bg / max( sqrt( sum( bg^2 ) ), .Machine$double.eps ), nrow = 1,
              dimnames = list( "BG", colnames( unstained.raw ) ) )
    },
    none = NULL
  )

  af.gate.spectra <- if ( !is.null( af.spectra ) ) af.spectra else background.basis

  fit <- function( y, spectra.fluor, af.library ) {

    fluorophores <- rownames( spectra.fluor )

    if ( bg.mode == "per.cell" ) {

      af.assign  <- .fix.assign.af( y, spectra.fluor, af.spectra,
                                    af.library$r, af.library$dots, n.threads )
      background <- af.assign$abundance * af.spectra[ af.assign$index, , drop = FALSE ]

      coefs <- unmix.ols.fast( y - background, spectra.fluor )
      colnames( coefs ) <- fluorophores

      return( list( abundance    = coefs,
                    residual     = y - background,
                    af.index     = af.assign$index,
                    af.abundance = af.assign$abundance ) )
    }

    design <- if ( is.null( background.basis ) ) spectra.fluor else
      rbind( background.basis, spectra.fluor )

    coefs <- unmix.ols.fast( y, design )
    colnames( coefs ) <- rownames( design )

    background <- 0
    if ( !is.null( background.basis ) ) {
      b <- coefs[ , rownames( background.basis ), drop = FALSE ]
      b[ b < 0 ] <- 0
      background <- b %*% background.basis
    }

    list( abundance = coefs[ , fluorophores, drop = FALSE ],
          residual  = y - background )
  }

  project <- function( spectra.fluor, keep = NULL ) {

    # How each library row reads in fluorophore space and the part of it the
    # panel cannot explain, for the per-event assignment. Both depend on the
    # spectra, so they follow the spectra each call is made with.
    af.library <- NULL

    if ( bg.mode == "per.cell" ) {
      unmixing <- solve.default( tcrossprod( spectra.fluor ), spectra.fluor )
      r <- t( af.spectra ) - t( spectra.fluor ) %*% ( unmixing %*% t( af.spectra ) )
      af.library <- list( r = r, dots = pmax( colSums( r^2 ), .Machine$double.eps ) )
    }

    y <- if ( is.null( keep ) ) stained.raw else stained.raw[ keep, , drop = FALSE ]

    list( unstained = fit( unstained.raw, spectra.fluor, af.library ),
          stained   = fit( y, spectra.fluor, af.library ) )
  }

  list( project          = project,
        background.basis = background.basis,
        af.gate.spectra  = af.gate.spectra )
}


#' Residual spillover slope for one source-target pair, by two estimators.
#'
#' `slope.truncated` is a robust line through the target-negative events, with
#' no binning and no zero-abundance anchor. Negativity is decided on the target
#' with this pair's current spillover estimate already removed, and the two are
#' iterated together. Deciding it on the raw target caps the fit's source range
#' at threshold over coefficient: the larger the residual spillover, the more of
#' the bright source population it pushes across the target's boundary, until
#' the events carrying the answer are the ones excluded and the pair reports
#' nothing at all. The estimator is then least able to measure exactly the
#' coefficients that matter most.
#'
#' `slope` is the binned envelope. Source-positive events are binned by source
#' abundance, everything below the source threshold is collapsed into a single
#' zero-abundance anchor, and in each bin the envelope is a fixed quantile of
#' the target-negative events only. It is retained as an alternative and for the
#' statistics the binning produces. When the population cannot support the
#' binning, the truncated estimate is still returned; the binned quantities come
#' back as `NA`.
#'
#' `coverage` is the fraction of target-negative events in the brightest third
#' of the source-positive population, measured on the corrected target, and is
#' the identifiability test for co-expression. A target positive only because
#' the source spills into it returns to the negative set once the estimate is
#' removed; a target the cell genuinely expresses does not, so when `coverage`
#' collapses the pair carries no information about spillover whichever slope is
#' taken. `slope.alt` is fitted to the plain median of the whole bin, so
#' disagreement between the two reports co-expression rather than instability.
#' @noRd
.fix.envelope.slope <- function( x.source, x.target,
                                 threshold.source, threshold.target,
                                 spread.var, neg.var,
                                 source.mask          = NULL,
                                 quantiles            = c( 0.05, 0.5 ),
                                 n.levels             = 10L,
                                 min.events           = 200L,
                                 min.bin.negative     = 25L,
                                 spread.addback       = TRUE,
                                 anchor.weight        = 1,
                                 max.truncated.events = 20000L,
                                 max.coefficient      = 0.2,
                                 max.mask.passes      = 3L,
                                 mask.tolerance       = 0.05,
                                 start.slope          = 0 ) {

  n <- length( x.source )
  if ( n < min.events ) return( NULL )

  if ( length( threshold.source ) == 1L )
    threshold.source <- rep( threshold.source, n )
  if ( length( threshold.target ) == 1L )
    threshold.target <- rep( threshold.target, n )

  # An event that clears the source threshold while another fluorophore is
  # brighter carries that other fluorophore's abundance read through a spectral
  # error, not this source's. It corrupts both estimators the same way, so the
  # filter is applied before either: for the robust line it removes exactly the
  # high-leverage points a wrong dye supplies, which is how a dim source
  # acquires a coefficient of order one.
  if ( !is.null( source.mask ) ) {

    keep.event <- source.mask | x.source <= threshold.source

    x.source         <- x.source[ keep.event ]
    x.target         <- x.target[ keep.event ]
    threshold.source <- threshold.source[ keep.event ]
    threshold.target <- threshold.target[ keep.event ]

    n <- length( x.source )
    if ( n < min.events ) return( NULL )
  }

  # Everything above the source threshold is kept and the negative bulk below it
  # is subsampled. The bulk sits at the origin and carries no leverage, so the
  # slope is unchanged while the robust fit's cost stops scaling with file size.
  select.negative <- function( slope ) {

    index <- which( x.target - slope * x.source < threshold.target )

    if ( length( index ) <= max.truncated.events ) return( index )

    bright <- index[ x.source[ index ] > threshold.source[ index ] ]
    bulk   <- setdiff( index, bright )
    n.bulk <- max( max.truncated.events - length( bright ), min.events )

    if ( length( bulk ) > n.bulk ) bulk <- sample( bulk, n.bulk )

    c( bright, bulk )
  }

  fit.truncated <- function( index, start = NULL ) {

    if ( length( index ) < min.events ) return( NA_real_ )

    .fix.huber.slope( x.source[ index ], x.target[ index ], start = start )[ 2 ]
  }

  start.slope <- if ( is.finite( start.slope ) &&
                      abs( start.slope ) <= max.coefficient ) start.slope else 0

  truncated.index <- select.negative( start.slope )
  slope.truncated <- fit.truncated( truncated.index, start = c( 0, start.slope ) )

  for ( pass in seq_len( as.integer( max.mask.passes ) ) ) {

    if ( !is.finite( slope.truncated ) ||
         abs( slope.truncated ) > max.coefficient ) break

    index.next <- select.negative( slope.truncated )
    slope.next <- fit.truncated( index.next, start = c( 0, slope.truncated ) )

    if ( !is.finite( slope.next ) ) break

    settled <- abs( slope.next - slope.truncated ) <=
      mask.tolerance * abs( slope.next )

    slope.truncated <- slope.next
    truncated.index <- index.next

    if ( settled ) break
  }

  slope.mask <- if ( is.finite( slope.truncated ) &&
                     abs( slope.truncated ) <= max.coefficient )
    slope.truncated else 0

  target.negative <- x.target - slope.mask * x.source < threshold.target

  source.positive <- x.source > threshold.source

  coverage <- 1

  if ( sum( source.positive ) >= min.bin.negative ) {

    bright.index <- which( source.positive )
    bright.cut   <- stats::quantile( x.source[ bright.index ], probs = 2 / 3,
                                     names = FALSE )
    top.index    <- bright.index[ x.source[ bright.index ] >= bright.cut ]

    if ( length( top.index ) >= min.bin.negative )
      coverage <- mean( target.negative[ top.index ] )
  }

  span.truncated <- if ( length( truncated.index ) > 1L )
    diff( range( x.source[ truncated.index ] ) ) else 0

  # The truncated estimator needs none of the binning, so a population that
  # cannot support the bins must not take the whole pair down with it.
  truncated.only <- function() {

    if ( !is.finite( slope.truncated ) ) return( NULL )

    list( slope           = NA_real_,
          slope.truncated = slope.truncated,
          se              = NA_real_,
          slope.alt       = NA_real_,
          disagreement    = NA_real_,
          coverage        = coverage,
          n               = n,
          span            = span.truncated,
          noise           = sqrt( max( neg.var, 0 ) ) )
  }

  if ( sum( source.positive ) < min.events ) return( truncated.only() )
  if ( sum( !source.positive ) < min.bin.negative ) return( truncated.only() )

  brk <- unique( stats::quantile(
    x.source[ source.positive ],
    probs = seq( 0, 1, length.out = n.levels + 1 ), names = FALSE ) )
  if ( length( brk ) < 4 ) return( truncated.only() )

  bin <- integer( n )
  bin[ source.positive ] <- as.integer( cut(
    x.source[ source.positive ], breaks = brk, include.lowest = TRUE ) )

  bins <- sort( unique( bin ) )
  if ( length( bins ) < 5 ) return( truncated.only() )

  # Each bin's row indices are split out once and reused for every per-bin
  # statistic below, rather than re-scanning `bin == b` over all n events
  # once per bin for each of centre / neg.n / bin.n / envelope.value /
  # median.value in turn.
  idx.by.bin <- split( seq_len( n ), bin )
  idx.by.bin <- idx.by.bin[ as.character( bins ) ]

  centre <- unname( vapply( idx.by.bin, function( idx )
    stats::median( x.source[ idx ] ), numeric( 1 ) ) )
  neg.n <- unname( vapply( idx.by.bin, function( idx )
    sum( target.negative[ idx ] ), integer( 1 ) ) )
  bin.n <- unname( lengths( idx.by.bin ) )

  usable <- neg.n >= min.bin.negative
  if ( sum( usable ) < 4L ) return( truncated.only() )

  sd.bin <- sqrt( pmax( neg.var + spread.var * pmax( centre, 0 ), 0 ) )
  z      <- if ( spread.addback ) -stats::qnorm( quantiles ) else c( 0, 0 )

  # Bins are weighted by the precision of their own quantile, events over
  # variance, and the zero-abundance anchor is capped at the best positive bin.
  # Uncapped, the anchor holds most of the file and sets the slope on its own,
  # which turns the fit into a comparison between two different populations
  # rather than a trend within one.
  anchor.cap <- function( w ) {
    if ( !any( bins > 0L ) || !any( bins == 0L ) ) return( w )
    w[ bins == 0L ] <- pmin( w[ bins == 0L ],
                             anchor.weight * max( w[ bins > 0L ] ) )
    w
  }

  weight.envelope <- anchor.cap( neg.n / pmax( sd.bin^2, .Machine$double.eps ) )
  weight.compare  <- anchor.cap( bin.n / pmax( sd.bin^2, .Machine$double.eps ) )

  fit.trace <- function( value, weight ) {

    keep <- usable & is.finite( value )
    if ( sum( keep ) < 4L ) return( list( slope = NA_real_, se = NA_real_ ) )

    x.bin <- centre[ keep ]
    y.bin <- value[ keep ]
    w.bin <- weight[ keep ]

    # Closed-form weighted least squares for the two-parameter fit, in place
    # of stats::lm.wfit(). At n.levels.pair bins a QR decomposition buys
    # nothing here, and dropping it removes this function's only call into
    # compiled linear algebra, which matters both for its own overhead, run
    # tens of thousands of times per correction, and for safety under
    # mclapply on macOS, where forking a process that has touched
    # Accelerate's threaded BLAS is unreliable.
    x.mean <- sum( w.bin * x.bin ) / sum( w.bin )
    y.mean <- sum( w.bin * y.bin ) / sum( w.bin )
    sxx    <- sum( w.bin * ( x.bin - x.mean )^2 )
    sxy    <- sum( w.bin * ( x.bin - x.mean ) * ( y.bin - y.mean ) )

    slope     <- if ( sxx > 0 ) sxy / sxx else NA_real_
    intercept <- y.mean - slope * x.mean
    residuals <- y.bin - ( intercept + slope * x.bin )

    dof <- sum( keep ) - 2L

    se <- if ( is.finite( slope ) && dof > 0 && sxx > 0 )
      sqrt( sum( w.bin * residuals^2 ) / dof / sxx ) else NA_real_

    list( slope = slope, se = se )
  }

  envelope.value <- unname( vapply( idx.by.bin, function( idx ) {
    v <- x.target[ idx ][ target.negative[ idx ] ]
    if ( length( v ) < min.bin.negative ) return( NA_real_ )
    stats::quantile( v, probs = quantiles[ 1 ], names = FALSE )
  }, numeric( 1 ) ) )

  median.value <- unname( vapply( idx.by.bin, function( idx )
    stats::quantile( x.target[ idx ], probs = quantiles[ 2 ], names = FALSE ),
    numeric( 1 ) ) )

  envelope <- fit.trace( envelope.value + z[ 1 ] * sd.bin, weight.envelope )
  compare  <- fit.trace( median.value   + z[ 2 ] * sd.bin, weight.compare )

  if ( !is.finite( envelope$slope ) && !is.finite( slope.truncated ) )
    return( NULL )

  disagreement <- abs( envelope$slope - compare$slope ) /
    ( abs( envelope$slope ) + abs( compare$slope ) + .Machine$double.eps )

  list( slope           = envelope$slope,
        slope.truncated = slope.truncated,
        se              = envelope$se,
        slope.alt       = compare$slope,
        disagreement    = disagreement,
        coverage        = coverage,
        n               = n,
        span            = max( centre ) - min( centre ),
        noise           = sqrt( max( neg.var, 0 ) ) )
}


#' An empty statistics block matching what `extract.raw.signature()` reports,
#' used for fluorophores whose population could not support a fit so that every
#' row of the signature log carries the same columns.
#' @noRd
.fix.empty.signature.stats <- function() {
  data.frame(
    n.bins           = NA_integer_,
    n.active         = NA_integer_,
    joint            = NA,
    x.span           = NA_real_,
    explained        = NA_real_,
    explained.total  = NA_real_,
    resid.rel        = NA_real_,
    intercept.rel    = NA_real_,
    bg.align         = NA_real_,
    clamp.frac       = NA_real_,
    anchor.rel       = NA_real_,
    vif.target       = NA_real_,
    vif.max          = NA_real_,
    vif.partner      = NA_character_,
    deg.change       = NA_real_,
    peak.curr        = NA_character_,
    peak.new         = NA_character_,
    peak.new.rel     = NA_real_,
    row.names        = NULL,
    stringsAsFactors = FALSE )
}


#' Median off-target abundance within a fluorophore's own population, as a
#' fraction of its own signal. Needs no ground truth, so it is the acceptance
#' test that survives into production.
#' @noRd
.fix.leakage <- function( raw.data, design, target, panel, max.events = 20000L ) {

  if ( nrow( raw.data ) > max.events )
    raw.data <- raw.data[ sample.int( nrow( raw.data ), max.events ), ,
                          drop = FALSE ]

  unmixed <- unmix.ols.fast( raw.data, design )
  colnames( unmixed ) <- rownames( design )

  on <- stats::median( unmixed[ , target ] )
  if ( !is.finite( on ) || on <= 0 ) return( NA_real_ )

  off <- setdiff( panel, target )
  if ( length( off ) == 0 ) return( NA_real_ )

  sum( abs( apply( unmixed[ , off, drop = FALSE ], 2, stats::median ) ) ) / on
}


#' Per-event autofluorescence library assignment for `bg.mode = "per.cell"`.
#'
#' Each event is assigned the single `af.spectra` row whose fitted removal
#' leaves the smallest total apparent fluorophore signal (L1) — `AutoSpectral`'s
#' `assign.af.fluorophores()`, or its Rcpp companion `assign.af.fluor.fast()`
#' when available. The assigned row's abundance is the scalar least-squares
#' fit of the event onto that row's component orthogonal to the panel,
#' `k = (y . r_j) / (r_j . r_j)`, clamped at zero so a poor match is never
#' added back as negative background.
#'
#' @param y Numeric matrix (events x detectors), raw data to assign.
#' @param spectra.fluor Numeric matrix (fluorophores x detectors), the panel,
#'   excluding any `af.name` row.
#' @param af.spectra Numeric matrix (n.af x detectors), the candidate library.
#' @param r.library Numeric matrix (detectors x n.af), the component of each
#'   library row the panel cannot explain.
#' @param r.dots Numeric vector, length n.af, `colSums(r.library^2)`.
#' @param n.threads Integer, OpenMP threads for
#'   `AutoSpectralRcpp::assign.af.fluor.fast()`.
#'
#' @return A list with `index` (integer, length `nrow(y)`, 1-based row of
#'   `af.spectra`) and `abundance` (numeric, length `nrow(y)`, clamped
#'   `>= 0`).
#'
#' @noRd
.fix.assign.af <- function( y, spectra.fluor, af.spectra, r.library, r.dots,
                            n.threads = 1L ) {

  k.raw <- sweep( y %*% r.library, 2, r.dots, "/" )

  af.idx <- if (
    requireNamespace( "AutoSpectralRcpp", quietly = TRUE ) &&
    "assign.af.fluor.fast" %in% ls( getNamespace( "AutoSpectralRcpp" ) )
  ) {
    AutoSpectralRcpp::assign.af.fluor.fast(
      raw.data   = y,
      spectra    = spectra.fluor,
      af.spectra = af.spectra,
      threads    = n.threads )
  } else {
    assign.af.fluorophores( raw.data = y, spectra = spectra.fluor,
                            af.spectra = af.spectra )
  }

  k.chosen <- k.raw[ cbind( seq_len( nrow( y ) ), af.idx ) ]
  k.chosen[ !is.finite( k.chosen ) | k.chosen < 0 ] <- 0

  list( index = af.idx, abundance = k.chosen )
}


#' Stratified subsample of raw events by dominant fluorophore, for
#' downsampling before the correction loops.
#'
#' A uniform subsample keeps every population in the same proportion the raw
#' file already had it in. In a fully stained sample that means the biggest
#' population, usually events with most of the panel silent, sets how hard
#' everything else gets thinned; a dim or rare fluorophore's own positive
#' population is cut at the same rate as a background bulk that was already
#' more than large enough to define its threshold, and it is exactly that
#' population the pair estimator and the signature fit need the most events
#' from.
#'
#' Membership is by positivity, not by a single winner-take-all dominance
#' label: an event above threshold for two fluorophores at once belongs to
#' both strata, not whichever scores higher. That matters specifically for a
#' fully stained sample, where co-expression is real and is exactly the
#' population downstream identifiability checks rely on; a mutually-exclusive
#' assignment would credit a co-positive event to only one of the two markers
#' and let a large stratum's subsampling discard it from the other's
#' population by chance. On a concatenated single-stain pool the distinction
#' is close to moot, since the spread-scaled boundary already keeps one
#' tube's events out of another's positive population except right at the
#' edge, so this reduces to the same partition either way there.
#'
#' Every fluorophore's stratum is sampled against its own quota: a floor that
#' keeps a dim dye's whole positive population intact whenever it is smaller
#' than the floor, and a proportional share of what is left of the budget for
#' strata that can use more. The chosen events are the union across
#' fluorophores, so a co-positive event drawn by more than one stratum is
#' kept once and effectively has a higher retention chance, which is the
#' right direction. Events positive for nothing are capped separately as a
#' flat fraction of the total, since additional events there mostly refine a
#' threshold that is already well determined.
#'
#' @param abundance Numeric matrix (events x fluorophores), a first-pass
#'   unmixed abundance under the spectra downsampling will feed into.
#' @param thresholds Named numeric vector, flat per-fluorophore positivity
#'   thresholds in the same units as `abundance`.
#' @param threshold.matrix Numeric matrix, same shape as `abundance`, the
#'   spread-scaled positivity boundary from `get.spread.thresholds()`.
#' @param n.total Integer, the total event budget.
#' @param background.frac Numeric in (0, 1), the share of `n.total` reserved
#'   for events positive for nothing. Default `0.3`.
#' @param min.stratum Integer, the floor below which a fluorophore's positive
#'   population is kept whole rather than thinned. Default `2000`.
#'
#' @return Integer vector of row indices to keep.
#'
#' @noRd
.fix.stratified.sample <- function( abundance, thresholds, threshold.matrix,
                                    n.total, background.frac = 0.3,
                                    min.stratum = 2000L ) {

  n <- nrow( abundance )
  if ( n <= n.total ) return( seq_len( n ) )

  positive <- ( abundance - threshold.matrix ) > 0

  background.idx <- which( rowSums( positive ) == 0L )
  positive.idx    <- lapply( seq_len( ncol( abundance ) ),
                             function( j ) which( positive[ , j ] ) )

  n.background <- min( length( background.idx ),
                       round( background.frac * n.total ) )
  budget.positive <- n.total - n.background

  sizes   <- vapply( positive.idx, length, integer( 1 ) )
  floor.n <- pmin( sizes, as.integer( min.stratum ) )

  # A panel with many fluorophores can demand more floor than the budget
  # holds; when it does, every stratum's floor is scaled down by the same
  # factor rather than starving whichever strata come later in column order.
  quota <- if ( sum( floor.n ) > budget.positive ) {
    pmin( sizes, floor( floor.n * budget.positive / max( sum( floor.n ), 1 ) ) )
  } else {
    spare           <- pmax( sizes - floor.n, 0 )
    leftover.budget <- budget.positive - sum( floor.n )
    extra <- if ( sum( spare ) > 0 )
      leftover.budget * spare / sum( spare ) else rep( 0, length( sizes ) )
    pmin( sizes, floor.n + round( extra ) )
  }

  chosen.positive <- unique( unlist( Map( function( idx, k )
    if ( length( idx ) <= k ) idx else sample( idx, k ),
    positive.idx, quota ), use.names = FALSE ) )

  chosen.background <- if ( n.background >= length( background.idx ) )
    background.idx else sample( background.idx, n.background )

  c( chosen.positive, chosen.background )
}

#' Closed-form Huber M-estimator for a two-parameter (intercept + slope) fit,
#' in place of fit.robust.linear.model()'s call to MASS::rlm(). Same default
#' behaviour: weights from psi.huber with tuning constant k = 1.345, scale
#' re-estimated each iteration from the MAD of the current residuals,
#' iterating to a coefficient-change tolerance, initial weights of 1 so the
#' first fit is OLS. On non-convergence it returns a zero slope rather than
#' falling back to OLS, matching fit.robust.linear.model()'s fix.unmix = TRUE
#' behaviour, since a coefficient this function cannot pin down is exactly
#' the case the caller needs to see as untrustworthy rather than as some
#' other estimator's best guess.
#'
#' Every step is the same closed-form weighted least squares fit.trace() in
#' `.fix.envelope.slope()` already uses, so this makes no call to
#' stats::lm.wfit(), solve(), or any other compiled linear-algebra routine,
#' which is the second and larger BLAS dependency this pair loop had.
#'
#' @param x,y Numeric vectors, the predictor and response.
#' @param k Numeric, the Huber tuning constant. Default `1.345`.
#' @param max.iter Integer, maximum IRLS iterations. Default `100`.
#' @param tol Numeric, relative coefficient-change convergence tolerance.
#'   Default `1e-4`.
#' @param start Optional numeric `c(intercept, slope)`, a warm start. When
#'   supplied, the first iteration's weights come from its residuals instead
#'   of OLS weights, so a caller iterating the same pair across outer passes
#'   can resume near its last estimate. Default `NULL`.
#'
#' @return Numeric vector `c(intercept, slope)`.
#'
#' @noRd
.fix.huber.slope <- function( x, y, k = 1.345, max.iter = 100L, tol = 1e-4,
                              start = NULL ) {

  if (
    requireNamespace( "AutoSpectralRcpp", quietly = TRUE ) &&
    "fix_huber_slope_rcpp" %in% ls( getNamespace( "AutoSpectralRcpp" ) )
  ) {
    AutoSpectralRcpp::fix_huber_slope_rcpp(
      x, y, k = k, max_iter = as.integer( max.iter ),
      tol = tol, start = start
    )
  } else {
    stop( "Install AutoSpectralRcpp" )
  }

}
