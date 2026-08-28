# Fix My Unmix

Corrects the component of a reference spectra error that lives inside
the row space of the spectra themselves, using a fully stained sample
and the fact that a marker-negative population must read zero in every
channel it is negative for.

This is the complement of
[`correct.unmixing.signatures()`](https://drcytometer.github.io/AutoSpectral/reference/correct.unmixing.signatures.md),
which corrects the shape of individual rows from detector-space
residuals and is structurally blind to in-span error. Run this second,
on the spectra that function returns, because in unmixed space a
row-shape error and a genuine spillover error are indistinguishable.

The correction runs in two phases. The first estimates the residual
spillover matrix pair by pair, from a robust line through the events
that are negative for the target channel: their apparent abundance there
is spillover from the source, so its slope against source abundance is
the residual coefficient. Truncating the target from above shrinks that
slope toward zero when the two markers are co-expressed, which is the
safe direction. The second phase does not back-solve a spectrum out of
the matrix. It uses the matrix only to define clean populations, then
re-measures each fluorophore's signature directly from background
subtracted raw data with
[`extract.raw.signature()`](https://drcytometer.github.io/AutoSpectral/reference/extract.raw.signature.md),
so the result is a physical spectrum rather than whichever combination
of the existing rows happens to reproduce the required compensation.

Both phases have a bias that can be measured without any ground truth.
Run on the control the reference spectra came from, this function should
return exactly what it was given; whatever it returns instead is its own
artefact. Supplying that run as `null.fit` subtracts the artefact from
the coefficients and from the signatures, which recovers substantially
more than refusing the same quantities on the same information.

Autofluorescence is fitted jointly with the panel from a multi-component
basis (`bg.mode = "af.deconv"`), not represented by a single row. In a
real sample the background a single row leaves behind tracks cell type,
cell type determines marker expression, and the leftover is then
indistinguishable from spillover by any channel-against-channel fit.

Every row update must pass an acceptance stack; a fluorophore that fails
any gate keeps its starting spectrum.

## Usage

``` r
fix.my.unmix(
  spectra,
  unstained.sample,
  fully.stained.sample,
  flow.control,
  asp,
  variants = NULL,
  af.name = "AF",
  af.basis = NULL,
  af.n.pc = "auto",
  bg.mode = c("af.deconv", "af.row", "global.mean", "none"),
  large.gate = TRUE,
  max.iter = 20L,
  downsample = 20000,
  downsample.background.frac = 0.3,
  downsample.min.stratum = 2000L,
  unstained.threshold = 0.99,
  unstained.margin = 1.3,
  spread.kappa = 2,
  envelope.quantiles = c(0.05, 0.5),
  min.negative.frac = 0.1,
  max.disagreement = 0.5,
  leakage.prior = TRUE,
  span.fraction = 0.6,
  min.negative.events = 200L,
  min.bin.negative = 25L,
  min.span = 5,
  min.rise = 1,
  n.levels = 60L,
  min.bin.events = 50L,
  n.levels.pair = 10L,
  multivariate = TRUE,
  ridge = 1e-06,
  max.truncated.events = 20000L,
  max.mask.passes = 3L,
  estimator = c("truncated", "envelope"),
  source.dominant = TRUE,
  spread.addback = FALSE,
  anchor.weight = 1,
  max.coefficient = 0.5,
  convergence.threshold = 0.01,
  convergence.quantile = 0.95,
  update.spectra = TRUE,
  step = 1,
  null.fit = NULL,
  max.resid.ratio = 3,
  max.intercept.ratio = 3,
  intercept = TRUE,
  min.explained = 0.8,
  max.explained = 1.2,
  max.resid = 0.03,
  max.intercept = 0.03,
  min.bg.align = -0.9,
  gate.on.bias = FALSE,
  min.impact.ratio = 2,
  max.angle = 15,
  max.clamp.frac = 0.15,
  max.anchor = 0.1,
  max.vif = 500,
  max.condition.increase = 1.05,
  peak.shift.min.rel = 0.7,
  max.hotspot = 5,
  leakage.margin = 0.05,
  n.threads = 1L,
  figures = TRUE,
  save = TRUE,
  verbose = TRUE
)
```

## Arguments

- spectra:

  The spectral matrix, fluorophores x detectors, L-infinity normalised.
  Ideally the output of
  [`correct.unmixing.signatures()`](https://drcytometer.github.io/AutoSpectral/reference/correct.unmixing.signatures.md).

- unstained.sample:

  File path and name for a raw unstained sample, acquired the same day
  and matching the autofluorescence of the fully stained sample.

- fully.stained.sample:

  File path and name for a raw fully stained sample.

- flow.control:

  The flow.control list.

- asp:

  The AutoSpectral parameter list.

- variants:

  The variant list returned by
  [`get.spectral.variants()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectral.variants.md).
  Only `variants$spillover.spread` is used, to set abundance-dependent
  positivity boundaries and to correct the envelope for the widening of
  the negative population. If `NULL`, boundaries are flat and the
  correction is biased at the bright end.

- af.name:

  Character or `NULL`, the name of an autofluorescence row in `spectra`.
  It is never treated as a panel fluorophore and never corrected.
  Default `"AF"`.

- af.basis:

  Optional matrix (components x detectors) from
  [`get.af.basis()`](https://drcytometer.github.io/AutoSpectral/reference/get.af.basis.md).
  When `NULL` and `bg.mode = "af.deconv"`, it is built from the
  unstained sample. Default `NULL`.

- af.n.pc:

  Integer or `"auto"`, passed to
  [`get.af.basis()`](https://drcytometer.github.io/AutoSpectral/reference/get.af.basis.md).
  Default `"auto"`.

- bg.mode:

  Character. `"af.deconv"` (default) fits a multi-component
  autofluorescence basis jointly with the panel; `"af.row"` uses the
  single `af.name` row already in `spectra`; `"global.mean"` uses the
  mean unstained spectrum as a single background row; `"none"` fits no
  background.

- large.gate:

  Logical, whether to use a large scatter gate. Default `TRUE`.

- max.iter:

  Integer, maximum spillover-matrix iterations. Default `20`.

- downsample:

  Logical or numeric. `FALSE` disables downsampling; a numeric gives the
  number of events to use. Values above the event count are reduced to
  it by a stratified sample over each event's dominant fluorophore under
  the starting spectra, so a dim or rare dye's own positive population
  is not thinned at the same rate as the background bulk. Default
  `20000`.

- downsample.background.frac:

  Numeric in (0, 1), the share of `downsample` reserved for events
  dominant for nothing. Default `0.3`.

- downsample.min.stratum:

  Integer, the floor below which a fluorophore's own positive population
  is kept whole by the stratified downsample rather than thinned
  further. Default `2000`.

- unstained.threshold:

  Numeric in (0, 1), the percentile of the unstained control defining
  positivity. Default `0.99`.

- unstained.margin:

  Numeric, multiplier applied to that threshold. Default `1.3`.

- spread.kappa:

  Numeric, how many spillover-spread standard deviations above the flat
  threshold still count as negative. Default `2`.

- envelope.quantiles:

  Numeric pair. The first is the envelope the slope is measured on; the
  second is a comparison quantile whose disagreement with the first
  reports co-expression. Default `c( 0.05, 0.5 )`.

- min.negative.frac:

  Numeric, the smallest fraction of target-negative events the brightest
  source bins may contain before the pair is declared unidentifiable and
  left alone. Default `0.10`.

- max.disagreement:

  Numeric, envelope-versus-median slope disagreement above which the
  coefficient's trust weight is reduced. Default `0.5`.

- leakage.prior:

  Logical, whether to weight coefficients by the spillover a plausible
  spectral error could produce, from
  [`get.variant.leakage.prior()`](https://drcytometer.github.io/AutoSpectral/reference/get.variant.leakage.prior.md).
  When `FALSE`, or when no variants are supplied, the estimate supplies
  its own prior variance, which lets a large spurious coefficient
  justify itself. Default `TRUE`.

- span.fraction:

  Numeric in (0, 1\], passed to
  [`get.variant.leakage.prior()`](https://drcytometer.github.io/AutoSpectral/reference/get.variant.leakage.prior.md).
  The fraction of real spectral error the variant family is expected to
  span; lower values give corrections more freedom outside the observed
  variant directions. Default `0.6`.

- min.negative.events:

  Integer, minimum events for a pair fit. Default `200`.

- min.bin.negative:

  Integer, minimum target-negative events an abundance bin must contain
  before its envelope quantile is used. Default `25`.

- min.span:

  Numeric, minimum source abundance span in units of that fluorophore's
  own threshold. Default `5`.

- min.rise:

  Numeric, the fitted rise across the source abundance span, in standard
  deviations of the target's negative population, below which a
  coefficient is treated as unresolvable and set to zero. Default `1`.

- n.levels:

  Integer, the maximum number of abundance bins for the signature fit.
  Must exceed the panel size by at least three when `multivariate` is
  `TRUE`. Default `60`.

- min.bin.events:

  Integer, the fewest events an abundance bin may contain, passed to
  [`extract.raw.signature()`](https://drcytometer.github.io/AutoSpectral/reference/extract.raw.signature.md).
  Default `50`.

- n.levels.pair:

  Integer, abundance bins for the pair estimator. Default `10`.

- multivariate:

  Logical, passed to
  [`extract.raw.signature()`](https://drcytometer.github.io/AutoSpectral/reference/extract.raw.signature.md),
  whether every active abundance is fitted jointly rather than the
  nuisance rows being subtracted beforehand. Default `TRUE`.

- ridge:

  Numeric, ridge penalty for that joint fit. Default `1e-6`.

- max.truncated.events:

  Integer, cap on the events used for the robust pair fit. Everything
  above the source threshold is kept and the negative bulk below it is
  subsampled, since the bulk sits at the origin and carries no leverage.
  Default `20000`.

- max.mask.passes:

  Integer, how many times the target-negative selection is recomputed
  with the current spillover estimate removed. Deciding negativity on
  the raw target caps the fit's source range at the target threshold
  divided by the coefficient, so the estimator loses its leverage
  precisely as the coefficient it is measuring grows. Default `3`.

- estimator:

  Character, which pair estimator supplies the coefficients.
  `"truncated"` (default) is a robust line through the target-negative
  events; `"envelope"` is the binned lower envelope.

- source.dominant:

  Logical, whether source-positive events dominated by another
  fluorophore are dropped from the pair fit. Their apparent source
  abundance is that other fluorophore's brightness read through a
  spectral error, so binning on it returns that error as this pair's
  spillover. Default `TRUE`.

- spread.addback:

  Logical, whether the measured spillover spread is added back to the
  envelope before its slope is taken. Default `FALSE`.

- anchor.weight:

  Numeric, the weight ceiling for the zero-abundance anchor bin, in
  multiples of the best positive bin. Default `1`.

- max.coefficient:

  Numeric, the largest residual spillover coefficient accepted. Rows are
  L-infinity normalised, so a coefficient this large means one row is
  wrong by that fraction of another row's whole spectrum, which violates
  the small-error premise the linearisation rests on. Default `0.2`.

- convergence.threshold:

  Numeric, residual spillover coefficient at which iteration stops.
  Default `0.01`.

- convergence.quantile:

  Numeric, the quantile of the off-diagonal coefficients the convergence
  test uses, so that one pathological pair cannot define convergence.
  Default `0.95`.

- update.spectra:

  Logical, whether to run the raw-space signature phase. When `FALSE`,
  only the spillover and compensation matrices are returned. Default
  `TRUE`.

- step:

  Numeric, the fraction of each accepted signature change applied.
  Default `1`.

- null.fit:

  Optional, the result of running this function on the control the
  reference spectra were extracted from, where it should return exactly
  what it was given. What it returns instead is its own artefact,
  measured with no ground truth, and it is subtracted from both the
  spillover matrix and the signatures. Its `signature.log` also sets
  `max.resid` and `max.intercept` from the panel's own distribution.
  Default `NULL`.

- max.resid.ratio:

  Numeric, multiple of the null run's median `resid.rel` above which a
  fit is refused. Ignored without `null.fit`. Default `3`.

- max.intercept.ratio:

  Numeric, the same for `intercept.rel`. Default `3`.

- intercept:

  Logical, whether the signature fit carries an intercept, which absorbs
  any constant offset the removal of the other panel rows left behind.
  Default `TRUE`.

- min.explained:

  Numeric, minimum `explained.total`, the fraction of the brightest
  bin's signal the whole signature fit reproduces. Default `0.8`.

- max.explained:

  Numeric, the maximum of the same quantity. A value above one means the
  nuisance removal over-subtracted. Default `1.2`.

- max.resid:

  Numeric, maximum relative fit residual. Overridden by `null.fit` when
  supplied. Default `0.03`.

- max.intercept:

  Numeric, maximum relative intercept. Overridden by `null.fit` when
  supplied. Default `0.03`.

- min.bg.align:

  Numeric, minimum cosine between the fitted intercept and the candidate
  signature. A value near minus one means the fit is trading a
  background floor against the slope, the event-level signature of an
  unmodelled background confound. Default `-0.9`.

- gate.on.bias:

  Logical, whether `min.impact.ratio` actually blocks a candidate,
  versus being computed and logged only. Default `FALSE`. Evidence from
  a counterfactual audit (score whether accepting each rejected
  candidate would have moved the row toward or away from ground truth):
  of every candidate this gate rejected, 79.9% would have helped if
  accepted (n=144, two substrates) - 90.7% on one substrate, a coin-flip
  47.2% on the other. The gate is net-harmful on the substrate where it
  discriminates worst and only weakly useful on the other, so it
  defaults off rather than removed - `bias.impact` and the ratio are
  still computed and logged either way, so the evidence for re-enabling
  it on a new substrate remains visible without code changes.

- min.impact.ratio:

  Numeric, how many times a proposed step's abundance effect must exceed
  the effect of the same phase's bias on its own control, measured
  through the current unmixing operator. Degrees weight every detector
  alike and every dye alike; this weights each row by the abundances it
  actually produces, so a large rotation of a dim row is cheap and a
  small rotation of a bright collinear row is not. Requires
  `null.spectra` and `gate.on.bias = TRUE`. Default `2`.

- max.angle:

  Numeric, maximum angular change of a row in degrees.

- max.clamp.frac:

  Numeric, maximum fraction of a candidate row's absolute mass that
  non-negativity clamping may remove. Default `0.15`.

- max.anchor:

  Numeric, maximum unremoved background relative to the brightest fitted
  signal. Default `0.10`.

- max.vif:

  Numeric, maximum variance inflation factor for the fluorophore's own
  abundance within its population. Default `500`.

- max.condition.increase:

  Numeric, the factor by which a single accepted row may increase the
  condition number of the unmixing design. Default `1.05`.

- peak.shift.min.rel:

  Numeric, a candidate row whose peak detector differs from the current
  row is only accepted if the current signature's own value at that new
  peak channel, relative to the current signature's own peak, is at
  least this large - a close secondary peak taking over is plausible, a
  shift to a channel the current signature barely uses is not. Use `Inf`
  to never allow a peak shift and `0` to always allow one. Default
  `0.7`.

- max.hotspot:

  Numeric, hotspot scale above which a fluorophore is considered
  inseparable from the autofluorescence basis and is frozen. Default
  `5`.

- leakage.margin:

  Numeric, the fraction by which held-out leakage must increase, not
  merely be numerically larger, before a candidate is refused for it.
  `.fix.leakage()` is a held-out statistic on a dominance population
  that can be small, and "any increase blocks" was itself producing
  false positives - candidates blocked on multiple consecutive fits that
  would demonstrably have helped if accepted. Default `0.05`.

- n.threads:

  Integer, OpenMP threads for the batched pair estimator
  (`fix_envelope_truncated_batch_rcpp()`). Keep at the default unless
  this call is not itself running inside another parallel context (e.g.
  one sample of several under `mclapply()`), since the two layers of
  parallelism would otherwise compete for the same cores. Default `1L`.

- figures:

  Logical, whether to write the spillover heatmap. Default `TRUE`.

- save:

  Logical, whether to write the csv outputs. Default `TRUE`.

- verbose:

  Logical, controls messaging. Default `TRUE`.

## Value

A named list:

- `spectra`:

  The corrected spectra, fluorophores x detectors. Rows that failed a
  gate are unchanged.

- `spectra.backsolved`:

  `spillover %*% spectra`, the algebraic back-solve. Diagnostic only: it
  reproduces the compensation exactly but is not constrained to be a
  physical spectrum.

- `spillover`:

  The estimated residual spillover matrix, fluorophores x fluorophores.

- `compensation`:

  Its inverse.

- `trust`:

  Per-coefficient trust weights used to damp the update.

- `coefficient.log`:

  Per-pair fit diagnostics from the final iteration.

- `signature.log`:

  Per-fluorophore signature statistics and gate outcomes.

- `convergence.log`:

  Per-iteration delta history.

- `af.basis`, `af.hotspot`, `af.frozen`:

  The autofluorescence basis, its coupling to the panel, and the
  fluorophores frozen because of it.
