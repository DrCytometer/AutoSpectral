# Correct Spectra by Graphical Lasso

Estimates the residual spillover matrix and re-measures fluorophore
signatures from a fully stained sample, the same task
[`fix.my.unmix()`](https://drcytometer.github.io/AutoSpectral/reference/fix.my.unmix.md)
performs, but replaces its pair-at-a-time negative-envelope estimator
with a joint, L1-penalised regression: for every target fluorophore, its
apparent (compensated) abundance among events negative for it is
regressed on every other fluorophore's abundance at once, with an L1
penalty that drives the coefficient of an uninvolved fluorophore to
exactly zero. This is neighbourhood selection (Meinshausen & Buhlmann
2006) applied one row at a time to the residual spillover matrix rather
than to a covariance matrix, so the returned coefficients are directly
the row of `spillover` that
[`fix.my.unmix()`](https://drcytometer.github.io/AutoSpectral/reference/fix.my.unmix.md)
estimates pair by pair.

The two functions share an identification strategy and diverge in how
they handle co-expression. Both restrict the fit to events negative for
the target, since a marker-negative population must read zero in every
channel it is negative for - whatever a negative population's apparent
abundance correlates with there is spillover, not real expression,
because the cells contributing it do not express the target at all.
[`fix.my.unmix()`](https://drcytometer.github.io/AutoSpectral/reference/fix.my.unmix.md)
then estimates one source at a time and relies on `source.dominant`
masking to stop a second bright fluorophore's own spillover from being
booked against the source under test. This function instead regresses
the target jointly on every other fluorophore's abundance in the same
fit, so a source that only looks correlated with the target because it
is itself correlated with the true culprit has that correlation absorbed
by the true culprit's own coefficient, and the L1 penalty turns the "is
this source really contributing" question into automatic variable
selection instead of a hand-tuned dominance heuristic. The trade-off is
that pairs which are both spectrally similar and genuinely co-expressed
remain as hard for this estimator to separate as for the pairwise one;
nothing about fitting them jointly changes what two collinear predictors
can and cannot identify.

As in
[`fix.my.unmix()`](https://drcytometer.github.io/AutoSpectral/reference/fix.my.unmix.md),
nothing is back-solved out of the spillover matrix. The matrix is used
only to refine which events count as negative for each target; the
spectra themselves are re-measured directly from background subtracted
raw data by
[`extract.raw.signature()`](https://drcytometer.github.io/AutoSpectral/reference/extract.raw.signature.md),
exactly as
[`fix.my.unmix()`](https://drcytometer.github.io/AutoSpectral/reference/fix.my.unmix.md)'s
second phase does, with one difference: the co-fluorophores carried into
that fit (`active`) are the ones this function's lasso selected for that
target, not the whole panel run through a ridge penalty. A pair the
lasso found no evidence of coupling for is left out of the signature fit
entirely, rather than being included and shrunk.

A coefficient's sign carries the identification, and it is not symmetric
in target and source. A marker-negative population cannot read below
zero from real emission - nothing biological subtracts signal - so where
a positive coefficient is confounded with genuine co-expression (see
`CONTEXT_information_sources.md`'s hypernegativity argument), a negative
one is unambiguous, and it identifies a different row than a positive
one would. Corrupting fluorophore j's reference row with eps of
fluorophore k's shape
(`S.wrong[j,] <- normalize(S.true[j,] + eps * S.true[k,])`) leaves the
pair fit the way one might expect - k's population regressed on j -
reading essentially zero, because unmixing is a change of basis and the
two rows actually coupled by that change are the other way round: j's
own population, regressed on k, reads close to `-eps`. The error is in
row j, but its fingerprint is a negative coefficient found while fitting
row k, not a positive one found while fitting row j's own. Every
target's per-row lasso fit already estimates every cell of this matrix
correctly, once every fluorophore has been run as a target in turn - the
fitted matrix itself is not wrong - but a target's active set, if it
were read only from its own row, would never see a negative coefficient
discovered while fitting someone else's row. `active.set` folds this in
explicitly: for every negative coefficient found anywhere in the fitted
matrix, it is the fluorophore whose row produced the fit (the source,
not the target of that particular fit) whose reference shape needs the
fitted target added to its own signature re-measurement.

## Usage

``` r
correct.spectra.glasso(
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
  downsample = 20000,
  downsample.background.frac = 0.3,
  downsample.min.stratum = 2000L,
  unstained.threshold = 0.99,
  unstained.margin = 1.3,
  spread.kappa = 2,
  min.negative.events = 200L,
  max.truncated.events = 20000L,
  max.mask.passes = 3L,
  mask.tolerance = 0.05,
  min.span = 5,
  min.rise = 1,
  max.coefficient = 0.5,
  coefficient.tolerance = 1e-04,
  n.lambda = 40L,
  lambda.min.ratio = 0.001,
  margin.frac = 0.05,
  max.iter = 5L,
  convergence.threshold = 0.01,
  convergence.quantile = 0.95,
  step.spillover = 1,
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
  max.clamp.frac = 0.15,
  max.anchor = 0.1,
  max.vif = 500,
  max.condition.increase = 1.05,
  peak.shift.min.rel = 0.7,
  max.angle = 15,
  max.hotspot = 5,
  leakage.margin = 0.05,
  n.levels = 60L,
  min.bin.events = 50L,
  multivariate = TRUE,
  ridge = 1e-06,
  output.suffix = "_glasso",
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
  positivity boundaries. If `NULL`, boundaries are flat. Default `NULL`.

- af.name:

  Character or `NULL`, the name of an autofluorescence row in `spectra`.
  Never treated as a panel fluorophore and never corrected. Default
  `"AF"`.

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

- downsample:

  Logical or numeric. `FALSE` disables downsampling; a numeric gives the
  number of events to use, stratified by dominant fluorophore under the
  starting spectra. Default `20000`.

- downsample.background.frac:

  Numeric in (0, 1), the share of `downsample` reserved for events
  dominant for nothing. Default `0.3`.

- downsample.min.stratum:

  Integer, the floor below which a fluorophore's own positive population
  is kept whole. Default `2000`.

- unstained.threshold:

  Numeric in (0, 1), the percentile of the unstained control defining
  positivity. Default `0.99`.

- unstained.margin:

  Numeric, multiplier applied to that threshold. Default `1.3`.

- spread.kappa:

  Numeric, how many spillover-spread standard deviations above the flat
  threshold still count as negative. Default `2`.

- min.negative.events:

  Integer, the fewest target-negative events a fluorophore must have
  before its row is fit at all; below this the fluorophore is left
  uncorrected for this iteration. Default `200`.

- max.truncated.events:

  Integer, cap on the target-negative events used per lasso fit; above
  the cap the population is subsampled, since the fit's cost scales with
  it directly. Default `20000`.

- max.mask.passes:

  Integer, how many times the target-negative selection is recomputed
  with the current fitted row's joint spillover contribution removed.
  Default `3`.

- mask.tolerance:

  Numeric, the fraction of the target-negative set that may change
  between mask passes before the mask is treated as settled. Default
  `0.05`.

- min.span:

  Numeric, minimum abundance span a selected source must have across the
  target-negative population, in units of the source's own flat
  threshold, or its coefficient is discarded. Default `5`.

- min.rise:

  Numeric, the fitted rise a selected source's coefficient must produce
  across its own span, in standard deviations of the target's negative
  population, or the coefficient is discarded. Default `1`.

- max.coefficient:

  Numeric, the largest residual spillover coefficient accepted; a larger
  fitted value is discarded rather than kept, on the same small-error
  premise
  [`fix.my.unmix()`](https://drcytometer.github.io/AutoSpectral/reference/fix.my.unmix.md)
  uses. Default `0.5`.

- coefficient.tolerance:

  Numeric, fitted coefficients smaller than this in absolute value are
  treated as exactly zero. Default `1e-4`.

- n.lambda:

  Integer, points on the lasso regularisation path. Default `40`.

- lambda.min.ratio:

  Numeric in (0, 1), the smallest path value as a fraction of
  `lambda.max`, the value at which every coefficient is already zero.
  Default `1e-3`.

- margin.frac:

  Numeric, the share of the two-way holdout error reduction a larger,
  sparser `lambda` may give up and still be preferred, measured from the
  intercept-only fit at `lambda.max` down to the best point on the grid.
  Scale free, so it keeps its meaning whatever the target-negative
  population's own variance. Default `0.05`.

- max.iter:

  Integer, maximum outer spillover-matrix iterations. Default `5`.

- convergence.threshold:

  Numeric, residual spillover coefficient at which iteration stops.
  Default `0.01`.

- convergence.quantile:

  Numeric, the quantile of the off-diagonal coefficients the convergence
  test uses. Default `0.95`.

- step.spillover:

  Numeric in (0, 1\], the fraction of each outer iteration's fitted
  spillover update applied. The lasso already shrinks and selects each
  row, so unlike
  [`fix.my.unmix()`](https://drcytometer.github.io/AutoSpectral/reference/fix.my.unmix.md)'s
  per-coefficient `trust` weighting this is a single damping factor for
  the whole update, only useful if the outer iteration oscillates.
  Default `1`.

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
  what it was given. Its bias is subtracted from both the spillover
  matrix and the accepted signatures. Default `NULL`.

- max.resid.ratio:

  Numeric, multiple of the null run's median `resid.rel` above which a
  fit is refused. Ignored without `null.fit`. Default `3`.

- max.intercept.ratio:

  Numeric, the same for `intercept.rel`. Default `3`.

- intercept:

  Logical, whether the signature fit carries an intercept. Default
  `TRUE`.

- min.explained:

  Numeric, minimum `explained.total`. Default `0.8`.

- max.explained:

  Numeric, maximum `explained.total`. Default `1.2`.

- max.resid:

  Numeric, maximum relative fit residual. Overridden by `null.fit`.
  Default `0.03`.

- max.intercept:

  Numeric, maximum relative intercept. Overridden by `null.fit`. Default
  `0.03`.

- min.bg.align:

  Numeric, minimum cosine between the fitted intercept and the candidate
  signature. Default `-0.9`.

- max.clamp.frac:

  Numeric, maximum fraction of a candidate row's absolute mass
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

  Numeric, see
  [`extract.raw.signature()`](https://drcytometer.github.io/AutoSpectral/reference/extract.raw.signature.md).
  Default `0.7`.

- max.angle:

  Numeric, maximum angular change of a row in degrees. Default `15`.

- max.hotspot:

  Numeric, hotspot scale above which a fluorophore is frozen as
  inseparable from the autofluorescence basis. Default `5`.

- leakage.margin:

  Numeric, the fraction by which held-out leakage must increase before a
  candidate is refused for it. Default `0.05`.

- n.levels:

  Integer, maximum abundance bins for the signature fit, passed to
  [`extract.raw.signature()`](https://drcytometer.github.io/AutoSpectral/reference/extract.raw.signature.md).
  Default `60`.

- min.bin.events:

  Integer, fewest events per abundance bin, passed to
  [`extract.raw.signature()`](https://drcytometer.github.io/AutoSpectral/reference/extract.raw.signature.md).
  Default `50`.

- multivariate:

  Logical, passed to
  [`extract.raw.signature()`](https://drcytometer.github.io/AutoSpectral/reference/extract.raw.signature.md).
  Default `TRUE`.

- ridge:

  Numeric, ridge penalty for that joint fit. Default `1e-6`.

- output.suffix:

  Character, appended to the csv and figure filenames so a run of this
  function does not overwrite
  [`fix.my.unmix()`](https://drcytometer.github.io/AutoSpectral/reference/fix.my.unmix.md)'s
  output in the same directory. Default `"_glasso"`.

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

  `spillover %*% spectra`, the algebraic back-solve. Diagnostic only.

- `spillover`:

  The estimated residual spillover matrix, fluorophores x fluorophores.
  Read `spillover[a, b]` as the estimated contribution of fluorophore
  a's true abundance to fluorophore b's apparent abundance under the
  starting spectra; it is filled cell by cell while each fluorophore in
  turn is the fit's target, and it is not symmetric.

- `compensation`:

  Its inverse.

- `coefficient.log`:

  Per-(source, target) pair, one row per pair, from the final iteration:
  the fitted coefficient (`0` if the lasso dropped it or a gate
  discarded it) and whether it was selected. Raw fit output; a negative
  `beta` here means the *source* row is the one implicated, not the
  target - see `active.set`.

- `lambda.log`:

  Per-target diagnostics from the final iteration: events used,
  candidate sources, sources selected, and the chosen `lambda`.

- `signature.log`:

  Per-fluorophore signature statistics and gate outcomes, from
  [`extract.raw.signature()`](https://drcytometer.github.io/AutoSpectral/reference/extract.raw.signature.md).

- `convergence.log`:

  Per-iteration delta history.

- `active.set`:

  Named list, the fluorophores carried into each target's signature fit:
  itself, whichever sources its own row selected with a positive
  coefficient (genuine spillover into it, which must be modelled to
  isolate its true signal), and whichever other fluorophores' own fits
  named it as a negative-coefficient source (evidence its own reference
  shape needs that fluorophore modelled out too). See `reversed.donors`
  for the second route on its own.

- `reversed.donors`:

  Named list, per fluorophore, the donors folded into its `active.set`
  by the negative-coefficient route: every other fluorophore whose own
  row-fit found this one as a negative- coefficient candidate source.

- `af.basis`, `af.hotspot`, `af.frozen`:

  The autofluorescence basis, its coupling to the panel, and the
  fluorophores frozen because of it.

## How the lasso row is estimated

For target \\j\\ and its target-negative event set \\N_j\\, let \\y\\ be
\\j\\'s compensated abundance over \\N_j\\ and \\X\\ be every other
fluorophore's abundance over the same events. The fit solves \$\$
\hat\beta_j = \arg\min\_\beta \tfrac{1}{2n}\lVert y - X\beta - \alpha
\rVert_2^2 + \lambda \sum\_{k \ne j} \|\beta_k\| \$\$ by cyclic
coordinate descent along a path of decreasing \\\lambda\\, with the
intercept \\\alpha\\ unpenalised so it can absorb a constant background
offset the way the pairwise estimator's anchor bin does. \\\lambda\\ is
chosen per target by a two-way split-half holdout: fit on one half,
score squared error on the other, and back; the two error curves are
summed over a shared \\\lambda\\ grid, and the largest \\\lambda\\
(sparsest fit) costing no more than `margin.frac` of the curve's total
error reduction is kept, then refit on the full target-negative set at
that \\\lambda\\. Preferring the sparser of two disagreeing halves is
the same safe-direction bias
[`fix.my.unmix()`](https://drcytometer.github.io/AutoSpectral/reference/fix.my.unmix.md)'s
target-negative truncation uses when co-expression would otherwise
shrink a slope toward zero.

In lay terms: instead of asking "does fluorophore A leak into
fluorophore B's channel, checked one pair at a time," this asks
fluorophore B's whole negative population "which of everyone else in the
panel actually explains the small amount of signal you have here," with
a rule that only keeps an answer if it clearly earns its place - a
fluorophore whose apparent contribution could just as well be noise, or
could just as well be somebody else's spillover showing up correlated
with it, gets a flat zero rather than a small uncertain number. And once
an answer does earn its place, which row it corrects depends on its
sign: a fluorophore that looks brighter than it should, in step with
something else in the panel, has real spillover to remove from its own
row; a fluorophore that looks dimmer than it should, in step with
something else, cannot really be losing signal to biology, so what is
wrong is not its own row but the other fluorophore's reference shape,
which has borrowed a slice of this one's.
