# Correct Unmixing Signatures

Detects and corrects shape errors in a reference `spectra` matrix
directly from stained sample data, for the situation where the available
reference spectra are subtly wrong for the sample being analysed -
typically because bead-derived spectra are applied to cells, or spectra
from another day, lot, or instrument state are in use.

Takes either `unstained.sample` and `fully.stained.sample` file paths,
reading and deriving `unmixed.thresholds`/`spillover.spread` itself (the
direct-use convention shared with
[`fix.my.unmix()`](https://drcytometer.github.io/AutoSpectral/reference/fix.my.unmix.md)
and
[`correct.spectra.glasso()`](https://drcytometer.github.io/AutoSpectral/reference/correct.spectra.glasso.md)),
or `raw.data` supplied directly as an already-read matrix, the
convention for a caller that already has it in memory. Supply one or the
other, not both.

Each event is assigned to the fluorophore it is most strongly positive
for, as a fraction of that fluorophore's own dynamic range above its
positivity threshold. Within each such dominance population, the data
are re-unmixed against a restricted design (the dominant fluorophore
plus any genuinely co-active ones), binned by abundance, anchored with
the background population, and the detector-space residual is regressed
on the restricted abundances. The dominant fluorophore's regression
slope estimates the error in its spectrum, and the row is updated by a
step chosen through a held-out search.

A correction is only applied when it passes every acceptance gate:

- the population spans enough abundance to identify a slope;

- the fluorophore's own term explains most of the background-subtracted
  signal at its brightest (`min.explained`), so a control too dim to
  characterise is left untouched;

- the held-out step search finds a step that reduces the restricted
  residual (`min.gain`);

- the correction is small relative to the row it corrects (`max.step`)
  and does not inflate the fluorophore's own apparent abundance
  (`max.span.drift`);

- the correction is not a background confound: when the apparent
  abundance trend is a common-mode background residual rather than a
  spectral error, the event-level regression splits that one physical
  direction into an intercept and a slope that are anti-collinear, and
  the correction is rejected (`max.bg.alignment`).

Rejected fluorophores keep their starting spectra; the method is
designed to leave a row untouched rather than risk making it worse.

Autofluorescence is handled by background removal before fitting rather
than as a design row, so a fluorophore's abundance and the background
amount are not estimated from the same spectral vector. For cells,
scatter-matched per-event subtraction (`bg.mode = "scatter.knn"`) uses
an unstained control matched on scatter, an independent measurement
channel. For beads and other uniform particles, whose scatter is
uninformative, the mean of the sample's own background events is
subtracted (`bg.mode = "global.mean"`).

When `scatter` is supplied and `gate.main` is `TRUE`, events are first
gated to the main scatter population by density, removing noise, debris
and most aggregates, whose distinct autofluorescence otherwise
contaminates the dominance populations and the background estimate.

## Usage

``` r
correct.unmixing.signatures(
  spectra,
  unstained.sample = NULL,
  fully.stained.sample = NULL,
  flow.control = NULL,
  asp = NULL,
  variants = NULL,
  af.name = "AF",
  raw.data = NULL,
  unmixed.thresholds = NULL,
  scatter = NULL,
  gate.main = FALSE,
  gate.level = 0.1,
  spillover.spread = NULL,
  spread.kappa = 2,
  bg.mode = c("global.mean", "scatter.knn", "none"),
  unstained = NULL,
  unstained.scatter = NULL,
  k.neighbors = 3L,
  unstained.threshold = 0.99,
  unstained.margin = 1.3,
  threshold.prior.weight = 0.3,
  threshold.max.ratio = 3,
  n.levels = 10L,
  n.iter = 6L,
  min.events = 200L,
  min.span = 5,
  min.explained = 0.5,
  min.gain = 0.002,
  step.grid = c(0, 0.03125, 0.0625, 0.125, 0.25, 0.5, 1),
  n.split.trials = 1L,
  min.split.frac = 0.6,
  max.step = 0.15,
  max.span.drift = 1.1,
  max.bg.alignment = -0.9,
  nuisance.frac = 0.5,
  footprint.frac = 0.02,
  footprint.min.channels = 3L,
  background.n = 5000L,
  true.spectra = NULL,
  min.deg.start = 0.1,
  verbose = TRUE
)
```

## Arguments

- spectra:

  Numeric matrix (fluorophores x detectors), the starting reference
  spectra to be corrected, L-infinity normalised.

- unstained.sample:

  File path and name for a raw unstained sample, used to derive
  `unmixed.thresholds` (any fluorophore not already covered by
  `variants$thresholds`) and, under `bg.mode = "scatter.knn"`, as the
  scatter-matched background reference. Ignored when `raw.data` is
  supplied directly. Required, together with `fully.stained.sample` and
  `flow.control`, whenever `raw.data` is not.

- fully.stained.sample:

  File path and name for a raw fully stained sample. Ignored when
  `raw.data` is supplied directly.

- flow.control:

  The flow.control list, used to select the scatter and spectral channel
  columns when reading `unstained.sample` and `fully.stained.sample`.
  Ignored when `raw.data` is supplied directly.

- asp:

  Optional AutoSpectral parameter list from
  [`get.autospectral.param()`](https://drcytometer.github.io/AutoSpectral/reference/get.autospectral.param.md).
  Used only to seed the random number generator (`asp$bird.seed`) for
  reproducible subsampling. Default `NULL`.

- variants:

  Optional variant list returned by
  [`get.spectral.variants()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectral.variants.md).
  Used as a shrinkage prior on top of the threshold this function
  derives from its own unmix (see `unmixed.thresholds`), not as the
  threshold itself:
  [`get.spectral.variants()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectral.variants.md)
  measures `variants$thresholds` from a per-cell-optimised unmix at the
  99.5th percentile, both of which make it systematically tighter than
  the bare-OLS, 99th-percentile estimate this function's own dominance
  assignment is judged against, especially for AF-collinear
  fluorophores. `variants$spillover.spread` supplies `spillover.spread`
  when that is not itself supplied. Default `NULL`.

- af.name:

  Character, the name of the autofluorescence row in `spectra`, or
  `NULL` if there is none. The AF row is never treated as a panel
  fluorophore and is never corrected. Default `"AF"`.

- raw.data:

  Optional numeric matrix (events x detectors), raw detector-space data.
  Pooled or concatenated single-stained controls, or a fully stained
  sample with well-separated populations. Columns must match the columns
  of `spectra`. Supply this directly to skip reading
  `fully.stained.sample` from disk – the internal calling convention,
  used when the raw matrix is already in memory. Default `NULL`, which
  requires `unstained.sample`, `fully.stained.sample` and `flow.control`
  instead.

- unmixed.thresholds:

  Optional named numeric vector covering every fluorophore in `spectra`
  (the autofluorescence row may be omitted), giving the positivity
  threshold in unmixed space. Default `NULL`, which computes the
  `unstained.margin`-scaled `unstained.threshold` percentile of
  `unstained`'s own bare unmix against `spectra` – the same
  pre-background-subtraction convention `dominant` below is judged
  against – then shrinks it toward `variants$thresholds` (see
  `threshold.prior.weight`, `threshold.max.ratio`).

- scatter:

  Optional numeric matrix (events x scatter parameters), row-matched to
  `raw.data`. Required for `gate.main` and for
  `bg.mode = "scatter.knn"`. Default `NULL`.

- gate.main:

  Logical, whether to gate events to the main scatter population before
  fitting. Requires `scatter`. Default `FALSE`.

- gate.level:

  Numeric in (0, 1). Events are kept when their 2D scatter density
  exceeds this fraction of the modal density. Default `0.1`.

- spillover.spread:

  Optional matrix from `get.spectral.variants()$spillover.spread`. When
  supplied, co-activity is judged against per-event thresholds widened
  by the spillover spread each bright fluorophore contributes, so
  spillover from the dominant dye is not mistaken for a co-active
  fluorophore. Default `NULL`, which falls back to
  `variants$spillover.spread` (flat thresholds if that is also
  unavailable).

- spread.kappa:

  Numeric, how many spillover-spread standard deviations above the flat
  threshold still count as negative. Default `2`.

- bg.mode:

  Character, background removal mode: `"global.mean"` (default),
  `"scatter.knn"`, or `"none"`. See Description.

- unstained:

  Optional numeric matrix (events x detectors), raw unstained control,
  required for `bg.mode = "scatter.knn"`.

- unstained.scatter:

  Optional numeric matrix (events x scatter parameters), row-matched to
  `unstained`, required for `bg.mode = "scatter.knn"`. When `gate.main`
  is `TRUE` the unstained control is gated with the same density rule.

- k.neighbors:

  Integer, neighbours for scatter-matched background subtraction.
  Default `3`.

- unstained.threshold:

  Numeric in (0, 1), the percentile of `unstained`'s own unmix defining
  positivity, used only when `unmixed.thresholds` is not supplied
  directly. Default `0.99` – less extreme than
  [`get.spectral.variants()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectral.variants.md)'s
  `0.995`, and therefore less sensitive to the handful of poorly
  AF-corrected events that dominate the very top of the tail for a
  collinear fluorophore.

- unstained.margin:

  Numeric, multiplier applied to that threshold. Default `1.3`.

- threshold.prior.weight:

  Numeric in `[0, 1]`, the log-space weight given to
  `variants$thresholds` when shrinking the internally-computed
  `unmixed.thresholds` toward it. `0` uses the internal estimate
  unchanged; `1` uses `variants$thresholds` outright. Ignored for any
  fluorophore `variants$thresholds` does not cover, or when `variants`
  is `NULL`. Default `0.3`.

- threshold.max.ratio:

  Numeric `> 1`, the magnitude ratio (either direction) between the
  internal estimate and `variants$thresholds` beyond which the internal
  side is treated as unreliable – most likely under-corrected
  autofluorescence – rather than partially trusted, and the threshold
  reverts to `variants$thresholds` outright instead of blending. Default
  `3`.

- n.levels:

  Integer, abundance bins per dominance population. Default `10`.

- n.iter:

  Integer, maximum correction iterations. Iteration stops early when no
  fluorophore accepts a step. Default `6`.

- min.events:

  Integer, minimum events in a dominance population. Default `200`.

- min.span:

  Numeric, minimum abundance span in units of the fluorophore's own
  threshold magnitude. Default `5`.

- min.explained:

  Numeric in (0, 1), minimum fraction of the background-subtracted
  signal the fluorophore's own term must account for at its brightest
  bin. Default `0.5`.

- min.gain:

  Numeric, minimum held-out relative residual reduction. Default
  `0.002`.

- step.grid:

  Numeric vector of candidate step sizes for the held-out search.
  Default `c( 0, 0.03125, 0.0625, 0.125, 0.25, 0.5, 1 )`; the small
  steps matter for tandem dyes, whose abundance-dependent variant
  mixture makes the residual objective a narrow valley.

- n.split.trials:

  Integer, number of independent random 50/50 splits used by the
  held-out step search. Default `1`, which reproduces the original
  fixed, un-reseeded split exactly. Raising this trades a single split's
  noise for a vote across `n.split.trials` splits (see
  `min.split.frac`); useful for fluorophores whose true correction is
  real but small relative to per-event noise, where a lone 50/50 split
  can land on the unlucky side.

- min.split.frac:

  Numeric in (0, 1\], the minimum fraction of `n.split.trials` splits
  that must independently find a beneficial step before one is accepted.
  Ignored when `n.split.trials = 1`. Default `0.6`.

- max.step:

  Numeric, maximum norm of the correction relative to the norm of the
  row it corrects. Larger proposed corrections are rejected outright
  rather than scaled down. Default `0.15`.

- max.span.drift:

  Numeric, maximum allowed growth of the fluorophore's apparent
  abundance span across iterations. Default `1.10`.

- max.bg.alignment:

  Numeric in `[-1, 0]`. The correction is rejected when the cosine
  between the event-level regression intercept and slope falls at or
  below this value. Default `-0.9`.

- nuisance.frac:

  Numeric in (0, 1), fraction of a dominance population that must be
  co-active for another fluorophore before it is carried as a nuisance
  term in the restricted design. Default `0.5`.

- footprint.frac:

  Numeric in `[0, 1)`. Restricts the slope fit and the held-out step
  search to detectors where the dominant dye's own current spectrum
  exceeds `footprint.frac` of its own peak. A dye cannot carry real
  shape-error signal in a detector it does not meaningfully emit into;
  for a narrow-emission dye read out on a wide detector array, those
  channels only dilute the held-out residual objective with noise from
  channels that are pure background for that dye. Abundance estimation
  and the `explained`/`bg.align` gates are unaffected; only the slope
  fit's response and the held-out objective's norm are restricted. `0`
  reproduces the previous, unrestricted behaviour exactly. Default
  `0.02`.

- footprint.min.channels:

  Integer, minimum detectors the `footprint.frac` mask must keep; below
  this the mask is dropped and every detector is used, so a
  pathologically narrow spectrum cannot leave too few channels to fit.
  Default `3L`.

- background.n:

  Integer, maximum background events used for the zero-abundance anchor.
  Default `5000`.

- true.spectra:

  Optional numeric matrix (fluorophores x detectors),
  independently-known ground truth with row names matching `spectra`.
  Purely diagnostic: when supplied, the returned `recovery` table
  reports the angular error before and after correction per fluorophore.

- min.deg.start:

  Numeric, degrees. Below this starting angular error, `recovered` is
  reported as `0` instead of `(deg.start - deg.after) / deg.start`,
  since a fluorophore that started (near) exactly correct makes that
  ratio blow up or divide by zero for a change of a fraction of a
  degree. Default `0.1`.

- verbose:

  Logical, controls messaging. Default `TRUE`.

## Value

A named list:

- `spectra`:

  The corrected reference spectra matrix. Rows that failed any
  acceptance gate are unchanged.

- `fit.log`:

  Data frame, one row per fluorophore per iteration, with the fit
  statistics and every gate quantity, including `bg.align`.

- `proposed.spectra`:

  Numeric matrix, one row per fluorophore per iteration
  (`"<fluorophore>.iter<n>"`), the L-infinity normalised candidate row
  that iteration would have produced, whether or not it was accepted.
  Row-matched to `fit.log` in the same order.

- `accepted`:

  Named logical vector, whether each panel fluorophore accepted at least
  one correction step.

- `dominant`:

  Integer vector over the (gated) events: the index into `panel` of each
  event's dominant fluorophore, `0` for background.

- `panel`:

  Character vector, the fluorophores eligible for correction (`spectra`
  rows minus `af.name`).

- `gate.keep`:

  Logical vector over the input rows of `raw.data`, `TRUE` for events
  inside the main-population gate, or `NULL` if no gating was applied.

- `recovery`:

  Data frame of angular errors against `true.spectra`. `recovered` is
  the fraction of the starting angular error removed,
  `(deg.start - deg.after) / deg.start` – `1` is fully recovered, `0` is
  no change, negative is worse; see `min.deg.start` for the
  near-zero-`deg.start` case. `NULL` if `true.spectra` was not supplied.
