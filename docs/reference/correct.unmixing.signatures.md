# Correct Unmixing Signatures

Detects and corrects shape errors in a reference `spectra` matrix
directly from stained sample data, for the situation where the available
reference spectra are subtly wrong for the sample being analysed -
typically because bead-derived spectra are applied to cells, or spectra
from another day, lot, or instrument state are in use.

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
  raw.data,
  spectra,
  unmixed.thresholds,
  asp = NULL,
  af.name = "AF",
  scatter = NULL,
  gate.main = TRUE,
  gate.level = 0.1,
  spillover.spread = NULL,
  spread.kappa = 2,
  bg.mode = c("global.mean", "scatter.knn", "none"),
  unstained = NULL,
  unstained.scatter = NULL,
  k.neighbors = 20L,
  n.levels = 10L,
  n.iter = 6L,
  min.events = 200L,
  min.span = 5,
  min.explained = 0.5,
  min.gain = 0.002,
  step.grid = c(0, 0.03125, 0.0625, 0.125, 0.25, 0.5, 1),
  max.step = 0.08,
  max.span.drift = 1.1,
  max.bg.alignment = -0.9,
  nuisance.frac = 0.5,
  background.n = 5000L,
  true.spectra = NULL,
  verbose = TRUE
)
```

## Arguments

- raw.data:

  Numeric matrix (events x detectors), raw detector-space data. Pooled
  or concatenated single-stained controls, or a fully stained sample
  with well-separated populations. Columns must match the columns of
  `spectra`.

- spectra:

  Numeric matrix (fluorophores x detectors), the starting reference
  spectra to be corrected, L-infinity normalised.

- unmixed.thresholds:

  Named numeric vector covering every fluorophore in `spectra` (the
  autofluorescence row may be omitted), giving the positivity threshold
  in unmixed space, typically the 99.5th percentile of an unstained
  control unmixed against `spectra`.

- asp:

  Optional AutoSpectral parameter list from
  [`get.autospectral.param()`](https://drcytometer.github.io/AutoSpectral/reference/get.autospectral.param.md).
  Used only to seed the random number generator (`asp$bird.seed`) for
  reproducible subsampling. Default `NULL`.

- af.name:

  Character, the name of the autofluorescence row in `spectra`, or
  `NULL` if there is none. The AF row is never treated as a panel
  fluorophore and is never corrected. Default `"AF"`.

- scatter:

  Optional numeric matrix (events x scatter parameters), row-matched to
  `raw.data`. Required for `gate.main` and for
  `bg.mode = "scatter.knn"`. Default `NULL`.

- gate.main:

  Logical, whether to gate events to the main scatter population before
  fitting. Requires `scatter`. Default `TRUE`.

- gate.level:

  Numeric in (0, 1). Events are kept when their 2D scatter density
  exceeds this fraction of the modal density. Default `0.1`.

- spillover.spread:

  Optional matrix from `get.spectral.variants()$spillover.spread`. When
  supplied, co-activity is judged against per-event thresholds widened
  by the spillover spread each bright fluorophore contributes, so
  spillover from the dominant dye is not mistaken for a co-active
  fluorophore. Default `NULL` (flat thresholds).

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
  Default `20`.

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

- max.step:

  Numeric, maximum norm of the correction relative to the norm of the
  row it corrects. Larger proposed corrections are rejected outright
  rather than scaled down. Default `0.08`.

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

- background.n:

  Integer, maximum background events used for the zero-abundance anchor.
  Default `5000`.

- true.spectra:

  Optional numeric matrix (fluorophores x detectors),
  independently-known ground truth with row names matching `spectra`.
  Purely diagnostic: when supplied, the returned `recovery` table
  reports the angular error before and after correction per fluorophore.

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

  Data frame of angular errors against `true.spectra`, or `NULL` if
  `true.spectra` was not supplied.
