# Get Spectral Variations for Fluorophores

Cycles through all fluorophores defined in `control.def.file`,
identifying variation in their spectral profiles via SOM clustering on
scatter-matched, per-event background-corrected data.

For each fluorophore the autofluorescence reference is derived **in
situ** from the paired universal-negative file (or internally from the
lower 25\\ used to project out autofluorescence and identify the
empirical peak detector. All positive events are scatter-matched to
unstained events and their per-event background is subtracted before SOM
clustering. This gives a comprehensive, population-level picture of true
fluorophore spectral variability without requiring a pre-computed
`af.spectra` matrix.

A cell-based unstained sample is still required at this stage, since it
anchors the positivity thresholds and per-node AF library. It is
normally read from the control file's `"AF"` row; if the control file
has no `"AF"` row (e.g. a bead-only or negative-free control setup),
supply one directly via `unstained.sample` instead.

The output is saved as an .rds file and per-fluorophore variant plots
are produced if requested.

Uses the Spillover Spreading Matrix built from
[`get.fluor.variants()`](https://drcytometer.github.io/AutoSpectral/reference/get.fluor.variants.md)'s
Residual Model regression rather than an empirical MAD ratio against a
separately-unmixed unstained baseline. There is no
`spread.denom.min.mad`/`"snr"` concept here: each control's own
regression already separates its baseline (intercept) from its
abundance-scaled spread (slope), estimated from its own dim-to-bright
event range rather than compared against a different sample's estimate.
A row is instead trusted once its regression used at least
`spread.min.events` events; below that, the hotspot-matrix fallback
(unchanged from `get.spectral.variants()`) fills the row when enough
trusted rows exist to calibrate against.

## Usage

``` r
get.spectral.variants(
  control.dir,
  control.def.file,
  asp,
  spectra,
  figures = TRUE,
  output.dir = NULL,
  parallel = TRUE,
  verbose = TRUE,
  threads = NULL,
  n.cells = 10000L,
  som.dim = 5L,
  k.neighbors = 3L,
  sim.threshold = 0.985,
  sim.threshold.floor = 0.9,
  af.collinear.threshold = 0.95,
  noise.floor.tail.fraction = 0.2,
  spread.min.events = 50L,
  spread.hotspot.fallback = TRUE,
  spread.hotspot.min.pairs = 20,
  huber.k = 1.345,
  huber.max.iter = 100L,
  variant.fill.color = "red",
  variant.fill.alpha = 0.7,
  median.line.color = "black",
  median.linewidth = 1,
  use.unmixed = TRUE,
  unstained.sample = NULL,
  stained.sample = NULL,
  optimize.necessity.threshold = 0.01,
  diagnostics = FALSE,
  ...
)
```

## Arguments

- control.dir:

  Character. Path to the single-stained control FCS files.

- control.def.file:

  Character. Path to the control definition CSV. Must pass
  [`check.control.file()`](https://drcytometer.github.io/AutoSpectral/reference/check.control.file.md).

- asp:

  The AutoSpectral parameter list from
  [`get.autospectral.param()`](https://drcytometer.github.io/AutoSpectral/reference/get.autospectral.param.md).

- spectra:

  Numeric matrix. Reference spectra; fluorophores in rows, detectors in
  columns.

- figures:

  Logical, default `TRUE`. Whether to save variant-spectrum plots.

- output.dir:

  Character or `NULL`. Directory for figures and the .rds output file.
  Defaults to `asp$variant.dir`.

- parallel:

  Logical, default `TRUE`. Enable parallel processing for SOM clustering
  (requires AutoSpectralRcpp).

- verbose:

  Logical, default `TRUE`. Set to `FALSE` to suppress messages.

- threads:

  Numeric or `NULL`. Number of parallel workers. Defaults to
  `asp$worker.process.n`.

- n.cells:

  Integer, default `10000`. Maximum positive events per fluorophore used
  for SOM clustering. Passed to `get.fluor.variants`.

- som.dim:

  Integer, default `5`. Side length of the square SOM grid; up to
  `som.dim^2` candidate variants per fluorophore before cosine QC.
  Passed to `get.fluor.variants`.

- k.neighbors:

  Integer, default `3`. Number of scatter-space nearest neighbours from
  the unstained pool used to estimate per-event background. Passed to
  `get.fluor.variants`.

- sim.threshold:

  Numeric, default `0.99`. Minimum cosine similarity to the reference
  spectrum for a SOM centroid to be retained as a variant. Passed to
  `get.fluor.variants`.

- sim.threshold.floor:

  Numeric, default `0.90`. Lower bound for adaptive relaxation of
  `sim.threshold` when the initial cutoff retains fewer than 20 events.

- af.collinear.threshold:

  Numeric, default `0.95`. Minimum cosine similarity between a
  fluorophore's reference spectrum and any of its paired unstained
  file's AF principal directions at or above which the AF-component
  projection step (and the low-rank fit feeding the Residual Model) is
  skipped.

- noise.floor.tail.fraction:

  Numeric in (0, 1), default `0.20`. Fraction of each detector's raw
  values (lowest end) used to estimate the per-control noise floor.
  Passed to `get.fluor.variants`.

- spread.min.events:

  Integer, default `50`. Minimum number of events behind a source
  fluorophore's Residual Model regression (`"spillover.spread.n"`)
  before its Spillover Spreading Matrix row is trusted. A 2-parameter
  regression needs more support than the old MAD point estimate did;
  rows below this are left blank unless filled by the hotspot-matrix
  fallback.

- spread.hotspot.fallback:

  Logical, default `TRUE`. When a source fluorophore's Spillover
  Spreading Matrix row fails the `spread.min.events` check (a weak or
  under-titrated control), fill that row from
  `calculate.hotspot.matrix(spectra)` instead of leaving it blank. The
  hotspot matrix is a purely geometric measure of pairwise spread
  susceptibility from the reference spectra alone; filling uses a single
  calibration constant (the median ratio of trusted `spillover.spread`
  entries to their hotspot-matrix counterparts), so it requires at least
  `spread.hotspot.min.pairs` trusted entries to calibrate against.
  Filled rows are tagged `"hotspot"` in the returned matrix's `"source"`
  attribute.

- spread.hotspot.min.pairs:

  Integer, default `20`. Minimum number of trusted (source, target)
  entries required to calibrate the hotspot-matrix fallback. Below this,
  weak controls are left blank as before and a message explains why.

- huber.k:

  Numeric, default `1.345`. Huber tuning constant passed to
  [`get.fluor.variants()`](https://drcytometer.github.io/AutoSpectral/reference/get.fluor.variants.md)'s
  spillover-spread regression.

- huber.max.iter:

  Integer, default `100L`. Maximum IRLS iterations passed to the same
  regression.

- variant.fill.color:

  Color for the shaded ribbon in variant plots. Default `"red"`.

- variant.fill.alpha:

  Alpha for `variant.fill.color`. Default `0.7`.

- median.line.color:

  Color for the reference-spectrum line. Default `"black"`.

- median.linewidth:

  Width of the reference-spectrum line. Default `1`.

- use.unmixed:

  Logical, default `TRUE`. Whether AF extraction and fluorophore variant
  assessment may use full-spectra OLS unmixing as part of their SOM
  clustering input, positivity selection, and Spillover Spreading Matrix
  construction. Set to `FALSE` when `spectra` contains several similar
  or collinear fluorophores (e.g. a bead-cell comparison panel). When
  `FALSE`, the returned `spillover.spread` is always `NULL`.

- unstained.sample:

  Optional file path to a cell-based unstained FCS file, used as the
  autofluorescence reference when the control file has no `"AF"` row.

- stained.sample:

  Optional file path to a representative stained FCS file, weighting the
  optimization necessity scores by fluorophore brightness. Pass `NULL`
  (default) to use purely geometric scores.

- optimize.necessity.threshold:

  Numeric in `[0, 1]`, default `0.01`. Passed to
  [`calculate.optimize.necessity()`](https://drcytometer.github.io/AutoSpectral/reference/calculate.optimize.necessity.md).

- diagnostics:

  Logical, default `FALSE`. When `TRUE`, prints additional messages on
  model characteristics to the console.

- ...:

  Ignored. Catches and warns on previously used deprecated arguments:
  `af.spectra`, `refine`, `problem.quantile`, `pos.quantile`.

## Value

A named list with elements:

- `thresholds`:

  Named numeric vector of positivity thresholds in the unmixed space,
  one per fluorophore.

- `neg.thresholds`:

  Named numeric vector, the 0.5th percentile of each fluorophore's
  unstained unmixed distribution.

- `variants`:

  Named list of variant-spectra matrices, one per fluorophore.

- `delta.list`:

  Named list of delta matrices (variant minus reference spectrum), one
  per fluorophore.

- `delta.norms`:

  Named list of Euclidean norms of the deltas.

- `noise.floor`:

  Named numeric vector, per-detector electronic noise floor in signal
  units (SD), pooled by minimum across controls.

- `spillover.spread`:

  Matrix (source fluorophore x target channel), the Residual Model
  Spillover Spreading Matrix: the Huber-robust slope of each target's
  squared residual-projection against the source's own recovered
  abundance – added unmixed variance per unit of the source's on-channel
  abundance. Diagonal entries are `NA`. Rows below `spread.min.events`
  are left `NA` across the row unless filled by the hotspot-matrix
  fallback. Carries `"n.events"` and `"source"` attributes (named
  integer/character vectors, one per source fluorophore): event count
  behind the regression, and whether the row came from the regression
  (`"residual"`) or the hotspot-matrix fallback (`"hotspot"`). `NULL` if
  no control supplied enough positive events. Saved as a heatmap when
  `figures = TRUE`.

- `spillover.spread.intercept`:

  Matrix, same shape as `spillover.spread`: each source's regression
  intercept per target channel – that source's own estimate of the
  target's baseline unmixed variance at zero abundance. Not filled by
  the hotspot fallback (the hotspot matrix has no calibrated intercept
  term); `NA` wherever `spillover.spread` came from that fallback or was
  left blank.

The list is also saved as an .rds file in `output.dir`.

## References

Cai X et al. (2026). Residual Model for unmixed spread prediction in
spectral flow cytometry. *bioRxiv* 2026.01.27.701929.
