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
  variant.fill.color = "red",
  variant.fill.alpha = 0.7,
  median.line.color = "black",
  median.linewidth = 1,
  use.unmixed = TRUE,
  unstained.sample = NULL,
  stained.sample = NULL,
  optimize.necessity.threshold = 0.01,
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
  for SOM clustering. Files with more events above threshold are
  randomly downsampled. Passed to `get.fluor.variants`.

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
  Relaxation is logged via
  [`warning()`](https://rdrr.io/r/base/warning.html) and the threshold
  actually used is returned as the `"cosine.threshold.used"` attribute.

- af.collinear.threshold:

  Numeric, default `0.95`. Minimum cosine similarity between `fluor`'s
  reference spectrum and any of its paired unstained file's AF principal
  directions (`af.pcs`) at or above which the AF-component projection
  step is skipped, since a joint OLS fit against near-collinear AF and
  fluorophore directions can push real fluorophore signal into the AF
  term. Recorded as the `"af.collinear"` attribute.

- noise.floor.tail.fraction:

  Numeric in (0, 1), default `0.20`. Fraction of each detector's raw
  values (lowest end) used to estimate the per-control noise floor.
  Passed to `get.fluor.variants`.

- variant.fill.color:

  Color for the shaded ribbon in variant plots. Default `"red"`.

- variant.fill.alpha:

  Alpha for `variant.fill.color`. Default `0.7`.

- median.line.color:

  Color for the reference-spectrum line. Default `"black"`.

- median.linewidth:

  Width of the reference-spectrum line. Default `1`.

- use.unmixed:

  Logical, default `TRUE`. Whether AF extraction
  ([`get.af.spectra()`](https://drcytometer.github.io/AutoSpectral/reference/get.af.spectra.md))
  and fluorophore variant assessment
  ([`get.fluor.variants()`](https://drcytometer.github.io/AutoSpectral/reference/get.fluor.variants.md))
  may use full-spectra OLS unmixing as part of their SOM clustering
  input, positivity selection, and Spillover Spreading Matrix
  construction. Set to `FALSE` when `spectra` contains several similar
  or collinear fluorophores (e.g. a bead-cell comparison panel), where a
  full-spectra unmix is itself unstable or unsolvable and would corrupt
  rather than inform those steps. When `FALSE`, clustering falls back to
  raw detector space only, the unstained-sample positivity thresholds
  used internally by
  [`get.fluor.variants()`](https://drcytometer.github.io/AutoSpectral/reference/get.fluor.variants.md)
  are not computed, and the returned `spillover.spread` is always
  `NULL`.

- unstained.sample:

  Optional file path to a cell-based unstained FCS file, used as the
  autofluorescence reference when the control file has no `"AF"` row.
  Required in that case; ignored (with a message) if the control file
  does have an `"AF"` row, since the in-situ unstained paired with the
  single-stained controls is used instead.

- stained.sample:

  Optional file path to a representative stained FCS file. When
  supplied, it is read and unmixed to obtain per-fluorophore median
  positive signal (MFI), which weights the optimization necessity scores
  by fluorophore brightness. Pass `NULL` (default) to use purely
  geometric scores.

- optimize.necessity.threshold:

  Numeric in `[0, 1]`, default `0.01`. Passed to
  [`calculate.optimize.necessity()`](https://drcytometer.github.io/AutoSpectral/reference/calculate.optimize.necessity.md).
  Fluorophores whose normalised leakage score falls below this value are
  flagged as not requiring per-cell spectral optimisation. The result is
  stored in `$optimize.recommended` in the returned list and used
  automatically by `unmix.autospectral.rcpp()` to skip unnecessary
  optimisation passes.

- ...:

  Ignored. Catches and warns on previously used deprecated arguments:
  `af.spectra`, `refine`, `problem.quantile`, `pos.quantile`.

## Value

A named list with elements:

- `thresholds`:

  Named numeric vector of positivity thresholds in the unmixed space,
  one per fluorophore.

- `variants`:

  Named list of variant-spectra matrices, one per fluorophore. Each
  matrix has variants in rows and detectors in columns.

- `delta.list`:

  Named list of delta matrices (variant minus reference spectrum), one
  per fluorophore.

- `delta.norms`:

  Named list of Euclidean norms of the deltas, one numeric vector per
  fluorophore.

- `noise.floor`:

  Named numeric vector, per-detector electronic noise floor in signal
  units (SD), pooled by minimum across controls. Matches the units of
  `noise.floor` elsewhere in the package
  ([`unmix.fcs()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.fcs.md),
  the C++ pipeline). Square it before passing to
  `estimate.noise.model(read.var.floor = ...)`, which expects a
  variance.

- `spillover.spread`:

  Matrix (source fluorophore x target channel), the Spillover Spreading
  Matrix: increase in unmixed variance a source fluorophore's positive
  population contributes to each other channel, per unit of its own
  on-channel signal. Diagonal entries are `NA`. `NULL` if no control
  supplied enough positive events. Saved as a heatmap when
  `figures = TRUE`.

The list is also saved as an .rds file in `output.dir`.
