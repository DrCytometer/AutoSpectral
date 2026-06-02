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

The output is saved as an .rds file and per-fluorophore variant plots
are produced if requested.

## Usage

``` r
get.spectral.variants.test(
  control.dir,
  control.def.file,
  asp,
  spectra,
  figures = TRUE,
  output.dir = NULL,
  parallel = FALSE,
  verbose = TRUE,
  threads = NULL,
  n.cells = 10000L,
  som.dim = 10L,
  k.neighbors = 3L,
  sim.threshold = 0.985,
  variant.fill.color = "red",
  variant.fill.alpha = 0.7,
  median.line.color = "black",
  median.linewidth = 1,
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

  Logical, default `FALSE`. Enable parallel processing across
  fluorophores.

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

  Integer, default `10`. Side length of the square SOM grid; up to
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

- variant.fill.color:

  Color for the shaded ribbon in variant plots. Default `"red"`.

- variant.fill.alpha:

  Alpha for `variant.fill.color`. Default `0.7`.

- median.line.color:

  Color for the reference-spectrum line. Default `"black"`.

- median.linewidth:

  Width of the reference-spectrum line. Default `1`.

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

The list is also saved as an .rds file in `output.dir`.
