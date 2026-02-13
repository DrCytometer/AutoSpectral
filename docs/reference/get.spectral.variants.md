# Get Spectral Variations for Fluorophores

This function cycles through all the fluorophores defined in
`control.def.file`, identifying variations in spectral profiles. It does
this by performing SOM clustering on the positive events in the cleaned
control data. The output is saved as an .rds file, and figures
summarizing the variation are saved, if desired. Note that the .rds file
contains all the needed information for downstream processing (per-cell
unmixing), so you can just load that using the `readRDS` function)
rather than re-running this process.

## Usage

``` r
get.spectral.variants(
  control.dir,
  control.def.file,
  asp,
  spectra,
  af.spectra,
  n.cells = 2000,
  som.dim = 10,
  figures = TRUE,
  output.dir = NULL,
  parallel = FALSE,
  verbose = TRUE,
  threads = NULL,
  refine = FALSE,
  problem.quantile = 0.95,
  ...
)
```

## Arguments

- control.dir:

  File path to the single-stained control FCS files.

- control.def.file:

  CSV file defining the single-color control file names, fluorophores
  they represent, marker names, peak channels, and gating requirements.

- asp:

  The AutoSpectral parameter list. Prepare using
  `get.autospectral.param`

- spectra:

  A matrix containing the spectral data. Fluorophores in rows, detectors
  in columns.

- af.spectra:

  Spectral signatures of autofluorescences, normalized between 0 and 1,
  with fluorophores in rows and detectors in columns. Prepare using
  `get.af.spectra`.

- n.cells:

  Numeric, default `2000`. Number of cells to use for defining the
  variation in spectra. Up to `n.cells` cells will be selected as
  positive events in the peak channel for each fluorophore, above the
  99.5th percentile level in the unstained sample.

- som.dim:

  Numeric, default `10`. Number of x and y dimensions to use in the SOM
  for clustering the spectral variation. The number of spectra returned
  for each fluorophore will increase with the quadratic of `som.dim`, so
  for 10, you will get up to 100 variants. Somewhere between 4 and 7
  appears to be sufficient, but with the pruning of variants implemented
  in
  [`unmix.autospectral()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.autospectral.md)
  in v1.0.0, this is less important.

- figures:

  Logical, controls whether the variation in spectra for each
  fluorophore is plotted in `output.dir`. Default is `TRUE`.

- output.dir:

  File path to whether the figures and .rds data file will be saved.
  Default is `NULL`, in which case `asp$variant.dir` will be used.

- parallel:

  Logical, default is `FALSE`, in which case sequential processing will
  be used. The new parallel processing should always be faster.

- verbose:

  Logical, default is `TRUE`. Set to `FALSE` to suppress messages.

- threads:

  Numeric, default is `NULL`, in which case `asp$worker.process.n` will
  be used. `asp$worker.process.n` is set by default to be one less than
  the available cores on the machine. Multi-threading is only used if
  `parallel` is `TRUE`.

- refine:

  Logical, default is `FALSE`. Controls whether to perform a second
  round of variation measurement on "problem cells", which are those
  with the highest spillover, as defined by `problem.quantile`. When
  `FALSE`, behavior is identical to versions of AutoSpectral prior to
  1.0.0. Setting to `TRUE` may help reduce spillover spread and unmixing
  errors. Using `refine=TRUE` does not impact subsequent unmixing
  calculation time in any significant way, unlike the same setting in
  [`get.af.spectra()`](https://drcytometer.github.io/AutoSpectral/reference/get.af.spectra.md).

- problem.quantile:

  Numeric, default `0.95`. The quantile for determining which cells will
  be considered "problematic" after unmixing with per-cell AF
  extraction. Cells in the `problem.quantile` or above with respect to
  total signal in the fluorophore (non-AF) channels after per-cell AF
  extraction will be used to determine additional autofluorescence
  spectra, using a second round of clustering and modulation of the
  previously selected autofluroescence spectra. A value of `0.95` means
  the top 5% of cells, those farthest from zero, will be selected for
  further investigation.

- ...:

  Ignored. Previously used for deprecated arguments such as
  `pos.quantile` and `sim.threshold`, which are now fixed internally and
  no longer user-settable.

## Value

A vector with the indexes of events inside the initial gate.
