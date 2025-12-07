# Unmix AutoSpectral

Unmix using the AutoSpectral method to extract autofluorescence and
optimize fluorophore signatures at the single cell level.

## Usage

``` r
unmix.autospectral(
  raw.data,
  spectra,
  af.spectra,
  asp,
  spectra.variants = NULL,
  weighted = FALSE,
  weights = NULL,
  calculate.error = FALSE,
  use.dist0 = TRUE,
  verbose = TRUE,
  parallel = TRUE,
  threads = NULL
)
```

## Arguments

- raw.data:

  Expression data from raw fcs files. Cells in rows and detectors in
  columns. Columns should be fluorescent data only and must match the
  columns in spectra.

- spectra:

  Spectral signatures of fluorophores, normalized between 0 and 1, with
  fluorophores in rows and detectors in columns.

- af.spectra:

  Spectral signatures of autofluorescences, normalized between 0 and 1,
  with fluorophores in rows and detectors in columns. Prepare using
  `get.af.spectra`.

- asp:

  The AutoSpectral parameter list.

- spectra.variants:

  Named list (names are fluorophores) carrying matrices of spectral
  signature variations for each fluorophore. Prepare using
  `get.spectral.variants`. Default is `NULL`.

- weighted:

  Logical, whether to use ordinary or weighted least squares unmixing as
  the base algorithm. Default is `FALSE` and will use OLS.

- weights:

  Optional numeric vector of weights (one per fluorescent detector).
  Default is `NULL`, in which case weighting will be done by channel
  means (Poisson variance). Only used if `weighted`.

- calculate.error:

  Logical, whether to calculate the RMSE unmixing model accuracy and
  include it as an output. Default is `FALSE`.

- use.dist0:

  Logical, controls whether the selection of the optimal AF signature
  for each cell is determined by which unmixing brings the fluorophore
  signals closest to 0 (`use.dist0` = `TRUE`) or by which unmixing
  minimizes the per-cell residual (`use.dist0` = `FALSE`). Default is
  `TRUE`.

- verbose:

  Logical, whether to send messages to the console. Default is `TRUE`.

- parallel:

  Logical, default is `FALSE`, in which case sequential processing will
  be used. The new parallel processing should always be faster.

- threads:

  Numeric, default is `NULL`, in which case `asp$worker.process.n` will
  be used. `asp$worker.process.n` is set by default to be one less than
  the available cores on the machine. Multi-threading is only used if
  `parallel` is `TRUE`.

## Value

Unmixed data with cells in rows and fluorophores in columns.
