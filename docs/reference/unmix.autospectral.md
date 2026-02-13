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
  use.dist0 = TRUE,
  verbose = TRUE,
  speed = c("fast", "medium", "slow"),
  parallel = TRUE,
  threads = NULL,
  n.variants = NULL,
  ...
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

- use.dist0:

  Logical, controls whether the selection of the optimal AF signature
  for each cell is determined by which unmixing brings the fluorophore
  signals closest to 0 (`use.dist0` = `TRUE`) or by which unmixing
  minimizes the per-cell residual (`use.dist0` = `FALSE`). Default is
  `TRUE`.

- verbose:

  Logical, default `TRUE`. Whether to send messages to the console.

- speed:

  Default is `fast`. Selector for the precision-speed trade-off in
  AutoSpectral per-cell fluorophore optimization. Options are `slow`,
  `medium` and `fast`. From v1.0.0, this controls the number of variants
  tested per cell (and per fluorophore). More variants takes longer, but
  gives better resolution in some unmixed data. When `speed = fast`, as
  single variant will be tested; for `medium`, three will be tested and
  for `slow`, 10 variants will be tested. From AutoSpectral v1.0.0, all
  options are available in the pure R version. Installation of
  `AutoSpectralRcpp` is strongly encouraged for speed, though.

- parallel:

  Logical, default is `TRUE`. The new parallel processing should always
  be faster.

- threads:

  Numeric, default is `NULL`, in which case `asp$worker.process.n` will
  be used. `asp$worker.process.n` is set by default to be one less than
  the available cores on the machine. Multi-threading is only used if
  `parallel` is `TRUE`.

- n.variants:

  Number of variants to test per cell. Allows explicit control over the
  number used, as opposed to `speed`, which selects from pre-defined
  choices. Providing a numeric value to `n.variants` will override
  `speed`, allowing up to `n.variants` (or the max available) variants
  to be tested. The default is `NULL`, in which case `n.variants` will
  be ignored.

- ...:

  Ignored. Previously used for deprecated arguments such as
  `calculate.error`.

## Value

Unmixed data with cells in rows and fluorophores in columns.
