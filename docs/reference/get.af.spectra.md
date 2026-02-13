# Get Autofluorescence Spectra

Extracts autofluorescence spectra from an unstained samples. Intended
for use with `unmix.autospectral`. Uses FlowSOM (EmbedSOM) clustering
for rapid identification of cells with similar AF profiles.

## Usage

``` r
get.af.spectra(
  unstained.sample,
  asp,
  spectra,
  som.dim = 10,
  figures = TRUE,
  plot.dir = NULL,
  table.dir = NULL,
  title = "Autofluorescence spectra",
  verbose = TRUE,
  refine = FALSE,
  problem.quantile = 0.99
)
```

## Arguments

- unstained.sample:

  Path and file name for a unstained sample FCS file. The sample type
  and processing (protocol) method should match the fully stained
  samples to which the AF will be applied, ideally.

- asp:

  The AutoSpectral parameter list.

- spectra:

  Spectral signatures of fluorophores, normalized between 0 and 1, with
  fluorophores in rows and detectors in columns.

- som.dim:

  Number of x and y dimensions for the SOM. Default is `10`.

- figures:

  Logical, whether to plot the spectral traces and heatmap for the AF
  signatures. Default is `TRUE`.

- plot.dir:

  Directory (folder) where the plots will be saved. Default is `NULL`,
  which inherits from `asp$figure.af.dir`.

- table.dir:

  Directory (folder) where the spectra csv file will be saved. Default
  is `NULL`, which inherits from `asp$table.af.dir`.

- title:

  Title for the output spectral plots and csv file. Default is
  `Autofluorescence spectra`.

- verbose:

  Logical, controls messaging. Default is `TRUE`.

- refine:

  Logical, default is `FALSE`. Controls whether to perform a second
  round of autofluorescence measurement on "problem cells", which are
  those with the highest spillover, as defined by `problem.quantile`.
  When `FALSE`, behavior is identical to versions of AutoSpectral prior
  to 1.0.0. If you are working with samples containing complex
  autofluorescence, e.g., tissues or tumors, using `refine=TRUE` will
  improve autofluorescence extraction in the unmixing at the cost of an
  increase in unmixing time. The increase in time will depend on the
  method used to assign autofluorescence spectra per cell (residual
  based assignment is very fast) and whether you have installed
  `AutoSpectralRcpp`, which will speed up assignment and unmixing.

- problem.quantile:

  Numeric, default `0.99`. The quantile for determining which cells will
  be considered "problematic" after unmixing with per-cell AF
  extraction. Cells in the `problem.quantile` or above with respect to
  total signal in the fluorophore (non-AF) channels after per-cell AF
  extraction will be used to determine additional autofluorescence
  spectra, using a second round of clustering and modulation of the
  previously selected autofluorescence spectra. A value of `0.99` means
  the top 1% of cells, those farthest from zero, will be selected for
  further investigation.

## Value

A matrix of autofluorescence spectra.
