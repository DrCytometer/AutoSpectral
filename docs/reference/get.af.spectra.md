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
  threads = NULL,
  som.dim = 10,
  figures = TRUE,
  plot.dir = NULL,
  table.dir = NULL,
  title = NULL
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

- threads:

  Numeric. Number of threads to use for parallel processing in the
  creation of the SOM. Default `NULL` reverts to `asp$worker.process.n`.

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

## Value

A matrix of autofluorescence spectra.
