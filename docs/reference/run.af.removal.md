# Run Autofluorescence Removal

This function runs the autofluorescence removal process on a list of
samples, using the specified parameters and settings.

## Usage

``` r
run.af.removal(
  clean.expr,
  af.removal.sample,
  spectral.channel,
  peak.channel,
  universal.negative,
  asp,
  scatter.param,
  negative.n = 500,
  positive.n = 1000,
  scatter.match = TRUE,
  intermediate.figures = FALSE,
  main.figures = TRUE,
  parallel = FALSE,
  threads = 1,
  verbose = TRUE
)
```

## Arguments

- clean.expr:

  List containing cleaned expression data.

- af.removal.sample:

  Vector of sample names for which autofluorescence removal is to be
  performed.

- spectral.channel:

  Vector of spectral channel names.

- peak.channel:

  Vector of peak detection channels for fluorophores.

- universal.negative:

  Name of the universal negative control.

- asp:

  The AutoSpectral parameter list. Prepare using
  `get.autospectral.param`

- scatter.param:

  Vector of scatter parameters.

- negative.n:

  Integer. Number of events to include in the downsampled negative
  population. Default is `500`.

- positive.n:

  Integer. Number of events to include in the downsampled positive
  population. Default is `1000`.

- scatter.match:

  Logical, default is `TRUE`. Whether to select negative events based on
  scatter profiles matching the positive events. Defines a region of FSC
  and SSC based on the distribution of selected positive events.

- intermediate.figures:

  Logical, if `TRUE` returns additional figures to show the inner
  workings of the cleaning, including definition of low-AF cell gates on
  the PCA-unmixed unstained and spectral ribbon plots of the AF
  exclusion from the unstained.

- main.figures:

  Logical, if `TRUE` creates the main figures to show the impact of
  intrusive autofluorescent event removal and scatter-matching for the
  negatives.

- parallel:

  Logical, default is `FALSE`, in which case parallel processing will
  not be used. Set to `TRUE` to run in parallel.

- threads:

  Number of cores to use for parallel processing, default is `1`.

- verbose:

  Logical, default is `TRUE`. Set to `FALSE` to suppress messages.

## Value

A list containing the expression data with autofluorescent events
removed for each sample.
