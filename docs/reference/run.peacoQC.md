# Run PeacoQC

This function runs PeacoQC to remove flow fluctuation errors from
expression data using parallel processing if specified.

## Usage

``` r
run.peacoQC(
  expr.data,
  spectral.channel,
  all.channels,
  asp,
  figures = TRUE,
  parallel = FALSE,
  threads = 1,
  verbose = TRUE
)
```

## Arguments

- expr.data:

  A list containing the expression data for each sample.

- spectral.channel:

  A character vector specifying the spectral channels.

- all.channels:

  A character vector specifying all channels.

- asp:

  The AutoSpectral parameter list. Prepare using
  `get.autospectral.param`

- figures:

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

A list containing the cleaned expression data for each sample.
