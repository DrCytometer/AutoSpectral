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
  not be used. Parallel processing will likely be faster when many small
  files are read in. If the data is larger, parallel processing may not
  accelerate the process much.

- verbose:

  Logical, default is `TRUE`. Set to `FALSE` to suppress messages.

## Value

A list containing the cleaned expression data for each sample.
