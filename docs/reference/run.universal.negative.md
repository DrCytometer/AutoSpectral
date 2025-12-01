# Run Universal Negative

This function processes universal negative samples to generate
expression data based on specified parameters.

## Usage

``` r
run.universal.negative(
  clean.expr,
  univ.sample,
  universal.negatives,
  scatter.param,
  peak.channels,
  downsample,
  negative.n,
  positive.n,
  spectral.channel,
  asp,
  control.type,
  scatter.match,
  intermediate.figures = FALSE,
  main.figures = TRUE,
  verbose = TRUE
)
```

## Arguments

- clean.expr:

  A matrix containing cleaned expression data.

- univ.sample:

  A character vector specifying the names of universal negative samples.

- universal.negatives:

  A list containing universal negative control parameters.

- scatter.param:

  A character vector specifying the scatter parameters.

- peak.channels:

  A character vector specifying the peak channels.

- downsample:

  A numeric value indicating the downsampling factor.

- negative.n:

  A numeric value indicating the number of negative events.

- positive.n:

  A numeric value indicating the number of positive events.

- spectral.channel:

  A character vector specifying the spectral channels.

- asp:

  The AutoSpectral parameter list. Prepare using
  `get.autospectral.param`

- control.type:

  A character string specifying the type of control: `beads` or `cells`

- scatter.match:

  A logical value indicating whether scatter matching is performed.

- intermediate.figures:

  Logical, if `TRUE` returns additional figures to show the inner
  workings of the cleaning, including definition of low-AF cell gates on
  the PCA-unmixed unstained and spectral ribbon plots of the AF
  exclusion from the unstained. Default is `FALSE` to speed up
  processing.

- main.figures:

  Logical, if `TRUE` creates the main figures to show the impact of
  intrusive autofluorescent event removal and scatter-matching for the
  negatives.

- verbose:

  Logical, default is `TRUE`. Set to `FALSE` to suppress messages.

## Value

A list containing the processed expression data for each universal
negative sample.
