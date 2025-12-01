# Get Universal Negative Control

This function identifies and processes the universal negative control
for a given sample, including scatter matching and plotting.

## Usage

``` r
get.universal.negative(
  clean.expr.data,
  samp,
  universal.negatives,
  scatter.param,
  peak.channels,
  downsample,
  negative.n,
  positive.n,
  spectral.channel,
  asp,
  control.type,
  scatter.match = TRUE,
  intermediate.figures = FALSE,
  main.figures = TRUE,
  verbose = TRUE
)
```

## Arguments

- clean.expr.data:

  List containing cleaned expression data.

- samp:

  Sample identifier.

- universal.negatives:

  Named vector mapping samples to their universal negatives.

- scatter.param:

  Vector of scatter parameters.

- peak.channels:

  Named vector mapping samples to their peak channels.

- downsample:

  Logical indicating whether to downsample the data.

- negative.n:

  Number of negative events to select.

- positive.n:

  Number of positive events to select.

- spectral.channel:

  Vector of spectral channel names.

- asp:

  The AutoSpectral parameter list. Prepare using
  `get.autospectral.param`

- control.type:

  Named vector mapping samples to their control types.

- scatter.match:

  Logical indicating whether to perform scatter matching. Default is
  `TRUE`.

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

A data frame containing the selected positive and scatter-matched
negative events.
