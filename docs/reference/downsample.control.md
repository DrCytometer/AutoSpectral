# Downsample Control Data

Downsamples control data by selecting a specified number of positive and
negative events based on peak channel values.

## Usage

``` r
downsample.control(
  clean.expr.data,
  samp,
  peak.channels,
  negative.n = 500,
  positive.n = 1000,
  verbose = TRUE
)
```

## Arguments

- clean.expr.data:

  A list containing cleaned expression data for each sample.

- samp:

  The sample identifier.

- peak.channels:

  A vector of peak channels for the samples.

- negative.n:

  Number of negative events to select, default `500`.

- positive.n:

  Number of positive events to select, default `1000`.

- verbose:

  Logical. Default is `TRUE`.

## Value

A matrix with the selected expression data.
