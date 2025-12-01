# Run Downsample

This function runs the downsampling process on a list of samples, using
the specified peak channels and parameters.

## Usage

``` r
run.downsample(
  clean.expr.data,
  downsample.sample,
  peak.channels,
  negative.n = 500,
  positive.n = 1000,
  verbose = TRUE
)
```

## Arguments

- clean.expr.data:

  List containing cleaned expression data.

- downsample.sample:

  Vector of sample names to be downsampled.

- peak.channels:

  Named vector mapping samples to their peak channels.

- negative.n:

  Number of negative events to select, default `500`.

- positive.n:

  Number of positive events to select, default `1000`.

- verbose:

  Logical, default is `TRUE`. Set to `FALSE` to suppress messages.

## Value

A list containing the downsampled expression data for each sample.
