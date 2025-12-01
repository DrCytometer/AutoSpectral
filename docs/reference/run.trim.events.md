# Run Trim Events

This function trims extreme events from multiple samples based on
specified peak channels and trim factors.

## Usage

``` r
run.trim.events(
  trim.sample.data,
  trim.sample,
  trim.peak.channels,
  trim.factor,
  asp
)
```

## Arguments

- trim.sample.data:

  A list containing the expression data for each sample.

- trim.sample:

  A character vector specifying the names of the samples.

- trim.peak.channels:

  A list containing the peak channels for each sample.

- trim.factor:

  A numeric value indicating the proportion of extreme events to trim
  from both ends of the peak channel distribution.

- asp:

  The AutoSpectral parameter list. Prepare using
  `get.autospectral.param`

## Value

A list containing the trimmed expression data for each sample.
