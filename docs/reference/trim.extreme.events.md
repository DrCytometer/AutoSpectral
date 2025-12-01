# Trim Extreme Events

This function trims extreme events from expression data based on a
specified peak channel and trim factor. Not recommended for most panels.

## Usage

``` r
trim.extreme.events(expr.data, peak.channel, trim.factor)
```

## Arguments

- expr.data:

  A matrix containing the expression data.

- peak.channel:

  A character string specifying the peak channel to be used for
  trimming.

- trim.factor:

  A numeric value indicating the proportion of extreme events to trim
  from both ends of the peak channel distribution.

## Value

A matrix with the extreme events trimmed.
