# Get Ungated Flow Expression Data

Retrieves flow cytometry expression data for specified samples, without
gating, removing out-of-range events.

## Usage

``` r
get.ungated.flow.expression.data(
  samp,
  file.name,
  control.dir,
  scatter.and.spectral.channel,
  spectral.channel,
  set.resolution
)
```

## Arguments

- samp:

  The sample identifier.

- file.name:

  A vector of file names for the samples.

- control.dir:

  The directory containing the control files.

- scatter.and.spectral.channel:

  A vector of scatter and spectral channels.

- spectral.channel:

  A vector of spectral channels.

- set.resolution:

  The resolution limit for the spectral channels.

## Value

A matrix with the flow expression data.
