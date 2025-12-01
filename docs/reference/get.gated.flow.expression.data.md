# Get Gated Flow Expression Data

Retrieves gated flow cytometry expression data for specified samples,
removing out-of-range events and applying gating boundaries.

## Usage

``` r
get.gated.flow.expression.data(
  samp,
  file.name,
  control.dir,
  scatter.and.spectral.channel,
  spectral.channel,
  set.resolution,
  flow.gate,
  gate.list,
  scatter.param,
  scatter.and.channel.label,
  asp
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

- flow.gate:

  A list of flow gates for the samples.

- gate.list:

  A list of gating boundaries.

- scatter.param:

  A vector of scatter parameters.

- scatter.and.channel.label:

  A label for scatter and channel.

- asp:

  The AutoSpectral parameter list. Prepare using
  `get.autospectral.param`

## Value

A matrix with the gated expression data.
