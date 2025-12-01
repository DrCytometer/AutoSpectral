# Do Gate

Performs gating on scatter parameters and returns a vector with the
indexes of events inside the initial gate.

The gating proceeds in three steps:

- Defines bounds by data trimming

- Defines a region around the target maximum found within the bounds

- Defines a gate around the target maximum, only within that region

The method uses numerical search of maxima over estimated densities and
Voronoi tessellations to improve density estimation around maxima.

## Usage

``` r
do.gate(
  gate.data,
  viability.gate,
  large.gate,
  samp,
  scatter.and.channel.label,
  control.type,
  asp
)
```

## Arguments

- gate.data:

  A data frame containing the gate data.

- viability.gate:

  A logical vector indicating the viability gate.

- large.gate:

  A logical vector indicating the large gate.

- samp:

  A sample identifier.

- scatter.and.channel.label:

  A label for scatter and channel.

- control.type:

  The type of control used, either "beads" or "cells".

- asp:

  The AutoSpectral parameter list, prepared using
  `get.autospectral.param`.

## Value

A vector with the indexes of events inside the initial gate.
