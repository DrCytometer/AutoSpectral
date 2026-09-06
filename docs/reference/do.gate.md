# Do Gate

Performs gating on scatter parameters and returns a vector with the
points describing the gate boundary.

The gating proceeds in three steps:

- Defines bounds from the region of scatter space the data actually
  occupies, then trims those bounds further

- Defines a region around the target maximum found within the bounds

- Defines a gate around the target maximum, only within that region

The method uses numerical search of maxima over estimated densities and
Voronoi tessellations to improve density estimation around maxima.

When `region.auto` is `TRUE` in the relevant `default.gate.param.cells`
/ `default.gate.param.beads` entry (the default), the bounds are
computed with
[`get.scatter.occupancy()`](https://drcytometer.github.io/AutoSpectral/reference/get.scatter.occupancy.md),
a density-relative method that is robust to a handful of extreme events
and to cytometers with a much larger theoretical scatter range than any
one sample occupies. `scatter.data.min.x` / `scatter.data.max.x` /
`scatter.data.min.y` / `scatter.data.max.y` still apply as an outer
sanity ceiling, and become the working bounds again if `region.auto` is
set to `FALSE`.

## Usage

``` r
do.gate(
  gate.data,
  viability.gate,
  large.gate,
  samp,
  scatter.and.channel.label,
  control.type,
  asp,
  color.palette = "plasma",
  max.points = 50000,
  gate.color = "darkgoldenrod1"
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

- color.palette:

  Optional character string defining the viridis color palette to be
  used for the fluorophore traces. Default is `plasma`. Use `rainbow` to
  be similar to FlowJo or SpectroFlo. Other options are the viridis
  color options: `magma`, `inferno`, `plasma`, `viridis`, `cividis`,
  `rocket`, `mako` and `turbo`.

- max.points:

  Number of points to plot (speeds up plotting). Default is `5e4`.

- gate.color:

  Color to plot the gate boundary line, default is `darkgoldenrod1`.

## Value

A set of points describing the gate boundary.

## See also

[`get.scatter.occupancy()`](https://drcytometer.github.io/AutoSpectral/reference/get.scatter.occupancy.md)
