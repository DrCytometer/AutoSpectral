# Plot Pre-defined Gate on Sample

This function plots a pre-defined gate on a sample, using ggplot2 and
other necessary packages.

## Usage

``` r
gate.sample.plot(
  samp,
  gate.data,
  gate.marker,
  gate.boundary,
  scatter.and.channel.label,
  control.type,
  asp,
  color.palette = "mako",
  max.points = 50000
)
```

## Arguments

- samp:

  Sample identifier.

- gate.data:

  Matrix containing gate data points.

- gate.marker:

  Vector containing gate marker names.

- gate.boundary:

  List containing gate boundary information.

- scatter.and.channel.label:

  Named vector mapping scatter and channel labels.

- control.type:

  Type of control: `beads` or `cells`. Deprecated.

- asp:

  The AutoSpectral parameter list. Prepare using
  `get.autospectral.param`

- color.palette:

  Optional character string defining the viridis color palette to be
  used for the fluorophore traces. Default is `rainbow`, which will be
  similar to FlowJo or SpectroFlo. Other pptions are the viridis color
  options: `magma`, `inferno`, `plasma`, `viridis`, `cividis`, `rocket`,
  `mako` and `turbo`.

- max.points:

  Number of points to plot (speeds up plotting). Default is `5e4`.

## Value

Saves the plot as a JPEG file in the specified directory.
