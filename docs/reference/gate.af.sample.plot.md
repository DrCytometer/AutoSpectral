# Plot Gate Autofluorescence Sample

This function plots the gate autofluorescence sample, including upper
and lower boundaries, using ggplot2 and other necessary packages.

## Usage

``` r
gate.af.sample.plot(
  plot.data,
  samp,
  af.boundary.upper,
  asp,
  max.points = 1e+05,
  color.palette = "viridis"
)
```

## Arguments

- plot.data:

  Matrix containing autofluorescence data points.

- samp:

  Sample identifier.

- af.boundary.upper:

  Matrix containing upper boundary information.

- asp:

  The AutoSpectral parameter list. Prepare using
  `get.autospectral.param`

- max.points:

  Number of points to plot (speeds up plotting). Default is `1e5`.

- color.palette:

  Optional character string defining the viridis color palette to be
  used for the fluorophore traces. Default is `viridis`. Options are the
  viridis color options: `magma`, `inferno`, `plasma`, `viridis`,
  `cividis`, `rocket`, `mako` and `turbo`.

## Value

Saves the plot as a JPEG file in the specified directory.
