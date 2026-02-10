# Plot Autofluorescence Identification Gate

This function plots the sample being used to identify intrusive
autofluorescence in the single-stained controls. The input data are
expected to be PCA projections of the unstained sample with an
accompanying region to identify the low-autofluorescence cell region.

## Usage

``` r
gate.af.identify.plot(
  gate.data,
  samp,
  gate.region,
  gate.bound.density,
  asp,
  color.palette = "rainbow",
  max.points = 1e+05
)
```

## Arguments

- gate.data:

  Matrix containing autofluorescence data points.

- samp:

  Sample identifier.

- gate.region:

  Dataframe containing region boundary information.

- gate.bound.density:

  Density (e.g., from MASS:kde2d) information

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

  Number of points to plot (speeds up plotting). Default is `1e5`.

## Value

Saves the plot as a JPEG file in the specified directory.
