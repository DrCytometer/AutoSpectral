# Gate Definition Plot

This function plots the gate during the definition step, including
boundaries and regions.

## Usage

``` r
gate.define.plot(
  samp,
  gate.data,
  gate.marker,
  gate.bound,
  gate.region,
  gate.population,
  scatter.and.channel.label,
  asp,
  color.palette = "plasma",
  max.points = 1e+05,
  gate.color = "black"
)
```

## Arguments

- samp:

  Sample identifier.

- gate.data:

  Matrix containing gate data points.

- gate.marker:

  Vector containing gate marker names.

- gate.bound:

  List containing gate boundary information.

- gate.region:

  List containing gate region information.

- gate.population:

  List containing gate population information.

- scatter.and.channel.label:

  Named vector mapping scatter and channel labels.

- asp:

  The AutoSpectral parameter list.

- color.palette:

  Optional character string defining the viridis color palette to be
  used for the fluorophore traces. Default is `plasma`. Use `rainbow` to
  be similar to FlowJo or SpectroFlo. Other options are the viridis
  color options: `magma`, `inferno`, `plasma`, `viridis`, `cividis`,
  `rocket`, `mako` and `turbo`.

- max.points:

  Number of points to plot. Default is `1e5`.

- gate.color:

  Color to plot the gate boundary line, default is `black`.

## Value

Saves the plot as a JPEG file in the specified directory.
