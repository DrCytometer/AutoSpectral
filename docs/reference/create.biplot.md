# Create Biplot

Create Biplot

## Usage

``` r
create.biplot(
  plot.data,
  x.dim,
  y.dim,
  asp,
  x.lab = NULL,
  y.lab = NULL,
  x.min = -5000,
  x.max = asp$expr.data.max,
  y.min = -5000,
  y.max = asp$expr.data.max,
  x.width.basis = -1000,
  y.width.basis = -1000,
  max.points = 5e+06,
  color.palette = "rainbow",
  save = TRUE,
  title = NULL,
  output.dir = NULL,
  width = 5,
  height = 5
)
```

## Arguments

- plot.data:

  A matrix or dataframe containing the flow cytometry data to be
  plotted. Column names should match the dimensions specified by `x.dim`
  and `y.dim`.

- x.dim:

  String specifying the column of `plot.data` for the x-axis of the
  plot.

- y.dim:

  String specifying the column of `plot.data` for the y-axis of the
  plot.

- asp:

  The AutoSpectral parameter list.

- x.lab:

  An optional label for the x-axis. If none is given (default `NULL`),
  the column name specified by `x.dim` will be used.

- y.lab:

  An optional label for the y-axis. If none is given (default `NULL`),
  the column name specified by `y.dim` will be used.

- x.min:

  Minimum value for the x-axis. Default is `-5000`.

- x.max:

  Maximum value for the x-axis. Default is the value specified by
  asp\$expr.data.max, which will be the maximum for the cytometer.

- y.min:

  Minimum value for the y-axis. Default is `-5000`.

- y.max:

  Maximum value for the y-axis. Default is the value specified by
  asp\$expr.data.max, which will be the maximum for the cytometer.

- x.width.basis:

  Width basis for the biexponential transform for the x-axis. Default is
  `-1000`.

- y.width.basis:

  Width basis for the biexponential transform for the x-axis. Default is
  `-1000`.

- max.points:

  Number of points to plot (speeds up plotting). Default is `5e6`.

- color.palette:

  Optional character string defining the viridis color palette to be
  used for the fluorophore traces. Default is `rainbow`, which will be
  similar to FlowJo or SpectroFlo. Other pptions are the viridis color
  options: `magma`, `inferno`, `plasma`, `viridis`, `cividis`, `rocket`,
  `mako` and `turbo`.

- save:

  Logical, if `TRUE`, saves a JPEG file to the `output.dir`. Otherwise,
  the plot will simply be created in the Viewer.

- title:

  Optional title for the plot filename. If `NULL`, defaults to `x.lab`
  vs. `y.lab`.

- output.dir:

  Optional output directory. Default is NULL, in which case the current
  working directory will be used.

- width:

  Numeric, width of the saved plot. Default is `5`.

- height:

  Numeric, height of the saved plot. Default is `5`.

## Value

Creates a biplot in the Viewer and optionally saves it as a JPEG file.
