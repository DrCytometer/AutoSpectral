# m \\\times\\ n Unmixed Biplot Grid

Produces a grid of scatter biplots comparing a set of m target
fluorophores (x-axis) against every other channel in the dataset
(y-axis), and saves the result as a multi-panel PDF. Each target
fluorophore defines one row of panels; the columns correspond to all
other selected channels.

Data are downsampled once to `max.points` events before any panel is
drawn, so all panels share the same cell population. Axis minima are
auto-scaled per channel using the 1st percentile of the data, floored at
the user-supplied `x.min` / `y.min`.

Panels can use either `geom_hex` (default, fastest) or the standard
`geom_scattermore` + `stat_density_2d` combination from
[`create.biplot()`](https://drcytometer.github.io/AutoSpectral/reference/create.biplot.md).

## Usage

``` r
unmixed.mxn.plot(
  unmixed.data,
  fluorophore,
  asp,
  channels = NULL,
  max.points = 50000,
  title = "mxn_biplot",
  biplot.size = 3,
  n.col = NULL,
  x.min = -1000,
  y.min = -1000,
  x.width.basis = -1000,
  y.width.basis = -1000,
  use.hex = TRUE,
  hex.bins = 64,
  color.palette = "viridis",
  output.dir = "."
)
```

## Arguments

- unmixed.data:

  Either a character string giving the path to an unmixed FCS file, or a
  numeric matrix / data.frame whose columns are fluorophore channels and
  whose rows are cells.

- fluorophore:

  Character vector of length \\m\\. Names of the target fluorophores to
  place on the x-axis (one row of panels per entry). Every element must
  be present in `colnames(unmixed.data)` (or the FCS columns).

- asp:

  The AutoSpectral parameter list, created by
  [`get.autospectral.param()`](https://drcytometer.github.io/AutoSpectral/reference/get.autospectral.param.md).

- channels:

  Optional character vector. Subset of channels to use as y-axis
  columns. When `NULL` (default) all columns except those in
  `fluorophore` are used. Must be present in the data if supplied.

- max.points:

  Integer. Total number of events to retain after a single random
  downsample applied before drawing any panel. Default `5e4`.

- title:

  Character string used as the PDF filename stem and as the overall plot
  title. Default `"mxn_biplot"`.

- biplot.size:

  Numeric. Edge length in inches of each individual panel. Default `3`.

- n.col:

  Integer. Number of panel columns per page. When `NULL` (default) this
  is set to the number of y-axis channels, giving one row of panels per
  target fluorophore on a single page.

- x.min:

  Numeric. Floor for the auto-scaled x-axis minimum (data units). The
  actual minimum used is `min(auto, x.min)`. Default `-1000`.

- y.min:

  Numeric. Floor for the auto-scaled y-axis minimum (data units).
  Default `-1000`.

- x.width.basis:

  Numeric. Width basis for the biexponential x-axis transform (passed to
  [`flowWorkspace::flowjo_biexp()`](https://rdrr.io/pkg/flowWorkspace/man/flowjo_biexp.html)).
  Default `-1000`.

- y.width.basis:

  Numeric. Width basis for the biexponential y-axis transform. Default
  `-1000`.

- use.hex:

  Logical. When `TRUE` (default), panels are rendered with
  [`ggplot2::geom_hex()`](https://ggplot2.tidyverse.org/reference/geom_hex.html)
  for speed. When `FALSE`, reverts to
  [`scattermore::geom_scattermore()`](https://rdrr.io/pkg/scattermore/man/geom_scattermore.html) +
  [`ggplot2::stat_density_2d()`](https://ggplot2.tidyverse.org/reference/geom_density_2d.html).

- hex.bins:

  Integer. Number of bins for `geom_hex` in each dimension. Ignored when
  `use.hex = FALSE`. Default `64`.

- color.palette:

  Character string. Viridis palette name (`"viridis"`, `"magma"`,
  `"inferno"`, `"plasma"`, `"cividis"`, `"rocket"`, `"mako"`, `"turbo"`)
  or `"rainbow"` to use the package default gradient. Default
  `"viridis"`.

- output.dir:

  Character string. Directory for the output PDF. Created automatically
  if absent. Default `"."`.

## Value

The combined
[`cowplot::plot_grid()`](https://wilkelab.org/cowplot/reference/plot_grid.html)
object is returned invisibly. The PDF is always written to `output.dir`.

## See also

[`unmixed.nxn.plot()`](https://drcytometer.github.io/AutoSpectral/reference/unmixed.nxn.plot.md),
[`create.biplot()`](https://drcytometer.github.io/AutoSpectral/reference/create.biplot.md)
