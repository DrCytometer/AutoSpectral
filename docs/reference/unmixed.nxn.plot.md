# n \\\times\\ n Unmixed Biplot Triangle

Produces a lower-triangular grid of scatter biplots for every pair of
channels in an unmixed dataset, with 1-D histograms on the diagonal, and
saves the result as a single large-format PDF.

For \\n\\ channels the grid is \\n \times n\\ panels. The upper triangle
and all cells above the diagonal are left blank, giving \\n(n-1)/2\\
scatter panels and \\n\\ histogram panels. Panel \\(i, j)\\ (row \\i\\,
column \\j\\) shows column \\j\\ on the x-axis and row \\i\\ on the
y-axis, consistent with the convention used by
[`unmixed.mxn.plot()`](https://drcytometer.github.io/AutoSpectral/reference/unmixed.mxn.plot.md):
the named fluorophore always appears on x.

Data are downsampled once to `max.points` events before any panel is
drawn so all panels share an identical cell population. Axis minima are
auto-scaled per channel from the 1st percentile of the data, floored at
the user-supplied `x.min` / `y.min`.

## Usage

``` r
unmixed.nxn.plot(
  unmixed.data,
  asp,
  channels = NULL,
  max.points = 50000,
  title = "nxn_biplot",
  biplot.size = 3,
  x.min = -1000,
  y.min = -1000,
  x.width.basis = -1000,
  y.width.basis = -1000,
  use.hex = TRUE,
  hex.bins = 64,
  color.palette = "viridis",
  hist.fill = "steelblue",
  output.dir = "."
)
```

## Arguments

- unmixed.data:

  Either a character string giving the path to an unmixed FCS file, or a
  numeric matrix / data.frame whose columns are fluorophore channels and
  whose rows are cells.

- asp:

  The AutoSpectral parameter list, created by
  [`get.autospectral.param()`](https://drcytometer.github.io/AutoSpectral/reference/get.autospectral.param.md).

- channels:

  Optional character vector. Subset of columns to include. When `NULL`
  (default) all columns are used. Large panels sets (\\n \> 20\\) will
  produce very large PDFs; consider subsetting.

- max.points:

  Integer. Total number of events retained after a single random
  downsample applied before any panel is drawn. Default `5e4`.

- title:

  Character string used as the PDF filename stem. Default
  `"nxn_biplot"`.

- biplot.size:

  Numeric. Edge length in inches of each individual panel (both scatter
  and histogram panels). The saved PDF will be `n * biplot.size`
  \\\times\\ `n * biplot.size` inches. Default `3`.

- x.min:

  Numeric. Floor for the auto-scaled x-axis minimum (data units). The
  actual minimum used is `min(auto, x.min)`. Default `-1000`.

- y.min:

  Numeric. Floor for the auto-scaled y-axis minimum (data units).
  Default `-1000`.

- x.width.basis:

  Numeric. Width basis for the biexponential x-axis transform. Default
  `-1000`.

- y.width.basis:

  Numeric. Width basis for the biexponential y-axis transform. Default
  `-1000`.

- use.hex:

  Logical. When `TRUE` (default), scatter panels are rendered with
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
  or `"rainbow"` to use the package default gradient. Applies to scatter
  panels. Default `"viridis"`.

- hist.fill:

  Character string. Fill color for the diagonal histogram bars. Default
  `"steelblue"`.

- output.dir:

  Character string. Directory for the output PDF. Created automatically
  if absent. Default `"."`.

## Value

The combined
[`cowplot::plot_grid()`](https://wilkelab.org/cowplot/reference/plot_grid.html)
object is returned invisibly. The PDF is always written to `output.dir`.

## See also

[`unmixed.mxn.plot()`](https://drcytometer.github.io/AutoSpectral/reference/unmixed.mxn.plot.md),
[`create.biplot()`](https://drcytometer.github.io/AutoSpectral/reference/create.biplot.md)
