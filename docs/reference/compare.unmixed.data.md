# Compare Two or More Sets of Already-Unmixed Data

Lays out biplots for two or more pre-unmixed data sets (matrices or data
frames of the same panel) as a grid, one row per data set, so that the
same channel pairs can be compared side-by-side across conditions (for
example, cells vs. beads, or before vs. after a correction step).

Channels are paired up sequentially rather than plotted N x N: for
`channels = c("APC", "FITC", "PE", "PerCP", "BV421")` the pairs are APC
x FITC, PE x PerCP, and BV421 x APC (the first channel is repeated to
complete the last pair when an odd number of channels is requested).
Each row uses
[`create.biplot()`](https://drcytometer.github.io/AutoSpectral/reference/create.biplot.md)
for every pair, and the panels are stitched together with
[`cowplot::plot_grid()`](https://wilkelab.org/cowplot/reference/plot_grid.html).

If `gate.boundary` is supplied and the `scatter.param` columns are
present in every data set,
[`apply.gate()`](https://drcytometer.github.io/AutoSpectral/reference/apply.gate.md)
is used to gate every data set to that boundary before plotting. If any
data set is missing a scatter column, the gate is not applied to any of
them (a warning names which data set(s) were missing it), so that every
row in the comparison reflects the same population.

If `data.list` is a named list, the names are used as row labels on the
right of the grid. If unnamed, and `data.list` was constructed inline as
`list(...)` at the call site, the deparsed expression for each element
(i.e. the variable name passed in) is used instead.

## Usage

``` r
compare.unmixed.data(
  data.list,
  channels,
  asp,
  gate.boundary = NULL,
  scatter.param = asp$default.scatter.parameter,
  variants = NULL,
  spread.kappa = 2,
  x.min = -5000,
  x.max = asp$expr.data.max,
  y.min = -5000,
  y.max = asp$expr.data.max,
  x.width.basis = -1000,
  y.width.basis = -1000,
  max.points = 5e+05,
  color.palette = "rainbow",
  panel.width = 3,
  row.height = 3,
  label.width = 1,
  label.size = 12,
  label.angle = -90,
  save = TRUE,
  title = "unmix_comparison",
  output.dir = NULL
)
```

## Arguments

- data.list:

  A list of two or more matrices or data frames of unmixed flow
  cytometry data, one per data set to compare. Every element must have
  named columns. Column names do not need to match exactly across data
  sets – see `channels` below.

- channels:

  Character vector of channel names to plot, in the order they should be
  paired (see Description). Must resolve to at least 2 channels present
  in every element of `data.list` after the intersect described below.

- asp:

  The AutoSpectral parameter list, prepared using
  [`get.autospectral.param()`](https://drcytometer.github.io/AutoSpectral/reference/get.autospectral.param.md).
  Passed through to
  [`create.biplot()`](https://drcytometer.github.io/AutoSpectral/reference/create.biplot.md)
  for every panel, and supplies the default for `scatter.param`.

- gate.boundary:

  Optional gate boundary, as returned by
  [`define.gate.landmarks()`](https://drcytometer.github.io/AutoSpectral/reference/define.gate.landmarks.md),
  [`define.gate.density()`](https://drcytometer.github.io/AutoSpectral/reference/define.gate.density.md),
  or
  [`do.gate()`](https://drcytometer.github.io/AutoSpectral/reference/do.gate.md)
  – a list containing at least numeric `x` and `y` components describing
  the polygon vertices. When supplied, every data set is gated with
  [`apply.gate()`](https://drcytometer.github.io/AutoSpectral/reference/apply.gate.md)
  before plotting, provided all data sets carry `scatter.param`. Default
  `NULL` (no gating).

- scatter.param:

  Character vector of length 2 giving the names of the two scatter
  columns to gate on, passed to
  [`apply.gate()`](https://drcytometer.github.io/AutoSpectral/reference/apply.gate.md).
  Only used when `gate.boundary` is supplied. Default
  `asp$default.scatter.parameter`.

- variants:

  The variant list returned by
  [`get.spectral.variants()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectral.variants.md),
  passed through to every
  [`create.biplot()`](https://drcytometer.github.io/AutoSpectral/reference/create.biplot.md)
  call for reference curves. Default `NULL` (no reference curves).

- spread.kappa:

  Numeric, passed through to
  [`create.biplot()`](https://drcytometer.github.io/AutoSpectral/reference/create.biplot.md).
  Default `2`.

- x.min, y.min:

  Numeric, axis minima (data units) passed through to
  [`create.biplot()`](https://drcytometer.github.io/AutoSpectral/reference/create.biplot.md).
  Default `-5000` for both.

- x.max, y.max:

  Numeric, axis maxima (data units) passed through to
  [`create.biplot()`](https://drcytometer.github.io/AutoSpectral/reference/create.biplot.md).
  Default `asp$expr.data.max` for both.

- x.width.basis, y.width.basis:

  Numeric, biexponential width bases passed through to
  [`create.biplot()`](https://drcytometer.github.io/AutoSpectral/reference/create.biplot.md).
  Default `-1000` for both.

- max.points:

  Numeric, per-panel point cap passed through to
  [`create.biplot()`](https://drcytometer.github.io/AutoSpectral/reference/create.biplot.md).
  Default `5e5`.

- color.palette:

  Character, passed through to
  [`create.biplot()`](https://drcytometer.github.io/AutoSpectral/reference/create.biplot.md).
  Default `"rainbow"`.

- panel.width, row.height:

  Numeric, width and height (inches) of a single biplot panel, used to
  size the saved figure. Defaults `3` and `3`.

- label.width:

  Numeric, width (inches) reserved for the row-label column on the right
  of the figure. Default `1`.

- label.size:

  Numeric, font size (points) for the row labels. Default `12`.

- label.angle:

  Numeric, rotation angle (degrees) for the row labels. Default `-90`
  (reads bottom-to-top, matching a facet_grid right strip).

- save:

  Logical, if `TRUE` (default), saves a JPEG file to `output.dir`.
  Otherwise the combined plot is only printed.

- title:

  Character, used as the JPEG filename stem. Default
  `"unmix_comparison"`.

- output.dir:

  Optional output directory. Default `NULL`, in which case the current
  working directory is used.

## Value

Invisibly, the combined `cowplot` object.

## See also

- [`create.biplot()`](https://drcytometer.github.io/AutoSpectral/reference/create.biplot.md)

- [`apply.gate()`](https://drcytometer.github.io/AutoSpectral/reference/apply.gate.md)

- [`compare.unmix()`](https://drcytometer.github.io/AutoSpectral/reference/compare.unmix.md)
