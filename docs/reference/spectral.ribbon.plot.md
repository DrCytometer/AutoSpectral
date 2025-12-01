# Spectral Ribbon Plot

This function generates spectral ribbon plots for positive and negative
expression data. The function gets called internally in
`clean.controls`. To use the function directly, pass data to `data.list`
in the form of a named list of matrices or data.frames.

## Usage

``` r
spectral.ribbon.plot(
  pos.expr.data = NULL,
  neg.expr.data = NULL,
  removed.data = NULL,
  spectral.channel,
  asp,
  fluor.name = NULL,
  title = NULL,
  figure.dir = NULL,
  factor.names = NULL,
  save = TRUE,
  color.palette = "rainbow",
  data.list = NULL,
  af = FALSE,
  plot.width = 15,
  plot.height = 10
)
```

## Arguments

- pos.expr.data:

  Internal argument for `clean.controls`. A matrix containing the
  positive expression data. Default is `NULL`.

- neg.expr.data:

  Internal argument for `clean.controls`. A matrix containing the
  negative expression data. Default is `NULL`.

- removed.data:

  Internal argument for `clean.controls`. A matrix containing the
  removed data, if applicable. Default is `NULL`. If omitted, only two
  groups (facets) are plotted.

- spectral.channel:

  A character vector specifying the spectral channels. Recommended: use
  `colnames(spectra)` or `flow.control$spectral.channel`.

- asp:

  The AutoSpectral parameter list. Prepare using get.autospectral.param.

- fluor.name:

  An optional character string specifying the fluorophore name for plot
  titles and filename. Default is `NULL`.

- title:

  An optional character string to prefix the plot file name. Default is
  `NULL`.

- figure.dir:

  Output folder where the figures will be created. Default is `NULL`,
  enabling automatic selection inside AutoSpectral. For user-supplied
  data, `asp$figure.spectral.ribbon.dir` will be used if `NULL`.

- factor.names:

  Optional titles for the facets on the plot. Default is `NULL`,
  enabling automatic selection inside AutoSpectral. If `data.list` is
  named, `factor.names` will be pulled from that.

- save:

  Logical, default is `TRUE`. If `TRUE`, the plot is saved to
  `figure.dir`. Otherwise, it is returned to the viewer only.

- color.palette:

  Optional character string defining the color palette to be used.
  Default is `rainbow`, mimicking a FlowJo scheme. Other choices are the
  viridis color options: `magma`, `inferno`, `plasma`, `viridis`,
  `cividis`, `rocket`, `mako` and `turbo`.

- data.list:

  Provide data here. A (named) list of matrices or dataframes for
  plotting. These should probably be flow expression data. Data provided
  will be subsetted to the columns in `spectral.channel`.

- af:

  Internal argument for `clean.controls`. A logical value indicating
  whether autofluorescence removal is being performed. Default is
  `FALSE`.

- plot.width:

  Width of the saved plot. Default is `15`.

- plot.height:

  Height of the saved plot. Default is `10`.

## Value

None. The function saves the generated spectral ribbon plot to a file.
