# Create Heatmap Plot

This function plots a matrix as a heatmap and saves it as a JPEG file.

## Usage

``` r
create.heatmap(
  matrix,
  number.labels = FALSE,
  title = NULL,
  legend.label = "heatmap",
  triangular = FALSE,
  plot.dir = NULL,
  fixed.scale = FALSE,
  scale.min = NULL,
  scale.max = NULL,
  color.palette = "viridis",
  show.legend = TRUE,
  figure.width = 8,
  figure.height = 6,
  save = TRUE
)
```

## Arguments

- matrix:

  Matrix or dataframe containing spectral data.

- number.labels:

  Logical indicating whether to add number labels to the heatmap.
  Default is `FALSE`.

- title:

  Optional prefix for the plot filename. Default is `NULL`, in which
  case the file will just be called `heatmap.jpg`

- legend.label:

  Character string that will appear on the heatmap legend. Default is
  `heatmap`

- triangular:

  Logical. Plot the lower triangle of the matrix only, diagonal
  included. Default is `FALSE`.

- plot.dir:

  Optional output directory. Default is `NULL`, in which case the
  working directory will be used.

- fixed.scale:

  Logical, determines whether to use an externally supplied fixed scale
  (min and max) for the heatmap color scale. Useful for putting multiple
  plots on the same scale. Default is `FALSE`

- scale.min:

  Optional numeric. Minimum for the fixed color scale. Default is
  `NULL`, for no fixed scaling.

- scale.max:

  Optional numeric. Maximum for the fixed color scale. Default is
  `NULL`, for no fixed scaling.

- color.palette:

  Optional character string defining the viridis color palette to be
  used for the fluorophore traces. Default is `viridis`. Options are the
  viridis color options: `magma`, `inferno`, `plasma`, `viridis`,
  `cividis`, `rocket`, `mako` and `turbo`.

- show.legend:

  Logical. If `TRUE`, figure legend will be included.

- figure.width:

  Numeric. Width of the heatmap figure. Default is `8`.

- figure.height:

  Numeric. Height of the heatmap figure. Default is `6`.

- save:

  Logical, if `TRUE`, saves a JPEG file to the `output.dir`. Otherwise,
  the plot will simply be created in the Viewer.

## Value

Saves the heatmap plot as a JPEG file and the SSM data as a CSV file in
the specified directory.
