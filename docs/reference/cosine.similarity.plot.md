# Cosine Similarity Plot

This function plots the matrix of cosine similarities (AKA "Similarity
Matrix") as a heatmap and saves it as a JPEG file. It also calculates
and displays the mixing matrix condition number (AKA "Complexity Index")
of the matrix.

## Usage

``` r
cosine.similarity.plot(
  spectra,
  filename = "autospectral_similarity_matrix",
  title = NULL,
  output.dir = "figure_similarity_heatmap",
  figure.width = 8,
  figure.height = 6,
  color.palette = "viridis",
  show.legend = TRUE,
  file.type = "jpg",
  save = TRUE
)
```

## Arguments

- spectra:

  Data frame or matrix containing spectral data.

- filename:

  Character string for the output file. Default is
  `autospectral_similarity_matrix`.

- title:

  Optional prefix for the plot filename. Default is `NULL`

- output.dir:

  File path where the plot will be created. Default is
  `figure_similarity_heatmap`. The directory will be created if it does
  not already exist.

- figure.width:

  Numeric. Width of output plot. Default is `8`.

- figure.height:

  Numeric. Height of output plot. Default is `6`.

- color.palette:

  Optional character string defining the viridis color palette to be
  used for the fluorophore traces. Default is `viridis`. Options are the
  viridis color options: `magma`, `inferno`, `plasma`, `viridis`,
  `cividis`, `rocket`, `mako` and `turbo`.

- show.legend:

  Logical. If `TRUE`, figure legend will be included.

- file.type:

  Character string specifying the output file format. One of `"jpg"`
  (default), `"jpeg"`, `"tiff"`, `"png"`, or `"pdf"`.

- save:

  Logical, if `TRUE`, saves a file to the `output.dir`. Otherwise, the
  plot will simply be created in the Viewer.

## Value

Saves the heatmap plot as an image file in the similarity directory.
