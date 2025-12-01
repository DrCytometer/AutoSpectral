# Spectral Heatmap

This function plots the spectral matrix as a heatmap using the specified
color palette, and saves it as a JPEG file.

## Usage

``` r
spectral.heatmap(
  spectra,
  title = NULL,
  plot.dir = NULL,
  legend.label = "Intensity",
  color.palette = "viridis",
  save = TRUE,
  show.legend = TRUE,
  plot.width = NULL,
  plot.height = NULL
)
```

## Arguments

- spectra:

  Matrix or dataframe containing spectral data format: fluorophores x
  detectors.

- title:

  Optional prefix for the plot filename.

- plot.dir:

  Optional output directory. Default is NULL, in which case the spectra
  figure folder will be used.

- legend.label:

  Character string that will appear on the heatmap legend.

- color.palette:

  Optional character string defining the viridis color palette to be
  used for the fluorophore traces. Default is `viridis`. Options are the
  viridis color options: `magma`, `inferno`, `plasma`, `viridis`,
  `cividis`, `rocket`, `mako` and `turbo`.

- save:

  Logical, if `TRUE`, saves a JPEG file to the `plot.dir`. Otherwise,
  the plot will simply be created in the Viewer.

- show.legend:

  Logical. If `TRUE`, figure legend will be included.

- plot.width:

  Width for the output plot. Default is `NULL`, in which case the width
  will be scaled automatically based on the number of detectors in
  `spectra`.

- plot.height:

  Height for the output plot. Default is `NULL`, in which case the
  height will be scaled automatically based on the number of rows in
  `spectra` plus a safety margin.

## Value

Saves the heatmap plot as a JPEG file in the specified directory.
