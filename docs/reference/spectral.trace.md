# Plot Fluorophore Spectra Traces

This function plots the fluorophore spectra as normalized traces,
optionally splitting by excitation lasers, and saves the plots as JPEG
files.

## Usage

``` r
spectral.trace(
  spectral.matrix,
  asp,
  title = "Fluorophore_Spectra",
  plot.dir = NULL,
  split.lasers = TRUE,
  figure.spectra.line.size = 1,
  figure.spectra.point.size = 1,
  color.palette = NULL,
  show.legend = TRUE,
  plot.width = NULL,
  plot.height = NULL,
  save = TRUE
)
```

## Arguments

- spectral.matrix:

  Matrix or dataframe containing spectral data. This should be in format
  fluorophores x detectors. Row names will be used as the fluorophore
  names. Column names will be used as the detectors (channels).

- asp:

  The AutoSpectral parameter list defined using
  `get.autospectral.param`.

- title:

  Title for the plot. Default is `Fluorophore_Spectra`

- plot.dir:

  Directory to save the plot files. Default is `NULL`, in which case the
  current working directory will be used.

- split.lasers:

  Logical indicating whether to create a second plot split by excitation
  lasers. Default is `TRUE`.

- figure.spectra.line.size:

  Numeric. Defines the width of the trace lines on the plot. Default is
  `1`.

- figure.spectra.point.size:

  Numeric. Defines the size of the points on the plot (one per
  detector). Default is `1`.

- color.palette:

  Optional character string defining the viridis color palette to be
  used for the fluorophore traces. Default is `NULL`, in which case
  default R Brewer colors will be assigned automatically. Options are
  the viridis color options: `magma`, `inferno`, `plasma`, `viridis`,
  `cividis`, `rocket`, `mako` and `turbo`.

- show.legend:

  Logical. If `TRUE`, figure legend will be included.

- plot.width:

  Optional numeric to manually set the plot width. Default is `NULL`.

- plot.height:

  Optional numeric to manually set the plot width. Default is `NULL`.

- save:

  Logical, if `TRUE`, saves a JPEG file to the `plot.dir`. Otherwise,
  the plot will simply be created in the Viewer.

## Value

Saves the plot(s) as JPEG files in the specified directory.
