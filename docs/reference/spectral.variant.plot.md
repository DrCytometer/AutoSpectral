# Spectral Variant Plot

This function plots a spectral ribbon showing the optimized fluorophore
signature with the range of variation from the spectral variants.

## Usage

``` r
spectral.variant.plot(
  spectra.variants,
  median.spectrum,
  title = "Spectral_variants",
  save = FALSE,
  plot.width = NULL,
  plot.height = 5,
  plot.dir = "./figure_spectral_variants",
  variant.fill.color = "red",
  variant.fill.alpha = 0.7,
  median.line.color = "black",
  median.linewidth = 1
)
```

## Arguments

- spectra.variants:

  A matrix of variations in the normalized spectra signature for a
  fluorophore, rows being the variants, columns being the detectors.

- median.spectrum:

  The optimized reference spectrum for the fluorophore, usually obtained
  from spectra prepared using AutoSpectral via `get.fluorophore.spectra`
  on cleaned control data.

- title:

  Optional title to pass to the plot and output filename. Default is
  `Spectral_variants`.

- save:

  Logical, controls whether the plot is displayed to the Viewer
  (`FALSE`) or also saved to disk (`TRUE`).

- plot.width:

  Optional numeric to control saved plot width. Default is `NULL`, in
  which case the plot width will be calculated based on the number of
  detectors in the data.

- plot.height:

  Optional numeric to control saved plot width. Default is `5`.

- plot.dir:

  Directory where the files will be saved if `save = TRUE`. Default is
  `./figure_spectral_variants`.

- variant.fill.color:

  Color for the shaded region indicating the range of variation in the
  spectra. Feeds to `fill` in `geom_ribbon`. Default is "red".

- variant.fill.alpha:

  Transparency (alpha) for the color in `variant.fill.color`. How
  intense the color of the variant spectra will be. Default is `0.7`

- median.line.color:

  Color for the line representing the median or optimized single
  spectrum. Default is "black".

- median.linewidth:

  Width of the line for the single optimized spectrum. Default is `1`.

## Value

Plots are displayed to the Viewer. Files are saved if saving enabled.
