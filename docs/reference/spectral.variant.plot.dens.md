# Spectral Variant Density Plot

Creates a static visualization showing multiple spectral variants
overlaid against a median spectral signature to illustrate consistency
and variance. Each variant is plotted as a semi-transparent line, and
the median spectrum is drawn on top as a solid reference trace.

## Usage

``` r
spectral.variant.plot.dens(
  spectra.variants,
  median.spectrum,
  title = "Spectral_variants_density",
  save = FALSE,
  plot.width = NULL,
  plot.height = 5,
  plot.dir = "./figure_spectral_variants",
  variant.color = "red",
  variant.alpha = 0.05,
  variant.linewidth = 1.5,
  median.line.color = "black",
  median.linewidth = 1
)
```

## Arguments

- spectra.variants:

  Matrix of spectral variants; rows are variants, columns are detectors.
  Column names are used as x-axis labels.

- median.spectrum:

  Numeric vector of the median signal intensity across detectors. Must
  have the same length as `ncol(spectra.variants)`.

- title:

  Plot title and output filename prefix. Default is
  `"Spectral_variants_density"`.

- save:

  Logical. If `TRUE`, saves the plot to `plot.dir`. Default is `FALSE`.

- plot.width:

  Width of the saved plot in inches. Default is `NULL`, which
  auto-scales based on the number of detectors.

- plot.height:

  Height of the saved plot in inches. Default is `5`.

- plot.dir:

  Output directory for the saved plot. Default is
  `"./figure_spectral_variants"`.

- variant.color:

  Line color for the individual variant traces. Default is `"red"`.

- variant.alpha:

  Alpha transparency for the variant lines. Default is `0.05`.

- variant.linewidth:

  Thickness of the variant lines. Default is `1.5`.

- median.line.color:

  Color of the median reference line. Default is `"black"`.

- median.linewidth:

  Thickness of the median reference line. Default is `1`.

## Value

The `ggplot` object representing the variant density plot. When
`save = TRUE`, a JPEG is also written to `plot.dir`.
