# Save a ggplot Directly to a Raster Device

Saves a ggplot object by opening the target device directly and printing
the plot, rather than going through
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)'s
own device- resolution and dimension-guessing logic. Produces the same
output as
`ggsave(filename, plot, device = device, width = width, height = height, dpi = dpi, ...)`
with `units = "in"`.

## Usage

``` r
.save.ggplot.fast(
  plot,
  filename,
  width,
  height,
  dpi = 300,
  device = ragg::agg_jpeg,
  ...
)
```

## Arguments

- plot:

  A ggplot object.

- filename:

  File path to save to.

- width:

  Numeric, plot width in inches.

- height:

  Numeric, plot height in inches.

- dpi:

  Numeric, resolution in dots per inch. Default `300`, matching
  [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)'s
  default.

- device:

  A ragg device-opening function, e.g.
  [`ragg::agg_jpeg`](https://ragg.r-lib.org/reference/agg_jpeg.html) or
  [`ragg::agg_png`](https://ragg.r-lib.org/reference/agg_png.html).
  Default
  [`ragg::agg_jpeg`](https://ragg.r-lib.org/reference/agg_jpeg.html).

- ...:

  Additional arguments passed to `device()`, e.g. `method = "fast"` or
  `quality = ` for
  [`ragg::agg_jpeg()`](https://ragg.r-lib.org/reference/agg_jpeg.html).

## Value

Invisibly, `filename`. Called for its side effect of writing the plot to
disk.
