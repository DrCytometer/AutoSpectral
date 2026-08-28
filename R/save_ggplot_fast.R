# save_ggplot_fast.r

#' @title Save a ggplot Directly to a Raster Device
#'
#' @description
#' Saves a ggplot object by opening the target device directly and printing
#' the plot, rather than going through `ggplot2::ggsave()`'s own device-
#' resolution and dimension-guessing logic. Produces the same output as
#' `ggsave(filename, plot, device = device, width = width, height = height,
#' dpi = dpi, ...)` with `units = "in"`.
#'
#' @param plot A ggplot object.
#' @param filename File path to save to.
#' @param width Numeric, plot width in inches.
#' @param height Numeric, plot height in inches.
#' @param dpi Numeric, resolution in dots per inch. Default `300`, matching
#' `ggplot2::ggsave()`'s default.
#' @param device A ragg device-opening function, e.g. `ragg::agg_jpeg` or
#' `ragg::agg_png`. Default `ragg::agg_jpeg`.
#' @param ... Additional arguments passed to `device()`, e.g. `method =
#' "fast"` or `quality = ` for `ragg::agg_jpeg()`.
#'
#' @return Invisibly, `filename`. Called for its side effect of writing the
#' plot to disk.

.save.ggplot.fast <- function(
    plot, filename, width, height, dpi = 300, device = ragg::agg_jpeg, ...
) {

  device( filename = filename, width = width, height = height,
          units = "in", res = dpi, ... )
  on.exit( grDevices::dev.off(), add = TRUE )

  print( plot )

  invisible( filename )
}
