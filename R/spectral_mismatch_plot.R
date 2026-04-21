# spectral_mismatch_plot.R

#' @title Spectral Mismatch Plot
#'
#' @description
#' Produces a two-panel diagnostic figure comparing two fluorophore spectra. The
#' upper panel overlays the normalised spectral traces for a reference and a test
#' spectrum across all detectors. The lower panel shows the per-detector spectral
#' difference (reference minus test) as a bar chart, making it easy to identify
#' channels where the spectra diverge.
#'
#' The combined figure is saved as a JPEG to `plot.dir` and returned invisibly as
#' a `ggplot` / cowplot object so it can be embedded in downstream reports or
#' further modified.
#'
#' @param ref.spectrum A named numeric vector (or single-row matrix) of
#' normalised spectral intensities for the reference spectrum. Names must
#' correspond to detector channel names.
#' @param test.spectrum A named numeric vector (or single-row matrix) of
#' normalised spectral intensities for the test spectrum. Must have the same
#' names as `ref.spectrum`.
#' @param ref.label Character string used as the legend label for the reference
#' spectrum. If `NULL` (default), the `rownames()` of `ref.spectrum` is used.
#' @param test.label Character string used as the legend label for the test
#' spectrum. If `NULL` (default), the `rownames()` of `test.spectrum` is used.
#' @param title Character string used as both the plot title and the JPEG
#' filename stem. Default is `"Spectral comparison"`.
#' @param bar.color Fill color for the spectral-difference bar chart.
#' Default is `"lightblue"`.
#' @param fluorophore.colors Character vector of length 2 giving the line colors
#' for the reference and test spectra respectively. Overrides `color.palette`
#' when non-`NULL`. Default is `c("black", "red")`.
#' @param plot.dir Directory to save the JPEG file. If `NULL` (default), the
#' current working directory is used. Created automatically if absent.
#' @param line.width Numeric. Line width for the spectral traces. Default `1`.
#' @param point.size Numeric. Point size for the spectral traces. Default `1`.
#' @param color.palette Optional viridis palette name (e.g., `"viridis"`,
#' `"plasma"`) used when `fluorophore.colors` is `NULL`. Ignored when
#' `fluorophore.colors` is set. Default is `NULL`.
#' @param show.legend Logical. Whether to display the color legend on the
#' upper panel. Default is `TRUE`.
#' @param plot.width Numeric. Width of the saved figure in inches. If `NULL`
#' (default), an automatic width is derived from the number of detectors
#' (approximately 12 inches per 64 channels, minimum 3 inches).
#' @param plot.height Numeric. Height of the saved figure in inches. If `NULL`
#' (default), an automatic height is derived from the number of spectra rows.
#'
#' @return The combined cowplot object is returned invisibly. The figure is
#' always saved to disk as `<plot.dir>/<title>.jpg`.
#'
#' @importFrom ggplot2 ggplot aes geom_path geom_point geom_bar labs
#' @importFrom ggplot2 theme_minimal theme element_text scale_color_manual
#' @importFrom ggplot2 scale_color_viridis_d
#' @importFrom cowplot plot_grid
#' @importFrom ragg agg_jpeg
#'
#' @export

spectral.mismatch.plot <- function(
    ref.spectrum,
    test.spectrum,
    ref.label          = NULL,
    test.label         = NULL,
    title              = "Spectral comparison",
    bar.color          = "lightblue",
    fluorophore.colors = c( "black", "red" ),
    plot.dir           = NULL,
    line.width         = 1,
    point.size         = 1,
    color.palette      = NULL,
    show.legend        = TRUE,
    plot.width         = NULL,
    plot.height        = NULL
  ) {

  # --- input validation ---
  if ( !is.numeric( ref.spectrum ) || !is.numeric( test.spectrum ) ) {
    stop( "`ref.spectrum` and `test.spectrum` must be numeric.", call. = FALSE )
  }

  ref.spectrum  <- as.matrix( ref.spectrum )
  test.spectrum <- as.matrix( test.spectrum )

  if ( !identical( colnames( ref.spectrum ), colnames( test.spectrum ) ) ) {
    stop(
      "`ref.spectrum` and `test.spectrum` must have identical column names.",
      call. = FALSE
    )
  }

  # --- assemble spectra data frame ---
  spectra           <- rbind( ref.spectrum, test.spectrum )
  spectra           <- data.frame( spectra, check.names = FALSE )
  spectra$Fluorophore <- rownames( spectra )

  if ( !is.null( ref.label ) )  spectra$Fluorophore[ 1 ] <- ref.label
  if ( !is.null( test.label ) ) spectra$Fluorophore[ 2 ] <- test.label

  pivot.cols <- setdiff( colnames( spectra ), "Fluorophore" )

  # --- output directory ---
  if ( is.null( plot.dir ) ) plot.dir <- getwd()
  if ( !dir.exists( plot.dir ) ) dir.create( plot.dir, recursive = TRUE )

  # --- automatic figure dimensions ---
  if ( is.null( plot.width ) )
    plot.width <- max( ( ( length( pivot.cols ) ) / 64 * 12 ), 3 )

  if ( is.null( plot.height ) )
    plot.height <- 5 + round( nrow( spectra ) / 8, 0 ) * 2

  # --- pivot to long format (base R, no tidyr dependency) ---
  spectra.long <- data.frame(
    Fluorophore      = rep( spectra$Fluorophore, times = length( pivot.cols ) ),
    Detector         = rep( pivot.cols, each = nrow( spectra ) ),
    Intensity        = as.vector( as.matrix( spectra[ , pivot.cols ] ) ),
    stringsAsFactors = FALSE
  )

  spectra.long$Detector <- factor(
    spectra.long$Detector,
    levels  = unique( spectra.long$Detector ),
    ordered = TRUE
  )

  # --- upper panel: spectral traces ---
  spectra.plot <- ggplot2::ggplot(
    spectra.long,
    ggplot2::aes(
      x     = Detector,
      y     = Intensity,
      group = Fluorophore,
      color = Fluorophore
    )
  ) +
    ggplot2::geom_path( linewidth = line.width ) +
    ggplot2::geom_point( size = point.size ) +
    ggplot2::labs(
      title = title,
      x     = "Detector",
      y     = "Normalized Intensity"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x    = ggplot2::element_text( angle = 45, hjust = 1 ),
      legend.position = "bottom"
    )

  # apply color palette
  if ( !is.null( fluorophore.colors ) ) {
    spectra.plot <- spectra.plot +
      ggplot2::scale_color_manual( values = fluorophore.colors )
  } else if ( !is.null( color.palette ) ) {
    spectra.plot <- spectra.plot +
      ggplot2::scale_color_viridis_d( option = color.palette )
  }

  if ( !show.legend ) {
    spectra.plot <- spectra.plot +
      ggplot2::theme( legend.position = "none" )
  }

  # --- lower panel: spectral difference ---
  # difference is always ref - test (row 1 - row 2 of original inputs)
  delta.vals <- as.numeric( ref.spectrum[ 1, ] ) - as.numeric( test.spectrum[ 1, ] )
  delta.long <- data.frame(
    Detector         = pivot.cols,
    Intensity        = delta.vals,
    stringsAsFactors = FALSE
  )
  delta.long$Detector <- factor(
    delta.long$Detector,
    levels  = pivot.cols,
    ordered = TRUE
  )

  delta.plot <- ggplot2::ggplot(
    delta.long,
    ggplot2::aes( x = Detector, y = Intensity )
  ) +
    ggplot2::geom_bar(
      stat  = "identity",
      fill  = bar.color,
      color = "black",
      alpha = 0.7
    ) +
    ggplot2::labs( x = "Detector", y = "Spectral Difference (ref \u2212 test)" ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text( angle = 45, hjust = 1 )
    )

  # --- combine panels ---
  combined.plot <- cowplot::plot_grid(
    spectra.plot,
    delta.plot,
    ncol        = 1,
    align       = "v",
    axis        = "lr",
    rel_heights = c( 2, 1 )
  )

  # --- save ---
  ggplot2::ggsave(
    file.path( plot.dir, sprintf( "%s.jpg", title ) ),
    combined.plot,
    device    = ragg::agg_jpeg,
    width     = plot.width,
    height    = plot.height,
    limitsize = FALSE
  )

  invisible( combined.plot )
}
