# spectral_variant_density_plot.R

utils::globalVariables( c(
  "detector_idx", "intensity", "id"
) )

#' @title Spectral Variant Density Plot
#'
#' @description
#' Creates a static visualization showing multiple spectral variants overlaid
#' against a median spectral signature to illustrate consistency and variance.
#' Each variant is plotted as a semi-transparent line, and the median spectrum
#' is drawn on top as a solid reference trace.
#'
#' @importFrom ggplot2 ggplot aes geom_line scale_x_continuous labs
#' @importFrom ggplot2 theme_minimal theme element_text ggsave coord_cartesian
#' @importFrom ragg agg_jpeg
#'
#' @param spectra.variants Matrix of spectral variants; rows are variants,
#'   columns are detectors. Column names are used as x-axis labels.
#' @param median.spectrum Numeric vector of the median signal intensity across
#'   detectors. Must have the same length as `ncol(spectra.variants)`.
#' @param title Plot title and output filename prefix. Default is
#'   `"Spectral_variants_density"`.
#' @param save Logical. If `TRUE`, saves the plot to `plot.dir`. Default is
#'   `FALSE`.
#' @param plot.width Width of the saved plot in inches. Default is `NULL`,
#'   which auto-scales based on the number of detectors.
#' @param plot.height Height of the saved plot in inches. Default is `5`.
#' @param plot.dir Output directory for the saved plot. Default is
#'   `"./figure_spectral_variants"`.
#' @param variant.color Line color for the individual variant traces. Default
#'   is `"red"`.
#' @param variant.alpha Alpha transparency for the variant lines. Default is
#'   `0.05`.
#' @param variant.linewidth Thickness of the variant lines. Default is `1.5`.
#' @param median.line.color Color of the median reference line. Default is
#'   `"black"`.
#' @param median.linewidth Thickness of the median reference line. Default is
#'   `1`.
#'
#' @param color.by.laser Logical, default `TRUE`. Colour the detector axis
#'   labels by excitation laser (from `cytometer_database.csv`). Falls back to
#'   black labels if the detectors do not all map to the instrument's detectors.
#' @param asp Optional AutoSpectral parameter list, used only to identify the
#'   cytometer for `color.by.laser`. If `NULL` (default), the cytometer is
#'   identified from the detector names.
#' @return The `ggplot` object representing the variant density plot. When
#'   `save = TRUE`, a JPEG is also written to `plot.dir`.
#'
#' @export

spectral.variant.plot.dens <- function(
    spectra.variants,
    median.spectrum,
    title             = "Spectral_variants_density",
    save              = FALSE,
    plot.width        = NULL,
    plot.height       = 5,
    plot.dir          = "./figure_spectral_variants",
    variant.color     = "red",
    variant.alpha     = 0.05,
    variant.linewidth = 1.5,
    median.line.color = "black",
    median.linewidth  = 1,
    color.by.laser    = TRUE,
    asp               = NULL
) {

  detector.names <- colnames( spectra.variants )
  num_detectors  <- ncol( spectra.variants )
  detector.n     <- seq_len( num_detectors )

  # long format, base R
  df_variants     <- data.frame( spectra.variants, check.names = FALSE )
  df_variants$id  <- seq_len( nrow( df_variants ) )

  long_variants <- stats::reshape(
    df_variants,
    varying   = detector.names,
    v.names   = "intensity",
    timevar   = "detector_idx",
    times     = detector.n,
    direction = "long"
  )

  # median reference line data
  median.df <- data.frame(
    detector_idx = detector.n,
    intensity    = as.numeric( median.spectrum )
  )

  y.upper <- max(
    1,
    max( c( as.matrix( spectra.variants ), median.spectrum ), na.rm = TRUE )
  )
  y.lower <- max(
    -1,
    min( c( as.matrix( spectra.variants ), median.spectrum ), na.rm = TRUE )
  )

  # labels need clip = "off"; only safe if no data extend beyond the y limits
  in.range <- all( c( as.matrix( spectra.variants ), median.spectrum ) >= y.lower &
                     c( as.matrix( spectra.variants ), median.spectrum ) <= y.upper,
                   na.rm = TRUE )
  if ( isTRUE( color.by.laser ) && !in.range )
    message( "Detector labels will be black: data extend beyond the y-axis limits." )
  laser.labels <- .laser.axis.labels(
    x = detector.n, labels = detector.names, asp = asp,
    color.by.laser = color.by.laser && in.range, angle = 45, hjust = 1, vjust = 1
  )

  variant.plot <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = long_variants,
      ggplot2::aes( x = detector_idx, y = intensity, group = id ),
      linewidth = variant.linewidth,
      color     = variant.color,
      alpha     = variant.alpha
    ) +
    ggplot2::geom_line(
      data = median.df,
      ggplot2::aes( x = detector_idx, y = intensity ),
      linewidth = median.linewidth,
      color     = median.line.color
    ) +
    ggplot2::scale_x_continuous(
      breaks = detector.n,
      labels = detector.names
    ) +
    ggplot2::labs(
      x     = "Detector",
      y     = "Intensity",
      title = title
    ) +
    ggplot2::coord_cartesian( ylim = c( y.lower, y.upper ),
                              clip = if ( laser.labels$mapped ) "off" else "on" ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text( angle = 45, hjust = 1 )
    ) +
    laser.labels$components

  if ( save ) {
    if ( is.null( plot.width ) )
      plot.width <- max( ( ( length( detector.names ) - 1 ) / 64 * 12 ), 3 )
    if ( !dir.exists( plot.dir ) ) dir.create( plot.dir )

    ggplot2::ggsave(
      file.path( plot.dir, paste0( title, ".jpg" ) ),
      plot      = variant.plot,
      device    = ragg::agg_jpeg,
      width     = plot.width,
      height    = plot.height,
      limitsize = FALSE
    )
  }

  return( variant.plot )
}
