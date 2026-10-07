# spectral_variant_plot.r

#' @title Spectral Variant Plot
#'
#' @description
#' This function plots a spectral ribbon showing the optimized fluorophore
#' signature with the range of variation from the spectral variants.
#'
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line scale_x_continuous labs
#' @importFrom ggplot2 theme_minimal theme ggsave coord_cartesian
#' @importFrom ragg agg_jpeg
#'
#' @param spectra.variants A matrix of variations in the normalized spectra
#' signature for a fluorophore, rows being the variants, columns being the
#' detectors.
#' @param median.spectrum The optimized reference spectrum for the fluorophore,
#' usually obtained from spectra prepared using AutoSpectral via
#' `get.fluorophore.spectra` on cleaned control data.
#' @param title Optional title to pass to the plot and output filename. Default is
#' `Spectral_variants`.
#' @param save Logical, controls whether the plot is displayed to the Viewer
#' (`FALSE`) or also saved to disk (`TRUE`).
#' @param plot.width Optional numeric to control saved plot width. Default is
#' `NULL`, in which case the plot width will be calculated based on the number
#' of detectors in the data.
#' @param plot.height Optional numeric to control saved plot width. Default is `5`.
#' @param plot.dir Directory where the files will be saved if `save = TRUE`.
#' Default is `./figure_spectral_variants`.
#' @param variant.fill.color Color for the shaded region indicating the range of
#' variation in the spectra. Feeds to `fill` in `geom_ribbon`. Default is "red".
#' @param variant.fill.alpha Transparency (alpha) for the color in
#' `variant.fill.color`. How intense the color of the variant spectra will be.
#' Default is `0.7`
#' @param median.line.color Color for the line representing the median or
#' optimized single spectrum. Default is "black".
#' @param median.linewidth Width of the line for the single optimized spectrum.
#' Default is `1`.
#'
#' @param color.by.laser Logical, default `TRUE`. Colour the detector axis
#'   labels by excitation laser (from `cytometer_database.csv`). Falls back to
#'   black labels if the detectors do not all map to the instrument's detectors.
#' @param asp Optional AutoSpectral parameter list, used only to identify the
#'   cytometer for `color.by.laser`. If `NULL` (default), the cytometer is
#'   identified from the detector names.
#' @return Plots are displayed to the Viewer. Files are saved if saving enabled.
#'
#' @export

spectral.variant.plot <- function(
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
    median.linewidth = 1,
    color.by.laser = TRUE,
    asp = NULL
  ) {

  detector.names <- colnames( spectra.variants )

  variant.data <- data.frame(
    detector.n = seq_len( ncol( spectra.variants ) ),
    detector = factor( detector.names, levels = detector.names ),
    median.spec = median.spectrum,
    min = apply( spectra.variants, 2, min, na.rm = TRUE ),
    max = apply( spectra.variants, 2, max, na.rm = TRUE )
  )

  y.upper <- 1
  y.lower <- max(
    -1,
    min( c( as.matrix( spectra.variants ), median.spectrum ), na.rm = TRUE )
  )

  # labels need clip = "off"; only safe if no data extend beyond the y limits
  in.range <- all( c( variant.data$min, variant.data$max, variant.data$median.spec ) >= y.lower &
                     c( variant.data$min, variant.data$max, variant.data$median.spec ) <= y.upper,
                   na.rm = TRUE )
  if ( isTRUE( color.by.laser ) && !in.range )
    message( "Detector labels will be black: data extend beyond the y-axis limits." )
  laser.labels <- .laser.axis.labels(
    x = variant.data$detector.n, labels = detector.names, asp = asp,
    color.by.laser = color.by.laser && in.range, angle = 45, hjust = 1, vjust = 1
  )

  variant.plot <- ggplot(
    data = variant.data,
    aes( x = detector.n )
    ) +
    geom_ribbon(
      aes( ymin = min, ymax = max ),
      fill = variant.fill.color,
      alpha = variant.fill.alpha
    ) +
    geom_line(
      aes( y = median.spec, group = 1 ),
      linewidth = median.linewidth,
      color = median.line.color
    ) +
    scale_x_continuous(
      breaks = variant.data$detector.n,
      labels = variant.data$detector
    ) +
    labs(
      x = "Detector",
      y = "Intensity",
      title = title
    ) +
    coord_cartesian( ylim = c( y.lower, y.upper ),
                     clip = if ( laser.labels$mapped ) "off" else "on" ) +
    theme_minimal() +
    theme( axis.text.x = element_text( angle = 45, hjust = 1 ) ) +
    laser.labels$components

  if ( save ) {
    if ( is.null( plot.width ) )
      plot.width <- max( ( ( length( detector.names ) - 1 ) / 64 * 12 ), 3 )

    if ( !dir.exists( plot.dir ) ) dir.create( plot.dir )

    ggsave(
      file.path( plot.dir, paste0( title, ".jpg" ) ),
      plot = variant.plot,
      device = ragg::agg_jpeg,
      width = plot.width,
      height = plot.height,
      limitsize = FALSE
    )
  } else {
    return( variant.plot )
  }

}
