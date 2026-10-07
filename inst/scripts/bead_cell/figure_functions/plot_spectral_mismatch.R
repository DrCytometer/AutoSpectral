# plot_spectral_mismatch.R
#
# Panel-ready version of AutoSpectral's spectral.mismatch.plot(): overlay of a
# reference (cell) and test (bead) spectrum above a bar chart of their
# difference, with
#   - detector labels coloured by laser (from cytometer_database.csv),
#   - the cosine similarity of the two spectra printed on the plot,
#   - optional exclusion of lasers (e.g. the ID7000 deep UV detectors) from
#     the display only (the cosine similarity is calculated on ALL detectors
#     unless `cosine.on.plotted = TRUE`),
#   - optional L2 (unit Euclidean norm) normalisation of both spectra before
#     the difference is taken, over ALL detectors (before display exclusion),
#   - no file saving (the caller saves the panel).

#' Laser colour key (single row) for one or more mismatch plots
#'
#' @param lasers Laser names to show, in display order.
#' @return A legend grob for cowplot::plot_grid().
make.laser.key <- function( lasers, laser.colors = default.laser.colors, base.size = 9 ) {
  lasers <- unique( lasers )
  key.df <- data.frame( Laser = factor( lasers, levels = lasers ), x = 1, y = 1 )
  key.plot <- ggplot2::ggplot( key.df, ggplot2::aes( x, y, colour = Laser ) ) +
    ggplot2::geom_point( size = 2.5 ) +
    ggplot2::scale_colour_manual(
      values = stats::setNames( laser.colors[ lasers ], lasers ), name = "Laser"
    ) +
    ggplot2::guides( colour = ggplot2::guide_legend( nrow = 1 ) ) +
    ggplot2::theme_void( base_size = base.size ) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.margin = ggplot2::margin( 0, 0, 0, 0 )
    )
  cowplot::get_legend( key.plot )
}

#' @param ref.spectrum,test.spectrum Named numeric vectors or 1-row matrices /
#'   data.frames with identical detector names.
#' @param detector.lasers data.frame from get.detector.lasers().
#' @param ref.label,test.label Legend labels.
#' @param normalization "max" (input as supplied, peak = 1), "L2" (each
#'   spectrum divided by its Euclidean norm over all detectors) or "none".
#' @param exclude.lasers Character vector of laser names to drop from the plot.
#' @param laser.colors Named colour vector (names = laser names).
#' @param title Optional plot title (e.g. cytometer name).
#' @param trace.colors Colours for reference and test spectrum.
#' @param bar.color Fill for the difference bars.
#' @param label.every Show every n-th detector label. `NULL` chooses
#'   automatically (about 70 labels maximum).
#' @param base.size Base font size.
#' @param axis.text.size Detector label font size.
#' @param cosine.on.plotted Logical. Calculate cosine similarity only on the
#'   plotted detectors (default `FALSE`: all detectors).
#' @param show.laser.legend Logical. Add a laser colour key below the plot.
#'   Set `FALSE` when several plots share one key (see make.laser.key()).
#' @param laser.key.height Relative height of the laser key.
#' @param heights Relative heights of overlay and difference panels.
#' @return A cowplot object. The cosine similarity is attached as
#'   `attr(x, "cosine.similarity")` and the displayed lasers as
#'   `attr(x, "lasers")`.
plot.spectral.mismatch <- function(
    ref.spectrum,
    test.spectrum,
    detector.lasers,
    ref.label = "Cells",
    test.label = "BD Beads",
    normalization = c( "max", "L2", "none" ),
    exclude.lasers = NULL,
    laser.colors = default.laser.colors,
    title = NULL,
    trace.colors = c( "black", "red" ),
    bar.color = "lightblue",
    line.width = 0.5,
    point.size = 0.7,
    label.every = NULL,
    base.size = 9,
    axis.text.size = 5.5,
    cosine.on.plotted = FALSE,
    cosine.digits = 3,
    show.laser.legend = TRUE,
    laser.key.height = 0.07,
    heights = c( 2, 1 )
) {

  normalization <- match.arg( normalization )

  as.named.vector <- function( x ) {
    if ( is.matrix( x ) || is.data.frame( x ) ) {
      v <- as.numeric( x[ 1, ] )
      names( v ) <- colnames( x )
      v
    } else x
  }
  ref  <- as.named.vector( ref.spectrum )
  test <- as.named.vector( test.spectrum )

  if ( !identical( names( ref ), names( test ) ) )
    stop( "`ref.spectrum` and `test.spectrum` must have identical detector names.",
          call. = FALSE )

  detectors <- names( ref )

  # normalise over all detectors, before any display exclusion
  ref  <- normalize.spectrum( ref, normalization )
  test <- normalize.spectrum( test, normalization )
  y.label <- switch( normalization,
    max = "Normalized intensity", L2 = "Intensity (L2-normalized)", none = "Intensity" )

  # laser for each detector; order detectors as in the cytometer database
  lasers <- detector.lasers$Laser[ match( detectors, detector.lasers$Detector ) ]
  lasers[ is.na( lasers ) ] <- "Unknown"
  db.order <- match( detectors, detector.lasers$Detector )
  ord <- order( ifelse( is.na( db.order ), Inf, db.order ) )
  detectors <- detectors[ ord ]
  lasers <- lasers[ ord ]

  # cosine similarity
  cos.idx <- if ( cosine.on.plotted ) !( lasers %in% exclude.lasers ) else
    rep( TRUE, length( detectors ) )
  cos.sim <- cosine.similarity( rbind( ref[ detectors ][ cos.idx ], test[ detectors ][ cos.idx ] ) )[ 1, 2 ]

  # drop excluded lasers from the display
  show <- !( lasers %in% exclude.lasers )
  detectors <- detectors[ show ]
  lasers <- lasers[ show ]
  n.det <- length( detectors )

  # axis label thinning and colouring
  if ( is.null( label.every ) ) label.every <- max( 1L, ceiling( n.det / 70 ) )
  label.idx <- seq( 1, n.det, by = label.every )
  label.breaks <- detectors[ label.idx ]
  label.cols <- unname( laser.colors[ lasers[ label.idx ] ] )
  label.cols[ is.na( label.cols ) ] <- "grey50"

  det.levels <- detectors

  long <- data.frame(
    Detector  = factor( rep( detectors, 2 ), levels = det.levels ),
    Spectrum  = factor( rep( c( ref.label, test.label ), each = n.det ),
                        levels = c( ref.label, test.label ) ),
    Intensity = c( ref[ detectors ], test[ detectors ] )
  )

  trace.plot <- ggplot2::ggplot(
    long,
    ggplot2::aes( x = Detector, y = Intensity, group = Spectrum, colour = Spectrum )
  ) +
    ggplot2::geom_line( linewidth = line.width ) +
    ggplot2::geom_point( size = point.size ) +
    ggplot2::scale_colour_manual( values = stats::setNames( trace.colors, c( ref.label, test.label ) ) ) +
    ggplot2::scale_x_discrete( breaks = label.breaks ) +
    ggplot2::annotate(
      "text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.6,
      size = base.size * 0.32,
      label = sprintf( "cosine similarity = %s", formatC( cos.sim, digits = cosine.digits, format = "f" ) )
    ) +
    ggplot2::labs( title = title, x = NULL, y = y.label ) +
    ggplot2::theme_minimal( base_size = base.size ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      legend.margin = ggplot2::margin( 0, 0, 0, 0 ),
      plot.title = ggplot2::element_text( face = "bold", size = base.size + 1 )
    )

  delta <- data.frame(
    Detector  = factor( detectors, levels = det.levels ),
    Intensity = ref[ detectors ] - test[ detectors ]
  )

  # detector labels drawn as a text layer below the panel (one colour per
  # label as a layer parameter), with space reserved via the axis-title margin
  label.df <- data.frame(
    Detector = factor( label.breaks, levels = det.levels ),
    label    = paste0( label.breaks, " " )
  )
  label.space <- max( nchar( label.breaks ) ) * axis.text.size * 0.6 + 4

  delta.plot <- ggplot2::ggplot( delta, ggplot2::aes( x = Detector, y = Intensity ) ) +
    ggplot2::geom_col( fill = bar.color, colour = "black", alpha = 0.7, linewidth = 0.15 ) +
    ggplot2::geom_text(
      data = label.df, ggplot2::aes( x = Detector, y = -Inf, label = label ),
      inherit.aes = FALSE, colour = label.cols, angle = 90, hjust = 1, vjust = 0.5,
      size = axis.text.size / ggplot2::.pt
    ) +
    ggplot2::scale_x_discrete( breaks = label.breaks ) +
    ggplot2::coord_cartesian( clip = "off" ) +
    ggplot2::labs( x = "Detector", y = paste0( ref.label, " − ", test.label ) ) +
    ggplot2::theme_minimal( base_size = base.size ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text( margin = ggplot2::margin( t = label.space, unit = "pt" ) )
    )

  combined <- cowplot::plot_grid(
    trace.plot, delta.plot,
    ncol = 1, align = "v", axis = "lr", rel_heights = heights
  )

  if ( show.laser.legend ) {
    laser.key <- make.laser.key( unique( lasers ), laser.colors, base.size )
    combined <- cowplot::plot_grid(
      combined, laser.key, ncol = 1, rel_heights = c( 1, laser.key.height )
    )
  }

  attr( combined, "cosine.similarity" ) <- cos.sim
  attr( combined, "lasers" ) <- unique( lasers )
  combined
}
