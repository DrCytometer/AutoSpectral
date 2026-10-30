# laser_space_shift.R
#
# Instrument-independent description of how a spectrum changes from cells to
# beads. Detectors differ between cytometers but the lasers are shared, so each
# laser's bank of detectors is treated as one axis ("laser space"):
#
#   1. restrict both spectra to detectors on the lasers common to all
#      instruments and normalise each to unit L2 norm over those detectors;
#   2. bank norm n_L = || s_L ||_2, the length of the spectrum's projection onto
#      the detectors of laser L. Because the spectrum has unit norm,
#      sum_L n_L^2 = 1, so (n_L) is a unit vector in laser space;
#   3. shift per laser, either
#        "clr"   : centred log2 ratio  log2( n_L,beads / n_L,cells ) - mean_L( . )
#                  A per-laser gain (laser power, detector gain applied to a whole
#                  bank) multiplies n_L for cells and beads alike, so it cancels
#                  in the ratio, and any common scaling cancels in the centring.
#                  This makes the shifts comparable between instruments.
#        "delta" : n_L,beads - n_L,cells (in units of the unit-normalised spectrum).
#
# Angles (degrees) decompose the cell-bead mismatch:
#   total   : spectral angle between the two unit spectra (common detectors);
#   between : angle between the two laser-space vectors (redistribution of
#             signal BETWEEN lasers). By Cauchy-Schwarz, between <= total;
#   within  : per laser, angle between the bank sub-spectra (change of shape
#             WITHIN a laser's emission profile).
#
# Concordance across instruments: mean pairwise cosine similarity of the
# instruments' shift vectors, tested against a laser-label permutation null
# (laser labels shuffled independently within each instrument).
#
# The calculations, laser.space.shift() and laser.shift.concordance(), are
# exported AutoSpectral functions (R/laser_space_shift.R). This file keeps the
# laser order, lasers.present() and the spider-chart plot.

laser.order <- c( "DeepUV", "UV", "Violet", "Blue", "YellowGreen", "Red", "IR" )

.as.named.vector <- function( x ) {
  if ( is.matrix( x ) || is.data.frame( x ) ) {
    v <- as.numeric( x[ 1, ] ); names( v ) <- colnames( x ); v
  } else x
}

#' Lasers present on a cytometer's spectrum (in laser.order)
lasers.present <- function( spectrum, detector.lasers ) {
  v <- .as.named.vector( spectrum )
  las <- unique( detector.lasers$Laser[ match( names( v ), detector.lasers$Detector ) ] )
  las <- las[ !is.na( las ) ]
  c( intersect( laser.order, las ), setdiff( las, laser.order ) )
}

# Polar coordinate system with straight edges between vertices (spider chart)
coord.radar <- function( start = 0, direction = 1, clip = "off" ) {
  ggplot2::ggproto(
    "CoordRadar", ggplot2::CoordPolar,
    theta = "x", r = "y", start = start, direction = sign( direction ), clip = clip,
    is_linear = function( coord ) TRUE
  )
}

#' Spider chart of per-laser shifts, one unfilled trace per instrument
#'
#' @param shift.matrix Instruments x lasers.
#' @param instrument.colors Named colour vector (names = instruments).
#' @param laser.colors Colours for the laser axis labels.
#' @param r.limits Radial limits; `NULL` = data range with padding (the centre
#'   is below the most negative value, so negative shifts plot inside the
#'   dashed zero ring).
#' @param r.label Radial axis title.
#' @return ggplot object.
plot.laser.radar <- function(
    shift.matrix,
    instrument.colors,
    laser.colors = default.laser.colors,
    r.limits = NULL,
    r.label = "Shift (beads vs cells)",
    title = NULL,
    subtitle = NULL,
    base.size = 9,
    line.width = 0.8,
    point.size = 1.6
) {

  m <- as.matrix( shift.matrix )
  lasers <- colnames( m )
  n <- length( lasers )

  df <- data.frame(
    Instrument = factor( rep( rownames( m ), each = n ), levels = rownames( m ) ),
    pos        = rep( seq_len( n ), nrow( m ) ),
    value      = as.vector( t( m ) )
  )

  if ( is.null( r.limits ) ) {
    rng <- range( c( df$value, 0 ) )
    span <- diff( rng )
    if ( span == 0 ) span <- 1
    r.limits <- c( rng[ 1 ] - 0.35 * span, rng[ 2 ] + 0.3 * span )
  }
  r.breaks <- pretty( c( min( df$value, 0 ), max( df$value, 0 ) ), n = 4 )
  r.breaks <- r.breaks[ r.breaks >= r.limits[ 1 ] & r.breaks <= r.limits[ 2 ] ]

  zero.ring <- data.frame( pos = seq( 0.5, n + 0.5, length.out = 361 ), value = 0 )

  label.cols <- unname( laser.colors[ lasers ] )
  label.cols[ is.na( label.cols ) ] <- "grey30"
  # laser names drawn as a text layer at the outer radius (one colour per label)
  label.df <- data.frame( pos = seq_len( n ), value = r.limits[ 2 ], label = lasers )

  ggplot2::ggplot( df, ggplot2::aes( x = pos, y = value ) ) +
    ggplot2::geom_path( data = zero.ring, ggplot2::aes( x = pos, y = value ),
                        inherit.aes = FALSE, linetype = "dashed", colour = "grey40",
                        linewidth = 0.4 ) +
    ggplot2::geom_polygon( ggplot2::aes( colour = Instrument, group = Instrument ),
                           fill = NA, linewidth = line.width, key_glyph = "path" ) +
    ggplot2::geom_point( ggplot2::aes( colour = Instrument ), size = point.size ) +
    ggplot2::geom_text( data = label.df, ggplot2::aes( x = pos, y = value, label = label ),
                        inherit.aes = FALSE, colour = label.cols, fontface = "bold",
                        size = base.size / ggplot2::.pt ) +
    ggplot2::scale_colour_manual( values = instrument.colors[ rownames( m ) ], name = NULL ) +
    ggplot2::scale_x_continuous( limits = c( 0.5, n + 0.5 ), breaks = seq_len( n ),
                                 labels = lasers, expand = c( 0, 0 ) ) +
    ggplot2::scale_y_continuous( limits = r.limits, breaks = r.breaks ) +
    coord.radar( start = -pi / n ) +
    ggplot2::labs( x = NULL, y = r.label, title = title, subtitle = subtitle ) +
    ggplot2::theme_minimal( base_size = base.size ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      plot.title = ggplot2::element_text( face = "bold", size = base.size + 1 ),
      plot.subtitle = ggplot2::element_text( size = base.size - 1 )
    )
}
