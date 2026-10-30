# plot_bead_angles.R
#
# Panels for Figure 2:
#   plot.bead.angles()  spectral angle (beads vs cells) per fluorophore, by bead
#                       type, one facet per cytometer; the same fluorophore is
#                       joined across bead types (paired design)
#   plot.bead.forest()  forest plot of the fold difference in spectral angle
#                       (ratio of geometric means) of each bead type versus a
#                       reference bead, per cytometer and pooled

# Defaults only (RColorBrewer "Set1", as in figure_palettes.R); the figure
# scripts pass manuscript.colors()$bead.
default.bead.colors <- c(
  BD                = "#E41A1C",
  BioLegend         = "#984EA3",
  SpectraComp       = "#FF7F00",
  UltraComp         = "#377EB8",
  UltraCompSpectral = "#4DAF4A"
)

#' Labels "1/4", "1/2", "1", "2", ... for a log2 ratio axis
.ratio.labels <- function( x ) ifelse( x < 1, paste0( "1/", round( 1 / x ) ), as.character( x ) )

#' Format exact P values for plot labels
#'
#' P >= 0.001 is printed with `digits` significant digits ("P = 0.0014",
#' "P = 0.19"); smaller values in scientific notation with a superscript
#' exponent ("P = 2.3 × 10⁻⁵"). No thresholds: every value is exact. For a
#' permutation test the smallest attainable value is 1 / (1 + permutations).
#' @param p Numeric vector of P values.
#' @param digits Significant digits.
#' @param prefix Text before the value.
.fmt.p <- function( p, digits = 2, prefix = "P = " ) {
  sup <- c( "-" = "\u207b", "0" = "\u2070", "1" = "\u00b9", "2" = "\u00b2", "3" = "\u00b3",
            "4" = "\u2074", "5" = "\u2075", "6" = "\u2076", "7" = "\u2077", "8" = "\u2078",
            "9" = "\u2079" )
  vapply( p, function( x ) {
    if ( is.na( x ) ) return( paste0( prefix, "NA" ) )
    if ( x >= 1e-3 )
      return( paste0( prefix, formatC( signif( x, digits ), digits = digits, format = "fg", flag = "#" ) ) )
    e <- floor( log10( x ) )
    m <- round( x / 10^e, digits - 1 )
    if ( m >= 10 ) { m <- m / 10; e <- e + 1 }
    paste0( prefix, formatC( m, digits = digits - 1, format = "f" ), " \u00d7 10",
            paste( sup[ strsplit( as.character( e ), "" )[[ 1 ]] ], collapse = "" ) )
  }, character( 1 ), USE.NAMES = FALSE )
}

#' Panel a: spectral angle by bead type, grouped by cytometer
#'
#' @param d Long table (Cytometer, Bead, Fluorophore, Angle), e.g. from
#'   prepare.angle.model.data().
#' @param bead.order,cytometer.order Display orders.
#' @param bead.colors Named colours.
#' @param bead.labels Named display labels for the x axis (names = bead types).
#' @param show.pairs Join each fluorophore across bead types (grey lines).
#' @param jitter.width Total horizontal spread of the points. Each fluorophore
#'   gets the same offset in every bead column (so the joining lines stay
#'   level); 0 = no spread.
#' @param y.breaks Breaks of the log10 y axis (degrees).
plot.bead.angles <- function(
    d,
    bead.order,
    cytometer.order,
    bead.colors = default.bead.colors,
    bead.labels = NULL,
    show.pairs = TRUE,
    jitter.width = 0.3,
    point.size = 0.9,
    point.alpha = 0.8,
    line.alpha = 0.35,
    y.breaks = c( 0.1, 0.3, 1, 3, 10, 30 ),
    y.label = "Spectral angle, beads vs cells (°)",
    title = NULL,
    base.size = 9
) {

  d <- as.data.frame( d )
  d <- d[ d$Bead %in% bead.order & d$Cytometer %in% cytometer.order, , drop = FALSE ]
  d$Cytometer <- factor( as.character( d$Cytometer ), levels = cytometer.order )
  d$Bead      <- factor( as.character( d$Bead ), levels = bead.order )

  # position of each bead within its own facet (free x scale drops absent beads)
  d$pos <- NA_real_
  for ( cy in levels( d$Cytometer ) ) {
    i <- d$Cytometer == cy
    d$pos[ i ] <- match( as.character( d$Bead[ i ] ),
                         intersect( bead.order, as.character( d$Bead[ i ] ) ) )
  }
  fl <- sort( unique( d$Fluorophore ) )
  off <- if ( length( fl ) > 1 ) ( match( d$Fluorophore, fl ) - 1 ) / ( length( fl ) - 1 ) - 0.5 else 0
  d$x.pt <- d$pos + jitter.width * off
  d <- d[ order( d$Cytometer, d$Fluorophore, d$pos ), ]

  if ( is.null( bead.labels ) ) bead.labels <- stats::setNames( bead.order, bead.order )
  y.lim <- range( d$Angle )
  y.breaks <- y.breaks[ y.breaks >= y.lim[ 1 ] / 3 & y.breaks <= y.lim[ 2 ] * 3 ]

  p <- ggplot2::ggplot( d, ggplot2::aes( x = Bead, y = Angle ) ) +
    ggplot2::geom_boxplot( ggplot2::aes( fill = Bead ), outlier.shape = NA, width = 0.6,
                           alpha = 0.3, colour = "grey25", linewidth = 0.3 )
  if ( show.pairs )
    p <- p + ggplot2::geom_line( ggplot2::aes( x = x.pt, group = Fluorophore ),
                                 colour = "grey55", linewidth = 0.2, alpha = line.alpha )
  p <- p +
    ggplot2::geom_point( ggplot2::aes( x = x.pt, colour = Bead ), size = point.size,
                         alpha = point.alpha, stroke = 0 ) +
    ggplot2::facet_grid( . ~ Cytometer, scales = "free_x", space = "free_x" ) +
    ggplot2::scale_y_log10( breaks = y.breaks, labels = as.character( y.breaks ) ) +
    ggplot2::scale_x_discrete( labels = bead.labels ) +
    ggplot2::scale_fill_manual( values = bead.colors, guide = "none" ) +
    ggplot2::scale_colour_manual( values = bead.colors, guide = "none" ) +
    ggplot2::labs( x = NULL, y = y.label, title = title ) +
    ggplot2::theme_bw( base_size = base.size ) +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text( angle = 45, hjust = 1, vjust = 1, colour = "black" ),
      axis.text.y      = ggplot2::element_text( colour = "black" ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect( fill = "grey92", colour = NA ),
      strip.text       = ggplot2::element_text( face = "bold" )
    )
  p
}

#' Panel b: forest plot of fold differences in spectral angle
#'
#' @param per.cyto Per-cytometer contrasts (fit.bead.angle.models()$per.cyto.vs.ref):
#'   columns Cytometer, Bead1, Bead2, Estimate (log2 ratio), Lower, Upper, P.holm.
#' @param pooled Pooled contrasts (fit.bead.angle.models()$pooled.vs.ref).
#' @param reference.bead Reference bead (denominator of the ratio).
#' @param pooled.label Facet label for the pooled estimates.
#' @param subtitle Optional subtitle (e.g. interaction test).
#' @param label.margin Right margin (pt) reserved for the estimate/P text.
#' Symbols are all filled; the Holm-adjusted P value is printed for every
#' estimate (no significance coding).
plot.bead.forest <- function(
    per.cyto,
    pooled,
    reference.bead,
    bead.order,
    cytometer.order,
    bead.colors = default.bead.colors,
    bead.labels = NULL,
    pooled.label = "Pooled",
    title = NULL,
    subtitle = NULL,
    point.size = 2.2,
    label.size = 2.6,
    label.margin = 140,
    base.size = 9
) {

  cols <- c( "Cytometer", "Bead1", "Estimate", "Lower", "Upper", "P.holm" )
  df <- rbind(
    data.frame( per.cyto[ , cols ], Pooled = FALSE, stringsAsFactors = FALSE ),
    data.frame( pooled[ , cols ], Pooled = TRUE, stringsAsFactors = FALSE )
  )
  df$Cytometer[ df$Pooled ] <- pooled.label
  df$Cytometer <- factor( df$Cytometer, levels = c( intersect( cytometer.order, df$Cytometer ), pooled.label ) )

  if ( is.null( bead.labels ) ) bead.labels <- stats::setNames( bead.order, bead.order )
  b.lev <- rev( setdiff( bead.order, reference.bead ) )
  df$Bead <- factor( df$Bead1, levels = b.lev )

  df$Ratio <- 2^df$Estimate
  df$Lo    <- 2^df$Lower
  df$Hi    <- 2^df$Upper
  df$Label <- sprintf( "%.2f (%.2f–%.2f)  %s", df$Ratio, df$Lo, df$Hi, .fmt.p( df$P.holm ) )

  x.lim <- c( min( df$Lo, 1 ) / 1.15, max( df$Hi, 1 ) * 1.15 )
  x.breaks <- 2^( -6:6 )
  x.breaks <- x.breaks[ x.breaks >= x.lim[ 1 ] & x.breaks <= x.lim[ 2 ] ]

  ggplot2::ggplot( df, ggplot2::aes( x = Ratio, y = Bead, colour = Bead1 ) ) +
    ggplot2::geom_vline( xintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 0.3 ) +
    ggplot2::geom_errorbar( ggplot2::aes( xmin = Lo, xmax = Hi ), width = 0.25, linewidth = 0.45,
                            orientation = "y" ) +
    ggplot2::geom_point( ggplot2::aes( shape = Pooled, fill = Bead1 ), size = point.size, stroke = 0.6 ) +
    ggplot2::geom_text( ggplot2::aes( x = Inf, label = Label ), hjust = -0.05, size = label.size,
                        colour = "black" ) +
    ggplot2::facet_grid( Cytometer ~ ., scales = "free_y", space = "free_y", switch = "y" ) +
    ggplot2::scale_x_continuous( trans = "log2", breaks = x.breaks, labels = .ratio.labels ) +
    ggplot2::coord_cartesian( xlim = x.lim, clip = "off" ) +
    ggplot2::scale_y_discrete( labels = bead.labels ) +
    ggplot2::scale_colour_manual( values = bead.colors, guide = "none" ) +
    ggplot2::scale_fill_manual( values = bead.colors, guide = "none" ) +
    ggplot2::scale_shape_manual( values = c( `FALSE` = 21, `TRUE` = 23 ), guide = "none" ) +
    ggplot2::labs(
      x = paste0( "Fold difference in spectral angle vs ", bead.labels[[ reference.bead ]],
                  "\n(ratio of geometric means, 95% CI)" ),
      y = NULL, title = title, subtitle = subtitle
    ) +
    ggplot2::theme_bw( base_size = base.size ) +
    ggplot2::theme(
      axis.text        = ggplot2::element_text( colour = "black" ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      strip.placement  = "outside",
      strip.background = ggplot2::element_rect( fill = "grey92", colour = NA ),
      strip.text.y.left = ggplot2::element_text( angle = 0, face = "bold", hjust = 1 ),
      plot.margin      = ggplot2::margin( 5.5, label.margin, 5.5, 5.5 )
    )
}
