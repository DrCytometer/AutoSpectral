# plot_fig4.R
#
# Panels for Figure 4 ("Is it the fluorophores?"):
#   plot.value.by.class()     a value (angle, variability) per fluorophore by
#                             dye class, one facet per cytometer; optional
#                             dodge by particle type (cells vs bead types)
#   plot.effect.forest()      forest plot of fold differences / fold changes
#                             (2^estimate) per cytometer and pooled; used for
#                             the class contrasts and for the slopes
#   plot.property.scatter()   angle vs a fluorophore property, per cytometer,
#                             with a descriptive least-squares line on the
#                             plotted (log) scales
#
# Shared helpers (.ratio.labels, .fmt.p, default.bead.colors) come from
# plot_bead_angles.R.

# Defaults only; Figure_4.R sets the colours from RColorBrewer (classes: Dark2;
# particle types: Set1 with cells in grey).
default.class.colors <- c(
  Small            = "#1B9E77",
  Brilliant        = "#D95F02",
  Brilliant_tandem = "#7570B3",
  Real             = "#E7298A",
  Protein          = "#66A61E",
  Protein_tandem   = "#E6AB02"
)
default.class.labels <- c(
  Small            = "Small molecule",
  Brilliant        = "Brilliant polymer",
  Brilliant_tandem = "Brilliant tandem",
  Real             = "Real",
  Protein          = "Protein",
  Protein_tandem   = "Protein tandem"
)

.fig4.theme <- function( base.size ) {
  ggplot2::theme_bw( base_size = base.size ) +
    ggplot2::theme(
      axis.text        = ggplot2::element_text( colour = "black" ),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect( fill = "grey92", colour = NA ),
      strip.text       = ggplot2::element_text( face = "bold" )
    )
}

.log.breaks <- function( lim, candidates = c( 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100 ) )
  candidates[ candidates >= lim[ 1 ] / 3 & candidates <= lim[ 2 ] * 3 ]

#' Value per fluorophore by dye class
#'
#' @param d data.frame with Cytometer, Fluorophore, Class, the value column and
#'   (if `dodge.by` is set) that column.
#' @param value Column to plot (log10 axis).
#' @param dodge.by NULL (one box per class, points coloured by class) or a
#'   column (e.g. "Particle") giving side-by-side boxes within each class.
#' @param dodge.colors Named colours for the dodge groups.
#' @param untested.classes Classes left out of the statistics; marked with a
#'   dagger in the axis labels.
plot.value.by.class <- function(
    d, value, y.label,
    class.order, class.labels = default.class.labels, class.colors = default.class.colors,
    cytometer.order,
    dodge.by = NULL, dodge.order = NULL, dodge.colors = NULL, dodge.labels = NULL,
    untested.classes = character( 0 ),
    point.size = 0.9, point.alpha = 0.8, jitter.width = 0.15,
    title = NULL, base.size = 9, seed = 42
) {
  d <- as.data.frame( d )
  d <- d[ is.finite( d[[ value ]] ) & d[[ value ]] > 0 & d$Class %in% class.order, , drop = FALSE ]
  d$Cytometer <- factor( d$Cytometer, levels = intersect( cytometer.order, unique( d$Cytometer ) ) )
  d$Class <- factor( d$Class, levels = intersect( class.order, unique( d$Class ) ) )
  d$Value <- d[[ value ]]

  x.labels <- class.labels[ levels( d$Class ) ]
  x.labels[ is.na( x.labels ) ] <- levels( d$Class )[ is.na( x.labels ) ]
  x.labels <- ifelse( levels( d$Class ) %in% untested.classes, paste0( x.labels, "†" ), x.labels )
  names( x.labels ) <- levels( d$Class )

  y.lim <- range( d$Value )
  y.breaks <- .log.breaks( y.lim )

  set.seed( seed )
  if ( is.null( dodge.by ) ) {
    p <- ggplot2::ggplot( d, ggplot2::aes( x = Class, y = Value ) ) +
      ggplot2::geom_boxplot( ggplot2::aes( fill = Class ), outlier.shape = NA, width = 0.6,
                             alpha = 0.3, colour = "grey25", linewidth = 0.3 ) +
      ggplot2::geom_point( ggplot2::aes( colour = Class ), size = point.size, alpha = point.alpha,
                           stroke = 0, position = ggplot2::position_jitter( width = jitter.width, height = 0 ) ) +
      ggplot2::scale_fill_manual( values = class.colors, guide = "none" ) +
      ggplot2::scale_colour_manual( values = class.colors, guide = "none" )
  } else {
    if ( is.null( dodge.order ) ) dodge.order <- sort( unique( d[[ dodge.by ]] ) )
    d$Dodge <- factor( d[[ dodge.by ]], levels = intersect( dodge.order, unique( d[[ dodge.by ]] ) ) )
    if ( is.null( dodge.labels ) ) dodge.labels <- stats::setNames( levels( d$Dodge ), levels( d$Dodge ) )

    # explicit x positions: within each cytometer, the particle types present
    # there are spread evenly across `dw` around the class position (facets
    # can hold different numbers of particle types)
    dw <- 0.8
    d$x.pos <- NA_real_
    n.max <- 1
    for ( cy in levels( d$Cytometer ) ) {
      i <- which( d$Cytometer == cy )
      present <- intersect( levels( d$Dodge ), as.character( d$Dodge[ i ] ) )
      n <- length( present )
      n.max <- max( n.max, n )
      off <- ( match( as.character( d$Dodge[ i ] ), present ) - ( n + 1 ) / 2 ) * dw / n
      d$x.pos[ i ] <- as.integer( d$Class[ i ] ) + off
    }
    box.w <- 0.85 * dw / n.max
    d$x.pt <- d$x.pos + stats::runif( nrow( d ), -0.3, 0.3 ) * box.w     # points over 60% of the box
    d$Box <- interaction( d$Cytometer, d$Class, d$Dodge, drop = TRUE )

    p <- ggplot2::ggplot( d, ggplot2::aes( y = Value, fill = Dodge, colour = Dodge ) ) +
      ggplot2::geom_boxplot( ggplot2::aes( x = x.pos, group = Box ), outlier.shape = NA, width = box.w,
                             alpha = 0.3, colour = "grey25", linewidth = 0.25 ) +
      ggplot2::geom_point( ggplot2::aes( x = x.pt ), size = point.size * 0.8, alpha = point.alpha,
                           stroke = 0 ) +
      ggplot2::scale_fill_manual( values = dodge.colors, labels = dodge.labels, name = NULL ) +
      ggplot2::scale_colour_manual( values = dodge.colors, labels = dodge.labels, name = NULL ) +
      ggplot2::guides( fill = ggplot2::guide_legend( nrow = 1, override.aes = list( alpha = 0.6 ) ),
                       colour = "none" )
  }

  x.scale <- if ( is.null( dodge.by ) ) {
    ggplot2::scale_x_discrete( labels = x.labels )
  } else {
    ggplot2::scale_x_continuous( breaks = seq_along( levels( d$Class ) ), labels = unname( x.labels ),
                                 expand = ggplot2::expansion( add = 0.5 ) )
  }

  p +
    ggplot2::facet_grid( . ~ Cytometer ) +
    ggplot2::scale_y_log10( breaks = y.breaks, labels = as.character( y.breaks ) ) +
    x.scale +
    ggplot2::labs( x = NULL, y = y.label, title = title ) +
    .fig4.theme( base.size ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text( angle = 45, hjust = 1, vjust = 1 ),
      panel.grid.major.x = ggplot2::element_blank(),
      legend.position = "top",
      legend.margin = ggplot2::margin( 0, 0, 0, 0 )
    )
}

#' Forest plot of 2^estimate with 95% CI, per cytometer and pooled
#'
#' @param est data.frame with Cytometer ("Pooled" for the pooled rows), Term
#'   (row label key), Estimate, Lower, Upper (log2 scale) and P.holm; optional
#'   column Group for side-by-side estimates within a term (e.g. cells vs beads).
#' @param term.order,term.labels Order and labels of the terms (y axis).
#' @param group.order,group.shapes Optional dodge groups (e.g. Cells, Beads);
#'   when given, the shape is set by group and the group is printed in the
#'   label. Symbols are all filled; the Holm-adjusted P value is printed for
#'   every estimate (no significance coding).
#' @param colour.by "none" (all symbols `point.colors[1]`, default black),
#'   "term" or "group": which variable `point.colors` (a named vector) maps.
#'   The terms are already labelled on the y axis, so colouring by term is
#'   optional; "none" avoids a second categorical palette next to the
#'   particle-type colours.
#' @param x.label Axis label.
plot.effect.forest <- function(
    est, term.order, term.labels = NULL,
    cytometer.order, pooled.label = "Pooled",
    group.order = NULL, group.shapes = NULL,
    colour.by = c( "none", "term", "group" ), point.colors = c( all = "black" ),
    x.label, title = NULL, subtitle = NULL,
    point.size = 2.2, label.size = 2.5, label.margin = 150, dodge = 0.32,
    base.size = 9
) {
  colour.by <- match.arg( colour.by )
  df <- as.data.frame( est )
  df <- df[ df$Term %in% term.order, , drop = FALSE ]
  if ( is.null( term.labels ) ) term.labels <- stats::setNames( term.order, term.order )

  df$Cytometer[ df$Cytometer == "Pooled" ] <- pooled.label
  df$Pooled <- df$Cytometer == pooled.label
  df$Cytometer <- factor( df$Cytometer, levels = c( intersect( cytometer.order, df$Cytometer ), pooled.label ) )

  terms <- rev( intersect( term.order, unique( df$Term ) ) )
  df$y <- match( df$Term, terms )
  has.group <- !is.null( group.order ) && "Group" %in% names( df )
  if ( has.group ) {
    g <- match( df$Group, group.order )
    ng <- length( group.order )
    df$y <- df$y + ( ( ng - 1 ) / 2 - ( g - 1 ) ) * dodge
  }

  df$Ratio <- 2^df$Estimate
  df$Lo <- 2^df$Lower
  df$Hi <- 2^df$Upper
  df$Label <- sprintf( "%s%.2f (%.2f–%.2f)  %s",
                       if ( has.group ) paste0( df$Group, ": " ) else "",
                       df$Ratio, df$Lo, df$Hi, .fmt.p( df$P.holm ) )
  df$Col <- switch( colour.by,
                    none  = rep( names( point.colors )[ 1 ], nrow( df ) ),
                    term  = as.character( df$Term ),
                    group = as.character( df$Group ) )
  df$Fill <- df$Col

  if ( has.group ) {
    if ( is.null( group.shapes ) ) group.shapes <- stats::setNames( c( 21, 24, 22, 25 )[ seq_along( group.order ) ], group.order )
    df$Shape <- factor( df$Group, levels = group.order )
  } else {
    df$Shape <- factor( ifelse( df$Pooled, "Pooled", "Cytometer" ), levels = c( "Cytometer", "Pooled" ) )
    group.shapes <- c( Cytometer = 21, Pooled = 23 )
  }

  x.lim <- c( min( df$Lo, 1, na.rm = TRUE ) / 1.15, max( df$Hi, 1, na.rm = TRUE ) * 1.15 )
  x.breaks <- 2^( -8:8 )
  x.breaks <- x.breaks[ x.breaks >= x.lim[ 1 ] & x.breaks <= x.lim[ 2 ] ]
  if ( length( x.breaks ) < 3 ) x.breaks <- signif( exp( seq( log( x.lim[ 1 ] ), log( x.lim[ 2 ] ), length.out = 4 ) ), 2 )

  lab.fun <- function( x ) vapply( x, function( v ) {
    if ( is.na( v ) ) return( "" )
    if ( v < 1 && abs( 1 / v - round( 1 / v ) ) < 1e-9 ) paste0( "1/", round( 1 / v ) ) else format( signif( v, 2 ) )
  }, character( 1 ) )

  ggplot2::ggplot( df, ggplot2::aes( x = Ratio, y = y, colour = Col ) ) +
    ggplot2::geom_vline( xintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 0.3 ) +
    ggplot2::geom_errorbar( ggplot2::aes( xmin = Lo, xmax = Hi ), width = 0.2, linewidth = 0.45,
                            orientation = "y" ) +
    ggplot2::geom_point( ggplot2::aes( shape = Shape, fill = Fill ), size = point.size, stroke = 0.6 ) +
    ggplot2::geom_text( ggplot2::aes( x = Inf, label = Label ), hjust = -0.05, size = label.size,
                        colour = "black" ) +
    ggplot2::facet_grid( Cytometer ~ ., switch = "y" ) +
    ggplot2::scale_x_continuous( trans = "log2", breaks = x.breaks, labels = lab.fun ) +
    ggplot2::scale_y_continuous( breaks = seq_along( terms ), labels = unname( term.labels[ terms ] ),
                                 limits = c( 0.4, length( terms ) + 0.6 ),
                                 expand = ggplot2::expansion( 0 ) ) +
    ggplot2::coord_cartesian( xlim = x.lim, clip = "off" ) +
    ggplot2::scale_colour_manual( values = point.colors, guide = "none" ) +
    ggplot2::scale_fill_manual( values = point.colors, guide = "none" ) +
    ggplot2::scale_shape_manual( values = group.shapes, name = NULL,
                                 guide = if ( has.group ) "legend" else "none" ) +
    ggplot2::labs( x = x.label, y = NULL, title = title, subtitle = subtitle ) +
    .fig4.theme( base.size ) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text( angle = 0, hjust = 1 ),
      legend.position = "top",
      plot.margin = ggplot2::margin( 5.5, label.margin, 5.5, 5.5 )
    )
}

#' Angle (or another log-scale response) vs a fluorophore property
#'
#' @param d data.frame with Cytometer, Panel (facet row label, optional), X, Y
#'   and the colour column.
#' @param colour.by Column used for the point colour.
#' @param x.log Log10 x axis.
#' @param x.breaks Optional x breaks for the log axis (default: decades/halves).
#' @param fit.line Draw a least-squares line per facet on the plotted scales
#'   (descriptive; the statistics come from the mixed models).
plot.property.scatter <- function(
    d, x.label, y.label, cytometer.order,
    colour.by, colour.values, colour.labels = NULL, colour.title = NULL,
    x.log = TRUE, y.log = TRUE, fit.line = TRUE, x.breaks = NULL,
    point.size = 1, point.alpha = 0.7, facet.rows = NULL,
    title = NULL, base.size = 9
) {
  d <- as.data.frame( d )
  ok <- is.finite( d$X ) & is.finite( d$Y )
  if ( x.log ) ok <- ok & d$X > 0
  if ( y.log ) ok <- ok & d$Y > 0
  d <- d[ ok, , drop = FALSE ]
  d$Cytometer <- factor( d$Cytometer, levels = intersect( cytometer.order, unique( d$Cytometer ) ) )
  d$Colour <- factor( d[[ colour.by ]], levels = intersect( names( colour.values ), unique( d[[ colour.by ]] ) ) )
  if ( is.null( colour.labels ) ) colour.labels <- stats::setNames( names( colour.values ), names( colour.values ) )

  p <- ggplot2::ggplot( d, ggplot2::aes( x = X, y = Y ) ) +
    ggplot2::geom_point( ggplot2::aes( colour = Colour ), size = point.size, alpha = point.alpha, stroke = 0 )
  if ( fit.line )
    p <- p + ggplot2::geom_smooth( method = "lm", formula = y ~ x, se = TRUE, colour = "black",
                                   linewidth = 0.4, fill = "grey70", alpha = 0.3 )
  if ( is.null( facet.rows ) ) {
    p <- p + ggplot2::facet_grid( . ~ Cytometer )
  } else {
    p <- p + ggplot2::facet_grid( stats::as.formula( paste( facet.rows, "~ Cytometer" ) ),
                                  scales = "free_x", switch = "y" )
  }
  if ( x.log ) {
    xb <- if ( is.null( x.breaks ) ) .log.breaks( range( d$X ) ) else x.breaks
    p <- p + ggplot2::scale_x_log10( breaks = xb, labels = as.character( xb ) )
  }
  if ( y.log ) {
    yb <- .log.breaks( range( d$Y ) )
    p <- p + ggplot2::scale_y_log10( breaks = yb, labels = as.character( yb ) )
  }
  p +
    ggplot2::scale_colour_manual( values = colour.values, labels = colour.labels, name = colour.title ) +
    ggplot2::guides( colour = ggplot2::guide_legend( nrow = 1, override.aes = list( size = 2, alpha = 1 ) ) ) +
    ggplot2::labs( x = x.label, y = y.label, title = title ) +
    .fig4.theme( base.size ) +
    ggplot2::theme( legend.position = "top", strip.placement = "outside",
                    legend.margin = ggplot2::margin( 0, 0, 0, 0 ) )
}
