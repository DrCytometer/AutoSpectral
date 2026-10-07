# create_brightness_biplot.R
#
# create.biplot()-style biexponential biplot for several populations at once,
# each coloured separately (e.g. one colour per antibody amount / brightness
# level), with spillover-spread boundary curves drawn exactly as
# create.biplot() draws them from `variants`:
#   y = flat.threshold +/- spread.kappa * sqrt( SS * |x| )
# where SS is the spillover spread of the x fluorophore into the y channel
# (see calculate.residual.spread()). Each curve stops where it reaches the
# axis limit (same rule as .spread.boundary() in AutoSpectral's
# create_biplot.R) instead of running along the plot edge.

# Boundary `flat + sign * kappa * sqrt( ss * |along| )` sampled along the
# other axis; the exact point where it reaches the axis limit is added and
# points beyond the limit are dropped. The curve is monotone in |along|, so the
# kept points are one contiguous run. No NA is returned: biexp.transform()
# cannot take NA (its convergence check `if ( residual > 1e-6 )` fails).
# Returns list( along, value ), possibly of length 0.
.fig.spread.boundary <- function( along, flat, ss, kappa, sign, lim.min, lim.max ) {
  if ( !is.finite( ss ) || ss < 0 ) ss <- 0
  lim   <- if ( sign > 0 ) lim.max else lim.min
  reach <- sign * ( lim - flat )
  if ( ss > 0 && kappa > 0 && reach > 0 ) {
    a.cross <- ( reach / kappa )^2 / ss
    a.cross <- c( -a.cross, a.cross )
    a.cross <- a.cross[ a.cross > min( along ) & a.cross < max( along ) ]
    along   <- sort( c( along, a.cross ) )
  }
  value <- flat + sign * kappa * sqrt( ss * abs( along ) )
  tol   <- 1e-9 * max( 1, abs( lim.min ), abs( lim.max ) )
  keep  <- value >= lim.min - tol & value <= lim.max + tol
  value <- pmin( pmax( value[ keep ], lim.min ), lim.max )
  list( along = along[ keep ], value = value )
}

#' @param populations Named list of matrices/data.frames of unmixed data (one
#'   per population, in plotting order: first is drawn first / at the back).
#' @param x.dim,y.dim Column names for the axes.
#' @param asp AutoSpectral parameter list.
#' @param population.colors Named colour vector (names = names(populations)).
#' @param legend.labels Optional named character vector of legend labels.
#' @param spread Optional list(ss = spillover spread x -> y, flat.pos = upper
#'   flat threshold in y, flat.neg = lower flat threshold in y). `NULL` draws
#'   no boundaries.
#' @param spread.kappa Number of spread SDs above the flat threshold (default 2,
#'   as in create.biplot()).
#' @param legend.inside Place the legend inside the panel (top left).
#' @return ggplot object.
create.brightness.biplot <- function(
    populations,
    x.dim,
    y.dim,
    asp,
    population.colors,
    legend.labels = NULL,
    spread = NULL,
    spread.kappa = 2,
    x.lab = NULL,
    y.lab = NULL,
    x.min = -5000,
    x.max = asp$expr.data.max,
    y.min = -5000,
    y.max = asp$expr.data.max,
    x.width.basis = -1000,
    y.width.basis = -1000,
    max.points.per.population = 5e4,
    point.size = asp$figure.gate.point.size,
    title = NULL,
    legend.inside = TRUE,
    legend.text.size = 7,
    spread.color = "red"
) {

  if ( is.null( x.lab ) ) x.lab <- x.dim
  if ( is.null( y.lab ) ) y.lab <- y.dim
  if ( is.null( legend.labels ) ) legend.labels <- stats::setNames( names( populations ), names( populations ) )

  # assemble plotting data (downsample large populations reproducibly)
  plot.data <- do.call( rbind, lapply( names( populations ), function( p ) {
    m <- populations[[ p ]]
    if ( nrow( m ) > max.points.per.population ) {
      set.seed( asp$bird.seed )
      m <- m[ sample( nrow( m ), max.points.per.population ), , drop = FALSE ]
    }
    data.frame( x = m[ , x.dim ], y = m[ , y.dim ], population = p )
  } ) )
  plot.data$population <- factor( plot.data$population, levels = names( populations ) )

  x.pos.log <- log10( x.max ) - 1
  y.pos.log <- log10( y.max ) - 1

  biexp.trans.x <- biexp.transform(
    channelRange = asp$default.transformation.param$length,
    maxValue = x.max, pos = x.pos.log,
    neg = asp$default.transformation.param$neg,
    widthBasis = x.width.basis, inverse = FALSE )
  biexp.trans.y <- biexp.transform(
    channelRange = asp$default.transformation.param$length,
    maxValue = y.max, pos = y.pos.log,
    neg = asp$default.transformation.param$neg,
    widthBasis = y.width.basis, inverse = FALSE )

  # axis breaks, extended to whole negative decades (as in create.biplot())
  extend.negative.breaks <- function( breaks, min.value ) {
    if ( !is.finite( min.value ) || min.value >= 0 ) return( breaks )
    neg.breaks <- breaks[ breaks < 0 ]
    max.decade.present <- if ( length( neg.breaks ) == 0 ) 0 else
      floor( log10( max( abs( neg.breaks ) ) ) )
    needed.decade <- floor( log10( abs( min.value ) ) )
    if ( needed.decade <= max.decade.present ) return( breaks )
    sort( c( breaks, -( 10 ^ seq( max.decade.present + 1, needed.decade ) ) ) )
  }
  x.breaks <- extend.negative.breaks( asp$ribbon.breaks, x.min )
  y.breaks <- extend.negative.breaks( asp$ribbon.breaks, y.min )
  x.breaks <- x.breaks[ x.breaks < x.max ]
  y.breaks <- y.breaks[ y.breaks < y.max ]
  make.labels <- function( b ) sapply( b, function( v )
    if ( v == 0 ) "0" else parse( text = paste0( "10^", log10( abs( v ) ) ) ) )

  plot.data$x.trans <- biexp.trans.x( plot.data$x )
  plot.data$y.trans <- biexp.trans.y( plot.data$y )

  biplot <- ggplot2::ggplot( plot.data, ggplot2::aes( x.trans, y.trans, colour = population ) ) +
    scattermore::geom_scattermore( pointsize = point.size, alpha = 1, na.rm = TRUE ) +
    ggplot2::scale_colour_manual(
      values = population.colors[ names( populations ) ],
      labels = legend.labels[ names( populations ) ],
      name = NULL
    ) +
    ggplot2::guides( colour = ggplot2::guide_legend( override.aes = list( size = 2.5 ) ) ) +
    ggplot2::scale_x_continuous(
      name = x.lab,
      breaks = biexp.trans.x( x.breaks ),
      limits = biexp.trans.x( c( x.min, x.max ) ),
      labels = make.labels( x.breaks )
    ) +
    ggplot2::scale_y_continuous(
      name = y.lab,
      breaks = biexp.trans.y( y.breaks ),
      limits = biexp.trans.y( c( y.min, y.max ) ),
      labels = make.labels( y.breaks )
    ) +
    ggplot2::labs( title = title ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(
        asp$figure.margin, asp$figure.margin, asp$figure.margin, asp$figure.margin ),
      axis.ticks = ggplot2::element_line( linewidth = asp$figure.panel.line.size ),
      axis.text = ggplot2::element_text( size = asp$figure.axis.text.size ),
      axis.title = ggplot2::element_text( size = asp$figure.axis.title.size ),
      panel.border = ggplot2::element_rect( fill = NA, linewidth = asp$figure.panel.line.size ),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text( size = asp$figure.axis.title.size, face = "bold" ),
      legend.text = ggplot2::element_text( size = legend.text.size ),
      legend.key.size = grid::unit( 0.35, "cm" ),
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank()
    )

  if ( legend.inside ) {
    biplot <- biplot + ggplot2::theme(
      legend.position = c( 0.02, 0.98 ),
      legend.justification = c( 0, 1 )
    )
  }

  # spillover-spread boundaries (same formula as create.biplot())
  if ( !is.null( spread ) ) {
    x.seq <- seq( x.min, x.max, length.out = 300 )
    upper <- .fig.spread.boundary( x.seq, spread$flat.pos, spread$ss, spread.kappa, 1, y.min, y.max )
    lower <- .fig.spread.boundary( x.seq, spread$flat.neg, spread$ss, spread.kappa, -1, y.min, y.max )
    for ( curve in list( upper, lower ) ) {
      if ( length( curve$along ) < 2 ) next
      curve.df <- data.frame( x = biexp.trans.x( curve$along ), y = biexp.trans.y( curve$value ) )
      biplot <- biplot +
        ggplot2::geom_line( data = curve.df, ggplot2::aes( x = x, y = y ),
                            colour = spread.color, inherit.aes = FALSE, na.rm = TRUE )
    }
  }

  biplot
}
