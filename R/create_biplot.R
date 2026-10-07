# create_biplot.r

#' @title Create Biplot
#'
#' @importFrom ggplot2 ggplot aes ggsave after_stat
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous theme_bw theme
#' @importFrom ggplot2 margin element_line element_text element_rect element_blank
#' @importFrom ggplot2 scale_fill_viridis_c scale_fill_gradientn stat_density_2d
#' @importFrom ggplot2 geom_line
#' @importFrom scattermore geom_scattermore
#' @importFrom ragg agg_jpeg
#'
#' @param plot.data A matrix or dataframe containing the flow cytometry data to
#' be plotted. Column names should match the dimensions specified by `x.dim` and
#' `y.dim`.
#' @param x.dim String specifying the column of `plot.data` for the x-axis of
#' the plot.
#' @param y.dim String specifying the column of `plot.data` for the y-axis of
#' the plot.
#' @param asp The AutoSpectral parameter list.
#' @param variants The variant list returned by `get.spectral.variants()`.
#'   When supplied, red curves are drawn for whichever of `x.dim`/`y.dim`
#'   `variants$thresholds` covers: each channel's flat threshold
#'   (`variants$thresholds`) plus `spread.kappa` spread standard deviations
#'   from the *other* axis's own spillover into it
#'   (`variants$spillover.spread`), and, below zero, the same spread widening
#'   subtracted from the directly measured negative-tail flat threshold
#'   (`variants$neg.thresholds`) rather than mirrored from the positive curve
#'   about zero. Falls back to the mirrored positive threshold, with a
#'   warning, when `variants$neg.thresholds` is absent (an older cached
#'   `variants` object). This is the one-source restriction of
#'   `get.spread.thresholds()`'s formula to whichever fluorophore the other
#'   axis actually shows - contributions from every other fluorophore in the
#'   panel are not visualisable on a 2D plot and are not included. `NULL`
#'   (default) draws no reference curves.
#' @param spread.kappa Numeric, spread standard deviations allowed above the
#'   flat threshold when `variants` is supplied - see `get.spread.thresholds()`.
#'   Default `2`.
#' @param x.lab An optional label for the x-axis. If none is given (default
#' `NULL`), the column name specified by `x.dim` will be used.
#' @param y.lab An optional label for the y-axis. If none is given (default
#' `NULL`), the column name specified by `y.dim` will be used.
#' @param x.min Minimum value for the x-axis. Default is `-5000`.
#' @param x.max Maximum value for the x-axis. Default is the value specified by
#' asp$expr.data.max, which will be the maximum for the cytometer.
#' @param y.min Minimum value for the y-axis. Default is `-5000`.
#' @param y.max Maximum value for the y-axis. Default is the value specified by
#' asp$expr.data.max, which will be the maximum for the cytometer.
#' @param x.width.basis Width basis for the biexponential transform for the
#' x-axis. Default is `-1000`.
#' @param y.width.basis Width basis for the biexponential transform for the
#' x-axis. Default is `-1000`.
#' @param max.points Number of points to plot (speeds up plotting). Default is
#' `5e6`.
#' @param color.palette Optional character string defining the viridis color
#' palette to be used for the fluorophore traces. Default is `rainbow`, which will
#' be similar to FlowJo or SpectroFlo. Other pptions are the viridis color
#' options: `magma`, `inferno`, `plasma`, `viridis`, `cividis`, `rocket`, `mako`
#' and `turbo`.
#' @param save Logical, if `TRUE`, saves a JPEG file to the `output.dir`.
#' Otherwise, the plot will simply be created in the Viewer.
#' @param title Optional title for the plot filename. If `NULL`, defaults to
#' `x.lab` vs. `y.lab`.
#' @param output.dir Optional output directory. Default is NULL, in which case
#' the current working directory will be used.
#' @param width Numeric, width of the saved plot. Default is `5`.
#' @param height Numeric, height of the saved plot. Default is `5`.
#'
#' @return Creates a biplot in the Viewer and optionally saves it as a JPEG file.
#'
#' @export

create.biplot <- function(
    plot.data,
    x.dim,
    y.dim,
    asp,
    variants = NULL,
    spread.kappa = 2,
    x.lab = NULL,
    y.lab = NULL,
    x.min = -5000,
    x.max = asp$expr.data.max,
    y.min = -5000,
    y.max = asp$expr.data.max,
    x.width.basis = -1000,
    y.width.basis = -1000,
    max.points = 5e5,
    color.palette = "rainbow",
    save = TRUE,
    title = NULL,
    output.dir = NULL,
    width = 5,
    height = 5
  ) {

  # check for x.dim, y.dim in colnames
  if ( !( x.dim %in% colnames( plot.data ) & y.dim %in% colnames( plot.data ) ) ) {
    print( colnames( plot.data ) )
    stop( "Either `xdim` or `y.dim` is not present in the data. See printed channels." )
  }

  # check inputs
  args <- list(
    x.min, x.max, y.min, y.max,
    x.width.basis, y.width.basis,
    max.points, width, height
  )

  arg.names <- c(
    "x.min", "x.max", "y.min", "y.max",
    "x.width.basis", "y.width.basis",
    "max.points", "width", "height"
  )

  for ( i in seq_along( args ) ) {
    if ( !is.numeric( args[[ i ]] ) ) {
      stop( paste( "Argument", arg.names[ i ], "must be numeric." ) )
    }
  }

  if ( is.null( x.lab ) )
    x.lab <- colnames( plot.data[ , x.dim, drop = FALSE ] )
  if ( is.null( y.lab ) )
    y.lab <- colnames( plot.data[ , y.dim, drop = FALSE ] )

  # downsample (faster plotting)
  if ( nrow( plot.data ) > max.points ) {
    # random sampling
    set.seed( asp$bird.seed )
    plot.data <- plot.data[ sample( seq_len( nrow( plot.data ) ), max.points ), ]
  }

  x.pos.log <- log10( x.max ) - 1
  y.pos.log <- log10( y.max ) - 1

  # set defaults
  if ( is.null( title ) )
    title <- paste( x.lab, "vs", y.lab )

  if ( is.null( output.dir ) )
    output.dir <- getwd()

  # If x.min/y.min (whether "auto"-computed or supplied directly) reaches
  # past whatever ribbon.breaks already covers on the negative side, add
  # whole-decade ticks (-1e4, -1e5, ...) down to the decade containing
  # x.min/y.min, so the axis isn't left with its deepest label sitting
  # well short of where the data (and the curve) actually extend.
  extend.negative.breaks <- function( breaks, min.value ) {

    if ( !is.finite( min.value ) || min.value >= 0 ) return( breaks )

    neg.breaks <- breaks[ breaks < 0 ]
    max.decade.present <- if ( length( neg.breaks ) == 0 ) 0 else
      floor( log10( max( abs( neg.breaks ) ) ) )

    needed.decade <- floor( log10( abs( min.value ) ) )

    if ( needed.decade <= max.decade.present ) return( breaks )

    new.decades <- seq( max.decade.present + 1, needed.decade )
    sort( c( breaks, -( 10 ^ new.decades ) ) )
  }

  # set plot limits
  x.breaks <- extend.negative.breaks( asp$ribbon.breaks, x.min )
  y.breaks <- extend.negative.breaks( asp$ribbon.breaks, y.min )
  x.breaks <- x.breaks[ x.breaks < x.max ]
  y.breaks <- y.breaks[ y.breaks < y.max ]
  x.axis.labels <- sapply( x.breaks, function( x ) {
    if ( x == 0 ) "0" else parse( text = paste0( "10^", log10( abs( x ) ) ) )
  } )
  y.axis.labels <- sapply( y.breaks, function( x ) {
    if ( x == 0 ) "0" else parse( text = paste0( "10^", log10( abs( x ) ) ) )
  } )
  x.limits <- c( x.min, x.max )
  y.limits <- c( y.min, y.max )

  # set transforms (one for x, one for y)
  biexp.trans.x <- biexp.transform(
    channelRange = asp$default.transformation.param$length,
    maxValue = x.max,
    pos = x.pos.log,
    neg = asp$default.transformation.param$neg,
    widthBasis = x.width.basis,
    inverse = FALSE )

  biexp.trans.y <- biexp.transform(
    channelRange = asp$default.transformation.param$length,
    maxValue = y.max,
    pos = y.pos.log,
    neg = asp$default.transformation.param$neg,
    widthBasis = y.width.basis,
    inverse = FALSE )

  # Spillover-spread reference curves. The boundary widens with the square
  # root of whichever fluorophore is spilling into the channel on the other
  # axis, so a flat line understates it at the bright end and overstates it
  # at the dim end - only a curve is honest about either. `abs()` on the
  # axis value keeps the curve real-valued across the biexponential display's
  # negative decades; `get.spread.thresholds()` itself is unaffected and
  # still clips negative abundance to zero contribution.
  spread.curve.data <- NULL

  if ( !is.null( variants ) ) {

    get.spillover.spread <- function( source, target ) {
      ss <- variants$spillover.spread
      if ( is.null( ss ) || !( source %in% rownames( ss ) ) ||
           !( target %in% colnames( ss ) ) ) return( 0 )
      val <- ss[ source, target ]
      if ( !is.finite( val ) || val < 0 ) val <- 0
      val
    }

    ss.xy <- get.spillover.spread( x.dim, y.dim )
    ss.yx <- get.spillover.spread( y.dim, x.dim )

    flat.y <- variants$thresholds[ y.dim ]
    flat.x <- variants$thresholds[ x.dim ]

    # Directly measured negative-tail flat component (see
    # `get.spectral.variants()$neg.thresholds`), falling back to the mirrored
    # positive threshold with a warning if an older cached `variants` object
    # does not carry it.
    flat.y.neg <- variants$neg.thresholds[ y.dim ]
    if ( length( flat.y.neg ) != 1 || !is.finite( flat.y.neg ) ) {
      flat.y.neg <- -flat.y
      if ( length( flat.y ) == 1 && is.finite( flat.y ) )
        warning( paste0( "No `neg.thresholds` in `variants` for ", y.dim,
                         "; mirroring the positive threshold about zero." ),
                 call. = FALSE )
    }

    flat.x.neg <- variants$neg.thresholds[ x.dim ]
    if ( length( flat.x.neg ) != 1 || !is.finite( flat.x.neg ) ) {
      flat.x.neg <- -flat.x
      if ( length( flat.x ) == 1 && is.finite( flat.x ) )
        warning( paste0( "No `neg.thresholds` in `variants` for ", x.dim,
                         "; mirroring the positive threshold about zero." ),
                 call. = FALSE )
    }

    x.seq <- seq( x.min, x.max, length.out = 300 )
    y.seq <- seq( y.min, y.max, length.out = 300 )

    spread.curve.data <- list()

    if ( length( flat.y ) == 1 && is.finite( flat.y ) ) {

      y.pos <- .spread.boundary( x.seq, flat.y, ss.xy, spread.kappa, 1, y.min, y.max )

      if ( length( y.pos$along ) > 1 )
        spread.curve.data$y.pos <- data.frame(
          x = biexp.trans.x( y.pos$along ), y = biexp.trans.y( y.pos$value ) )
    }

    if ( length( flat.y.neg ) == 1 && is.finite( flat.y.neg ) ) {

      y.neg <- .spread.boundary( x.seq, flat.y.neg, ss.xy, spread.kappa, -1, y.min, y.max )

      if ( length( y.neg$along ) > 1 )
        spread.curve.data$y.neg <- data.frame(
          x = biexp.trans.x( y.neg$along ), y = biexp.trans.y( y.neg$value ) )
    }

    if ( length( flat.x ) == 1 && is.finite( flat.x ) ) {

      x.pos <- .spread.boundary( y.seq, flat.x, ss.yx, spread.kappa, 1, x.min, x.max )

      if ( length( x.pos$along ) > 1 )
        spread.curve.data$x.pos <- data.frame(
          x = biexp.trans.x( x.pos$value ), y = biexp.trans.y( x.pos$along ) )
    }

    if ( length( flat.x.neg ) == 1 && is.finite( flat.x.neg ) ) {

      x.neg <- .spread.boundary( y.seq, flat.x.neg, ss.yx, spread.kappa, -1, x.min, x.max )

      if ( length( x.neg$along ) > 1 )
        spread.curve.data$x.neg <- data.frame(
          x = biexp.trans.x( x.neg$value ), y = biexp.trans.y( x.neg$along ) )
    }
  }

  # convert to data frame for plotting
  plot.data <- data.frame(
    x = plot.data[ , x.dim ],
    y = plot.data[ , y.dim ] )

  # apply transformation
  plot.data$x.trans <- biexp.trans.x( plot.data$x )
  plot.data$y.trans <- biexp.trans.y( plot.data$y )

  # set up the plot
  biplot <- ggplot( plot.data, aes( x.trans, y.trans ) ) +
    geom_scattermore(
      pointsize = asp$figure.gate.point.size,
      color = "black",
      alpha = 1,
      na.rm = TRUE
    ) +
    stat_density_2d(
      aes( fill = after_stat( level ) ),
      geom = "polygon",
      na.rm = TRUE
    ) +
    scale_x_continuous(
      name = x.lab,
      breaks = biexp.trans.x( x.breaks ),
      limits = biexp.trans.x( x.limits ),
      labels = x.axis.labels
    ) +
    scale_y_continuous(
      name = y.lab,
      breaks = biexp.trans.y( y.breaks ),
      limits = biexp.trans.y( y.limits ),
      labels = y.axis.labels
    ) +
    theme_bw() +
    theme(
      plot.margin = margin(
        asp$figure.margin, asp$figure.margin, asp$figure.margin, asp$figure.margin
      ),
      legend.position = "none",
      axis.ticks = element_line( linewidth = asp$figure.panel.line.size ),
      axis.text = element_text( size = asp$figure.axis.text.size ),
      axis.title = element_text( size = asp$figure.axis.title.size ),
      panel.border = element_rect( fill = NA, linewidth = asp$figure.panel.line.size ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

  if ( !is.null( spread.curve.data ) ) {
    if ( !is.null( spread.curve.data$y.pos ) )
      biplot <- biplot +
        geom_line( data = spread.curve.data$y.pos, aes( x = x, y = y ),
                   color = "red", inherit.aes = FALSE, na.rm = TRUE )
    if ( !is.null( spread.curve.data$y.neg ) )
      biplot <- biplot +
        geom_line( data = spread.curve.data$y.neg, aes( x = x, y = y ),
                   color = "red", inherit.aes = FALSE, na.rm = TRUE )
    # x-direction curves are functions of |y|: draw them in data (y) order
    # with geom_path(); geom_line() would sort by x and interleave the
    # positive- and negative-y branches into a zigzag
    if ( !is.null( spread.curve.data$x.pos ) )
      biplot <- biplot +
        geom_path( data = spread.curve.data$x.pos, aes( x = x, y = y ),
                   color = "red", inherit.aes = FALSE, na.rm = TRUE )
    if ( !is.null( spread.curve.data$x.neg ) )
      biplot <- biplot +
        geom_path( data = spread.curve.data$x.neg, aes( x = x, y = y ),
                   color = "red", inherit.aes = FALSE, na.rm = TRUE )
  }

  # color options
  viridis.colors <- c(
    "magma", "inferno", "plasma", "viridis",
    "cividis", "rocket", "mako", "turbo"
  )

  # add fill layer for color palette
  if ( color.palette %in% viridis.colors ) {
    biplot <- biplot + scale_fill_viridis_c( option = color.palette )
  } else {
    biplot <- biplot +
      scale_fill_gradientn( colors = asp$density.palette.base.color )
  }

  # save or return the plot
  if ( save )
    ggsave(
      file.path( output.dir, sprintf( "%s.jpg", title ) ),
      plot = biplot,
      device = ragg::agg_jpeg,
      width = width,
      height = height,
      limitsize = FALSE
    )
  else
    return( biplot )

  print( biplot )
}

# Spillover-spread boundary `flat + sign * kappa * sqrt( ss * |along| )`,
# sampled along the other axis (`along`). The exact `along` value at which
# the curve reaches the axis limit it is heading towards is added to the
# samples, and points beyond that limit are dropped, so the drawn curve stops
# at the plot edge instead of running along it. The curve is monotone in
# |along|, so the kept points are one contiguous run. No NA is returned:
# biexp.transform() cannot take NA (its `if ( residual > 1e-6 )` check fails).
#
# @return list( along, value ), same length, possibly 0.
# @noRd
.spread.boundary <- function( along, flat, ss, kappa, sign, lim.min, lim.max ) {

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
