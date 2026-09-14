# create_biplot.r

#' @title Create Biplot
#'
#' @description
#' Creates a biexponential-transformed 2D density biplot for two channels
#' of flow cytometry data. Every argument beyond `plot.data` is optional:
#' if `x.dim`/`y.dim` are omitted, the first two columns of `plot.data` are
#' used; if `asp` is omitted, every cytometer-specific setting it would
#' otherwise supply falls back to a literal default matching the Aurora
#' cytometer profile (`get.autospectral.param.aurora()`). This means
#' `create.biplot( plot.data = flow.data[ , c( 3, 6 ) ] )` works standalone.
#'
#' @details
#' For every argument documented below as "falls back to `asp$...`", the
#' precedence is: the explicit argument (if supplied) beats `asp`'s value
#' (if `asp` is supplied) beats the literal default. `asp`, when supplied,
#' is therefore a convenience for setting several of these at once for a
#' non-Aurora cytometer -- it does not override an explicitly-passed
#' argument.
#'
#' `x.min`/`y.min` and `x.width.basis`/`y.width.basis` each accept the
#' string `"auto"` (the default) in place of a number. In auto mode, the
#' width basis is calculated per-channel from the plotted data itself,
#' following the data-driven approach of Parks, Roederer & Moore (2006):
#' \deqn{w = \frac{M - \log_{10}(T / |r|)}{2}}
#' where \code{M} is \code{pos}, \code{T} is \code{x.max}/\code{y.max},
#' and \code{r} is the \code{logicle.q}-th quantile (default the 5th
#' percentile) of that channel's own *negative-valued* raw data only --
#' "the most negative value to be included in the display" in the
#' terminology of the paper. Restricting the quantile to values below
#' zero (rather than the whole channel) is deliberate: a channel where
#' fewer than `logicle.q` of events are negative would otherwise return
#' a positive `r`, which the formula does not define sensibly. This `w`
#' is converted back to this package's `widthBasis` convention via
#' `widthBasis = -10^(2w)` (see `biexp.transform()`). In lay terms: the
#' transform looks at how far the negative tail of this specific channel
#' actually extends and picks just enough curvature near zero to show
#' it, rather than using one fixed guess for every channel. `x.min` and
#' `x.width.basis` resolve to the *same* value when either is `"auto"` --
#' auto mode does not give the axis floor separate headroom beyond the
#' width basis. `logicle.w.min` (default `0.5`, matching
#' `biexp.transform()`'s own default `widthBasis = -10`) is the floor
#' applied when a channel has no visible negative population to measure.
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
#' `y.dim`. Must have at least two columns.
#' @param x.dim String specifying the column of `plot.data` for the x-axis of
#' the plot. Default `NULL`, in which case the first column of `plot.data`
#' is used.
#' @param y.dim String specifying the column of `plot.data` for the y-axis of
#' the plot. Default `NULL`, in which case the second column of `plot.data`
#' is used.
#' @param asp Optional AutoSpectral parameter list. When supplied, its
#' fields are used for any of `x.max`, `y.max`, `pos`, `channel.range`,
#' `neg`, `bird.seed`, `ribbon.breaks`, `figure.gate.point.size`,
#' `figure.margin`, `figure.panel.line.size`, `figure.axis.text.size`,
#' `figure.axis.title.size` and `density.palette.base.color` that were not
#' passed explicitly. Default `NULL` -- see Details for the literal
#' (Aurora-profile) fallback used for each.
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
#' @param x.min Minimum value for the x-axis, or `"auto"` (default) to
#' calculate it from the data -- see Details.
#' @param x.max Maximum value for the x-axis. Default `NULL`, falling back
#' to `asp$expr.data.max`, then to `4194304` (Aurora).
#' @param y.min Minimum value for the y-axis, or `"auto"` (default) to
#' calculate it from the data -- see Details.
#' @param y.max Maximum value for the y-axis. Default `NULL`, falling back
#' to `asp$expr.data.max`, then to `4194304` (Aurora).
#' @param x.width.basis Width basis for the biexponential transform for the
#' x-axis, or `"auto"` (default) to calculate it from the data -- see
#' Details.
#' @param y.width.basis Width basis for the biexponential transform for the
#' y-axis, or `"auto"` (default) to calculate it from the data -- see
#' Details.
#' @param pos Number of positive decades (`M`) spanned by the
#' biexponential transform, applied to both axes. Default `NULL`,
#' falling back to `asp$default.transformation.param$pos` (the
#' canonical, package-wide value), then to `log10(x.max) - 1` only when
#' `asp` itself is not supplied.
#' @param channel.range Display channel range (`channelRange`) passed to
#' `biexp.transform()`. Default `NULL`, falling back to
#' `asp$default.transformation.param$length`, then to `256`.
#' @param neg Must be `0` -- see `biexp.transform()`. Default `NULL`,
#' falling back to `asp$default.transformation.param$neg`, then to `0`.
#' @param logicle.q Quantile of each channel's raw data used as `r` (the
#' most negative value to bring into the display) in the `"auto"`
#' width-basis calculation -- see Details. Default `0.05`, matching
#' Parks, Roederer & Moore (2006) / `flowCore::estimateLogicle()`.
#' @param logicle.w.min Floor, in logicle decades, on the `"auto"`-computed
#' width when a channel has no visible negative population -- see Details.
#' Default `0.5`.
#' @param max.points Number of points to plot (speeds up plotting). Default is
#' `5e5`.
#' @param color.palette Optional character string defining the viridis color
#' palette to be used for the fluorophore traces. Default is `rainbow`, which will
#' be similar to FlowJo or SpectroFlo. Other pptions are the viridis color
#' options: `magma`, `inferno`, `plasma`, `viridis`, `cividis`, `rocket`, `mako`
#' and `turbo`.
#' @param bird.seed Integer random seed used when downsampling to
#' `max.points`. Default `NULL`, falling back to `asp$bird.seed`, then to
#' the package's usual seed.
#' @param ribbon.breaks Numeric vector of raw-scale axis break positions.
#' Default `NULL`, falling back to `asp$ribbon.breaks`, then to
#' `c( -1e3, 0, 1e3, 1e4, 1e5, 1e6 )` (Aurora).
#' @param figure.gate.point.size Point size for the scattermore layer.
#' Default `NULL`, falling back to `asp$figure.gate.point.size`, then to
#' `0.8`.
#' @param figure.margin Plot margin (all four sides). Default `NULL`,
#' falling back to `asp$figure.margin`, then to `4.0`.
#' @param figure.panel.line.size Line width for axis ticks and the panel
#' border. Default `NULL`, falling back to `asp$figure.panel.line.size`,
#' then to `0.5`.
#' @param figure.axis.text.size Axis text size. Default `NULL`, falling
#' back to `asp$figure.axis.text.size`, then to `12.0`.
#' @param figure.axis.title.size Axis title size. Default `NULL`, falling
#' back to `asp$figure.axis.title.size`, then to `12.0`.
#' @param density.palette.base.color Fallback fill gradient colors used
#' when `color.palette` is not one of the viridis options. Default `NULL`,
#' falling back to `asp$density.palette.base.color`, then to
#' `c( "blue", "cyan", "green", "yellow", "red" )`.
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
#' @references
#' Parks DR, Roederer M, Moore WA (2006). A new "Logicle" display method
#' avoids deceptive effects of logarithmic scaling for low signals and
#' compensated data. Cytometry A, 69(6):541-551.
#'
#' @export

create.biplot.test <- function(
    plot.data,
    x.dim = NULL,
    y.dim = NULL,
    asp = NULL,
    variants = NULL,
    spread.kappa = 2,
    x.lab = NULL,
    y.lab = NULL,
    x.min = "auto",
    x.max = NULL,
    y.min = "auto",
    y.max = NULL,
    x.width.basis = "auto",
    y.width.basis = "auto",
    pos = NULL,
    channel.range = NULL,
    neg = NULL,
    logicle.q = 0.05,
    logicle.w.min = 0.5,
    max.points = 5e5,
    color.palette = "rainbow",
    bird.seed = NULL,
    ribbon.breaks = NULL,
    figure.gate.point.size = NULL,
    figure.margin = NULL,
    figure.panel.line.size = NULL,
    figure.axis.text.size = NULL,
    figure.axis.title.size = NULL,
    density.palette.base.color = NULL,
    save = TRUE,
    title = NULL,
    output.dir = NULL,
    width = 5,
    height = 5
) {

  # fall back to the first two columns of plot.data when x.dim/y.dim are
  # not supplied, so create.biplot( plot.data ) works on its own
  if ( is.null( x.dim ) || is.null( y.dim ) ) {
    if ( ncol( plot.data ) < 2 ) {
      stop( "`plot.data` must have at least two columns when `x.dim`/`y.dim` are not supplied." )
    }
    if ( is.null( x.dim ) ) x.dim <- colnames( plot.data )[ 1 ]
    if ( is.null( y.dim ) ) y.dim <- colnames( plot.data )[ 2 ]
  }

  # check for x.dim, y.dim in colnames
  if ( !( x.dim %in% colnames( plot.data ) & y.dim %in% colnames( plot.data ) ) ) {
    print( colnames( plot.data ) )
    stop( "Either `xdim` or `y.dim` is not present in the data. See printed channels." )
  }

  # explicit argument > asp field (if asp supplied) > literal Aurora-profile
  # default, so create.biplot() works with no `asp` at all while `asp`
  # still overrides the literal defaults for other cytometers
  resolve.param <- function( val, asp.val, default ) {
    if ( !is.null( val ) ) return( val )
    if ( !is.null( asp.val ) ) return( asp.val )
    default
  }

  x.max <- resolve.param( x.max, asp$expr.data.max, 4194304 )
  y.max <- resolve.param( y.max, asp$expr.data.max, 4194304 )
  pos <- resolve.param( pos, asp$default.transformation.param$pos, log10( x.max ) )
  channel.range <- resolve.param( channel.range, asp$default.transformation.param$length, 256 )
  neg <- resolve.param( neg, asp$default.transformation.param$neg, 0 )
  bird.seed <- resolve.param(
    bird.seed, asp$bird.seed,
    as.integer( prod( which( letters %in% strsplit( "hummingbird", "" )[[ 1 ]] ) ) )
  )
  ribbon.breaks <- resolve.param( ribbon.breaks, asp$ribbon.breaks, c( -1e3, 0, 1e3, 1e4, 1e5, 1e6 ) )
  figure.gate.point.size <- resolve.param( figure.gate.point.size, asp$figure.gate.point.size, 0.8 )
  figure.margin <- resolve.param( figure.margin, asp$figure.margin, 4.0 )
  figure.panel.line.size <- resolve.param( figure.panel.line.size, asp$figure.panel.line.size, 0.5 )
  figure.axis.text.size <- resolve.param( figure.axis.text.size, asp$figure.axis.text.size, 12.0 )
  figure.axis.title.size <- resolve.param( figure.axis.title.size, asp$figure.axis.title.size, 12.0 )
  density.palette.base.color <- resolve.param(
    density.palette.base.color, asp$density.palette.base.color,
    c( "blue", "cyan", "green", "yellow", "red" )
  )

  # check inputs that must be plain numbers
  numeric.args <- list(
    x.max = x.max, y.max = y.max,
    pos = pos, channel.range = channel.range, neg = neg,
    bird.seed = bird.seed, ribbon.breaks = ribbon.breaks,
    figure.gate.point.size = figure.gate.point.size,
    figure.margin = figure.margin,
    figure.panel.line.size = figure.panel.line.size,
    figure.axis.text.size = figure.axis.text.size,
    figure.axis.title.size = figure.axis.title.size,
    max.points = max.points, width = width, height = height
  )

  for ( nm in names( numeric.args ) ) {
    if ( !is.numeric( numeric.args[[ nm ]] ) ) {
      stop( paste( "Argument", nm, "must be numeric." ) )
    }
  }

  # x.min/y.min/x.width.basis/y.width.basis may additionally be "auto"
  auto.args <- list(
    x.min = x.min, y.min = y.min,
    x.width.basis = x.width.basis, y.width.basis = y.width.basis
  )

  for ( nm in names( auto.args ) ) {
    val <- auto.args[[ nm ]]
    if ( !( is.numeric( val ) || identical( val, "auto" ) ) ) {
      stop( paste( "Argument", nm, "must be numeric or \"auto\"." ) )
    }
  }

  if ( is.null( x.lab ) )
    x.lab <- colnames( plot.data[ , x.dim, drop = FALSE ] )
  if ( is.null( y.lab ) )
    y.lab <- colnames( plot.data[ , y.dim, drop = FALSE ] )

  # downsample (faster plotting)
  if ( nrow( plot.data ) > max.points ) {
    # random sampling
    set.seed( bird.seed )
    plot.data <- plot.data[ sample( seq_len( nrow( plot.data ) ), max.points ), ]
  }

  # Data-driven width basis (Parks, Roederer & Moore 2006)
  get.auto.width.basis <- function( x, pos, max.value, q, w.min ) {

    x.neg <- x[ is.finite( x ) & x < 0 ]
    if ( length( x.neg ) == 0 ) {
      return( list( r = -1, width.basis = -10 ^ ( 2 * w.min ) ) )
    }

    r <- stats::quantile( x.neg, probs = q, names = FALSE, na.rm = TRUE )
    r <- max( min( r, -1 ), -max.value )

    w <- ( pos - log10( max.value / abs( r ) ) ) / 2
    w <- min( max( w, w.min ), pos / 2 )

    list( r = r, width.basis = -10 ^ ( 2 * w ) )
  }

  if ( identical( x.width.basis, "auto" ) || identical( x.min, "auto" ) ) {
    x.auto <- get.auto.width.basis(
      x = plot.data[ , x.dim ], pos = pos, max.value = x.max,
      q = logicle.q, w.min = logicle.w.min
    )
    if ( identical( x.width.basis, "auto" ) ) x.width.basis <- x.auto$width.basis
    if ( identical( x.min, "auto" ) ) x.min <- x.width.basis * 2
  }

  if ( identical( y.width.basis, "auto" ) || identical( y.min, "auto" ) ) {
    y.auto <- get.auto.width.basis(
      x = plot.data[ , y.dim ], pos = pos, max.value = y.max,
      q = logicle.q, w.min = logicle.w.min
    )
    if ( identical( y.width.basis, "auto" ) ) y.width.basis <- y.auto$width.basis
    if ( identical( y.min, "auto" ) ) y.min <- y.width.basis * 2
  }

  # set defaults
  if ( is.null( title ) ) title <- paste( x.lab, "vs", y.lab )

  if ( is.null( output.dir ) ) output.dir <- getwd()

  # add more negative axes ticks if needed
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
  x.breaks <- extend.negative.breaks( ribbon.breaks, x.min )
  y.breaks <- extend.negative.breaks( ribbon.breaks, y.min )
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
    channelRange = channel.range,
    maxValue = x.max,
    pos = pos - 1,
    neg = neg,
    widthBasis = x.width.basis,
    inverse = FALSE )

  biexp.trans.y <- biexp.transform(
    channelRange = channel.range,
    maxValue = y.max,
    pos = pos - 1,
    neg = neg,
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

      y.pos.raw <- flat.y + spread.kappa * sqrt( ss.xy * abs( x.seq ) )
      y.pos.raw <- pmin( pmax( y.pos.raw, y.min ), y.max )

      spread.curve.data$y.pos <- data.frame(
        x = biexp.trans.x( x.seq ), y = biexp.trans.y( y.pos.raw ) )
    }

    if ( length( flat.y.neg ) == 1 && is.finite( flat.y.neg ) ) {

      y.neg.raw <- flat.y.neg - spread.kappa * sqrt( ss.xy * abs( x.seq ) )
      y.neg.raw <- pmin( pmax( y.neg.raw, y.min ), y.max )

      spread.curve.data$y.neg <- data.frame(
        x = biexp.trans.x( x.seq ), y = biexp.trans.y( y.neg.raw ) )
    }

    if ( length( flat.x ) == 1 && is.finite( flat.x ) ) {

      x.pos.raw <- flat.x + spread.kappa * sqrt( ss.yx * abs( y.seq ) )
      x.pos.raw <- pmin( pmax( x.pos.raw, x.min ), x.max )

      spread.curve.data$x.pos <- data.frame(
        x = biexp.trans.x( x.pos.raw ), y = biexp.trans.y( y.seq ) )
    }

    if ( length( flat.x.neg ) == 1 && is.finite( flat.x.neg ) ) {

      x.neg.raw <- flat.x.neg - spread.kappa * sqrt( ss.yx * abs( y.seq ) )
      x.neg.raw <- pmin( pmax( x.neg.raw, x.min ), x.max )

      spread.curve.data$x.neg <- data.frame(
        x = biexp.trans.x( x.neg.raw ), y = biexp.trans.y( y.seq ) )
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
      pointsize = figure.gate.point.size,
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
        figure.margin, figure.margin, figure.margin, figure.margin
      ),
      legend.position = "none",
      axis.ticks = element_line( linewidth = figure.panel.line.size ),
      axis.text = element_text( size = figure.axis.text.size ),
      axis.title = element_text( size = figure.axis.title.size ),
      panel.border = element_rect( fill = NA, linewidth = figure.panel.line.size ),
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
    if ( !is.null( spread.curve.data$x.pos ) )
      biplot <- biplot +
        geom_line( data = spread.curve.data$x.pos, aes( x = x, y = y ),
                   color = "red", inherit.aes = FALSE, na.rm = TRUE )
    if ( !is.null( spread.curve.data$x.neg ) )
      biplot <- biplot +
        geom_line( data = spread.curve.data$x.neg, aes( x = x, y = y ),
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
      scale_fill_gradientn( colors = density.palette.base.color )
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
