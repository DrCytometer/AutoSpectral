# create_biplot.r

#' @title Create Biplot
#'
#' @importFrom ggplot2 ggplot aes ggsave after_stat
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous theme_bw theme
#' @importFrom ggplot2 margin element_line element_text element_rect element_blank
#' @importFrom ggplot2 scale_fill_viridis_c scale_fill_gradientn stat_density_2d
#' @importFrom flowWorkspace flowjo_biexp
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
    x.lab = NULL,
    y.lab = NULL,
    x.min = -5000,
    x.max = asp$expr.data.max,
    y.min = -5000,
    y.max = asp$expr.data.max,
    x.width.basis = -1000,
    y.width.basis = -1000,
    max.points = 5e6,
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
    set.seed( 42 )
    plot.data <- plot.data[ sample( seq_len( nrow( plot.data ) ), max.points ), ]
  }

  ### The FlowJo Biexponential transformation caps out at -1000 for the width
  # parameter in both FlowJo and in the R implementation in flowWorkspace. This
  # is inadequate for a lot of spectral flow data. As a crude work-around, we can
  # use any "excess" width basis to feed into modifying the number of positive log
  # decades used in the transformation, which has a (sort of) similar effect.
  # At some point, this needs to be implemented in C or C++ with a full range of
  # width bases, but that's beyond me at the moment.
  if ( x.width.basis < -1000 ) {
    x.excess.width.basis <- x.width.basis + 1000
    x.pos.log.delta <- log10( abs( x.excess.width.basis ) )
    x.pos.log <- log10( x.max ) - 1 - x.pos.log.delta
    x.pos.log <- pmax( x.pos.log, 2 )
  } else {
    x.pos.log <- log10( x.max ) - 1
  }

  if ( y.width.basis < -1000 ) {
    y.excess.width.basis <- y.width.basis + 1000
    y.pos.log.delta <- log10( abs( y.excess.width.basis ) )
    y.pos.log <- log10( y.max ) - 1 - y.pos.log.delta
    y.pos.log <- pmax( y.pos.log, 2 )
  } else {
    y.pos.log <- log10( y.max ) - 1
  }

  # set defaults
  if ( is.null( title ) )
    title <- paste( x.lab, "vs", y.lab )

  if ( is.null( output.dir ) )
    output.dir <- getwd()

  # set plot limits
  x.breaks <- asp$ribbon.breaks[ asp$ribbon.breaks < x.max ]
  y.breaks <- asp$ribbon.breaks[ asp$ribbon.breaks < y.max ]
  x.axis.labels <- sapply( x.breaks, function( x ) {
    if ( x == 0 ) "0" else parse( text = paste0( "10^", log10( abs( x ) ) ) )
  } )
  y.axis.labels <- sapply( y.breaks, function( x ) {
    if ( x == 0 ) "0" else parse( text = paste0( "10^", log10( abs( x ) ) ) )
  } )
  x.limits <- c( x.min, x.max )
  y.limits <- c( y.min, y.max )

  # set transforms (one for x, one for y)
  biexp.transform.x <- flowWorkspace::flowjo_biexp(
    channelRange = asp$default.transformation.param$length,
    maxValue = x.max,
    pos = x.pos.log,
    neg = asp$default.transformation.param$neg,
    widthBasis = x.width.basis,
    inverse = FALSE )

  biexp.transform.y <- flowWorkspace::flowjo_biexp(
    channelRange = asp$default.transformation.param$length,
    maxValue = y.max,
    pos = y.pos.log,
    neg = asp$default.transformation.param$neg,
    widthBasis = y.width.basis,
    inverse = FALSE )

  # convert to data frame for plotting
  plot.data <- data.frame(
    x = plot.data[ , x.dim ],
    y = plot.data[ , y.dim ] )

  # apply transformation
  plot.data$x.trans <- biexp.transform.x( plot.data$x )
  plot.data$y.trans <- biexp.transform.y( plot.data$y )

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
      breaks = biexp.transform.x( x.breaks ),
      limits = biexp.transform.x( x.limits ),
      labels = x.axis.labels
    ) +
    scale_y_continuous(
      name = y.lab,
      breaks = biexp.transform.y( y.breaks ),
      limits = biexp.transform.y( y.limits ),
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


