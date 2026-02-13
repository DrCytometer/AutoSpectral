# gate_af_sample_plot.r

#' @title Plot Autofluorescence Gates on Samples
#'
#' @description
#' This function plots the autofluorescence exclusion gate on the sample(s).
#'
#' @importFrom ggplot2 ggplot aes scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 theme_bw theme element_line geom_path after_stat
#' @importFrom ggplot2 element_text element_rect margin ggsave
#' @importFrom ggplot2 scale_fill_viridis_c scale_fill_gradientn stat_density_2d
#' @importFrom scattermore geom_scattermore
#' @importFrom ragg agg_jpeg
#'
#' @param plot.data Matrix containing autofluorescence data points.
#' @param samp Sample identifier.
#' @param af.boundary.upper Matrix containing upper boundary information.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#' @param max.points Number of points to plot (speeds up plotting). Default is
#' `5e4`.
#' @param color.palette Optional character string defining the viridis color
#' palette to be used for the fluorophore traces. Default is `viridis`. Options
#' are the viridis color options: `magma`, `inferno`, `plasma`, `viridis`,
#' `cividis`, `rocket`, `mako` and `turbo`.
#'
#' @return Saves the plot as a JPEG file in the specified directory.

gate.af.sample.plot <- function(
    plot.data,
    samp,
    af.boundary.upper,
    asp,
    max.points = 5e4,
    color.palette = "viridis"
  ) {

  valid.idx <- which( !is.na( plot.data[ , 1 ] ) & !is.na( plot.data[ , 2 ] ) )
  if ( length( valid.idx ) == 0 ) {
    warning( "AF plot.data has no valid x/y values; skipping density plot." )
    return( invisible( NULL ) )
  }
  plot.data <- plot.data[ valid.idx, ]

  # downsample (faster plotting)
  if ( nrow( plot.data ) > max.points ) {
    # random sampling
    set.seed( 42 )
    plot.data <- plot.data[ sample( seq_len( nrow( plot.data ) ), max.points ), ]
  }

  # get axis labels
  x.lab <- colnames( plot.data )[ 1 ]
  y.lab <- colnames( plot.data )[ 2 ]

  # convert to data frame for plotting
  plot.data <- data.frame(
    x = plot.data[ , 1 ],
    y = plot.data[ , 2 ] )

  # set transform
  biexp.transform <- flowWorkspace::flowjo_biexp(
    channelRange = asp$default.transformation.param$length,
    maxValue = asp$default.transformation.param$max.range,
    pos = asp$default.transformation.param$pos,
    neg = asp$default.transformation.param$neg,
    widthBasis = asp$default.transformation.param$width,
    inverse = FALSE )

  plot.data.ggp <- data.frame(
    x.trans = biexp.transform( plot.data[ , 1] ),
    y.trans = biexp.transform( plot.data[ , 2] )
  )

  # set plot limits
  breaks <- asp$ribbon.breaks
  limits <- c( asp$ribbon.plot.min, asp$expr.data.max )
  axis.labels <- sapply( breaks, function( x ) {
    if ( x == 0 ) "0" else parse( text = paste0( "10^", log10( abs( x ) ) ) )
  } )

  # set up the base plot
  gate.plot <- ggplot( plot.data.ggp, aes( x.trans, y.trans ) ) +
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
      breaks = biexp.transform( breaks ),
      limits = biexp.transform( limits ),
      labels = axis.labels
    ) +
    scale_y_continuous(
      name = y.lab,
      breaks = biexp.transform( breaks ),
      limits = biexp.transform( limits ),
      labels = axis.labels
    ) +
    theme_bw() +
    theme(
      plot.margin = margin(
        asp$figure.margin, asp$figure.margin, asp$figure.margin, asp$figure.margin
      ),
      legend.position = "none",
      axis.ticks = element_line( linewidth = asp$figure.panel.line.size ),
      axis.text = element_text( size = asp$figure.axis.text.size ),
      axis.text.x = element_text( angle = 45, hjust = 1 ),
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
    gate.plot <- gate.plot + scale_fill_viridis_c( option = color.palette )
  } else {
    gate.plot <- gate.plot +
      scale_fill_gradientn( colors = asp$density.palette.base.color )
  }

  # add AF gate boundary
  if ( !is.null( af.boundary.upper ) ) {
    af.boundary.upper.ggp <- data.frame(
      x = c( af.boundary.upper$x,
             af.boundary.upper$x[ 1 ] ),
      y = c( af.boundary.upper$y,
             af.boundary.upper$y[ 1 ] )
    )

    af.boundary.upper.ggp$x <- biexp.transform( af.boundary.upper.ggp$x )
    af.boundary.upper.ggp$y <- biexp.transform( af.boundary.upper.ggp$y )

    gate.plot <- gate.plot +
      geom_path(
        data = af.boundary.upper.ggp,
        aes( x, y ),
        color = "black",
        linewidth = asp$figure.gate.line.size
      )
  }

  # save the final plot
  ggsave(
    file.path(
      asp$figure.clean.control.dir,
      paste( asp$af.plot.filename, samp, ".jpg", sep = "_" ) ),
    plot = gate.plot,
    device = ragg::agg_jpeg,
    width = asp$figure.width,
    height = asp$figure.height,
    limitsize = FALSE
  )

}
