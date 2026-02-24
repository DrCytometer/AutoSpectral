# gate_af_sample_plot.r

#' @title Plot Autofluorescence Gates on Samples
#'
#' @description
#' This function plots the autofluorescence exclusion gate on the sample(s).
#'
#' @importFrom ggplot2 ggplot aes scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 theme_bw theme element_line geom_path after_stat
#' @importFrom ggplot2 element_text element_rect margin ggsave
#' @importFrom ggplot2 scale_fill_viridis_d scale_fill_manual coord_cartesian
#' @importFrom scattermore geom_scattermore
#' @importFrom ragg agg_jpeg
#' @importFrom MASS kde2d
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
#' `cividis`, `rocket`, `mako` and `turbo`. Use `rainbow`to be similar to FlowJo
#' or SpectroFlo.
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

  # get density for plotting
  bw <- apply( plot.data.ggp, 2, bandwidth.nrd )

  if ( requireNamespace("AutoSpectralRcpp", quietly = TRUE ) &&
       "fast_kde2d_cpp" %in% ls( getNamespace( "AutoSpectralRcpp" ) ) &&
       nrow( plot.data.ggp ) > 10000 ) {
    # use C++ function to get density
    gate.bound.density <- AutoSpectralRcpp::fast_kde2d_cpp(
      x = plot.data.ggp[ , 1 ],
      y = plot.data.ggp[ , 2 ],
      n = 100,
      h = bw * 0.1,
      x_limits = range( plot.data.ggp[ , 1 ] ),
      y_limits = range( plot.data.ggp[ , 2 ] )
    )
  } else {
    # use slower MASS call
    gate.bound.density <- MASS::kde2d(
      x = plot.data.ggp[ , 1 ],
      y = plot.data.ggp[ , 2 ],
      n = 60,
      h = bw * 0.8
    )
  }

  # format the density for plotting
  density.df <- expand.grid(
    x = gate.bound.density$x,
    y = gate.bound.density$y
  )
  density.df$z <- as.vector( gate.bound.density$z )
  density.df$z[ is.na( density.df$z ) ] <- 0
  max.z <- max( density.df$z, na.rm = TRUE )
  density.breaks <- seq( 0.05 * max.z, max.z, length.out = 11 )
  if( diff( range( density.breaks ) ) == 0) {
    density.breaks <- seq( 0, max.z + 0.1, length.out = 11)
  }

  # set up the base plot
  gate.plot <- ggplot( plot.data.ggp, aes( x.trans, y.trans ) ) +
    geom_scattermore(
      pointsize = asp$figure.gate.point.size,
      color = "black",
      alpha = 1,
      na.rm = TRUE
    ) +
    geom_contour_filled(
      data = density.df,
      aes( x = x, y = y, z = z ),
      breaks = density.breaks,
      alpha = 1,
      inherit.aes = FALSE,
      na.rm = TRUE
    ) +
    scale_x_continuous(
      name = x.lab,
      breaks = biexp.transform( breaks ),
      #limits = biexp.transform( limits ),
      labels = axis.labels
    ) +
    scale_y_continuous(
      name = y.lab,
      breaks = biexp.transform( breaks ),
      #limits = biexp.transform( limits ),
      labels = axis.labels
    ) +
    coord_cartesian(
      xlim = biexp.transform( limits ),
      ylim = biexp.transform( limits )
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
    gate.plot <- gate.plot + scale_fill_viridis_d( option = color.palette )
  } else {
    n.bins <- max( 1, length( density.breaks ) - 1 )
    rainbow.palette <- grDevices::colorRampPalette( asp$density.palette.base.color )( n.bins )

    gate.plot <- gate.plot +
      scale_fill_manual( values = rainbow.palette )
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
