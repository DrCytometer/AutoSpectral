# gate_sample_plot.r

#' @title Plot Pre-defined Gate on Sample
#'
#' @description
#' This function plots a pre-defined gate on a sample, using ggplot2 and other
#' necessary packages.
#'
#' @importFrom ggplot2 ggplot aes scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 theme_bw theme element_line geom_path after_stat
#' @importFrom ggplot2 element_text element_rect margin expansion ggsave
#' @importFrom ggplot2 scale_fill_viridis_d scale_fill_gradientn geom_contour_filled
#' @importFrom ggplot2 scale_fill_viridis_c stat_density_2d
#' @importFrom scattermore geom_scattermore
#' @importFrom ragg agg_jpeg
#' @importFrom MASS kde2d
#'
#' @param samp Sample identifier.
#' @param gate.data Matrix containing gate data points.
#' @param gate.marker Vector containing gate marker names.
#' @param gate.boundary List containing gate boundary information.
#' @param scatter.and.channel.label Named vector mapping scatter and
#' channel labels.
#' @param control.type Type of control: `beads` or `cells`. Deprecated.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#' @param color.palette Optional character string defining the viridis color
#' palette to be used for the fluorophore traces. Default is `mako`. Use `rainbow`
#' to be similar to FlowJo or SpectroFlo. Other options are the viridis color
#' options: `magma`, `inferno`, `plasma`, `viridis`, `cividis`, `rocket`, `mako`
#' and `turbo`.
#' @param max.points Number of points to plot (speeds up plotting). Default is
#' `5e4`.
#' @param switch.n Number of points to trigger the switch to using slower but
#' more robust density plotting. Default is `1e4`.
#'
#' @return Saves the plot as a JPEG file in the specified directory.

gate.sample.plot <- function(
    samp,
    gate.data,
    gate.marker,
    gate.boundary,
    scatter.and.channel.label,
    control.type,
    asp,
    color.palette = "mako",
    max.points = 5e4,
    switch.n = 2e4
  ) {

  # downsample (faster plotting)
  n.points <- nrow( gate.data )
  if ( n.points > max.points ) {
    # random sampling
    set.seed( 42 )
    gate.data <- gate.data[ sample( seq_len( nrow( gate.data ) ), max.points ), ]
    n.points <- max.points
  }

  # restrict to preset limits for plotting
  gate.data[ , 1 ] <- pmin( gate.data[ , 1 ], asp$scatter.data.max.x )
  gate.data[ , 2 ] <- pmin( gate.data[ , 2 ], asp$scatter.data.max.x )

  # get density for plotting
  bw <- apply( gate.data, 2, bandwidth.nrd )

  if ( requireNamespace("AutoSpectralRcpp", quietly = TRUE ) &&
       "fast_kde2d_cpp" %in% ls( getNamespace( "AutoSpectralRcpp" ) ) &&
       n.points > switch.n ) {
    # use C++ function to get density
    gate.bound.density <- AutoSpectralRcpp::fast_kde2d_cpp(
      x = gate.data[ , 1 ],
      y = gate.data[ , 2 ],
      n = 100,
      h = bw * 0.1,
      x_limits = range( gate.data[ , 1 ] ),
      y_limits = range( gate.data[ , 2 ] )
    )
  } else {
    # use slower MASS call
    gate.bound.density <- MASS::kde2d(
      x = gate.data[ , 1 ],
      y = gate.data[ , 2 ],
      n = 60,
      h = bw * 0.8
    )
  }

  # format the density for plotting
  density.df <- expand.grid(
    x = gate.bound.density$x,
    y = gate.bound.density$y
  )
  density.df <- expand.grid( x = gate.bound.density$x, y = gate.bound.density$y )
  density.df$z <- as.vector( gate.bound.density$z )
  density.df <- density.df[ !is.na( density.df$z ), ]
  density.df <- density.df[ !duplicated( density.df[ , c( "x", "y" ) ] ), ]
  max.z <- max( density.df$z, na.rm = TRUE )
  density.breaks <- seq( 0.05 * max.z, max.z, length.out = 11 )

  # convert points to data frame for plotting
  gate.data.ggp <- data.frame(
    x = gate.data[ , 1 ],
    y = gate.data[ , 2 ] )

  # ensure gate limits are drawn onscale
  gate.boundary.ggp <- data.frame(
    x = c( gate.boundary$x, gate.boundary$x[ 1 ] ),
    y = c( gate.boundary$y, gate.boundary$y[ 1 ] )
  )
  gate.boundary.ggp$x <- pmin( gate.boundary.ggp$x, asp$scatter.data.max.x )
  gate.boundary.ggp$y <- pmin( gate.boundary.ggp$y, asp$scatter.data.max.y )

  # get axis labels
  x.lab <- names( which( scatter.and.channel.label == gate.marker[ 1 ] ) )
  y.lab <- names( which( scatter.and.channel.label == gate.marker[ 2 ] ) )

  # create axes labels
  x.limits <- c( asp$scatter.data.min.x, asp$scatter.data.max.x )
  x.breaks <- seq( asp$scatter.data.min.x, asp$scatter.data.max.x, asp$data.step )
  x.labels <- paste0( round( x.breaks / 1e6, 1 ), "e6" )
  y.limits <- c( asp$scatter.data.min.y, asp$scatter.data.max.y )
  y.breaks <- seq( asp$scatter.data.min.y, asp$scatter.data.max.y, asp$data.step )
  y.labels <- paste0( round( y.breaks / 1e6, 1 ), "e6" )

  # set up the plot
  gate.plot <- ggplot( gate.data.ggp, aes( x, y ) ) +
    geom_scattermore(
      pointsize = asp$figure.gate.point.size,
      color = "black",
      alpha = 1,
      na.rm = TRUE
    )

  # use fast contours for many events, slow ggplot for few events, just dots for minimal events
  if ( n.points > switch.n ) {
    gate.plot <- gate.plot +
      geom_contour_filled(
        data = density.df,
        aes( x = x, y = y, z = z ),
        breaks = density.breaks,
        alpha = 1,
        inherit.aes = FALSE
      )
  } else if ( n.points > 1000 ) {
    gate.plot <- gate.plot +
      stat_density_2d(
        aes( fill = after_stat( level ) ),
        geom = "polygon",
        n = 60,
        na.rm = TRUE
      )
  }

  gate.plot <- gate.plot +
    geom_path(
      data = gate.boundary.ggp,
      aes( x, y ),
      color = "darkgoldenrod1 ",
      linewidth = asp$figure.gate.line.size
    ) +
    scale_x_continuous(
      name = x.lab,
      breaks = x.breaks,
      labels = x.labels,
      limits = x.limits,
      expand = expansion( asp$figure.gate.scale.expand )
    ) +
    scale_y_continuous(
      name = y.lab,
      breaks = y.breaks,
      labels = y.labels,
      limits = y.limits,
      expand = expansion( asp$figure.gate.scale.expand )
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
    if ( n.points > switch.n ) {
      gate.plot <- gate.plot + scale_fill_viridis_d( option = color.palette )
    } else {
      gate.plot <- gate.plot + scale_fill_viridis_c( option = color.palette )
    }

  } else {
    gate.plot <- gate.plot +
      scale_fill_gradientn( colors = asp$density.palette.base.color )
  }

  ggsave(
    file.path( asp$figure.gate.dir, sprintf( "%s.jpg", samp ) ),
    plot = gate.plot,
    device = ragg::agg_jpeg,
    width = asp$figure.width,
    height = asp$figure.height,
    limitsize = FALSE
  )

}
