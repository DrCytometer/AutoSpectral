# gate_define_plot.r

#' @title Gate Definition Plot
#'
#' @description
#' This function plots the gate during the definition step,
#' including boundaries and regions.
#'
#' @importFrom ggplot2 ggplot aes scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 theme_bw theme element_line after_stat
#' @importFrom ggplot2 element_text element_rect margin expansion ggsave
#' @importFrom ggplot2 scale_fill_viridis_c scale_fill_gradientn stat_density_2d
#' @importFrom ggplot2 geom_text geom_point geom_path
#' @importFrom scattermore geom_scattermore
#' @importFrom ragg agg_jpeg
#'
#' @param samp Sample identifier.
#' @param gate.data Matrix containing gate data points.
#' @param gate.marker Vector containing gate marker names.
#' @param gate.bound List containing gate boundary information.
#' @param gate.region List containing gate region information.
#' @param gate.population List containing gate population information.
#' @param scatter.and.channel.label Named vector mapping scatter and channel labels.
#' @param asp The AutoSpectral parameter list.
#' @param color.palette Optional character string defining the viridis color
#' palette to be used for the fluorophore traces. Default is `rainbow`, which will
#' be similar to FlowJo or SpectroFlo. Other pptions are the viridis color
#' options: `magma`, `inferno`, `plasma`, `viridis`, `cividis`, `rocket`, `mako`
#' and `turbo`.
#' @param max.points Number of points to plot. Default is `1e5`.
#'
#' @return Saves the plot as a JPEG file in the specified directory.

gate.define.plot <- function(
    samp,
    gate.data,
    gate.marker,
    gate.bound,
    gate.region,
    gate.population,
    scatter.and.channel.label,
    asp,
    color.palette = "plasma",
    max.points = 1e5
  ) {

  # downsample (faster plotting)
  if ( nrow( gate.data ) > max.points ) {
    # random sampling
    set.seed( 42 )
    gate.data <- gate.data[ sample( seq_len( nrow( gate.data ) ), max.points ), ]
  }

  # convert to data frame for plotting
  gate.data.ggp <- data.frame(
    x = gate.data[ , 1 ],
    y = gate.data[ , 2 ] )

  # define search regions and gate boundary asdata frames
  gate.bound.ggp <- data.frame(
    x = c(
      gate.bound$x.low,
      gate.bound$x.high,
      gate.bound$x.high,
      gate.bound$x.low,
      gate.bound$x.low
    ),
    y = c(
      gate.bound$y.low,
      gate.bound$y.low,
      gate.bound$y.high,
      gate.bound$y.high,
      gate.bound$y.low
    )
  )

  gate.region.ggp <- data.frame(
    x = c(
      gate.region$x.low,
      gate.region$x.high,
      gate.region$x.high,
      gate.region$x.low,
      gate.region$x.low
    ),
    y = c(
      gate.region$y.low,
      gate.region$y.low,
      gate.region$y.high,
      gate.region$y.high,
      gate.region$y.low )
  )

  gate.population$boundary$x <- pmin( gate.population$boundary$x, asp$scatter.data.max.x )
  gate.population$boundary$y <- pmin( gate.population$boundary$y, asp$scatter.data.max.y )

  gate.boundary.ggp <- data.frame(
    x = c( gate.population$boundary$x, gate.population$boundary$x[ 1 ] ),
    y = c( gate.population$boundary$y, gate.population$boundary$y[ 1 ] )
  )

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
    ) +
    stat_density_2d(
      aes( fill = after_stat( level ) ),
      geom = "polygon",
      na.rm = TRUE
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
    # search boundary
    geom_path(
      data = gate.bound.ggp,
      aes( x, y ),
      color = "gray30",
      linewidth = asp$figure.gate.line.size,
      linetype = "dashed"
      ) +
    # valid search region
    geom_path(
      data = gate.region.ggp,
      aes( x, y ),
      color = "black",
      linewidth = asp$figure.gate.line.size
      ) +
    # final gate
    geom_path(
      data = gate.boundary.ggp,
      aes( x, y ),
      color = "black",
      linewidth = asp$figure.gate.line.size
      ) +
    # max density points (see asp$density.max.target)
    geom_point(
      data = gate.bound$density.max,
      size = 1.9 * asp$figure.gate.point.size,
      stroke = 0.1 * asp$figure.gate.point.size,
      color = asp$gate.tesselation.color
      ) +
    # label the density points, nudging the label off the point
    geom_text(
      data = gate.bound$density.max,
      aes( label = num.label ),
      hjust = 0, vjust = 0, size = asp$figure.axis.text.size / 2.5,
      color = asp$gate.tesselation.color
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

  ggsave(
    file.path( asp$figure.gate.dir, sprintf( "%s.jpg", samp ) ),
    plot = gate.plot,
    device = ragg::agg_jpeg,
    width = asp$figure.width,
    height = asp$figure.height,
    limitsize = FALSE
  )

}

