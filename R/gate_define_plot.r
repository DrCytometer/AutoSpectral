# gate_define_plot.r

#' @title Gate Definition Plot
#'
#' @description
#' This function plots the gate during the definition step,
#' including boundaries and regions.
#'
#' @importFrom ggplot2 ggplot aes scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 scale_color_gradientn theme_bw theme element_line
#' @importFrom ggplot2 element_text element_rect margin expansion ggsave
#' @importFrom ggplot2 guide_colorbar geom_text geom_point
#' @importFrom scattermore geom_scattermore
#' @importFrom fields interp.surface
#' @importFrom rlang .data
#'
#' @param samp Sample identifier.
#' @param gate.data Matrix containing gate data points.
#' @param gate.marker Vector containing gate marker names.
#' @param gate.bound List containing gate boundary information.
#' @param gate.region List containing gate region information.
#' @param gate.population List containing gate population information.
#' @param scatter.and.channel.label Named vector mapping scatter and
#' channel labels.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#' @param color.palette Optional character string defining the viridis color
#' palette to be used for the fluorophore traces. Default is `rainbow`, which will
#' be similar to FlowJo or SpectroFlo. Other pptions are the viridis color
#' options: `magma`, `inferno`, `plasma`, `viridis`, `cividis`, `rocket`, `mako`
#' and `turbo`.
#'
#' @return Saves the plot as a JPEG file in the specified directory.

gate.define.plot <- function( samp, gate.data, gate.marker, gate.bound,
                              gate.region, gate.population,
                              scatter.and.channel.label, asp,
                              color.palette = "rainbow" )
{

  gate.data.ggp <- data.frame(
    x = gate.data[ , 1 ],
    y = gate.data[ , 2 ] )

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

  gate.population$boundary$x[
    gate.population$boundary$x > asp$scatter.data.max.x ] <- asp$scatter.data.max.x
  gate.population$boundary$y[
    gate.population$boundary$y > asp$scatter.data.max.y ] <- asp$scatter.data.max.y

  gate.boundary.ggp <- data.frame(
    x = c( gate.population$boundary$x,
           gate.population$boundary$x[ 1 ] ),
    y = c( gate.population$boundary$y,
           gate.population$boundary$y[ 1 ] )
  )

  x.lab.idx <- which( scatter.and.channel.label == gate.marker[ 1 ] )
  x.lab <- names( scatter.and.channel.label[ x.lab.idx ] )
  y.lab.idx <- which( scatter.and.channel.label == gate.marker[ 2 ] )
  y.lab <- names( scatter.and.channel.label[ y.lab.idx ] )

  gate.plot <- ggplot( gate.data.ggp, aes( .data$x, .data$y ) ) +
    geom_scattermore(
      pointsize = asp$figure.gate.point.size,
      alpha = 1, na.rm = TRUE ) +
    stat_density_2d(
      aes( fill = after_stat( level ) ),
      geom = "polygon",
      contour = TRUE,
      na.rm = TRUE ) +
    scale_x_continuous(
      name = x.lab,
      breaks = seq(
        asp$scatter.data.min.x, asp$scatter.data.max.x, asp$data.step
        ),
      labels = paste0(
        round( seq( asp$scatter.data.min.x, asp$scatter.data.max.x,
                    asp$data.step ) / 1e6, 1 ), "e6"
        ),
      limits = c(
        asp$scatter.data.min.x, asp$scatter.data.max.x
        ),
      expand = expansion( asp$figure.gate.scale.expand ) ) +
    scale_y_continuous(
      name = y.lab,
      breaks = seq(
        asp$scatter.data.min.y, asp$scatter.data.max.y, asp$data.step
        ),
      labels = paste0(
        round(
          seq(
            asp$scatter.data.min.y, asp$scatter.data.max.y, asp$data.step
            ) / 1e6, 1 ), "e6"
        ),
      limits = c( asp$scatter.data.min.y, asp$scatter.data.max.y ),
      expand = expansion( asp$figure.gate.scale.expand ) ) +
    geom_path(
      aes( .data$x, .data$y, color = NULL ),
      data = gate.bound.ggp, linewidth = asp$figure.gate.line.size,
      linetype = "dashed"
      ) +
    geom_path(
      aes( .data$x, .data$y, color = NULL ),
      data = gate.region.ggp, linewidth = asp$figure.gate.line.size
      ) +
    geom_path(
      aes( .data$x, .data$y, color = NULL ),
      data = gate.boundary.ggp, linewidth = asp$figure.gate.line.size
      ) +
    geom_point(
      data = gate.bound$density.max,
      size = 1.9 * asp$figure.gate.point.size,
      stroke = 0.1 * asp$figure.gate.point.size,
      color = asp$gate.tesselation.color
      ) +
    geom_text(
      data = gate.bound$density.max,
      aes( label = .data$num.label ),
      hjust = 0, vjust = 0, size = asp$figure.axis.text.size / 2.5,
      color = asp$gate.tesselation.color
      ) +
    theme_bw() +
    theme(
      plot.margin = margin( asp$figure.margin, asp$figure.margin,
                            asp$figure.margin, asp$figure.margin ),
      legend.position = "none",
      axis.ticks = element_line( linewidth = asp$figure.panel.line.size ),
      axis.text = element_text( size = asp$figure.axis.text.size ),
      axis.text.x = element_text( angle = 45, hjust = 1 ),
      axis.title = element_text( size = asp$figure.axis.title.size ),
      panel.border = element_rect( linewidth = asp$figure.panel.line.size ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
      )

  # color options
  virids.colors <- c( "magma", "inferno", "plasma", "viridis", "cividis",
                      "rocket", "mako", "turbo" )
  if ( color.palette %in% virids.colors ) {
    gate.plot <- gate.plot +
      scale_fill_viridis_c( option = color.palette )
  } else {
    gate.plot <- gate.plot +
      scale_fill_gradientn(
        colours = asp$density.palette.base.color, values = asp$ribbon.scale.values )
  }

  ggsave(
    file.path( asp$figure.gate.dir, sprintf( "%s.jpg", samp ) ),
    plot = gate.plot,
    width = asp$figure.width,
    height = asp$figure.height
    )

}

