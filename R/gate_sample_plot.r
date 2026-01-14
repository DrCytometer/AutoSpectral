# gate_sample_plot.r

#' @title Plot Pre-defined Gate on Sample
#'
#' @description
#' This function plots a pre-defined gate on a sample, using ggplot2 and other
#' necessary packages.
#'
#' @importFrom ggplot2 ggplot aes scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 scale_fill_gradientn theme_bw theme element_line
#' @importFrom ggplot2 element_text element_rect margin expansion ggsave
#' @importFrom ggplot2 stat_density_2d after_stat geom_path scale_fill_viridis_c
#' @importFrom scattermore geom_scattermore
#' @importFrom rlang .data
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
#' palette to be used for the fluorophore traces. Default is `rainbow`, which will
#' be similar to FlowJo or SpectroFlo. Other pptions are the viridis color
#' options: `magma`, `inferno`, `plasma`, `viridis`, `cividis`, `rocket`, `mako`
#' and `turbo`.
#'
#' @return Saves the plot as a JPEG file in the specified directory.

gate.sample.plot <- function( samp, gate.data, gate.marker, gate.boundary,
                              scatter.and.channel.label, control.type, asp,
                              color.palette = "rainbow" )
{

  gate.data.ggp <- data.frame(
    x = gate.data[ , 1 ],
    y = gate.data[ , 2 ] )

  gate.boundary$x[ gate.boundary$x > asp$scatter.data.max.x ] <- asp$scatter.data.max.x
  gate.boundary$y[ gate.boundary$y > asp$scatter.data.max.y ] <- asp$scatter.data.max.y

  gate.boundary.ggp <- data.frame(
    x = c( gate.boundary$x,
           gate.boundary$x[ 1 ] ),
    y = c( gate.boundary$y,
           gate.boundary$y[ 1 ] )
  )

  x.lab.idx <- which( scatter.and.channel.label == gate.marker[ 1 ] )
  x.lab <- names( scatter.and.channel.label[ x.lab.idx ] )
  y.lab.idx <- which( scatter.and.channel.label == gate.marker[ 2 ] )
  y.lab <- names( scatter.and.channel.label[ y.lab.idx ] )

  gate.plot <- ggplot( gate.data.ggp, aes( x = .data$x, y = .data$y ) ) +
    geom_scattermore(
      pointsize = asp$figure.gate.point.size,
      alpha = 1, na.rm = TRUE
    ) +
    stat_density_2d(
      aes( fill = after_stat( level ) ),
      geom = "polygon",
      contour = TRUE,
      na.rm = TRUE ) +
    geom_path(
      aes( x = .data$x, y = .data$y ),
      data = gate.boundary.ggp,
      color = "black",
      linewidth = asp$figure.gate.line.size
    ) +
    scale_x_continuous(
      name = x.lab,
      breaks = seq(
        asp$scatter.data.min.x, asp$scatter.data.max.x, asp$data.step
        ),
      labels = paste0(
        round(
          seq(
            asp$scatter.data.min.x, asp$scatter.data.max.x, asp$data.step
            ) / 1e6, 1 ), "e6" ),
      limits = c( asp$scatter.data.min.x, asp$scatter.data.max.x ),
      expand = expansion( asp$figure.gate.scale.expand )
    ) +
    scale_y_continuous(
      name = y.lab,
      breaks = seq( asp$scatter.data.min.y, asp$scatter.data.max.y, asp$data.step ),
      labels = paste0(
        round(
          seq(
            asp$scatter.data.min.y, asp$scatter.data.max.y, asp$data.step
            ) / 1e6, 1 ), "e6" ),
      limits = c( asp$scatter.data.min.y, asp$scatter.data.max.y ),
      expand = expansion( asp$figure.gate.scale.expand )
    ) +
    theme_bw() +
    theme(
      plot.margin = margin(
        asp$figure.margin, asp$figure.margin,
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
