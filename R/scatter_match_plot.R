# scatter_match_plot.r

#' @title Plot Scatter-Matching of Universal Negative
#'
#' @description
#' This function generates scatter plots for matching positive and negative
#' expression data, selected based on scatter parameters gates.
#'
#' @importFrom ggplot2 ggplot aes facet_wrap after_stat
#' @importFrom ggplot2 scale_fill_gradientn scale_fill_viridis_c stat_density_2d
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous theme_bw
#' @importFrom ggplot2 theme element_rect element_text margin ggsave
#' @importFrom scattermore geom_scattermore
#' @importFrom ragg agg_jpeg
#'
#' @param pos.expr.data A matrix containing the positive expression data.
#' @param neg.expr.data A matrix containing the negative expression data.
#' @param fluor.name A character string specifying the fluorophore name.
#' @param scatter.param A character vector specifying the scatter parameters.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#' @param color.palette Optional character string defining the viridis color
#' palette to be used for the fluorophore traces. Default is `rainbow`, which will
#' be similar to FlowJo or SpectroFlo. Other pptions are the viridis color
#' options: `magma`, `inferno`, `plasma`, `viridis`, `cividis`, `rocket`, `mako`
#' and `turbo`.
#'
#' @return None. The function saves the generated scatter plot to a file.
#'
#' @export

scatter.match.plot <- function(
    pos.expr.data,
    neg.expr.data,
    fluor.name,
    scatter.param,
    asp,
    color.palette = "rainbow"
  ) {

  # convert to data frames for plotting
  pos.scatter.plot <- data.frame(
    pos.expr.data[ , scatter.param ],
    check.names = FALSE
  )
  neg.scatter.plot <- data.frame(
    neg.expr.data[ , scatter.param ],
    check.names = FALSE
  )

  # set group names
  pos.scatter.plot$group <- fluor.name
  neg.scatter.plot$group <- "Negative"

  # aggregate and set groups as factors
  scatter.plot.data <- rbind( pos.scatter.plot, neg.scatter.plot )
  scatter.plot.data$group <- factor(
    scatter.plot.data$group,
    levels = c( "Negative", fluor.name )
  )
  data.ggp <- data.frame(
    x = scatter.plot.data[ , 1 ],
    y = scatter.plot.data[ , 2 ],
    group = scatter.plot.data$group
  )

  # create axes labels
  x.limits <- c( asp$scatter.data.min.x, asp$scatter.data.max.x )
  x.breaks <- seq( asp$scatter.data.min.x, asp$scatter.data.max.x, asp$data.step )
  x.labels <- paste0( round( x.breaks / 1e6, 1), "e6" )
  y.limits <- c( asp$scatter.data.min.y, asp$scatter.data.max.y )
  y.breaks <- seq( asp$scatter.data.min.y, asp$scatter.data.max.y, asp$data.step )
  y.labels <- paste0( round( y.breaks / 1e6, 1), "e6" )

  # plot the data
  scatter.plot <- suppressWarnings(
    ggplot( data.ggp, aes( x, y, group ) ) +
      geom_scattermore(
        pointsize = asp$figure.gate.point.size * 3,
        alpha = 1,
        na.rm = TRUE
      ) +
      stat_density_2d(
        aes( fill = after_stat( level ) ),
        geom = "polygon",
        contour = TRUE,
        na.rm = TRUE
      ) +
      facet_wrap( ~ group, ncol = 2 ) +
      scale_x_continuous(
        name = scatter.param[ 1 ],
        breaks = x.breaks,
        labels = x.labels,
        limits = x.limits,
        expand = expansion( asp$figure.gate.scale.expand )
      ) +
      scale_y_continuous(
        name = scatter.param[ 2 ],
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
        strip.background = element_rect( fill = "white" ),
        strip.text = element_text(
          size = asp$scatter.match.plot.text.size,
          face = asp$scatter.match.plot.text.face
        ),
        axis.ticks = element_line( linewidth = asp$figure.panel.line.size ),
        axis.text = element_text( size = asp$figure.axis.text.size ),
        axis.title = element_text( size = asp$figure.axis.title.size ),
        panel.border = element_rect( fill = NA, linewidth = asp$figure.panel.line.size ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
  )

  # color options
  viridis.colors <- c(
    "magma", "inferno", "plasma", "viridis",
    "cividis", "rocket", "mako", "turbo"
  )


  # set the color palette on the plot
  if ( color.palette %in% viridis.colors ) {
    scatter.plot <- scatter.plot +
      scale_fill_viridis_c( option = color.palette )
  } else {
    scatter.plot <- scatter.plot +
      scale_fill_gradientn(
        colours = asp$density.palette.base.color,
        values = asp$ribbon.scale.values
      )
  }

  scatter.plot.filename <- paste(
    fluor.name,
    asp$scatter.match.plot.filename,
    sep = "_"
  )

  suppressWarnings(
    ggsave(
      scatter.plot.filename,
      path = asp$figure.scatter.dir.base,
      plot = scatter.plot,
      device = ragg::agg_jpeg,
      width = asp$scatter.match.plot.width,
      height = asp$scatter.match.plot.height,
      limitsize = FALSE
    )
  )

}
