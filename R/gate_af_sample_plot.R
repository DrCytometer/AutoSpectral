# gate_af_sample_plot.r


#' @title Plot Gate Autofluorescence Sample
#'
#' @description
#' This function plots the gate autofluorescence sample, including upper and
#' lower boundaries, using ggplot2 and other necessary packages.
#'
#' @importFrom ggplot2 ggplot aes scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 stat_density_2d theme_bw theme element_line
#' @importFrom ggplot2 element_text element_rect margin ggsave
#' @importFrom ggplot2 scale_fill_viridis_c geom_path
#' @importFrom scattermore geom_scattermore
#'

#' @param plot.data Matrix containing autofluorescence data points.
#' @param samp Sample identifier.
#' @param af.boundary.upper Matrix containing upper boundary information.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#' @param max.points Number of points to plot (speeds up plotting). Default is
#' `1e5`.
#' @param color.palette Optional character string defining the viridis color
#' palette to be used for the fluorophore traces. Default is `viridis`. Options
#' are the viridis color options: `magma`, `inferno`, `plasma`, `viridis`,
#' `cividis`, `rocket`, `mako` and `turbo`.
#'
#' @return Saves the plot as a JPEG file in the specified directory.

gate.af.sample.plot <- function( plot.data, samp, af.boundary.upper, asp,
                                 max.points = 1e5,
                                 color.palette = "viridis" ) {

  # set plot limits
  breaks <- asp$ribbon.breaks

  axis.labels <- sapply( breaks, function( x ) {
    if ( x == 0 ) "0" else parse( text = paste0( "10^", log10( abs( x ) ) ) )
  } )

  limits <- c( asp$ribbon.plot.min, asp$expr.data.max )

  # get axis labels
  x.lab <- colnames( plot.data )[ 1 ]
  y.lab <- colnames( plot.data )[ 2 ]

  # set transform
  biexp.transform <- flowWorkspace::flowjo_biexp(
    channelRange = asp$default.transformation.param$length,
    maxValue = asp$default.transformation.param$max.range,
    pos = asp$default.transformation.param$pos,
    neg = asp$default.transformation.param$neg,
    widthBasis = asp$default.transformation.param$width,
    inverse = FALSE )
  biexp.inverse <- flowWorkspace::flowjo_biexp(
    channelRange = asp$default.transformation.param$length,
    maxValue = asp$default.transformation.param$max.range,
    pos = asp$default.transformation.param$pos,
    neg = asp$default.transformation.param$neg,
    widthBasis = asp$default.transformation.param$width,
    inverse = TRUE )

  plot.biexp.transform <- scales::trans_new(
    name = "biexp",
    transform = biexp.transform,
    inverse = biexp.inverse
  )

  # downsample data
  if ( nrow( plot.data ) < max.points )
    max.points <- nrow( plot.data )

  plot.data <- data.frame(
    x = plot.data[ 1:max.points, 1 ],
    y = plot.data[ 1:max.points, 2 ]
  )

  # transform data
  plot.data$x.trans <- plot.biexp.transform$transform( plot.data$x )
  plot.data$y.trans <- plot.biexp.transform$transform( plot.data$y )

  gate.plot <- suppressWarnings(
    ggplot( plot.data, aes( x = x.trans, y = y.trans ) ) +
      geom_scattermore(
        # aes( x = x.trans, y = y.trans ),
        pointsize = asp$figure.gate.point.size,
        alpha = 1,
        na.rm = TRUE
      ) +
      stat_density_2d(
        aes( fill = after_stat( level ) ),
        geom = "polygon",
        contour = TRUE,
        na.rm = TRUE ) +
      scale_x_continuous(
        name = x.lab,
        breaks = plot.biexp.transform$transform( breaks ),
        limits = plot.biexp.transform$transform( limits ),
        labels = axis.labels
      ) +
      scale_y_continuous(
        name = y.lab,
        breaks = plot.biexp.transform$transform( breaks ),
        limits = plot.biexp.transform$transform( limits ),
        labels = axis.labels
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
  )

  if ( !is.null( af.boundary.upper ) ){
    af.boundary.upper.ggp <- data.frame(
      x = c( af.boundary.upper$x,
             af.boundary.upper$x[ 1 ] ),
      y = c( af.boundary.upper$y,
             af.boundary.upper$y[ 1 ] )
    )

    af.boundary.upper.ggp$x <- plot.biexp.transform$transform( af.boundary.upper.ggp$x )
    af.boundary.upper.ggp$y <- plot.biexp.transform$transform( af.boundary.upper.ggp$y )

    gate.plot <- gate.plot +
      geom_path( aes( .data$x, .data$y, color = NULL ),
               data = af.boundary.upper.ggp, linewidth = asp$figure.gate.line.size )
  }

  # color options
  virids.colors <- c( "magma", "inferno", "plasma", "viridis", "cividis",
                      "rocket", "mako", "turbo" )
  if ( color.palette %in% virids.colors ) {
    gate.plot <- gate.plot +
      scale_fill_viridis_c( option = color.palette )
  } else {
    gate.plot <- gate.plot +
      scale_fill_gradientn( colours = asp$density.palette.base.color,
                            values = asp$ribbon.scale.values )
  }

  suppressWarnings(
    ggsave(
      file.path(
        asp$figure.clean.control.dir,
        paste( asp$af.plot.filename, samp, ".jpg", sep = "_" ) ),
      plot = gate.plot, width = asp$figure.width,
      height = asp$figure.height
    )
  )

}
