# gate_af_identify_plot.r

#' @title Plot Autofluorescence Identification Gate
#'
#' @description
#' This function plots the sample being used to identify intrusive
#' autofluorescence in the single-stained controls. The input data are expected
#' to be PCA projections of the unstained sample with an accompanying region to
#' identify the low-autofluorescence cell region.
#'
#' @importFrom ggplot2 ggplot aes scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 scale_color_gradientn theme_bw theme element_line
#' @importFrom ggplot2 element_text element_rect margin expansion ggsave
#' @importFrom ggplot2 guide_colorbar
#' @importFrom scattermore geom_scattermore
#' @importFrom fields interp.surface
#' @importFrom rlang .data
#'
#' @param gate.data Matrix containing autofluorescence data points.
#' @param samp Sample identifier.
#' @param gate.region Dataframe containing region boundary information.
#' @param gate.bound.density Density (e.g., from MASS:kde2d) information
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#' @param color.palette Optional character string defining the viridis color
#' palette to be used for the fluorophore traces. Default is `rainbow`, which will
#' be similar to FlowJo or SpectroFlo. Other pptions are the viridis color
#' options: `magma`, `inferno`, `plasma`, `viridis`, `cividis`, `rocket`, `mako`
#' and `turbo`.
#'
#' @return Saves the plot as a JPEG file in the specified directory.

gate.af.identify.plot <- function( gate.data, samp, gate.region,
                                   gate.bound.density, asp,
                                   color.palette = "rainbow" ) {

  if ( nrow( gate.data ) == 0 ||
      all( is.na( gate.data[ , 1 ] ) ) ||
      all( is.na( gate.data[ , 2 ] ) ) ) {

    warning( "gate.data has no valid x/y values; skipping density plot." )

    return( invisible( NULL ) )
  }

  gate.data.ggp <- data.frame(
    x = gate.data[ , 1 ],
    y = gate.data[ , 2 ]
    )

  # get axis labels
  axes.labels <- colnames( gate.data )

  # get data range & step
  x.min <- min( gate.data[ , 1 ] )
  x.max <- max( gate.data[ , 1 ] )

  y.min <- min( gate.data[ , 2 ] )
  y.max <- max( gate.data[ , 2 ] )

  x.breaks <- round( seq( x.min, x.max, length.out = 10 ) )
  y.breaks <- round( seq( y.min, y.max, length.out = 10 ) )

  gate.plot <- suppressWarnings(
    ggplot( gate.data.ggp, aes( .data$x, .data$y ) ) +
      geom_scattermore(
        pointsize = 1.2 * asp$figure.gate.point.size,
        alpha = 1, na.rm = TRUE ) +
      stat_density_2d(
        aes( fill = after_stat( level ) ),
        geom = "polygon",
        contour = TRUE,
        na.rm = TRUE ) +
      scale_x_continuous(
        name = axes.labels[ 1 ],
        breaks = x.breaks,
        labels = x.breaks,
        limits = c( x.min, x.max ),
        expand = expansion( asp$af.figure.gate.scale.expand ) ) +
      scale_y_continuous(
        name = axes.labels[ 2 ],
        breaks = y.breaks,
        labels = y.breaks,
        limits = c( y.min, y.max ),
        expand = expansion( asp$af.figure.gate.scale.expand ) ) +
      geom_path( aes( .data$x, .data$y, color = NULL ),
                 data = gate.region, linewidth = asp$figure.gate.line.size ) +
      theme_bw() +
      theme( plot.margin = margin(
        asp$figure.margin, asp$figure.margin,
        asp$figure.margin, asp$figure.margin ),
        legend.position = "none",
        axis.ticks = element_line( linewidth = asp$figure.panel.line.size ),
        axis.text = element_text( size = asp$figure.axis.text.size ),
        axis.text.x = element_text( angle = 45, hjust = 1 ),
        axis.title = element_text( size = asp$figure.axis.title.size ),
        panel.border = element_rect( linewidth = asp$figure.panel.line.size ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() )
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

  suppressWarnings(
    ggsave(
      file.path(
        asp$figure.clean.control.dir,
        paste( asp$af.plot.define.filename, samp, ".jpg", sep = "_" ) ),
      plot = gate.plot, width = asp$figure.width,
      height = asp$figure.height
    )
  )

}
