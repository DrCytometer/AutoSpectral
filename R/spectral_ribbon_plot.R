# spectral_ribbon_plot.r

#' @title Spectral Ribbon Plot
#'
#' @description
#' This function generates spectral ribbon plots for positive and negative
#' expression data. The function gets called internally in `clean.controls`. To
#' use the function directly, pass data to `data.list` in the form of a named
#' list of matrices or data.frames.
#'
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes scale_y_continuous geom_bin2d facet_wrap xlab
#' @importFrom ggplot2 ylab scale_fill_gradientn theme_minimal theme element_text
#' @importFrom ggplot2 element_blank ggsave scale_fill_viridis_c
#' @importFrom flowWorkspace flowjo_biexp
#' @importFrom scales trans_new
#'
#' @param pos.expr.data Internal argument for `clean.controls`. A matrix
#' containing the positive expression data. Default is `NULL`.
#' @param neg.expr.data Internal argument for `clean.controls`. A matrix
#' containing the negative expression data. Default is `NULL`.
#' @param removed.data Internal argument for `clean.controls`. A matrix
#' containing the removed data, if applicable. Default is `NULL`. If omitted,
#' only two groups (facets) are plotted.
#' @param spectral.channel A character vector specifying the spectral channels.
#' Recommended: use `colnames(spectra)` or `flow.control$spectral.channel`.
#' @param asp The AutoSpectral parameter list. Prepare using get.autospectral.param.
#' @param fluor.name An optional character string specifying the fluorophore
#' name for plot titles and filename. Default is `NULL`.
#' @param title An optional character string to prefix the plot file name.
#' Default is `NULL`.
#' @param figure.dir Output folder where the figures will be created. Default is
#' `NULL`, enabling automatic selection inside AutoSpectral. For user-supplied
#' data, `asp$figure.spectral.ribbon.dir` will be used if `NULL`.
#' @param factor.names Optional titles for the facets on the plot. Default is
#' `NULL`, enabling automatic selection inside AutoSpectral. If `data.list` is
#' named, `factor.names` will be pulled from that.
#' @param save Logical, default is `TRUE`. If `TRUE`, the plot is saved to
#' `figure.dir`. Otherwise, it is returned to the viewer only.
#' @param color.palette Optional character string defining the color palette to
#' be used. Default is `rainbow`, mimicking a FlowJo scheme. Other choices
#' are the viridis color options: `magma`, `inferno`, `plasma`, `viridis`,
#' `cividis`, `rocket`, `mako` and `turbo`.
#' @param data.list Provide data here. A (named) list of matrices or dataframes
#' for plotting. These should probably be flow expression data. Data provided
#' will be subsetted to the columns in `spectral.channel`.
#' @param af Internal argument for `clean.controls`. A logical value indicating
#' whether autofluorescence removal is being performed. Default is `FALSE`.
#' @param plot.width Width of the saved plot. Default is `15`.
#' @param plot.height Height of the saved plot. Default is `10`.
#'
#' @return None. The function saves the generated spectral ribbon plot to
#' a file.
#'
#' @export

spectral.ribbon.plot <- function(
    pos.expr.data = NULL,
    neg.expr.data = NULL,
    removed.data = NULL,
    spectral.channel,
    asp,
    fluor.name = NULL,
    title = NULL,
    figure.dir = NULL,
    factor.names = NULL,
    save = TRUE,
    color.palette = "rainbow",
    data.list = NULL,
    af = FALSE,
    plot.width = 15,
    plot.height = 10 ) {


  # for user-provided data: list of data.frames/matrices to `data.list`
  if ( !is.null( data.list ) ) {
    data.frames <- lapply( data.list, function( df ) {
      df[ , spectral.channel, drop = FALSE ]
    } )

    if ( is.null( factor.names ) ) {
      if ( is.null( names( data.list ) ) ) {
        factor.names <- paste0( "Fluorophore ", seq_along( data.frames ) )
      } else {
        factor.names <- names( data.list )
      }
    }

    if ( length( factor.names) != length( data.frames ) )
      stop( "Length of factor.names must match number of datasets in data.list." )

    if ( is.null( figure.dir ) )
      figure.dir <- asp$figure.spectral.ribbon.dir

    if ( is.null( title ) & is.null( fluor.name ) )
      title <- paste( factor.names, collapse = "_" )

  } else {
    # internal AutoSpectral control-cleaning data wrangling
    data.frames <- list()

    if ( !af ) {
      # scatter matching with background subtraction
      pos.mfi <- apply(
        neg.expr.data[ , spectral.channel, drop = FALSE ], 2, median
      )
      pos.minus.bg <- sweep(
        pos.expr.data[ , spectral.channel, drop = FALSE ],
        2,
        pos.mfi,
        FUN = "-"
      )
      data.frames <- list(
        pos.minus.bg,
        pos.expr.data[ , spectral.channel, drop = FALSE ],
        neg.expr.data[ , spectral.channel, drop = FALSE ]
      )

      if ( is.null( factor.names ) )
        factor.names <- c( fluor.name, paste( "Raw", fluor.name ), "Negative" )

      if ( is.null( title ) )
        title <- "Scatter match"

      if ( is.null( figure.dir ) )
        figure.dir <- asp$figure.spectral.ribbon.dir

    } else {
      # Intrusive AF event cleaning
      data.frames <- list(
        pos.expr.data[ , spectral.channel, drop = FALSE ],
        neg.expr.data[ , spectral.channel, drop = FALSE ],
        removed.data[ , spectral.channel, drop = FALSE ]
      )

      if ( is.null( factor.names ) )
        factor.names <- c(
          paste( "Original", fluor.name ),
          paste( "Cleaned", fluor.name ),
          "Removed events"
        )

      if ( is.null( title ) )
        title <- "AF removal"

      if ( is.null( figure.dir ) )
        figure.dir <- asp$figure.clean.control.dir

    }
  }

  # labels
  names( data.frames ) <- factor.names
  data.frames <- Map( function( df, nm ) {
    df <- data.frame( df, check.names = FALSE )
    df$group <- nm
    df
  }, data.frames, factor.names )

  # rearrange
  ribbon.plot.data <- do.call( rbind, data.frames )
  ribbon.plot.data$group <- factor( ribbon.plot.data$group, levels = factor.names )

  ribbon.plot.long <- tidyr::pivot_longer(
    ribbon.plot.data,
    cols = -group,
    names_to = "channel",
    values_to = "value"
  )
  ribbon.plot.long$channel <- factor( ribbon.plot.long$channel,
                                      levels = unique( ribbon.plot.long$channel ) )

  # setting scales
  ribbon.breaks <- asp$ribbon.breaks
  ribbon.labels <- sapply( ribbon.breaks, function( x ) {
    if ( x == 0 ) "0" else parse( text = paste0( "10^", log10( abs( x ) ) ) )
  } )
  ribbon.limits <- c( asp$ribbon.plot.min, asp$expr.data.max )

  biexp.transform <- flowWorkspace::flowjo_biexp(
    channelRange = asp$default.transformation.param$length,
    maxValue     = asp$default.transformation.param$max.range,
    pos          = asp$default.transformation.param$pos,
    neg          = asp$default.transformation.param$neg,
    widthBasis   = asp$default.transformation.param$width,
    inverse      = FALSE
  )
  biexp.inverse <- flowWorkspace::flowjo_biexp(
    channelRange = asp$default.transformation.param$length,
    maxValue     = asp$default.transformation.param$max.range,
    pos          = asp$default.transformation.param$pos,
    neg          = asp$default.transformation.param$neg,
    widthBasis   = asp$default.transformation.param$width,
    inverse      = TRUE
  )
  plot.biexp.transform <- scales::trans_new(
    name = "biexp",
    transform = biexp.transform,
    inverse = biexp.inverse
  )

  # create plot
  ribbon.plot <- suppressWarnings(
    ggplot( ribbon.plot.long, aes( channel, value ) ) +
      scale_y_continuous(
        trans = plot.biexp.transform,
        breaks = ribbon.breaks,
        limits = ribbon.limits,
        labels = ribbon.labels
      ) +
      geom_bin2d(
        bins = c(
          length(
            unique(
              ribbon.plot.long$channel
            )
          ),
          asp$ribbon.bins
        ),
        boundary = 0.5
      ) +
      facet_wrap( ~group, ncol = 1 ) +
      xlab( "Detector" ) +
      ylab( "Intensity" ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(
          angle = asp$ribbon.plot.axis.text.angle,
          vjust = 1,
          hjust = 1
        ),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.text = element_text(
          size = asp$ribbon.plot.strip.text.size,
          face = asp$ribbon.plot.strip.text.face
        )
      )
  )

  # color options
  virids.colors <- c( "magma", "inferno", "plasma", "viridis", "cividis",
                      "rocket", "mako", "turbo" )
  if ( color.palette %in% virids.colors ) {
    ribbon.plot <- ribbon.plot +
      scale_fill_viridis_c( option = color.palette )
  } else {
    ribbon.plot <- ribbon.plot +
      scale_fill_gradientn(
        colours = asp$density.palette.base.color,
        values = asp$ribbon.scale.values
      )
  }

  # save or return
  if ( save ) {
    ribbon.plot.filename <- paste(
      title, fluor.name, asp$ribbon.plot.filename,
      sep = "_"
    )

    suppressWarnings(
      ggsave(
        ribbon.plot.filename,
        plot = ribbon.plot,
        path = figure.dir,
        width = plot.width,
        height = plot.height,
        create.dir = TRUE
      )
    )

  } else {
    suppressWarnings( print( ribbon.plot ) )
  }
}
