# spectral_ribbon_plot.r

#' @title Spectral Ribbon Plot
#'
#' @description
#' This function generates spectral ribbon plots for positive and negative
#' expression data. The function gets called internally in `clean.controls`. To
#' use the function directly, pass data to `data.list` in the form of a named
#' list of matrices or data.frames.
#'
#' @importFrom ggplot2 ggplot aes scale_y_continuous scale_x_continuous
#' @importFrom ggplot2 annotation_raster annotate ggtitle xlab ylab
#' @importFrom ggplot2 theme_minimal theme element_text element_blank
#' @importFrom ggplot2 ggsave expansion
#' @importFrom patchwork wrap_plots
#' @importFrom flowWorkspace flowjo_biexp
#' @importFrom ragg agg_jpeg
#' @importFrom viridis viridis
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
#' @param max.points Number of events per panel to plot. Default is `5e4`.
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
    plot.height = 10,
    max.points = 5e4
) {

  # ---------------------------------------------------------------------------
  # 1. Assemble list of raw matrices (one per panel), spectral columns only
  # ---------------------------------------------------------------------------

  if ( !is.null( data.list ) ) {
    # user-supplied path
    data.frames <- lapply( data.list, function( df ) {
      as.matrix( df[ , spectral.channel, drop = FALSE ] )
    } )

    if ( is.null( factor.names ) ) {
      if ( is.null( names( data.list ) ) ) {
        factor.names <- paste0( "Fluorophore ", seq_along( data.frames ) )
      } else {
        factor.names <- names( data.list )
      }
    }

    if ( length( factor.names ) != length( data.frames ) )
      stop( "Length of factor.names must match number of datasets in data.list." )

    if ( is.null( figure.dir ) ) figure.dir <- asp$figure.spectral.ribbon.dir
    if ( is.null( title ) & is.null( fluor.name ) )
      title <- paste( factor.names, collapse = "_" )

  } else {
    # internal AutoSpectral paths
    if ( !af ) {
      # downsample inputs first
      pos.ds <- if ( nrow( pos.expr.data ) > max.points ) {
        pos.expr.data[ sample( nrow( pos.expr.data ), max.points ), , drop = FALSE ]
      } else pos.expr.data

      neg.ds <- if ( nrow( neg.expr.data ) > max.points ) {
        neg.expr.data[ sample( nrow( neg.expr.data ), max.points ), , drop = FALSE ]
      } else neg.expr.data

      # now compute background subtraction on downsampled data
      pos.mfi <- apply(
        neg.ds[ , spectral.channel, drop = FALSE ], 2, stats::median
      )
      pos.minus.bg <- sweep(
        pos.ds[ , spectral.channel, drop = FALSE ],
        2, pos.mfi, FUN = "-"
      )
      data.frames <- list(
        pos.minus.bg,
        pos.ds[ , spectral.channel, drop = FALSE ],
        neg.ds[ , spectral.channel, drop = FALSE ]
      )

      if ( is.null( factor.names ) )
        factor.names <- c( fluor.name, paste( "Raw", fluor.name ), "Negative" )
      if ( is.null( title ) )      title      <- "Scatter match"
      if ( is.null( figure.dir ) ) figure.dir <- asp$figure.spectral.ribbon.dir

    } else {
      # intrusive AF event cleaning — downsample each matrix
      downsample <- function( m ) {
        if ( nrow( m ) > max.points )
          m[ sample( nrow( m ), max.points ), , drop = FALSE ]
        else m
      }
      data.frames <- list(
        downsample( pos.expr.data[ , spectral.channel, drop = FALSE ] ),
        downsample( neg.expr.data[ , spectral.channel, drop = FALSE ] ),
        downsample( removed.data[ , spectral.channel, drop = FALSE ] )
      )

      if ( is.null( factor.names ) )
        factor.names <- c(
          paste( "Original", fluor.name ),
          paste( "Cleaned",  fluor.name ),
          "Removed events"
        )
      if ( is.null( title ) )      title      <- "AF removal"
      if ( is.null( figure.dir ) ) figure.dir <- asp$figure.clean.control.dir
    }
  }

  if ( !dir.exists( figure.dir ) ) dir.create( figure.dir )

  # ---------------------------------------------------------------------------
  # 3. Build biexp transform and axis parameters
  # ---------------------------------------------------------------------------

  biexp.transform <- flowWorkspace::flowjo_biexp(
    channelRange = asp$default.transformation.param$length,
    maxValue     = asp$default.transformation.param$max.range,
    pos          = asp$default.transformation.param$pos,
    neg          = asp$default.transformation.param$neg,
    widthBasis   = asp$default.transformation.param$width,
    inverse      = FALSE
  )

  ribbon.breaks <- asp$ribbon.breaks
  ribbon.labels <- sapply( ribbon.breaks, function( x ) {
    if ( x == 0 ) "0" else parse( text = paste0( "10^", log10( abs( x ) ) ) )
  } )
  ribbon.limits <- c( asp$ribbon.plot.min, asp$expr.data.max )

  n_x    <- length( spectral.channel )
  n_y    <- asp$ribbon.bins
  y_min  <- biexp.transform( ribbon.limits[1] )
  y_max  <- biexp.transform( ribbon.limits[2] )
  y_breaks_r <- seq( y_min, y_max, length.out = n_y + 1 )

  # ---------------------------------------------------------------------------
  # 4. Build colour palette function
  # ---------------------------------------------------------------------------

  viridis.colors <- c(
    "magma", "inferno", "plasma", "viridis",
    "cividis", "rocket", "mako", "turbo"
  )
  if ( color.palette %in% viridis.colors ) {
    pal_fn <- function( n ) viridis::viridis( n, option = color.palette )
  } else {
    pal_fn <- grDevices::colorRampPalette( asp$density.palette.base.color )
  }
  # prepend NA so index 0 (empty bins) maps to transparent
  palette.256 <- c( NA_character_, pal_fn( 256 ) )

  # ---------------------------------------------------------------------------
  # 5. make.raster: transform -> bin via C++ -> colour-map
  # ---------------------------------------------------------------------------

  make.raster <- function( m ) {
    # apply biexp transform to the full matrix in one call
    trans <- matrix(
      biexp.transform( as.vector( m ) ),
      nrow = nrow( m ),
      ncol = n_x
    )

    # bin into count matrix via C++ (column-major, cache-friendly)
    if ( requireNamespace( "AutoSpectralRcpp", quietly = TRUE ) &&
         "bin_matrix_cpp" %in% ls( getNamespace( "AutoSpectralRcpp" ) ) ) {
      mat <- AutoSpectralRcpp::bin_matrix_cpp( trans, y_breaks_r, n_y )
    } else {
      # pure-R fallback
      mat <- matrix( 0L, nrow = n_y, ncol = n_x )
      for ( col in seq_len( n_x ) ) {
        yi <- findInterval( trans[ , col ], y_breaks_r, rightmost.closed = TRUE )
        yi <- pmax( 1L, pmin( n_y, yi ) )
        mat[ , col ] <- tabulate( yi, nbins = n_y )
      }
    }

    # log-scale counts and map to palette
    mat_log  <- log1p( mat )
    mat_max  <- max( mat_log, na.rm = TRUE )
    if ( mat_max == 0 ) return( matrix( NA_character_, nrow = n_y, ncol = n_x ) )
    mat_norm <- mat_log / mat_max

    # col_idx: 0 for empty bins, 1-256 for data; shift by +1 to index palette.256
    col_idx <- as.integer( ceiling( mat_norm * 256 ) ) + 1L
    col_mat <- matrix( palette.256[ col_idx ], nrow = n_y, ncol = n_x )

    # flip vertically: annotation_raster row 1 is the top of the panel
    col_mat[ n_y:1, , drop = FALSE ]
  }

  raster.list <- lapply( data.frames, make.raster )

  # ---------------------------------------------------------------------------
  # 6. Shared theme (matches original facet_wrap version exactly)
  # ---------------------------------------------------------------------------

  panel.theme <- theme_minimal() +
    theme(
      axis.text.x      = element_text(
        angle = asp$ribbon.plot.axis.text.angle,
        vjust = 1,
        hjust = 1
      ),
      panel.grid.minor = element_blank(),
      legend.position  = "none",
      strip.text       = element_text(
        size = asp$ribbon.plot.strip.text.size,
        face = asp$ribbon.plot.strip.text.face
      )
    )

  # ---------------------------------------------------------------------------
  # 7. Build one ggplot panel per group, combine with patchwork
  # ---------------------------------------------------------------------------

  panels <- lapply( seq_along( factor.names ), function( i ) {
    # suppress x-axis text on all but the bottom panel
    x.text <- if ( i == length( factor.names ) ) {
      element_text(
        angle = asp$ribbon.plot.axis.text.angle,
        vjust = 1,
        hjust = 1
      )
    } else {
      element_blank()
    }

    ggplot() +
      # white background so NA (empty) cells show as white
      annotate(
        "rect",
        xmin = 0.5, xmax = n_x + 0.5,
        ymin = y_min, ymax = y_max,
        fill = "white", colour = NA
      ) +
      annotation_raster(
        raster.list[[ i ]],
        xmin = 0.5,   xmax = n_x + 0.5,
        ymin = y_min, ymax = y_max,
        interpolate = FALSE
      ) +
      scale_x_continuous(
        name   = if ( i == length( factor.names ) ) "Detector" else NULL,
        breaks = seq_len( n_x ),
        labels = spectral.channel,
        limits = c( 0.5, n_x + 0.5 ),
        expand = expansion( 0 )
      ) +
      scale_y_continuous(
        name   = "Intensity",
        limits = c( y_min, y_max ),
        breaks = biexp.transform( ribbon.breaks ),
        labels = ribbon.labels,
        expand = expansion( 0 )
      ) +
      ggtitle( factor.names[[ i ]] ) +
      panel.theme +
      theme( axis.text.x = x.text )
  } )

  ribbon.plot <- patchwork::wrap_plots( panels, ncol = 1 )

  # ---------------------------------------------------------------------------
  # 8. Save or return
  # ---------------------------------------------------------------------------

  if ( save ) {
    ribbon.plot.filename <- paste(
      title, fluor.name, asp$ribbon.plot.filename,
      sep = "_"
    )

    ggsave(
      ribbon.plot.filename,
      plot   = ribbon.plot,
      device = ragg::agg_jpeg,
      path   = figure.dir,
      width  = plot.width,
      height = plot.height,
      limitsize = FALSE,
      create.dir = TRUE
    )

    if ( !is.null( data.list ) ) return( ribbon.plot )

  } else {
    return( ribbon.plot )
  }
}
