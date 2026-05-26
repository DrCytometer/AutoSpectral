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
#' @importFrom ggplot2 geom_contour_filled annotation_raster annotate
#' @importFrom ragg agg_jpeg
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
#' `cividis`, `rocket`, `mako` and `turbo`. Use `rainbow` to be similar to FlowJo
#' or SpectroFlo.
#' @param raster.bins Number of pixels on each axis of the rasterised scatter
#' image. Default is `256L`.
#'
#' @return Saves the plot as a JPEG file in the specified directory.

gate.af.sample.plot <- function(
    plot.data,
    samp,
    af.boundary.upper,
    asp,
    max.points = 5e4,
    color.palette = "viridis",
    raster.bins = 256L
  ) {

  # ---------------------------------------------------------------------------
  # 1. Validate, filter NAs, downsample
  # ---------------------------------------------------------------------------

  valid.idx <- which( !is.na( plot.data[ , 1 ] ) & !is.na( plot.data[ , 2 ] ) )
  if ( length( valid.idx ) == 0 ) {
    warning( "AF plot.data has no valid x/y values; skipping density plot." )
    return( invisible( NULL ) )
  }
  plot.data <- plot.data[ valid.idx, , drop = FALSE ]

  if ( nrow( plot.data ) > max.points ) {
    set.seed( 42 )
    plot.data <- plot.data[ sample( seq_len( nrow( plot.data ) ), max.points ), , drop = FALSE ]
  }

  # ---------------------------------------------------------------------------
  # 2. Biexp transform
  # ---------------------------------------------------------------------------

  biexp.transform <- flowWorkspace::flowjo_biexp(
    channelRange = asp$default.transformation.param$length,
    maxValue     = asp$default.transformation.param$max.range,
    pos          = asp$default.transformation.param$pos,
    neg          = asp$default.transformation.param$neg,
    widthBasis   = asp$default.transformation.param$width,
    inverse      = FALSE
  )

  x.lab <- colnames( plot.data )[ 1 ]
  y.lab <- colnames( plot.data )[ 2 ]

  # transform in a single vectorised call
  trans.mat <- matrix(
    biexp.transform( as.vector( plot.data ) ),
    nrow = nrow( plot.data ),
    ncol = 2L
  )

  # ---------------------------------------------------------------------------
  # 3. Axis geometry (transformed space)
  # ---------------------------------------------------------------------------

  breaks <- asp$ribbon.breaks
  limits <- c( asp$ribbon.plot.min, asp$expr.data.max )
  axis.labels <- sapply( breaks, function( x ) {
    if ( x == 0 ) "0" else parse( text = paste0( "10^", log10( abs( x ) ) ) )
  } )

  x.lim.trans <- biexp.transform( limits )
  y.lim.trans <- biexp.transform( limits )
  x.brk.trans <- biexp.transform( breaks )
  y.brk.trans <- biexp.transform( breaks )

  # ---------------------------------------------------------------------------
  # 4. Rasterised scatter layer
  # ---------------------------------------------------------------------------

  x_breaks_r <- seq( x.lim.trans[1], x.lim.trans[2], length.out = raster.bins + 1L )
  y_breaks_r <- seq( y.lim.trans[1], y.lim.trans[2], length.out = raster.bins + 1L )

  xi  <- pmax( 1L, pmin( raster.bins, findInterval( trans.mat[ , 1 ], x_breaks_r ) ) )
  yi  <- pmax( 1L, pmin( raster.bins, findInterval( trans.mat[ , 2 ], y_breaks_r ) ) )
  idx <- ( xi - 1L ) * raster.bins + yi
  counts <- tabulate( idx, nbins = raster.bins * raster.bins )
  scatter.mat <- matrix( counts, nrow = raster.bins, ncol = raster.bins )

  scatter.log  <- log1p( scatter.mat )
  scatter.max  <- max( scatter.log )
  if ( scatter.max > 0 ) {
    scatter.norm <- scatter.log / scatter.max
    grey.vals    <- grDevices::grey( 1 - scatter.norm * 0.85 )
    grey.vals[ scatter.mat == 0L ] <- NA_character_
  } else {
    grey.vals <- matrix( NA_character_, nrow = raster.bins, ncol = raster.bins )
  }
  scatter.raster <- matrix( grey.vals, nrow = raster.bins, ncol = raster.bins )
  scatter.raster <- scatter.raster[ raster.bins:1, , drop = FALSE ]   # flip Y

  # ---------------------------------------------------------------------------
  # 5. KDE density contours (always fast_kde2d_cpp when available)
  # ---------------------------------------------------------------------------

  plot.data.ggp <- data.frame( x.trans = trans.mat[ , 1 ], y.trans = trans.mat[ , 2 ] )
  bw <- apply( plot.data.ggp, 2, MASS::bandwidth.nrd )

  gate.bound.density <- if (
    requireNamespace( "AutoSpectralRcpp", quietly = TRUE ) &&
    "fast_kde2d_cpp" %in% ls( getNamespace( "AutoSpectralRcpp" ) )
  ) {
    AutoSpectralRcpp::fast_kde2d_cpp(
      x        = plot.data.ggp[ , 1 ],
      y        = plot.data.ggp[ , 2 ],
      n        = 100,
      h        = bw * 0.1,
      x_limits = x.lim.trans,
      y_limits = y.lim.trans
    )
  } else {
    MASS::kde2d(
      x = plot.data.ggp[ , 1 ],
      y = plot.data.ggp[ , 2 ],
      n = 60,
      h = bw * 0.8
    )
  }

  # direct data.frame construction — avoids expand.grid allocation
  n.grid     <- length( gate.bound.density$x )
  density.df <- data.frame(
    x = rep( gate.bound.density$x, times = n.grid ),
    y = rep( gate.bound.density$y, each  = n.grid ),
    z = as.vector( gate.bound.density$z )
  )
  density.df$z[ is.na( density.df$z ) ] <- 0

  max.z          <- max( density.df$z )
  density.breaks <- seq( 0.05 * max.z, max.z, length.out = 11 )
  if ( diff( range( density.breaks ) ) == 0 ) {
    density.breaks <- seq( 0, max.z + 0.1, length.out = 11 )
  }

  # ---------------------------------------------------------------------------
  # 6. Build plot
  # ---------------------------------------------------------------------------

  gate.plot <- ggplot() +
    annotate(
      "rect",
      xmin = x.lim.trans[1], xmax = x.lim.trans[2],
      ymin = y.lim.trans[1], ymax = y.lim.trans[2],
      fill = "white", colour = NA
    ) +
    annotation_raster(
      scatter.raster,
      xmin = x.lim.trans[1], xmax = x.lim.trans[2],
      ymin = y.lim.trans[1], ymax = y.lim.trans[2],
      interpolate = FALSE
    ) +
    geom_contour_filled(
      data        = density.df,
      aes( x = x, y = y, z = z ),
      breaks      = density.breaks,
      alpha       = 1,
      inherit.aes = FALSE,
      na.rm       = TRUE
    ) +
    scale_x_continuous(
      name   = x.lab,
      breaks = x.brk.trans,
      labels = axis.labels
    ) +
    scale_y_continuous(
      name   = y.lab,
      breaks = y.brk.trans,
      labels = axis.labels
    ) +
    coord_cartesian(
      xlim = x.lim.trans,
      ylim = y.lim.trans
    ) +
    theme_bw() +
    theme(
      plot.margin      = margin(
        asp$figure.margin, asp$figure.margin, asp$figure.margin, asp$figure.margin
      ),
      legend.position  = "none",
      axis.ticks       = element_line( linewidth = asp$figure.panel.line.size ),
      axis.text        = element_text( size = asp$figure.axis.text.size ),
      axis.text.x      = element_text( angle = 45, hjust = 1 ),
      axis.title       = element_text( size = asp$figure.axis.title.size ),
      panel.border     = element_rect( fill = NA, linewidth = asp$figure.panel.line.size ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

  # ---------------------------------------------------------------------------
  # 7. Colour scale for density contours
  # ---------------------------------------------------------------------------

  viridis.colors <- c(
    "magma", "inferno", "plasma", "viridis",
    "cividis", "rocket", "mako", "turbo"
  )

  if ( color.palette %in% viridis.colors ) {
    gate.plot <- gate.plot + scale_fill_viridis_d( option = color.palette )
  } else {
    n.bins          <- max( 1L, length( density.breaks ) - 1L )
    rainbow.palette <- grDevices::colorRampPalette( asp$density.palette.base.color )( n.bins )
    gate.plot       <- gate.plot + scale_fill_manual( values = rainbow.palette )
  }

  # ---------------------------------------------------------------------------
  # 8. AF gate boundary (already in raw scatter space — transform before drawing)
  # ---------------------------------------------------------------------------

  if ( !is.null( af.boundary.upper ) ) {
    af.boundary.upper.ggp <- data.frame(
      x = biexp.transform( c( af.boundary.upper$x, af.boundary.upper$x[ 1 ] ) ),
      y = biexp.transform( c( af.boundary.upper$y, af.boundary.upper$y[ 1 ] ) )
    )

    gate.plot <- gate.plot +
      geom_path(
        data      = af.boundary.upper.ggp,
        aes( x, y ),
        color     = "black",
        linewidth = asp$figure.gate.line.size
      )
  }

  # ---------------------------------------------------------------------------
  # 9. Save
  # ---------------------------------------------------------------------------

  ggsave(
    file.path(
      asp$figure.clean.control.dir,
      paste( asp$af.plot.filename, samp, ".jpg", sep = "_" )
    ),
    plot      = gate.plot,
    device    = ragg::agg_jpeg,
    width     = asp$figure.width,
    height    = asp$figure.height,
    limitsize = FALSE
  )

}
