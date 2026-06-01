# gate_sample_plot.r

#' @title Plot Pre-defined Gate on Sample
#'
#' @description
#' This function plots a pre-defined gate on a sample, using ggplot2 and other
#' necessary packages.
#'
#' @importFrom ggplot2 ggplot aes scale_x_continuous scale_y_continuous element_blank
#' @importFrom ggplot2 theme_bw theme element_line geom_path after_stat coord_cartesian
#' @importFrom ggplot2 element_text element_rect margin expansion ggsave
#' @importFrom ggplot2 scale_fill_viridis_d geom_contour_filled annotation_raster annotate
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ragg agg_jpeg
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
#' @param gate.color Color to plot the gate boundary line, default is `darkgoldenrod1`.
#' @param switch.n Minimum number of points required to render density contours.
#' Below this threshold only the rasterised scatter layer is shown. Default is `1e4`.
#' @param raster.bins Number of pixels on each axis of the rasterised scatter
#' image. Higher values give finer detail at the cost of a little extra memory.
#' Default is `256L`.
#'
#' @return Saves the plot as a JPEG file in the specified directory.
#'
#' @export

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
    gate.color = "darkgoldenrod1",
    switch.n = 1e4,
    raster.bins = 256L
  ) {

  # ---------------------------------------------------------------------------
  # 1. Downsample and clip to axis limits
  # ---------------------------------------------------------------------------

  n.points <- nrow( gate.data )
  if ( n.points > max.points ) {
    set.seed( asp$bird.seed )
    gate.data <- gate.data[ sample( seq_len( n.points ), max.points ), , drop = FALSE ]
    n.points  <- max.points
  }

  # clip to preset scatter limits before binning
  gate.data[ , 1 ] <- pmin( gate.data[ , 1 ], asp$scatter.data.max.x )
  gate.data[ , 2 ] <- pmin( gate.data[ , 2 ], asp$scatter.data.max.y )

  # ---------------------------------------------------------------------------
  # 2. Axis geometry (computed once, reused by raster + contours + scales)
  # ---------------------------------------------------------------------------

  x.limits <- c( asp$scatter.data.min.x, asp$scatter.data.max.x )
  y.limits <- c( asp$scatter.data.min.y, asp$scatter.data.max.y )
  x.breaks <- seq( asp$scatter.data.min.x, asp$scatter.data.max.x, asp$data.step )
  y.breaks <- seq( asp$scatter.data.min.y, asp$scatter.data.max.y, asp$data.step )
  x.labels <- paste0( round( x.breaks / 1e6, 1 ), "e6" )
  y.labels <- paste0( round( y.breaks / 1e6, 1 ), "e6" )

  x.lab <- names( which( scatter.and.channel.label == gate.marker[ 1 ] ) )
  y.lab <- names( which( scatter.and.channel.label == gate.marker[ 2 ] ) )

  # ---------------------------------------------------------------------------
  # 3. Rasterised scatter layer  (replaces geom_scattermore)
  # ---------------------------------------------------------------------------

  # pixel break vectors spanning the axis limits
  x_breaks_r <- seq( x.limits[1], x.limits[2], length.out = raster.bins + 1L )
  y_breaks_r <- seq( y.limits[1], y.limits[2], length.out = raster.bins + 1L )

  # bin points into a count matrix —— use C++ path when available
  scatter.mat <- if (
    requireNamespace( "AutoSpectralRcpp", quietly = TRUE ) &&
    "bin_matrix_cpp" %in% ls( getNamespace( "AutoSpectralRcpp" ) )
  ) {
    xi <- pmax( 1L, pmin( raster.bins, findInterval( gate.data[, 1], x_breaks_r ) ) )
    yi <- pmax( 1L, pmin( raster.bins, findInterval( gate.data[, 2], y_breaks_r ) ) )
    m  <- matrix( 0L, nrow = raster.bins, ncol = raster.bins )
    idx <- (xi - 1L) * raster.bins + yi
    counts <- tabulate( idx, nbins = raster.bins * raster.bins )
    matrix( counts, nrow = raster.bins, ncol = raster.bins )
  } else {
    xi <- pmax( 1L, pmin( raster.bins, findInterval( gate.data[, 1], x_breaks_r ) ) )
    yi <- pmax( 1L, pmin( raster.bins, findInterval( gate.data[, 2], y_breaks_r ) ) )
    m  <- matrix( 0L, nrow = raster.bins, ncol = raster.bins )
    idx <- (xi - 1L) * raster.bins + yi
    counts <- tabulate( idx, nbins = raster.bins * raster.bins )
    matrix( counts, nrow = raster.bins, ncol = raster.bins )
  }

  # map log-counts to greyscale: empty = NA (transparent), data = grey ramp
  scatter.log  <- log1p( scatter.mat )
  scatter.max  <- max( scatter.log )
  if ( scatter.max > 0 ) {
    scatter.norm <- scatter.log / scatter.max
    grey.vals    <- grDevices::grey( 1 - scatter.norm * 0.85 )  # light → dark
    grey.vals[ scatter.mat == 0L ] <- NA_character_
  } else {
    grey.vals <- matrix( NA_character_, nrow = raster.bins, ncol = raster.bins )
  }
  # annotation_raster: rows = y (bottom to top after flip), cols = x
  scatter.raster <- matrix( grey.vals, nrow = raster.bins, ncol = raster.bins )
  scatter.raster <- scatter.raster[ raster.bins:1, , drop = FALSE ]   # flip Y

  # ---------------------------------------------------------------------------
  # 4. KDE density contours (always fast_kde2d_cpp if available)
  # ---------------------------------------------------------------------------

  density.df     <- NULL
  density.breaks <- NULL

  if ( n.points >= switch.n ) {
    bw <- apply( gate.data, 2, MASS::bandwidth.nrd )

    gate.bound.density <- if (
      requireNamespace( "AutoSpectralRcpp", quietly = TRUE ) &&
      "fast_kde2d_cpp" %in% ls( getNamespace( "AutoSpectralRcpp" ) )
    ) {
      AutoSpectralRcpp::fast_kde2d_cpp(
        x        = gate.data[ , 1 ],
        y        = gate.data[ , 2 ],
        n        = 100,
        h        = bw * 0.1,
        x_limits = x.limits,
        y_limits = y.limits
      )
    } else {
      MASS::kde2d(
        x = gate.data[ , 1 ],
        y = gate.data[ , 2 ],
        n = 60,
        h = bw * 0.8
      )
    }

    # format density — direct construction avoids expand.grid overhead
    n.grid     <- length( gate.bound.density$x )
    density.df <- data.frame(
      x = rep( gate.bound.density$x, times = n.grid ),
      y = rep( gate.bound.density$y, each  = n.grid ),
      z = as.vector( gate.bound.density$z )
    )
    density.df$z[ is.na( density.df$z ) ] <- 0

    max.z          <- max( density.df$z )
    density.breaks <- seq( 0.05 * max.z, max.z, length.out = 11 )
  }

  # ---------------------------------------------------------------------------
  # 5. Gate boundary data frame
  # ---------------------------------------------------------------------------

  gate.boundary.ggp <- data.frame(
    x = c( gate.boundary$x, gate.boundary$x[ 1 ] ),
    y = c( gate.boundary$y, gate.boundary$y[ 1 ] )
  )
  gate.boundary.ggp$x <- pmin( gate.boundary.ggp$x, asp$scatter.data.max.x )
  gate.boundary.ggp$y <- pmin( gate.boundary.ggp$y, asp$scatter.data.max.y )

  # ---------------------------------------------------------------------------
  # 6. Build the plot
  # ---------------------------------------------------------------------------

  gate.plot <- ggplot() +
    # white panel background so NA raster pixels show as white
    annotate(
      "rect",
      xmin = x.limits[1], xmax = x.limits[2],
      ymin = y.limits[1], ymax = y.limits[2],
      fill = "white", colour = NA
    ) +
    # rasterised scatter
    annotation_raster(
      scatter.raster,
      xmin = x.limits[1], xmax = x.limits[2],
      ymin = y.limits[1], ymax = y.limits[2],
      interpolate = FALSE
    )

  # density contours on top of scatter (only when enough points)
  if ( !is.null( density.df ) ) {
    gate.plot <- gate.plot +
      geom_contour_filled(
        data        = density.df,
        aes( x = x, y = y, z = z ),
        breaks      = density.breaks,
        alpha       = 1,
        inherit.aes = FALSE,
        na.rm       = TRUE
      )
  }

  # gate boundary line
  gate.plot <- gate.plot +
    geom_path(
      data      = gate.boundary.ggp,
      aes( x, y ),
      color     = gate.color,
      linewidth = asp$figure.gate.line.size
    ) +
    scale_x_continuous(
      name   = x.lab,
      breaks = x.breaks,
      labels = x.labels,
      expand = expansion( asp$figure.gate.scale.expand )
    ) +
    scale_y_continuous(
      name   = y.lab,
      breaks = y.breaks,
      labels = y.labels,
      expand = expansion( asp$figure.gate.scale.expand )
    ) +
    coord_cartesian( xlim = x.limits, ylim = y.limits ) +
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

  if ( !is.null( density.df ) ) {
    viridis.colors <- c(
      "magma", "inferno", "plasma", "viridis",
      "cividis", "rocket", "mako", "turbo"
    )

    if ( color.palette %in% viridis.colors ) {
      gate.plot <- gate.plot + scale_fill_viridis_d( option = color.palette )
    } else {
      n.bins         <- max( 1L, length( density.breaks ) - 1L )
      rainbow.palette <- grDevices::colorRampPalette( asp$density.palette.base.color )( n.bins )
      gate.plot      <- gate.plot + scale_fill_manual( values = rainbow.palette )
    }
  }

  # ---------------------------------------------------------------------------
  # 8. Save
  # ---------------------------------------------------------------------------

  ggsave(
    file.path( asp$figure.gate.dir, sprintf( "%s.jpg", samp ) ),
    plot      = gate.plot,
    device    = ragg::agg_jpeg,
    width     = asp$figure.width,
    height    = asp$figure.height,
    limitsize = FALSE
  )

}
