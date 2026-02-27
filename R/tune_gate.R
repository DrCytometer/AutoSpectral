
#' @title Define Gate by Density
#'
#' @description
#' A short description...
#'
#' @importFrom cowplot save_plot plot_grid
#' @importFrom ggplot2 ggplot aes scale_x_continuous scale_y_continuous element_blank
#' @importFrom ggplot2 theme_bw theme element_line geom_path after_stat coord_cartesian
#' @importFrom ggplot2 element_text element_rect margin expansion ggsave labs
#' @importFrom ggplot2 scale_fill_viridis_d geom_contour_filled scale_fill_manual
#' @importFrom scattermore geom_scattermore
#' @importFrom MASS kde2d
#'
#' @param control.file File path and name for the CSV file defining the single-
#' color control file names, fluorophores they represent, marker names, peak
#' channels, and gating requirements.
#' @param control.dir File path to the single-stained control FCS files.
#' @param asp The AutoSpectral parameter list defined using
#' `get.autospectral.param`.
#' @param n.cells The number of cells to use for defining the gate boundary. The
#' default is `c(100, 500, 2000)` and will test each of those values. This many
#' cells will be selected from the peak channel (brightest first) in the single-
#' color controls. For example, if you set `200` and marked files such as
#' `CD3-PE.fcs` and `CD19-FITC.fcs` as `gate.define=TRUE` in the control file,
#' the brightest 200 events in the YG1 channel from the CD3-PE file and the
#' brightest 200 events in the B1 channel for the CD19-FITC file would be used
#' to define the gate.
#' @param percentiles Numeric 1 - 100, default `c(30, 50, 70)`. The percentile
#' cutoffs to test for density in the scatter to use for defining the gate. For
#' example, a value of `50` would take the 50% of cells closest to the density
#' peak, more or less. Smaller numbers will define a tighter gate.
#' @param grid.n Numeric, default `100`. The binning grid for the kernel density
#' estimation. If `n.cells` is very low, you may wish to lower this number to
#' compress the search space.
#' @param bandwidth.factor Numeric, default `1`. A multiplier for the bandwidth
#' for the kernel density estimation. Larger numbers will smooth the density,
#' reducing discrimination between peaks in the density (such as between live
#' and dead cells).
#' @param fsc.channel Channel to use for Forward Scatter. Default `NULL` will
#' use the `asp$default.scatter.parameter[1]`, which is appropriate for your
#' machine.
#' @param ssc.channel Channel to use for Side Scatter. Default `NULL` will
#' use the `asp$default.scatter.parameter[2]`, which is appropriate for your
#' machine. For machines with multiple side scatter measurements, you can change
#' this.
#' @param fsc.lims Numeric vector. Limits for plotting the FSC. The default `NULL`
#' uses c`c(asp$scatter.data.min.x, asp$scatter.data.max.x)`.
#' @param ssc.lims Numeric vector. Limits for plotting the SSC. The default `NULL`
#' uses c`c(asp$scatter.data.min.y, asp$scatter.data.max.y)`.
#' @param output.dir File path where you want to save the results. Default is
#' `./figure_gate_tuning`.
#' @param gate.name Character, name for the gate. Useful for distinguishing gates
#' when you have multiple types. Default is `cell_gate`.
#' @param filename Character, name for the output files. Default is `gate_tuning`.
#' @param color.palette Optional character string defining the viridis color
#' palette to be used for the fluorophore traces. Default is `plasma`. Use `rainbow`
#' to be similar to FlowJo or SpectroFlo. Other options are the viridis color
#' options: `magma`, `inferno`, `plasma`, `viridis`, `cividis`, `rocket`, `mako`
#' and `turbo`.
#' @param boundary.color Color for the gate boundary line on the plot. Default
#' is `darkgoldenrod1`.
#' @param points.to.plot Numeric, default `1e5`. Maximum number of points to show
#' on the plot. More points will take longer, but you really shouldn't have even
#' close to this number when defining the gate with landmarks.
#' @param width Numeric, default `4`. Width of the saved plot.
#' @param height Numeric, default `4`. Height of the saved plot.
#'
#' @seealso#'
#' * [define.gate.landmarks()]
#' * [define.gate.density()]
#' * [gate.sample.plot()]
#'
#' @return A combined plot of all gates. Saves plots of the gates defined using
#' combinations of the specified parameters.
#'
#' @export
#'
#' @references
#' Laniewski, Nathan. \emph{flowstate}.
#' \url{https://github.com/nlaniewski/flowstate}
#'

tune.gate <- function(
    control.file,
    control.dir,
    asp,
    n.cells = c(100, 500, 2000),
    percentiles = c(30, 50, 70),
    grid.n = 100,
    bandwidth.factor = 1,
    fsc.channel = NULL,
    ssc.channel = NULL,
    fsc.lims = NULL,
    ssc.lims = NULL,
    output.dir = "./figure_gate_tuning",
    gate.name = "cell_gate",
    filename = "gate_tuning",
    color.palette = "mako",
    boundary.color = "red",
    points.to.plot = 1e5,
    width = 9,
    height = 10
) {
  if ( is.null( fsc.channel ) ) fsc.channel <- asp$default.scatter.parameter[ 1 ]
  if ( is.null( ssc.channel ) ) ssc.channel <- asp$default.scatter.parameter[ 2 ]

  # check that the inputs have the right structure and type
  stopifnot( is.numeric( n.cells ) )
  stopifnot( is.numeric( percentiles ) )

  # screen out duplicates
  n.cells <- unique( n.cells )
  percentiles <- unique( percentiles )

  # read control info
  control.table <- utils::read.csv(
    control.file,
    stringsAsFactors = FALSE,
    strip.white = TRUE
  )

  # trim white space, convert blanks to NAs
  control.table[] <- lapply( control.table, function( x ) {
    if ( is.character( x ) ) {
      x <- trimws( x )
      x[ x == "" ] <- NA
      x
    } else x
  } )

  # identify which files are requested to be used to draw the gate
  file.idx <- which( control.table$gate.name == gate.name & control.table$gate.define == TRUE )
  files.to.use <- control.table$filename[ file.idx ]

  # check that these are all the same type of sample
  sample.type <- control.table$control.type[ file.idx ]

  # check for unstained samples, remove if present and warn
  sample.fluors <- control.table$fluorophore[ file.idx ]
  unstained.samples <- grepl( "AF|negative", sample.fluors, ignore.case = TRUE )

  if ( any( unstained.samples ) ) {
    file.idx <- file.idx[ !unstained.samples ]
    unstained.samples <- control.table$filename[ file.idx[ unstained.samples ] ]
    files.to.use <- control.table$filename[ file.idx ]

    message(
      paste(
        "Unstained samples were included in the files to be used for gating.",
        "These have been removed:",
        unstained.samples
      )
    )
  }

  # re-check to be sure we still have samples
  sample.n <- length( file.idx )

  if ( sample.n < 1 ) {
    stop( "No samples left after excluding unstained samples. Please check control file and try again.",
          call. = FALSE )
  }

  # check that everything matches (same type of particles, same gating definitions)
  check.consistency <- function( vec, var.name ) {
    u <- unique( vec )
    if ( length( u ) > 1 ) {
      stop( paste( "Inconsistent values for", var.name, "found in selected samples." ), call. = FALSE )
    }
    return( u )
  }
  sample.type <- check.consistency( control.table$control.type[ file.idx ], "control.type" )
  large.gate  <- check.consistency( control.table$large.gate[ file.idx ], "large.gate" )
  viability   <- check.consistency( control.table$is.viability[ file.idx ], "is.viability" )

  # read in up to 100k events, split between files
  message( "Combining scatter data from all files for plotting" )
  max.event.n <- ceiling( 1e5 / sample.n )

  plot.data <- do.call( rbind, lapply( files.to.use, function( f ) {
    # read in scatter data
    scatter.data <- readFCS(
      file.path( control.dir, f ),
      return.keywords = FALSE
    )[ , c( fsc.channel, ssc.channel ) ]

    # set seed and sample
    if ( nrow( scatter.data ) > max.event.n ) {
      set.seed( 42 )
      scatter.data <- scatter.data[ sample( nrow( scatter.data ), max.event.n ), ]
    }
    return( scatter.data )
  } ) )

  # match file names and peak channels
  files.channels <- control.table$channel[ file.idx ]
  names( files.channels ) <- control.table$filename[ file.idx ]
  # exclude unstained samples--shouldn't be any at this point
  files.channels <- files.channels[
    !grepl( "AF|negative", control.table$fluorophore, ignore.case = TRUE ) ]
  files.channels <- files.channels[ !is.na( files.channels ) ]

  # read in up to maximum requested cell number of pos events for each control
  max.cells <- max( n.cells )

  pooled.scatter.list <- lapply( names( files.channels ), function( f ) {
    get.top.events(
      file.path( control.dir, f ),
      files.channels[[ f ]],
      max.cells,
      scatter.param = c( fsc.channel, ssc.channel )
    )
  } )
  names( pooled.scatter.list ) <- names( files.channels )

  # set up for plotting
  message( "Generating scatter plot" )

  if ( !is.null( fsc.lims ) ) {
    asp$scatter.data.min.x <- fsc.lims[ 1 ]
    asp$scatter.data.max.x <- fsc.lims[ 2 ]
  }
  if ( !is.null( ssc.lims ) ) {
    asp$scatter.data.min.y <- ssc.lims[ 1 ]
    asp$scatter.data.max.y <- ssc.lims[ 2 ]
  }

  # restrict data to preset limits for plotting
  plot.data[ , 1 ] <- pmin( plot.data[ , 1 ], asp$scatter.data.max.x )
  plot.data[ , 2 ] <- pmin( plot.data[ , 2 ], asp$scatter.data.max.y )

  # calculate density once for plotting
  bw <- apply( plot.data, 2, bandwidth.nrd ) * bandwidth.factor

  if ( requireNamespace( "AutoSpectralRcpp", quietly = TRUE ) &&
       "fast_kde2d_cpp" %in% ls( getNamespace( "AutoSpectralRcpp" ) ) &&
       nrow( plot.data ) > 1e4 ) {
    # use C++ function to get density
    gate.bound.density <- AutoSpectralRcpp::fast_kde2d_cpp(
      x = plot.data[ , 1 ],
      y = plot.data[ , 2 ],
      n = grid.n,
      h = bw * 0.1,
      x_limits = range( plot.data[ , 1 ] ),
      y_limits = range( plot.data[ , 2 ] )
    )
  } else {
    # use slower MASS call
    gate.bound.density <- MASS::kde2d(
      x = plot.data[ , 1 ],
      y = plot.data[ , 2 ],
      n = grid.n,
      h = bw
    )
  }

  # format the density for plotting
  density.df <- data.frame( expand.grid(
    x = gate.bound.density$x,
    y = gate.bound.density$y
  ) )
  density.df$z <- as.vector( gate.bound.density$z )
  density.df <- density.df[ !is.na( density.df$z ), ]
  density.df <- density.df[ !duplicated( density.df[ , c( "x", "y" ) ] ), ]
  max.z <- max( density.df$z, na.rm = TRUE )
  density.breaks <- seq( 0.05 * max.z, max.z, length.out = 11 )

  # convert points to data frame for plotting
  gate.data.ggp <- data.frame(
    x = plot.data[ , 1 ],
    y = plot.data[ , 2 ] )

  # get axis labels
  x.lab <- fsc.channel
  y.lab <- ssc.channel

  # create axes labels
  x.limits <- c( asp$scatter.data.min.x, asp$scatter.data.max.x )
  x.breaks <- seq( asp$scatter.data.min.x, asp$scatter.data.max.x, asp$data.step )
  x.labels <- paste0( round( x.breaks / 1e6, 1 ), "e6" )
  y.limits <- c( asp$scatter.data.min.y, asp$scatter.data.max.y )
  y.breaks <- seq( asp$scatter.data.min.y, asp$scatter.data.max.y, asp$data.step )
  y.labels <- paste0( round( y.breaks / 1e6, 1 ), "e6" )

  # set up base plot once
  base.plot <- ggplot( gate.data.ggp, aes( x, y ) ) +
    geom_scattermore(
      pointsize = asp$figure.gate.point.size,
      color = "black",
      alpha = 1,
      na.rm = TRUE
    ) +
    geom_contour_filled(
      data = density.df,
      aes( x = x, y = y, z = z ),
      breaks = density.breaks,
      alpha = 1,
      inherit.aes = FALSE,
      na.rm = TRUE
    ) +
    scale_x_continuous(
      name = x.lab,
      breaks = x.breaks,
      labels = x.labels,
      #limits = x.limits,
      expand = expansion( asp$figure.gate.scale.expand )
    ) +
    scale_y_continuous(
      name = y.lab,
      breaks = y.breaks,
      labels = y.labels,
      #limits = y.limits,
      expand = expansion( asp$figure.gate.scale.expand )
    ) +
    coord_cartesian( xlim = x.limits, ylim = y.limits ) +
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
    base.plot <- base.plot + scale_fill_viridis_d( option = color.palette )
  } else {
    n.bins <- max( 1, length( density.breaks ) - 1 )
    rainbow.palette <- grDevices::colorRampPalette( asp$density.palette.base.color )( n.bins )

    base.plot <- base.plot +
      scale_fill_manual( values = rainbow.palette )
  }

  print( base.plot )

  # store plots in a list
  plot.list <- list()
  counter <- 1

  # now, for loops for both cell.n and percentiles
  message( "Generating contour gates for each parameter combination" )

  for ( n in n.cells ) {
    # subset the landmark cells
    landmark.data <- do.call( rbind, lapply( pooled.scatter.list, function( f ) {
      f[ 1:min( n, nrow( f ) ), ]
    } ) )

    for ( p in percentiles ) {
      # kernel density estimation
      dens <- MASS::kde2d(
        landmark.data[, fsc.channel ],
        landmark.data[, ssc.channel ],
        h = bw,
        n = 60
      )

      # find the density threshold
      z.sort <- sort( dens$z, decreasing = TRUE)
      cumulative.dens <- cumsum( z.sort ) / sum( z.sort )
      threshold <- z.sort[ which.min( abs( cumulative.dens - p/100 ) ) ]

      # extract the contour line at that threshold
      contour.lines <- grDevices::contourLines( dens$x, dens$y, dens$z, levels = threshold )

      # take the longest contour (in case there are multiple)
      main.contour <- contour.lines[[ which.max( sapply( contour.lines, function(l) length( l$x ) ) ) ]]
      gate.coords <- cbind( main.contour$x, main.contour$y )

      # define a smooth convex hull around those coordinates
      gate.hull <- gate.coords[ grDevices::chull( gate.coords ), ]

      # convert boundary to data.frame for plotting
      gate.boundary.ggp <- as.data.frame( gate.hull )
      colnames( gate.boundary.ggp ) <- c( "x", "y" )
      gate.boundary.ggp <- rbind( gate.boundary.ggp, gate.boundary.ggp[ 1, ] )

      # ensure gate limits are drawn onscale
      gate.boundary.ggp$x <- pmin( gate.boundary.ggp$x, asp$scatter.data.max.x )
      gate.boundary.ggp$y <- pmin( gate.boundary.ggp$y, asp$scatter.data.max.y )

      # create plot of the data
      plot.title <- paste0( "n=", n, ", p=", p, "%" )
      gp <- base.plot +
        geom_path(
          data = gate.boundary.ggp,
          aes( x, y ),
          color = boundary.color,
          linewidth = asp$figure.gate.line.size
        ) +
        labs( title = plot.title )

      plot.list[[ counter ]] <- gp
      counter <- counter + 1
    }
  }

  # save all the plots together
  if ( !dir.exists( output.dir ) ) dir.create( output.dir, recursive = TRUE )
  ts <- format( Sys.time(), "%Y%m%d_%H%M%S" )
  combined.plot <- cowplot::plot_grid(
    plotlist = plot.list,
    ncol = length( percentiles )
  )

  final.filename <- file.path(
    output.dir,
    paste0( filename, "_", gate.name, "_", ts, ".pdf" )
  )
  cowplot::save_plot(
    final.filename,
    combined.plot,
    base_width = width,
    base_height = height
  )

  message( paste( "Tuning grid saved to:", final.filename ) )
  return( combined.plot )
}
