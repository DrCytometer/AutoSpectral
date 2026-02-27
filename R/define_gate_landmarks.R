# define_gate_landmarks.r

#' @title Define Gate by Landmarks
#'
#' @description
#' Performs gating on scatter parameters using events selected on the basis of
#' high fluorescence intensity in the expected peak channel in the single-stained
#' control samples. This approach allows for identification of cells fairly
#' reliably by "landmarks", especially when using well-defined, abundant
#' populations such a T cells (e.g., CD3, CD4, CD8 single-stained samples),
#' monocytes (CD14-stained sample) or neutrophils (CD66, CD33).
#'
#' @importFrom MASS kde2d bandwidth.nrd
#'
#' @param control.file File path and name for the CSV file defining the single-
#' color control file names, fluorophores they represent, marker names, peak
#' channels, and gating requirements.
#' @param control.dir File path to the single-stained control FCS files.
#' @param asp The AutoSpectral parameter list defined using
#' `get.autospectral.param`.
#' @param gating.params Previously saved gating parameters. Load in the .rds file
#' using `readRDS` and pass the result here if you wish to replicate a previous
#' run. This is essentially just an updated version of `asp` containing any
#' modifications based on information passed as arguments to `define.gate.landmarks`.
#' @param n.cells The number of cells to use for defining the gate boundary. This
#' many cells will be selected from the peak channel (brightest first) in the
#' single-color controls. For example, if you set `200` and marked files such as
#' `CD3-PE.fcs` and `CD19-FITC.fcs` as `gate.define=TRUE` in the control file,
#' the brightest 200 events in the YG1 channel from the CD3-PE file and the
#' brightest 200 events in the B1 channel for the CD19-FITC file would be used
#' to define the gate.
#' @param percentile Numeric 1 - 100, default `70`. The percentile cutoff for
#' density in the scatter to use for defining the gate. For example, a value of
#' `50` would take the 50% of cells closest to the density peak, more or less.
#' Smaller numbers will define a tighter gate.
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
#' @param fsc.lims Numeric vector. Limits for searching and plotting the FSC. The
#' default `NULL` uses c`c(asp$scatter.data.min.x, asp$scatter.data.max.x)`.
#' @param ssc.lims Numeric vector. Limits for searching and plotting the SSC. The
#' default `NULL` uses c`c(asp$scatter.data.min.y, asp$scatter.data.max.y)`.
#' @param output.dir File path where you want to save the results. Default is
#' `./figure_gate`.
#' @param gate.name Character, name for the gate. Useful for distinguishing gates
#' when you have multiple types. Default is `cell_gate`.
#' @param filename Character, name for the output files. Default is `gate_definition`.
#' @param color.palette Optional character string defining the viridis color
#' palette to be used for the fluorophore traces. Default is `plasma`. Use `rainbow`
#' to be similar to FlowJo or SpectroFlo. Other options are the viridis color
#' options: `magma`, `inferno`, `plasma`, `viridis`, `cividis`, `rocket`, `mako`
#' and `turbo`.
#' @param boundary.color Color for the gate boundary line on the plot. Default
#' is `black`.
#' @param points.to.plot Numeric, default `1e5`. Maximum number of points to show
#' on the plot. More points will take longer, but you really shouldn't have even
#' close to this number when defining the gate with landmarks.
#' @param width Numeric, default `5`. Width of the saved plot.
#' @param height Numeric, default `5`. Height of the saved plot.
#'
#' @seealso
#' * [tune.gate()]
#' * [define.gate.density()]
#' * [gate.define.plot()]
#'
#' @return A set of points describing the gate boundary.
#'
#' @export
#'
#' @references
#' Laniewski, Nathan. \emph{flowstate}.
#' \url{https://github.com/nlaniewski/flowstate}
#'

define.gate.landmarks <- function(
    control.file,
    control.dir,
    asp,
    gating.params = NULL,
    n.cells = 2000,
    percentile = 70,
    grid.n = 100,
    bandwidth.factor = 1,
    fsc.channel = NULL,
    ssc.channel = NULL,
    fsc.lims = NULL,
    ssc.lims = NULL,
    output.dir = "./figure_gate",
    gate.name = "cell_gate",
    filename = "landmark_gate_definition_",
    color.palette = "plasma",
    boundary.color = "black",
    points.to.plot = 1e5,
    width = 5,
    height = 5
) {

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

  # fill in blanks
  control.table$is.viability[ is.na( control.table$is.viability ) ] <- FALSE
  control.table$large.gate[ is.na( control.table$large.gate ) ] <- FALSE

  # identify which files have been requested for drawing the gate
  file.idx <- which( control.table$gate.name == gate.name & control.table$gate.define == TRUE )
  files.to.use <- control.table$filename[ file.idx ]

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

  # if the user has provided gating.params (asp), replace asp and start here as base
  if ( !is.null( gating.params ) ) asp <- gating.params

  if ( !is.null( fsc.lims ) ) {
    asp$scatter.data.min.x <- fsc.lims[ 1 ]
    asp$scatter.data.max.x <- fsc.lims[ 2 ]
  }
  if ( !is.null( ssc.lims ) ) {
    asp$scatter.data.min.y <- ssc.lims[ 1 ]
    asp$scatter.data.max.y <- ssc.lims[ 2 ]
  }

  # create output folder
  asp$figure.gate.dir <- output.dir
  if ( !dir.exists( asp$figure.gate.dir ) ) dir.create( asp$figure.gate.dir )

  # update figure dimensions
  asp$figure.width <- width
  asp$figure.height <- height

  # set scatter channels
  if ( is.null( fsc.channel ) ) fsc.channel <- asp$default.scatter.parameter[ 1 ]
  if ( is.null( ssc.channel ) ) ssc.channel <- asp$default.scatter.parameter[ 2 ]

  percentile <- percentile * 0.01 # convert to percentage (quantile)

  message( paste( "Defining", gate.name, "gate using", length( file.idx ), "FCS files." ) )

  files.channels <- control.table$channel[ file.idx ]
  names( files.channels ) <- control.table$filename[ file.idx ]
  # exclude unstained samples--shouldn't be any at this point
  files.channels <- files.channels[
    !grepl( "AF|negative", control.table$fluorophore, ignore.case = TRUE ) ]
  files.channels <- files.channels[ !is.na( files.channels ) ]

  pooled.scatter.data <- do.call( rbind, lapply( names( files.channels ), function( f ) {
    get.top.events(
      file.path( control.dir, f ),
      files.channels[[ f ]],
      n.cells,
      scatter.param = c( fsc.channel, ssc.channel )
    )
  } ) )

  # kernel density estimation
  bw <- apply( pooled.scatter.data, 2, bandwidth.nrd )
  dens <- MASS::kde2d(
    pooled.scatter.data[, fsc.channel ],
    pooled.scatter.data[, ssc.channel ],
    h = bandwidth.factor * bw,
    n = grid.n
  )

  # can add namespace switch to use kde2d from AutoSpectralRcpp if available
  # and if this proves slow

  # find the density threshold
  z.sort <- sort( dens$z, decreasing = TRUE)
  cumulative.dens <- cumsum( z.sort ) / sum( z.sort )
  threshold <- z.sort[ which.min( abs( cumulative.dens - percentile ) ) ]

  # extract the contour line at that threshold
  contour.lines <- grDevices::contourLines( dens$x, dens$y, dens$z, levels = threshold )

  # take the longest contour (in case there are multiple)
  main.contour <- contour.lines[[ which.max( sapply( contour.lines, function(l) length( l$x ) ) ) ]]
  gate.coords <- cbind( main.contour$x, main.contour$y )

  # define a smooth convex hull around those coordinates
  gate.hull <- gate.coords[ grDevices::chull( gate.coords ), ]
  gate.hull <- rbind( gate.hull, gate.hull[ 1, ] )

  gate.boundary <- list(
    x = gate.hull[ , 1 ],
    y = gate.hull[ , 2 ]
  )
  gate.population <- list( boundary = gate.boundary )
  names( gate.boundary ) <- gate.name

  # format regions for plotting
  gate.bound <- list(
    density = dens,
    voronoi = NULL,
    density.max = NULL,
    x.low = asp$scatter.data.min.x,
    x.high = asp$scatter.data.max.x,
    y.low = asp$scatter.data.min.y,
    y.high = asp$scatter.data.max.y
  )

  # plotting
  tryCatch(
    expr = {
      gate.define.plot(
        samp = paste0( filename, gate.name ),
        gate.data = pooled.scatter.data,
        gate.marker = c( fsc.channel, ssc.channel ),
        gate.bound = gate.bound,
        gate.region = NULL,
        gate.population = gate.population,
        scatter.and.channel.label = c( fsc.channel, ssc.channel ),
        asp = asp,
        color.palette = color.palette,
        max.points = points.to.plot,
        gate.color = boundary.color
      )
    },
    error = function( e ) {
      message( "Error in plotting gate: ", e$message )
      return( NULL )
    }
  )

  # generate timestamp
  ts <- format( Sys.time(), "%Y%m%d_%H%M%S" )

  # save the gate as an R object for later use
  saveRDS(
    gate.boundary,
    file = file.path( output.dir, paste0( filename, "_", gate.name, "_", ts, ".rds" ) )
  )
  saveRDS(
    asp,
    file = file.path( output.dir, paste0( filename, "_params_", gate.name, "_", ts, ".rds" ) )
  )

  message( paste( "Files saved with timestamp suffix:", ts ) )

  return( gate.boundary )

}
