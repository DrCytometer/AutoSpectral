# define_gate_density.r

#' @title Define Gate by Density
#'
#' @description
#' Performs gating on scatter parameters based on a density search.
#'
#' The gating proceeds in three steps:
#'
#' - Defines bounds by data trimming
#' - Defines a region around the target maximum found within the bounds
#' - Defines a gate around the target maximum, only within that region
#'
#' The method uses numerical search of maxima over estimated densities and
#' Voronoi tessellations to improve density estimation around maxima.
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
#' many cells will be selected from the specified FCS files. The default, `NULL`,
#' will revert to `asp$gate.downsample.n`. Around 1e5 is recommended.
#' @param grid.n Numeric, default `100`. The binning grid for the kernel density
#' estimation. If `n.cells` is very low, you may wish to lower this number to
#' compress the search space.
#' @param bandwidth.factor Numeric, default `NULL` will revert to
#' `asp$gate.bound.density.bw.factor` (usually 1). A multiplier for the
#' bandwidth for the kernel density estimation. Larger numbers will smooth the
#' density, reducing discrimination between peaks in the density (such as between
#' live and dead cells).
#' @param target.pop Numeric 0-n, default `NULL` will revert to
#' `asp$gate.bound.density.max.target` (usually 1). In the output plots, the
#' numbers correspond to density peaks. Whichever number is set here will select
#' the population 1 less than that number as the target for defining the gate
#' region.
#' @param neighbors Number of neighboring events to consider when finding local
#' density maxima in the kernel density estimation. Default `NULL` will revert
#' to `asp$gate.bound.density.neigh.size`, which is usually 3.
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
#' @param fsc.search.min The minimum bound for the search space for dense event
#' populations for defining the gate on the x-axis. The default `NULL` reverts to
#' `asp$default.gate.param$region.factor.x.low`, which is usually 0.05 or 5% of
#' the FSC axis.
#' @param fsc.search.max The maximum bound for the search space for dense event
#' populations for defining the gate on the x-axis. The default `NULL` reverts to
#' `asp$default.gate.param$region.factor.x.high`, which is usually 0.8 or 80% of
#' the FSC axis.
#' @param ssc.search.min The minimum bound for the search space for dense event
#' populations for defining the gate on the y-axis. The default `NULL` reverts to
#' `asp$default.gate.param$region.factor.y.low`, which is usually 0.05 or 5% of
#' the SSC axis.
#' @param ssc.search.max The maximum bound for the search space for dense event
#' populations for defining the gate on the y-axis. The default `NULL` reverts to
#' `asp$default.gate.param$region.factor.y.high`, which is usually 0.8 or 80% of
#' the SSC axis.#'
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
#' * [define.gate.landmarks()]
#' * [do.gate()]
#' * [gate.define.plot()]
#'
#' @return A set of points describing the gate boundary.
#'
#' @export
#'
#' @references Roca, Carlos P et al. "AutoSpill is a principled framework that
#' simplifies the analysis of multichromatic flow cytometry data" \emph{Nature
#' Communications} 12 (2890) 2021.


define.gate.density <- function(
    control.file,
    control.dir,
    asp,
    gating.params = NULL,
    n.cells = NULL,
    grid.n = NULL,
    bandwidth.factor = NULL,
    target.pop = NULL,
    neighbors = NULL,
    fsc.channel = NULL,
    ssc.channel = NULL,
    fsc.lims = NULL,
    ssc.lims = NULL,
    fsc.search.min = NULL,
    fsc.search.max = NULL,
    ssc.search.min = NULL,
    ssc.search.max = NULL,
    output.dir = "./figure_gate",
    gate.name = "cell_gate",
    filename = "density_gate_definition_",
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

  # create output folder
  asp$figure.gate.dir <- output.dir
  if ( !dir.exists( asp$figure.gate.dir ) ) dir.create( asp$figure.gate.dir )

  # update figure dimensions
  asp$figure.width <- width
  asp$figure.height <- height

  # set scatter channels
  if ( is.null( fsc.channel ) ) fsc.channel <- asp$default.scatter.parameter[ 1 ]
  if ( is.null( ssc.channel ) ) ssc.channel <- asp$default.scatter.parameter[ 2 ]

  ### update preferences with user-provided data
  # determine suffix based on sample type
  suffix <- ifelse( sample.type == "beads", ".beads", ".cells" )

  # create a named list of parameters to update
  updates <- list(
    n.cells = "gate.downsample.n",
    target.pop = "gate.bound.density.max.target",
    bandwidth.factor = "gate.bound.density.bw.factor",
    neighbors = "gate.bound.density.neigh.size",
    grid.n = "gate.bound.density.grid.n"
  )

  # apply updates
  for ( arg.name in names( updates ) ) {
    val <- get( arg.name )
    if ( !is.null( val ) ) {
      asp.key <- paste0( updates[[ arg.name ]], suffix )
      asp[[ asp.key ]] <- val
    }
  }

  # update gate search region if needed
  target <- paste0( "default.gate.param", suffix )
  if ( !is.null( fsc.lims ) ) {
    asp[[ target ]]$region.factor.x.low  <- fsc.lims[ 1 ]
    asp[[ target ]]$region.factor.x.high <- fsc.lims[ 2 ]
  }
  if ( !is.null( ssc.lims ) ) {
    asp[[ target ]]$region.factor.y.low  <- ssc.lims[ 1 ]
    asp[[ target ]]$region.factor.y.high <- ssc.lims[ 2 ]
  }
  if ( !is.null( fsc.search.min ) ) asp[[ target ]]$region.factor.x.low  <- fsc.search.min
  if ( !is.null( fsc.search.max ) ) asp[[ target ]]$region.factor.x.high <- fsc.search.max
  if ( !is.null( ssc.search.min ) ) asp[[ target ]]$region.factor.y.low  <- ssc.search.min
  if ( !is.null( ssc.search.max ) ) asp[[ target ]]$region.factor.y.high <- ssc.search.max

  message( paste( "Defining", gate.name, "gate using", length( file.idx ), "FCS files." ) )

  # read in up to 100k events, split between files
  max.event.n <- ceiling( 1e5 / length( files.to.use ) )
  gate.data <- lapply( files.to.use, function( f ) {
    # read in scatter data
    scatter.data <- readFCS(
      file.path( control.dir, f ),
      return.keywords = FALSE
    )[ , c( fsc.channel, ssc.channel ) ]

    event.n <- nrow( scatter.data )

    # set seed and sample
    if ( event.n > max.event.n ) {
      set.seed( 42 )
      scatter.data <- scatter.data[ sample( event.n, max.event.n ), ]
    }

    return( scatter.data )
  } )

  gate.data <- do.call( rbind, gate.data )

  gate <- do.gate.test(
    gate.data = gate.data,
    viability.gate = viability,
    large.gate = large.gate,
    samp = paste0( filename, gate.name ),
    scatter.and.channel.label = c( fsc.channel, ssc.channel ),
    control.type = sample.type,
    asp = asp,
    color.palette = color.palette,
    max.points = points.to.plot,
    gate.color = boundary.color
  )
  names( gate ) <- gate.name

  # generate timestamp
  ts <- format( Sys.time(), "%Y%m%d_%H%M%S" )

  # save the gate as an R object for later use
  saveRDS(
    gate,
    file = file.path( output.dir, paste0( filename, "_", gate.name, "_", ts, ".rds" ) )
  )
  saveRDS(
    asp,
    file = file.path( output.dir, paste0( filename, "_params_", gate.name, "_", ts, ".rds" ) )
  )

  message( paste( "Files saved with timestamp suffix:", ts ) )

  return( gate )
}
