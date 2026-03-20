# define_flow_control.r

##' @title Define Flow Control
#'
#' @description
#' A complex function designed to convert the single-stained control FCS files
#' and the metadata `control.def.file` into a data structure for AutoSpectral.
#' It reads the metadata file, locates the corresponding FCS files, determines
#' the gates needed based on variables such as beads or cells, `large.gate`, and
#' `viability.gate`, and then creates gates for each combination. It imports
#' and gates the data in the FCS files and assigns various factors to track
#' the data. Parallel processing `parallel=TRUE` will likely speed up the run
#' considerably on Mac and Linux systems supporting forking, but will likely
#' not help much on Windows unless >10 cores are available.
#'
#' @param control.dir File path to the single-stained control FCS files.
#' @param control.def.file CSV file defining the single-color control file names,
#' fluorophores they represent, marker names, peak channels, and gating requirements.
#' @param asp The AutoSpectral parameter list defined using
#' `get.autospectral.param`.
#' @param gate Logical, default is `TRUE`, in which case, automated gating will
#' be performed. If `FALSE`, the FCS files will be imported without automatically
#' generated gates applied. That is, all data in the files will be used. This is
#' intended to allow the user to pre-gate the files in commercial software.
#' @param gating.system Character string selecting the automated gating system to
#' employ in defining the initial scatter gates for identifying cells in the FCS
#' files. Options are `landmarks` and `density`. The `density` option uses the
#' original gating system from AutoSpill, picking cell populations based on dense
#' regions on FSC/SSC. The default `landmarks` system picks out the cell regions
#' using the brightest events in the peak channel for the single-stained control(s).
#' This approach is generally more robust. A good way to use this is to utilize
#' `landmarks` in combination with specifying which single-stained control files
#' should be used for defining the gates. For example, abundant, bright cell
#' markers such as CD4 will reliably identify the lymphocyte region, as CD14 will
#' identify the monocyte region, etc. For more instructions, see the help pages
#' on GitHub.
#' @param gate.list Optional named list of gates. To use this, pre-define the
#' gates using `define.gate.landmarks()` and/or `define.gate.density()`, ensure
#' that the names of the gates correspond to the names in the `control.def.file`,
#' and ensure that the `gate.name` column has been filled in for the
#' `control.def.file`. Default `NULL` will revert to creating new gates.
#' @param parallel Logical, default is `FALSE`, in which case parallel processing
#' will not be used. Parallel processing will likely be faster when many small
#' files are read in. If the data is larger, parallel processing may not
#' accelerate the process much.
#' @param verbose Logical, default is `TRUE`. Set to `FALSE` to suppress messages.
#' @param threads Numeric, number of threads to use for parallel processing.
#' Default is `NULL` which will revert to `asp$worker.process.n` if
#' `parallel=TRUE`.
#' @param color.palette Optional character string defining the viridis color
#' palette to be used for the fluorophore traces. Use `rainbow`
#' to be similar to FlowJo or SpectroFlo. Other options are the viridis color
#' options: `magma`, `inferno`, `plasma`, `viridis`, `cividis`, `rocket`, `mako`
#' and `turbo`.
#'
#' @return A list (`flow.control`) with the following components:
#' - `filename`: Names of the single-color control files.
#' - `fluorophore`: Corresponding fluorophores used in the experiment.
#' - `antigen`: Corresponding markers used in the experiment.
#' - `control.type`: Type of control used (beads or cells).
#' - `universal.negative`: Corresponding universal negative for each control.
#' - `viability`: Logical factor; whether a control represents a viability marker.
#' - `large.gate`: Logical factor; large gate setting.
#' - `autof.channel.idx`: Index of the autofluorescence marker channel.
#' - `event.number.width`: Width of the event number.
#' - `expr.data.max`: Maximum expression data value.
#' - `expr.data.max.ceil`: Ceiling of the maximum expression data value.
#' - `expr.data.min`: Minimum expression data value.
#' - `channel`: Preliminary peak channels for the fluorophores.
#' - `channel.n`: Number of channels.
#' - `spectral.channel`: Spectral channel information.
#' - `spectral.channel.n`: Number of spectral channels.
#' - `sample`: Sample names (fluorophores).
#' - `scatter.and.channel`: FSC, SSC, and peak channel information.
#' - `scatter.and.channel.label`: Labels for scatter and channel.
#' - `scatter.and.channel.spectral`: FSC, SSC, and spectral channels.
#' - `scatter.parameter`: Scatter parameters used for gating.
#' - `event`: Event factor.
#' - `event.n`: Number of events.
#' - `event.sample`: Sample information for events. Links events to samples.
#' - `event.type`: Type of events.
#' - `expr.data`: Expression data used for extracting spectra.
#'
#' @seealso
#' * [tune.gate()]
#' * [define.gate.landmarks()]
#' * [define.gate.density()]
#' * [do.gate()]
#' * [gate.define.plot()]
#'
#' @export
#'
#' @references Roca, Carlos P et al. "AutoSpill is a principled framework that
#' simplifies the analysis of multichromatic flow cytometry data" \emph{Nature
#' Communications} 12 (2890) 2021.

define.flow.control.test <- function(
    control.dir,
    control.def.file,
    asp,
    gate = TRUE,
    gating.system = c("landmarks", "density"),
    gate.list = NULL,
    parallel = FALSE,
    verbose = TRUE,
    threads = NULL,
    color.palette = NULL
  ) {

  if ( verbose ) message( "\033[34mChecking control file for errors \033[0m" )
  check.control.file( control.dir, control.def.file, asp, strict = TRUE )

  # read control info
  if ( verbose ) message( "\033[34mReading control information \033[0m" )
  control.table <- utils::read.csv(
    control.def.file,
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

  # check the user-supplied gates for consistency and structure
  if ( !is.null( gate.list ) ) check.gates( gate.list, control.table, asp )

  # read channels from an FCS file
  all.channels <- colnames(
    readFCS( file.path( control.dir, control.table$filename[ 1 ] ) )
  )

  # remove unnecessary channels
  non.spectral.pattern <- paste0( asp$non.spectral.channel, collapse = "|" )

  if ( grepl( "Discover", asp$cytometer ) ) {
    spec.idx <- grep( asp$spectral.channel, all.channels )
  } else {
    spec.idx <- grep( non.spectral.pattern, all.channels, invert = TRUE )
  }

  spectral.channel <- all.channels[ spec.idx ]

  # reorganize channels if necessary
  spectral.channel <- check.channels( spectral.channel, asp )
  spectral.channel.n <- length( spectral.channel )

  # record and store voltages for checks during unmixing
  header <- readFCSheader( file.path( control.dir, control.table$filename[ 1 ] ) )[[ 1 ]]
  # extract the $PnV values
  spectral.voltages <- vapply( spec.idx, function( idx ) {
    v <- header[[ paste0( "$P", idx, "V" ) ]]
    if ( is.null( v ) ) return( NA_character_ ) else return( as.character( v ) )
  }, character( 1 ) )
  names( spectral.voltages ) <- all.channels[ spec.idx ]

  # set samples and gate combos
  control.table$sample <- control.table$fluorophore
  gating.system <- match.arg( gating.system )
  control.table <- assign.gates( control.table, gating.system, gate, verbose )

  # set factors and fill in missing data
  flow.fluorophore <- control.table$fluorophore
  flow.fluorophore[ is.na( flow.fluorophore ) ] <- "Negative"
  flow.antigen <- control.table$marker

  flow.control.type <- control.table$control.type

  flow.antigen[ flow.fluorophore == "AF" ] <- "AF"
  flow.antigen[ is.na( flow.antigen ) ] <- "other"

  flow.viability <- control.table$is.viability
  flow.viability[ is.na( flow.viability ) ] <- FALSE
  names( flow.viability ) <- control.table$sample

  flow.universal.negative <- control.table$universal.negative
  flow.universal.negative[ is.na( flow.universal.negative ) ] <- FALSE
  names( flow.universal.negative ) <- control.table$sample

  flow.large.gate <- control.table$large.gate
  flow.large.gate[ is.na( flow.large.gate ) ] <- FALSE
  names( flow.large.gate ) <- control.table$sample

  # set default AF channel if none has been provided
  flow.channel <- control.table$channel
  if ( any( flow.fluorophore == "AF" ) ) {
    idx <- which( flow.fluorophore == "AF" )
    if ( length( flow.channel[ idx ] ) == 0 || all( is.na( flow.channel[ idx ] ) ) ) {
      flow.channel[ idx ] <- asp$af.channel
    }
  }
  flow.channel[ is.na( flow.channel ) ] <- "other"

  flow.autof.marker.idx <- which( flow.antigen == "AF" )
  if ( length( flow.autof.marker.idx ) != 1 ) flow.autof.marker.idx <- NULL

  # read scatter parameters
  if ( verbose ) message( "\033[34mDetermining channels to be used \033[0m" )

  flow.scatter.parameter <- read.scatter.parameter( asp )

  # set scatter parameters and channels
  flow.scatter.and.channel <- c(
    asp$default.time.parameter,
    flow.scatter.parameter,
    flow.channel )
  flow.scatter.and.channel.spectral <- c(
    asp$default.time.parameter,
    flow.scatter.parameter,
    spectral.channel )

  flow.scatter.and.channel.matched.bool <-
    flow.scatter.and.channel.spectral %in% all.channels

  if ( !all( flow.scatter.and.channel.matched.bool ) ) {
    channel.matched <-
      flow.scatter.and.channel.spectral[ flow.scatter.and.channel.matched.bool ]
    flow.scatter.and.channel.unmatched <- paste0(
      sort( setdiff( flow.scatter.and.channel.spectral, channel.matched ) ),
      collapse = ", " )
    flow.set.unmatched <- paste0(
      sort( setdiff( all.channels, channel.matched ) ),
      collapse = ", " )
    error.msg <- sprintf(
      "wrong channel name, not found in fcs data\n\texpected: %s\n\tfound: %s",
      flow.scatter.and.channel.unmatched, flow.set.unmatched )
    stop( error.msg, call. = FALSE )
  }

  names( flow.channel ) <- control.table$sample

  if ( anyDuplicated( flow.scatter.and.channel.spectral ) != 0 )
    stop( "Names for channels overlap", call. = FALSE )

  # set labels for time, scatter parameters and channels
  flow.scatter.and.channel.label <- c(
    "Time", flow.scatter.parameter,
    ifelse( ! is.na( flow.antigen ),
            paste0( flow.antigen, " - ", flow.fluorophore ), flow.channel )
    )
  names( flow.scatter.and.channel.label ) <- flow.scatter.and.channel

  # get range of fcs data
  flow.set.resolution <- asp$expr.data.max
  flow.expr.data.min <- asp$expr.data.min
  flow.expr.data.max <- asp$expr.data.max
  flow.expr.data.max.ceil <- ceiling( flow.expr.data.max / asp$data.step ) *
    asp$data.step

  # create figure and table directories
  if ( verbose ) message( "\033[34mCreating output folders \033[0m" )
  create.directory( asp )

  # Ensure flow.file.name and other vectors match the expanded control.table
  flow.file.name <- control.table$filename
  names( flow.file.name ) <- control.table$sample

  final.gate.list <- list()

  if ( gate ) {
    if ( verbose ) message( "\033[34mDefining gates \033[0m" )
    unique.names <- unique( control.table$gate.name )

    for ( g.name in unique.names ) {
      is.orphan <- grepl( "density_orphan", g.name )
      # prioritize pre-defined gate from user
      if ( !is.null( gate.list ) && g.name %in% names( gate.list ) ) {
        final.gate.list[[ g.name ]] <- gate.list[[ g.name ]]
      } else {
        # define gate using the chosen system (landmarks/density)
        if ( is.orphan || gating.system == "density" ) {
          final.gate.list[[ g.name ]] <- define.gate.density(
            control.file = NULL,
            control.dir = control.dir,
            asp = asp,
            output.dir = asp$figure.gate.dir,
            gate.name = g.name,
            verbose = FALSE,
            control.table = control.table,
            color.palette = if ( is.null(color.palette) ) "plasma" else color.palette
          )
        } else {
          final.gate.list[[ g.name ]] <- define.gate.landmarks(
            control.file = NULL,
            control.dir = control.dir,
            asp = asp,
            output.dir = asp$figure.gate.dir,
            gate.name = g.name,
            verbose = FALSE,
            control.table = control.table,
            check = FALSE,
            color.palette = if ( is.null(color.palette) ) "plasma" else color.palette
          )
        }
      }
    }
  }

  # map sample names to gate names
  flow.gate <- control.table$gate.name
  names( flow.gate ) <- control.table$sample

  # read in FCS files
  if ( verbose ) message( "\033[34mReading FCS files \033[0m" )

  args.list <- list(
    file.name = flow.file.name,
    control.dir = control.dir,
    scatter.and.spectral.channel = flow.scatter.and.channel.spectral,
    spectral.channel = spectral.channel,
    set.resolution = flow.set.resolution,
    flow.gate = flow.gate,
    gate.list = final.gate.list,
    scatter.param = flow.scatter.parameter,
    scatter.and.channel.label = flow.scatter.and.channel.label,
    asp = asp,
    apply.gate = gate,
    color.palette = if ( is.null(color.palette) ) "mako" else color.palette
  )

  # set up parallel processing
  if ( parallel ) {
    if ( is.null( threads ) ) threads <- asp$worker.process.n

    if ( verbose & gate ) message( "\033[34mPlotting gates... \033[0m" )

    exports <- c( "control.table$sample", "args.list", "gate.sample.plot",
                  "get.gated.flow.expression.data", "readFCS" )
    result <- create.parallel.lapply(
      asp,
      exports,
      parallel = parallel,
      threads = threads,
      export.env = environment(),
      allow.mclapply.mac = TRUE
    )
    lapply.function <- result$lapply
  } else {
    lapply.function <- lapply
    result <- list( cleanup = NULL )
  }

  # main call to read in flow data
  flow.expr.data <- tryCatch( {
    lapply.function( control.table$sample, function( f ) {
      do.call( get.gated.flow.expression.data, c( list( f ), args.list ) )
    } )
  }, finally = {
    # clean up cluster when done if applicable
    if ( !is.null( result$cleanup ) ) result$cleanup()
  } )

  names( flow.expr.data ) <- control.table$sample

  # organize data
  if ( verbose ) message( "\033[34mOrganizing control info \033[0m" )

  flow.sample.n <- length( control.table$sample )
  flow.sample.event.number.max <- 0

  for ( fs.idx in 1 : flow.sample.n ) {
    flow.sample.event.number <- nrow( flow.expr.data[[ fs.idx ]]  )

    rownames( flow.expr.data[[ fs.idx ]] ) <- paste(
      control.table$sample[ fs.idx ], seq_len( flow.sample.event.number ), sep = "_" )

    if ( flow.sample.event.number < 500 ) {
      warning(
        paste0(
          "\033[31m",
          "Warning! Fewer than 500 gated events in ",
          flow.file.name[ fs.idx ],
          "\033[0m",
          "\n"
        )
      )
    }

    if ( flow.sample.event.number > flow.sample.event.number.max )
      flow.sample.event.number.max <- flow.sample.event.number
  }

  flow.event.number.width <-
    floor( log10( flow.sample.event.number.max ) ) + 1
  flow.event.regexp <- sprintf( "\\.[0-9]{%d}$", flow.event.number.width )

  # set rownames
  for ( fs.idx in 1 : flow.sample.n ) {
    flow.sample.event.number <- nrow( flow.expr.data[[ fs.idx ]]  )
    flow.the.sample <- control.table$sample[ fs.idx ]
    flow.the.event <- sprintf(
      "%s.%0*d", flow.the.sample,
      flow.event.number.width, 1 : flow.sample.event.number
      )
    rownames( flow.expr.data[[ fs.idx ]] ) <- flow.the.event
  }

  flow.expr.data <- do.call( rbind, flow.expr.data )

  # set events
  flow.event <- rownames( flow.expr.data )
  flow.event.n <- length( flow.event )
  flow.event.sample <- sub( flow.event.regexp, "", flow.event )
  flow.event.sample <- factor( flow.event.sample, levels = control.table$sample )
  event.type.factor <- control.table$sample
  names( event.type.factor ) <- flow.control.type
  flow.event.type <- factor(
    flow.event.sample,
    levels = event.type.factor,
    labels = names( event.type.factor ) )
  names( flow.control.type ) <- flow.fluorophore

  # quickly re-determine peak AF channel empirically
  if ( any( flow.fluorophore == "AF" ) ) {
    idx <- which( flow.fluorophore == "AF" )
    af.data <- flow.expr.data[ which( flow.event.sample == "AF" ), ]
    af.max <- which.max( colMeans( af.data[ , spectral.channel ] ) )
    flow.channel[ idx ] <- spectral.channel[ af.max ]
  }

  # make control info
  flow.control <- list(
    filename = flow.file.name,
    fluorophore = flow.fluorophore,
    antigen = flow.antigen,
    control.type = flow.control.type,
    universal.negative = flow.universal.negative,
    viability = flow.viability,
    large.gate = flow.large.gate,
    autof.channel.idx = flow.autof.marker.idx,
    event.number.width = flow.event.number.width,
    expr.data.max = flow.expr.data.max,
    expr.data.max.ceil = flow.expr.data.max.ceil,
    expr.data.min = flow.expr.data.min,
    channel = flow.channel,
    channel.n = length( flow.channel ),
    spectral.channel = spectral.channel,
    spectral.channel.n = spectral.channel.n,
    voltages = spectral.voltages,
    sample = control.table$sample,
    scatter.and.channel = flow.scatter.and.channel,
    scatter.and.channel.label = flow.scatter.and.channel.label,
    scatter.and.channel.spectral = flow.scatter.and.channel.spectral,
    scatter.parameter = flow.scatter.parameter,
    event = flow.event,
    event.n = flow.event.n,
    event.sample = flow.event.sample,
    event.type = flow.event.type,
    expr.data = flow.expr.data
  )

  if ( verbose && gate ) {
    message(
      paste0(
        "\033[32m",
        "Control setup complete!",
        "\n",
        "Review gates in figure_gate.",
        "\033[0m"
      )
    )
  } else if ( verbose ) {
    message(
      paste0(
        "\033[32m",
        "Control setup complete!",
        "\033[0m"
      )
    )
  }

  return( flow.control )
}
