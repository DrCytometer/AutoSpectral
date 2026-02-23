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
#' @param parallel Logical, default is `FALSE`, in which case parallel processing
#' will not be used. Parallel processing will likely be faster when many small
#' files are read in. If the data is larger, parallel processing may not
#' accelerate the process much.
#' @param verbose Logical, default is `TRUE`. Set to `FALSE` to suppress messages.
#' @param threads Numeric, number of threads to use for parallel processing.
#' Default is `NULL` which will revert to `asp$worker.process.n` if
#' `parallel=TRUE`.
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
#' @export

define.flow.control <- function(
    control.dir,
    control.def.file,
    asp,
    gate = TRUE,
    parallel = FALSE,
    verbose = TRUE,
    threads = NULL
  ) {

  if ( verbose ) message( "\033[34mChecking control file for errors \033[0m" )

  check.control.file( control.dir, control.def.file, asp, strict = TRUE )

  if ( parallel && is.null( threads ) ) threads <- asp$worker.process.n

  # read channels from controls
  if ( verbose ) message( "\033[34mReading control information \033[0m" )

  # read control info
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

  # read channels from an FCS file
  flow.set.channel <- colnames(
    readFCS( file.path( control.dir, control.table$filename[ 1 ] ) )
  )

  # remove unnecessary channels
  non.spectral.channel <- asp$non.spectral.channel
  non.spectral.channel <- paste0( non.spectral.channel, collapse = "|" )
  flow.spectral.channel <- flow.set.channel[
    !grepl( non.spectral.channel, flow.set.channel ) ]

  if ( grepl( "Discover", asp$cytometer ) )
    flow.spectral.channel <- flow.spectral.channel[
      grep( asp$spectral.channel, flow.spectral.channel ) ]

  # reorganize channels if necessary
  flow.spectral.channel <- check.channels( flow.spectral.channel, asp )
  flow.spectral.channel.n <- length( flow.spectral.channel )

  # get fluorophores and markers
  control.table$sample <- control.table$fluorophore
  control.table$is.viability[ is.na( control.table$is.viability ) ] <- FALSE
  control.table$large.gate[ is.na( control.table$large.gate ) ] <- FALSE

  # universal.negative must be strictly filename or NA
  control.table$universal.negative[
    control.table$universal.negative == "" |
      is.na(control.table$universal.negative)
  ] <- NA

  # identify universal negative types
  negative.types <- data.frame(
    negative = control.table$universal.negative,
    large.gate = control.table$large.gate,
    viability = control.table$is.viability,
    stringsAsFactors = FALSE
  )

  # identify unique ones
  unique.neg <- unique( negative.types[ !is.na( negative.types$negative ), ] )

  # define samples (replicate variously gated negatives)
  if ( verbose ) message( "\033[34mMatching negatives for the controls \033[0m" )

  if ( nrow( unique.neg ) > 0 ) {
    for ( i in seq_len( nrow( unique.neg ) ) ) {

      neg.file <- unique.neg$negative[ i ]
      lg       <- unique.neg$large.gate[ i ]
      vi       <- unique.neg$viability[ i ]

      # does a negative with this gate combination already exist?
      match.found <- any(
        control.table$filename == neg.file &
          control.table$large.gate == lg &
          control.table$is.viability == vi
      )

      if ( !match.found ) {
        # create a new row copying from the source negative
        matching.row <- control.table[control.table$filename == neg.file, ]
        new.row <- matching.row[ 1, ]
        new.row$large.gate        <- lg
        new.row$is.viability      <- vi
        new.row$fluorophore       <- paste(new.row$control.type, "Negative", i)
        new.row$sample            <- new.row$fluorophore
        new.row$universal.negative <- NA   # new row is a replicate, not a source
        control.table <- rbind(control.table, new.row)
      }
    }
  }

  # define gates needed
  if ( verbose ) message( "\033[34mDetermining the gates that will be needed \033[0m" )

  gate.list <- list()
  gate.types <- data.frame(
    negative        = control.table$universal.negative,
    type        = control.table$control.type,
    viability   = control.table$is.viability,
    large.gate  = control.table$large.gate,
    stringsAsFactors = FALSE
  )

  # replace FALSE with NA so missing universal negatives will be treated identically
  unique.gates <- unique( gate.types )
  match.gate <- function( row, unique.gates ) {
    matches <- which(
      ( ( row$negative %in% unique.gates$negative ) |
          ( is.na( row$negative ) & is.na( unique.gates$negative ) ) ) &
        ( row$type       == unique.gates$type ) &
        ( row$viability  == unique.gates$viability ) &
        ( row$large.gate == unique.gates$large.gate )
    )

    if ( length( matches ) == 0 )
      return( NA_integer_ )  # no match found (should not happen)

    matches[ 1 ]
  }

  control.table$gate <- vapply(
    seq_len( nrow( gate.types ) ),
    function( i ) match.gate( gate.types[ i, ], unique.gates ),
    FUN.VALUE = integer( 1 )
  )
  flow.gate <- control.table$gate

  # rewrite universal.negative as sample names rather than filenames
  control.table$universal.negative <- sapply(
    seq_len( nrow( control.table ) ), function( i ) {
    neg.file <- control.table$universal.negative[ i ]

    # If no assigned universal negative, keep NA
    if ( is.na( neg.file ) ) return( NA )

    # Try matching within the same gate
    match.row <- subset(
      control.table,
      filename == neg.file &
        gate == control.table$gate[ i ]
    )

    if ( nrow( match.row ) > 0 )
      return( match.row$sample[ 1 ] )

    # Else fallback to any matching filename
    match.row <- subset( control.table, filename == neg.file )
    if ( nrow( match.row ) > 0 )
      return( match.row$sample[ 1 ] )

    return( NA )
  } )

  # set samples and gate combos
  flow.sample <- control.table$sample
  flow.sample.n <- length( flow.sample )

  if ( gate )
    names( flow.gate ) <- flow.sample

  flow.file.name <- control.table$filename
  names( flow.file.name ) <- flow.sample
  flow.fluorophore <- control.table$fluorophore
  flow.fluorophore[ is.na( flow.fluorophore ) ] <- "Negative"
  flow.antigen <- control.table$marker
  flow.channel <- control.table$channel
  flow.control.type <- control.table$control.type
  flow.universal.negative <- control.table$universal.negative
  flow.viability <- control.table$is.viability
  flow.large.gate <- control.table$large.gate
  flow.antigen[ flow.fluorophore == "AF" ] <- "AF"
  flow.antigen[ is.na( flow.antigen ) ] <- "other"

  # set default AF channel if none has been provided
  if ( any( flow.fluorophore == "AF" ) ) {
    idx <- which( flow.fluorophore == "AF" )
    if ( length( flow.channel[ idx ] ) == 0 || all( is.na( flow.channel[ idx ] ) ) ) {
      flow.channel[ idx ] <- asp$af.channel
    }
  }

  flow.channel[ is.na( flow.channel ) ] <- "other"
  flow.universal.negative[ is.na( flow.universal.negative ) ] <- FALSE
  names( flow.universal.negative ) <- flow.sample
  flow.viability[ is.na( flow.viability ) ] <- FALSE
  flow.large.gate[ is.na( flow.large.gate ) ] <- FALSE
  flow.channel.n <- length( flow.channel )
  flow.autof.marker.idx <- which( flow.antigen == "AF" )
  if ( length( flow.autof.marker.idx ) != 1 )
    flow.autof.marker.idx <- NULL

  names( flow.viability ) <- flow.sample
  names( flow.large.gate ) <- flow.sample

  # read scatter parameters
  if ( verbose ) message( "\033[34mDetermining channels to be used \033[0m" )

  flow.scatter.parameter <- read.scatter.parameter( asp )

  # set scatter parameters and channels
  flow.scatter.and.channel <- c(
    asp$default.time.parameter,
    flow.scatter.parameter, flow.channel )
  flow.scatter.and.channel.spectral <- c(
    asp$default.time.parameter,
    flow.scatter.parameter,
    flow.spectral.channel )

  flow.scatter.and.channel.matched.bool <-
    flow.scatter.and.channel.spectral %in% flow.set.channel

  if ( ! all( flow.scatter.and.channel.matched.bool ) ) {
    channel.matched <-
      flow.scatter.and.channel.spectral[ flow.scatter.and.channel.matched.bool ]
    flow.scatter.and.channel.unmatched <- paste0(
      sort( setdiff( flow.scatter.and.channel.spectral, channel.matched ) ),
      collapse = ", " )
    flow.set.unmatched <- paste0(
      sort( setdiff( flow.set.channel, channel.matched ) ),
      collapse = ", " )
    error.msg <- sprintf(
      "wrong channel name, not found in fcs data\n\texpected: %s\n\tfound: %s",
      flow.scatter.and.channel.unmatched, flow.set.unmatched )
    stop( error.msg, call. = FALSE )
  }

  names( flow.channel ) <- flow.sample

  if ( anyDuplicated( flow.scatter.and.channel.spectral ) != 0 )
    stop( "internal error: names for channels overlap", call. = FALSE )

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

  if ( gate ) {
    # define gates on downsampled pooled fcs by type
    if ( verbose ) message( "\033[34mDefining the gates \033[0m" )

    # Remove rows with NA filenames before making gate list
    valid.rows <- !is.na( control.table$filename )
    control.table.valid <- control.table[valid.rows, ]

    for( gate.id in unique( control.table.valid$gate ) ) {
      files.to.gate <- unique(
        control.table.valid[ control.table.valid$gate == gate.id, ]$filename
      )

      files.n <- length( files.to.gate )
      downsample.n <- ceiling( asp$gate.downsample.n.cells / files.n )
      scatter.coords <- lapply( files.to.gate, function( f ) {
        sample.fcs.file(
          f,
          control.dir = control.dir,
          downsample.n = downsample.n,
          asp = asp
        )
      } )

      scatter.coords <- do.call( rbind, scatter.coords )

      viability.gate <- unique(
        control.table.valid[ control.table.valid$gate == gate.id, ]$is.viability
      )
      is.viability <- ifelse(
        viability.gate,
        "viabilityMarker",
        "nonViability"
      )
      large.gate <- unique(
        control.table.valid[ control.table.valid$gate == gate.id, ]$large.gate
      )
      is.large.gate <- ifelse( large.gate, "largeGate", "smallGate" )
      control.type <- unique(
        control.table.valid[ control.table.valid$gate == gate.id, ]$control.type
      )
      samp <- paste( control.type, is.viability, is.large.gate, sep = "_" )

      stopifnot(
        length( viability.gate ) == 1,
        length( large.gate ) == 1,
        length( control.type ) == 1
      )

      gate.boundary <- do.gate(
        scatter.coords,
        viability.gate,
        large.gate,
        samp,
        flow.scatter.and.channel.label,
        control.type,
        asp
      )

      gate.list[[ gate.id ]] <- gate.boundary
    }

  }

  # read in FCS files
  if ( verbose ) message( "\033[34mReading FCS files \033[0m" )

  args.list <- list(
    file.name = flow.file.name,
    control.dir = control.dir,
    scatter.and.spectral.channel = flow.scatter.and.channel.spectral,
    spectral.channel = flow.spectral.channel,
    set.resolution = flow.set.resolution,
    flow.gate = flow.gate,
    gate.list = gate.list,
    scatter.param = flow.scatter.parameter,
    scatter.and.channel.label = flow.scatter.and.channel.label,
    asp = asp,
    apply.gate = gate
  )

  # set up parallel processing
  if ( parallel ) {
    if ( verbose & gate ) message( "\033[34mPlotting gates... \033[0m" )

    exports <- c( "flow.sample", "args.list", "gate.sample.plot",
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

  FUN_local <- get.gated.flow.expression.data

  # main call to read in flow data
  flow.expr.data <- tryCatch( {
    lapply.function( flow.sample, function( f ) {
      do.call( FUN_local, c( list( f ), args.list ) )
    } )
  }, finally = {
    # clean up cluster when done
    if ( !is.null( result$cleanup ) ) result$cleanup()
  } )

  names( flow.expr.data ) <- flow.sample

  # organize data
  if ( verbose ) message( "\033[34mOrganizing control info \033[0m" )

  flow.sample.event.number.max <- 0

  for ( fs.idx in 1 : flow.sample.n )
  {
    flow.sample.event.number <- nrow( flow.expr.data[[ fs.idx ]]  )

    rownames( flow.expr.data[[ fs.idx ]] ) <- paste(
      flow.sample[ fs.idx ], seq_len( flow.sample.event.number ), sep = "_" )

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
    flow.the.sample <- flow.sample[ fs.idx ]
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
  flow.event.sample <- factor( flow.event.sample, levels = flow.sample )
  event.type.factor <- flow.sample
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
    af.max <- which.max( colMeans( af.data[ , flow.spectral.channel ] ) )
    flow.channel[ idx ] <- flow.spectral.channel[ af.max ]
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
    channel.n = flow.channel.n,
    spectral.channel = flow.spectral.channel,
    spectral.channel.n = flow.spectral.channel.n,
    sample = flow.sample,
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
