# reload_flow_control.r

#' @title Reload Flow Control Information
#'
#' @description
#' This function reloads essential information from control files to permit
#' rapid unmixing at a later date, without recalculating spectra or gates.
#'
#' @param control.dir file path to the single stained control FCS files
#' @param control.def.file csv file defining the single color control file
#' names, fluorophores they represent, marker names, peak channels and
#' gating requirements.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#'
#' @return A list containing the reloaded flow control information.
#'
#' @export

reload.flow.control <- function(
    control.dir,
    control.def.file,
    asp
  ) {

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

  if ( anyDuplicated( control.table$filename ) != 0 )
    stop( "duplicated filenames in fcs data", call. = FALSE )

  # read channels from an FCS file
  all.channels <- colnames(
    readFCS( file.path( control.dir, control.table$filename[ 1 ] ) )
  )

  # remove unnecessary channels
  non.spectral.pattern <- paste0( asp$non.spectral.channel, collapse = "|" )
  spec.idx <- grep( non.spectral.pattern, all.channels, invert = TRUE )
  spectral.channel <- all.channels[ spec.idx ]

  if ( grepl( "Discover", asp$cytometer ) ) {
    spec.idx <- grep( asp$spectral.channel, spectral.channel )
    spectral.channel <- spectral.channel[ spec.idx ]
  }

  # reorganize channels if necessary
  spectral.channel <- check.channels( spectral.channel, asp )
  spectral.channel.n <- length( spectral.channel )

  ## record and store voltages for checks during unmixing
  # read header from the first file
  header <- readFCSheader( file.path( control.dir, control.table$filename[ 1 ] ) )[[ 1 ]]

  # get parameter names
  p.names <- unlist( header[ grep( "^\\$P\\d+N$", names( header ) ) ] )

  # initialize with fallback
  spectral.voltages <- stats::setNames(
    rep( NA_character_, length( spectral.channel ) ), spectral.channel
  )

  # ID7000 doesn't store voltage/gain info
  if ( !grepl( "ID7000", asp$cytometer, ignore.case = TRUE ) ) {
    spectral.voltages <- tryCatch({
      vapply( spectral.channel, function( ch ) {
        # match channel name to header index key (e.g., "$P10N")
        p.idx.key <- names( p.names )[ which( p.names == ch ) ]
        if ( length( p.idx.key ) == 0 ) return( NA_character_ )

        # extract the numeric part (the 'n')
        n <- gsub( "[^0-9]", "", p.idx.key )

        # Mosaic uses $PnG; others $PnV
        pnv.id <- ifelse( grepl( "Mosaic", asp$cytometer, ignore.case = TRUE ), "G", "V" )
        val <- header[[ paste0( "$P", n, pnv.id ) ]]

        if ( is.null( val ) ) return( NA_character_ ) else return( as.character( val ) )

      }, character( 1 ) )
    }, error = function( e ) {
      warning( "Failed to extract spectral voltages/gains: ", e$message,
               call. = FALSE )
      return( spectral.voltages )
    })
  }

  names( spectral.voltages ) <- spectral.channel

  # get fluorophores and markers
  flow.fluorophore <- control.table$fluorophore
  flow.fluorophore[ is.na( flow.fluorophore ) ] <- "Negative"

  flow.control.type <- control.table$control.type
  names( flow.control.type ) <- flow.fluorophore

  flow.antigen <- control.table$marker
  flow.channel <- control.table$channel

  # read scatter parameters
  flow.scatter.parameter <- read.scatter.parameter( asp )

  # set labels for time, scatter parameters and channels
  flow.scatter.and.channel <- c(
    asp$default.time.parameter, flow.scatter.parameter, flow.channel
    )
  flow.scatter.and.channel.spectral <- c(
    asp$default.time.parameter,
    flow.scatter.parameter,
    spectral.channel
  )

  flow.scatter.and.channel.label <- c(
    "Time", flow.scatter.parameter,
    ifelse(
      ! is.na( flow.antigen ),
      paste0( flow.antigen, " - ", flow.fluorophore ),
      flow.channel
    )
  )
  names( flow.scatter.and.channel.label ) <- flow.scatter.and.channel

  # make control info
  flow.control <- list(
    filename = control.table$filename,
    fluorophore = flow.fluorophore,
    control.type = flow.control.type,
    antigen = flow.antigen,
    channel = flow.channel,
    expr.data.max = asp$expr.data.max,
    expr.data.min = asp$expr.data.min,
    spectral.channel = spectral.channel,
    voltages = spectral.voltages,
    sample = control.table$sample,
    scatter.and.channel.label = flow.scatter.and.channel.label,
    scatter.and.channel.spectral = flow.scatter.and.channel.spectral,
    scatter.parameter = flow.scatter.parameter
  )

  return( flow.control )

}
