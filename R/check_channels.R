# check_channels.r

#' @title Check Channels
#'
#' @description
#' A helper function to reorganize the spectral channels in a nice order for
#' plotting. Puts them in excitation/emission order.
#'
#' @param spectral.channels Vector of initial spectral channel names.
#' @param asp The AutoSpectral parameter list. Generate using
#' `get.autospectral.param`
#'
#' @return Returns the vector of spectral channels re-organized in excitation-
#' emission order. That is, narrowest to longest excitation laser, wtih narrowest
#' to longest emission wavelength inside each laser group.
#'
#' @export

check.channels <- function( spectral.channels, asp ) {

  # check for `FJ-Comp`
  fj.comp <- grepl( "FJ-Comp", spectral.channels )

  if ( any( fj.comp ) )
    warning( "FJ-Comp channels have been detected as the `spectral.channels`.
             This may cause errors in reading FCS files or unmixing due to channel mismatches." )

  # match against reference channels for this cytometer
  database.path <- system.file(
    "extdata", "cytometer_database.csv",
    package = "AutoSpectral"
  )

  cytometers <- utils::read.csv( database.path )

  if ( asp$cytometer == "Aurora" ) {
    detectors <- cytometers$Aurora
    if ( asp$cytometer.version == "NL" ) {
      detectors <- cytometers$NorthernLights
    } else {
      detectors <- cytometers$Aurora
    }
  } else if ( asp$cytometer == "ID7000" ) {
    detectors <- cytometers$ID7000
  } else if ( grepl( "Discover", asp$cytometer ) ) {
    detectors <- cytometers$Discover
  } else if ( asp$cytometer == "Opteon" ) {
    detectors <- cytometers$Opteon
  } else if ( asp$cytometer == "Mosaic" ) {
    detectors <- cytometers$Mosaic
  } else if ( asp$cytometer == "Xenith" ) {
    detectors <- cytometers$Xenith
  } else if ( asp$cytometer == "Symphony" ) {
    detectors <- cytometers$A5SE
  } else {
    stop( "Unsupported cytometer" )
  }

  # re-arrange in excitation/emission order
  common.channels <- intersect( detectors, spectral.channels )
  extra.channels  <- setdiff( spectral.channels, detectors )
  spectral.channels <- c( common.channels, extra.channels )

  return( spectral.channels )
}

#' @title Spectral Voltage/Gain Keyword Suffix
#' @description
#' Internal helper. Determines which FCS keyword suffix ($PnV, $PnG, or
#' $PnR) should be compared as the "voltage" for a given cytometer.
#' - Mosaic uses $PnG (gain) instead of $PnV.
#' - ID7000 has no per-channel PMT voltage.
#' - All other supported cytometers use $PnV.
.spectral.voltage.suffix <- function( asp, header ) {
  if ( grepl( "ID7000", asp$cytometer, ignore.case = TRUE ) ) return( NA_character_ )
  if ( grepl( "Mosaic", asp$cytometer, ignore.case = TRUE ) ) return( "G" )
  "V"
}

#' @title Extract Spectral Voltages/Gains
#' @description
#' Internal helper. Given an FCS header and a keyword suffix from
#' `.spectral.voltage.suffix()`, extracts the voltage/gain value for each
#' named spectral channel (matched by channel name, not position).
.extract.spectral.voltages <- function( header, spectral.channel, suffix ) {
  if ( is.na( suffix ) ) {
    return( stats::setNames(
      rep( NA_character_, length( spectral.channel ) ), spectral.channel
    ) )
  }
  p.names <- unlist( header[ grep( "^\\$P\\d+N$", names( header ) ) ] )
  stats::setNames(
    vapply( spectral.channel, function( ch ) {
      p.idx.key <- names( p.names )[ which( p.names == ch ) ]
      if ( length( p.idx.key ) == 0 ) return( NA_character_ )
      n <- gsub( "[^0-9]", "", p.idx.key )
      val <- header[[ paste0( "$P", n, suffix ) ]]
      if ( is.null( val ) ) NA_character_ else as.character( val )
    }, character( 1 ) ),
    spectral.channel
  )
}

#' @title Check Control Voltage Consistency
#' @description
#' Internal helper. Compares detector voltage/gain settings for each spectral
#' channel across a set of single-stained control FCS files. The first file
#' in `filenames` is the reference; every other file is compared against it,
#' channel by channel. Returns a data frame of mismatches (zero rows if
#' everything is consistent, or if the cytometer/file exposes no usable
#' voltage/gain keyword).
.check.control.voltages <- function( control.dir, filenames, spectral.channel, asp ) {

  empty <- data.frame(
    file = character(), channel = character(),
    ref.voltage = character(), current.voltage = character(),
    stringsAsFactors = FALSE
  )

  if ( length( filenames ) < 2 ) return( empty )

  ref.header <- readFCSheader( file.path( control.dir, filenames[ 1 ] ) )[[ 1 ]]
  suffix     <- .spectral.voltage.suffix( asp, ref.header )

  # no usable voltage/gain keyword for this cytometer/file -> nothing to check
  if ( is.na( suffix ) ) return( empty )

  ref.voltages <- .extract.spectral.voltages( ref.header, spectral.channel, suffix )

  mismatches <- list()

  for ( f in filenames[ -1 ] ) {
    curr.header   <- readFCSheader( file.path( control.dir, f ) )[[ 1 ]]
    curr.voltages <- .extract.spectral.voltages( curr.header, spectral.channel, suffix )

    diff.idx <- which( !mapply( identical, ref.voltages, curr.voltages ) )

    if ( length( diff.idx ) > 0 ) {
      mismatches[[ length( mismatches ) + 1 ]] <- data.frame(
        file            = f,
        channel         = spectral.channel[ diff.idx ],
        ref.voltage     = ref.voltages[ diff.idx ],
        current.voltage = curr.voltages[ diff.idx ],
        stringsAsFactors = FALSE
      )
    }
  }

  if ( length( mismatches ) == 0 ) return( empty )
  do.call( rbind, mismatches )
}
