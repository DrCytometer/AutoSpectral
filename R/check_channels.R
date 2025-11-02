# check_channels.r

#' @title Check Channels
#'
#' @description
#' A helper function to reorganize the spectral channels in a nice order for
#' plotting. Puts them in excitation/emission order.
#'
#' @importFrom utils read.csv
#'
#' @param spectral.channels Vector of initial spectral channel names.
#' @param cytometer The cytometer. Use `asp$cytometer`.
#'
#' @return Returns the vector of spectral channels re-organized in excitation-
#' emission order. That is, narrowest to longest excitation laser, wtih narrowest
#' to longest emission wavelength inside each laser group.
#'
#' @export

check.channels <- function( spectral.channels, cytometer ) {
  # check for `FJ-Comp`
  fj.comp <- grepl( "FJ-Comp", spectral.channels )

  if ( any( fj.comp ) )
    warning( "FJ-Comp channels have been detected as the `spectral.channels`.
             This may cause errors in reading FCS files or unmixing due to channel mismatches." )

  # match against reference channels for this cytometer
  database.path <- system.file( "extdata", "cytometer_database.csv",
                                     package = "AutoSpectral" )
  cytometers <- read.csv( database.path )

  if ( cytometer == "Aurora" ) {
    detectors <- cytometers$Aurora
  } else if ( cytometer == "ID7000" ) {
    detectors <- cytometers$ID7000
  } else if ( grepl( "Discover", cytometer ) ) {
    detectors <- cytometers$Discover
  } else if ( cytometer == "Opteon" ) {
    detectors <- cytometers$Opteon
  } else {
    stop( "Unsupported cytometer" )
  }

  # re-arrange in excitation/emission order
  common.channels <- intersect( detectors, spectral.channels )
  extra.channels  <- setdiff( spectral.channels, detectors )
  spectral.channels <- c( common.channels, extra.channels )

  return( spectral.channels )
}
