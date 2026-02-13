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
  database.path <- system.file( "extdata", "cytometer_database.csv",
                                package = "AutoSpectral" )

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
