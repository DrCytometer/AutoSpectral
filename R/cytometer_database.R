# cytometer_database.r

#' @title Cytometer Database
#'
#' @description
#' Reference names for fluorescence spectral detectors for the supported cytometers.
#'
#' @format A data frame with the following columns:
#' - `Aurora`: Cytek Aurora channels
#' - `Aurora_laser`: Laser associated with each channel on the Aurora
#' - `NorthernLights`: Cytek Northern Lights channels
#' - `NothernLights_laser`: Laser associated with each channel on the NorthernLights
#' - `ID7000`: Sony ID7000 channels (7-laser is used as archetype)
#' - `ID7000_laser`: Laser associated with each channel on the ID7000
#' - `Discover`: BD FACSDiscover (A8 and S8) channels
#' - `Discover_laser`: Laser associated with each channel on the Discover
#' - `Opteon`: Agilent Novocyte Opteon channels
#' - `Opteon_laser`: Laser associated with each channel on the Opteon
#' - `Mosaic`: Beckman Coulter CytoFLEX Mosaic channels
#' - `Mosaic_laser`: Laser associated with each channel on the Mosaic
#' - `Xenith`: ThermoFisher Attune Xenith channels
#' - `Xenith_description`: Filters for the Xenith
#' - `Xenith_laser`: Laser associated with each channel on the Xenith
#' - `A5SE`: BD FACSymphony A5 SE channels
#' - `A5SE_laser`: Laser associated with each channel on the A5 SE
#' - `CytoStellar_laser`: Laser associated with each channel on the CytoStellar
