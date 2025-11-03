# fluorophore_database.r

#' @title Fluorophore Database
#'
#' @description
#' Information about fluorophores and their detection on various cytometers.
#'
#' @format A data frame with the following columns:
#' - `fluorophore`: Fluorophore name
#' - `synonym1`: Fluorophore synonym 1
#' - `synonym2`: Fluorophore synonym 2
#' - `synonym3`: Fluorophore synonym 3
#' - `channel.Aurora`: Peak channels for fluorophores on the 5-laser Aurora
#' - `channel.ID7000`: Peak channels for fluorophores on the 5-laser ID7000
#' - `channel.s8`: Peak channels for fluorophores on the FACSDiscover S8
#' - `channel.opteon`: Peak channels for fluorophores on the Opteon
#' - `channel.mosaic`: Peak channels for fluorophores on the Mosaic#'
#' - `channel.xenith`: Peak channels for fluorophores on the Xenith
#' - `excitation.laser`: Excitation laser
#' - `nominal.wavelength`: Numeric; nominal peak emission wavelength
#' - `is.viability`: Logical; whether the fluorophore is a viability dye
