# read_bd_spectra.r

#' @title Extract Spectra From BD FCS File
#'
#' @description
#' Extracts the spectra (spillover) values from BD FCS files (e.g., A8, S8)
#'
#' @importFrom flowCore read.FCSheader
#'
#' @param fcs.file Path and filename to the FCS file to be read.
#'
#' @return The spillover matrix as fluorophores x detectors, with fluorophore
#' and detector names applied.
#' @export


read.bd.spectra <- function( fcs.file ) {

  spill.vector <- suppressWarnings(
    flowCore::read.FCSheader(
      fcs.file,
      keyword = "SPILL"
    )[[ 1 ]]
  )

  spill.vector <- strsplit( spill.vector, "," )[[ 1 ]]
  n.detectors <- as.numeric( spill.vector[ 1 ] )
  detector.names <- spill.vector[ 2:( n.detectors + 1 ) ]

  spill.matrix <- matrix(
    as.numeric(
      spill.vector[ ( n.detectors + 2 ) : length( spill.vector ) ]
      ),
    nrow = n.detectors, ncol = n.detectors
    )
  spill.matrix <- t( spill.matrix )

  marker.names <- suppressWarnings(
    flowCore::read.FCSheader(
      fcs.file,
      keyword = "BDSPECTRAL UNMIXED"
    )[[ 1 ]]
  )

  marker.names <- strsplit( marker.names, "," )[[ 1 ]]
  marker.n <- length( marker.names )

  spill.matrix <- spill.matrix[ 1:marker.n, ]
  colnames( spill.matrix ) <- detector.names
  rownames( spill.matrix ) <- marker.names

  return( spill.matrix )
}
