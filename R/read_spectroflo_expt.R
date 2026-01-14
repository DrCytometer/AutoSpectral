# read_spectroflo_expt.r

#' @title Extract Spectra From SpectroFlo Expt File
#'
#' @description
#' Reads an Experiment (.Expt) file from SpectroFlo and extracts the spillover
#' matrix that was used for the unmixing.
#'
#' @importFrom xml2 read_xml xml_find_all xml_text
#' @importFrom utils write.csv
#' @importFrom flowCore read.FCSheader
#'
#' @param expt.file File name and path to the .Expt file to be read.
#' @param output.dir Directory where the spillover .csv file will be written.
#' @param output.filename Name for the output spillover file. Default is
#' `"SpectroFlo_spillover_matrix.csv"`.
#' @param fcs.file File name and path to an FCS file to be read. Optional, but
#' inclusion of an FCS file will allow for the detector names (channels) to be
#' added to the returned spillover matrix.
#'
#' @return Spillover matrix (fluorophores x detectors). Fluorophore names will
#' be automatically extracted from the control .fcs files linked to the .Expt
#' file, with an attempt to clean them up. Detector names (columns) are not
#' extracted.
#'
#' @export

read.spectroflo.expt <- function( expt.file, output.dir,
                                  output.filename = "SpectroFlo_spillover_matrix.csv",
                                  fcs.file = NULL ) {

  # read SpectroFlo .Expt file (XML format)
  expt <- xml2::read_xml( expt.file )

  # find the spillover value nodes
  vector.nodes <- xml_find_all( expt, ".//*[local-name() = '_SpilloverVectorArea']" )

  spillover.list <- lapply( vector.nodes, function( node ) {
    float.nodes <- xml_find_all( node, ".//*[local-name()='float']" )
    as.numeric( xml_text( float.nodes ) )
  } )

  # convert to a matrix
  spillover.matrix <- do.call( rbind, spillover.list )

  # spillover.matrix <- t( apply( spillover.matrix, 1, function( x ) x / max( x ) ) )

  # find the associated fluorophore names
  url.nodes <- xml_find_all( expt, ".//*[local-name() = '_Url']" )
  url.paths <- xml_text( url.nodes )

  # extract fluorophore names from file names (before "_Controls.fcs")
  fluor.names <- gsub( ".*\\\\|_Controls\\.fcs.*", "", url.paths )

  # remove leading plate/rack code (e.g., "B3 ") if present
  fluor.names <- gsub( "^[A-Z]\\d+\\s+", "", fluor.names )

  # remove trailing " (Cells)" or similar
  fluor.names <- gsub(" \\(.*\\)", "", fluor.names )

  # extract bit before underscore if present (when library controls are used)
  fluor.names <- gsub( "_.*", "", fluor.names )

  rownames( spillover.matrix ) <- fluor.names

  # if an FCS file is provided, extract detector names and add those
  if ( !is.null( fcs.file ) ) {
    spill.vector <- flowCore::read.FCSheader(
      fcs.file,
      keyword = "$SPILLOVER"
    )[[ 1 ]]

    spill.vector <- strsplit( spill.vector, "," )[[ 1 ]]
    n.detectors <- as.numeric( spill.vector[ 1 ] )
    detector.names <- spill.vector[ 2:( n.detectors + 1 ) ]
    colnames( spillover.matrix ) <- detector.names
  }

  write.csv(
    spillover.matrix,
    file.path( output.dir, output.filename )
  )

  return( spillover.matrix )
}
