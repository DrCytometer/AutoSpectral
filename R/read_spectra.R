# read_spectra.r

#' @title Read In Saved Spectra
#'
#' @description
#' Reads in CSV files created by `AutoSpectral` or in the same format
#' (fluorophores in rows, detectors in columns, first row is detector names,
#' first column contains fluorophore names).
#'
#' @importFrom utils read.csv
#'
#' @param spectra.file File name for the spectra CSV file to be read.
#' @param spectra.dir File path to the folder containing `spectra.file`. Default
#' is `table_spectra`, where `AutoSpectral` saves the spectra files.
#' @param remove.af Logical, default is `FALSE`. If `TRUE`, returns the spectral
#' matrix without the default autofluorescence spectrum.
#' @param af.param Name of the autofluorescence parameter. Default is `AF`. Note
#' that any fluorophores can be removed from the matrix by supplying a character
#' vector, e.g., `c("BUV395", "PE")`, if desired.
#'
#' @return A matrix containing the fluorophore spectra (fluorophore x detectors).
#'
#' @export

read.spectra <- function( spectra.file,
                          spectra.dir = "./table_spectra",
                          remove.af = FALSE,
                          af.param = "AF" ) {

  spectral.matrix <- as.matrix(
    read.csv(
      file.path( spectra.dir, spectra.file ),
      check.names = FALSE,
      row.names = 1
    )
  )

  # normalize
  spectral.matrix <- t(
    apply(
      spectral.matrix, 1, function( x ) x / max( x ) )
  )

  if ( remove.af )
    spectral.matrix <- spectral.matrix[ !rownames( spectral.matrix ) %in% af.param, ]

  return( spectral.matrix )
}
