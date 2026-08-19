# read_spectra.r

#' @title Read In Saved Spectra
#'
#' @description
#' Reads in CSV files created by `AutoSpectral` or in the same format
#' (fluorophores in rows, detectors in columns, first row is detector names,
#' first column contains fluorophore names).
#'
#' @param spectra.file File name for the spectra CSV file to be read.
#' @param spectra.dir File path to the folder containing `spectra.file`. Default
#' is `table_spectra`, where `AutoSpectral` saves the spectra files.
#' @param remove.af Logical, default is `FALSE`. If `TRUE`, returns the spectral
#' matrix without the default autofluorescence spectrum.
#' @param af.param Name of the autofluorescence parameter. Default is `AF`. Note
#' that any fluorophores can be removed from the matrix by supplying a character
#' vector, e.g., `c("BUV395", "PE")`, if desired.
#' @param check.collinearity Logical, default `TRUE`. CSV files cannot carry
#' the `"fluorophore"` identity attribute that `check.spectra.duplicates()`
#' checks exactly elsewhere, and a `sample`-disambiguated rowname (e.g.
#' `"PE (cells)"` / `"PE (cells) (CD4)"`) is unique by construction, so a
#' rownames-only check would never catch two controls for the same dye. When
#' `TRUE`, screens for suspiciously similar rows as an early, load-time
#' warning (the same screen `check.spectra.duplicates()` falls back to later).
#' @param collinearity.threshold Numeric, default `0.95`. Cosine similarity
#' above which two rows trigger the warning. Matches the package default for
#' `asp$similarity.warning.n`.
#'
#' @return A matrix containing the fluorophore spectra (fluorophore x detectors).
#'
#' @export

read.spectra <- function(
    spectra.file,
    spectra.dir = "./table_spectra",
    remove.af = FALSE,
    af.param = "AF",
    check.collinearity = TRUE,
    collinearity.threshold = 0.95
) {

  spectral.matrix <- as.matrix(
    utils::read.csv(
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

  # turn off collinearity check for AF files
  check.collinearity <- !grepl( "Autofluorescence", spectra.file )

  # see check.spectra.duplicates() -- this is the same cosine-similarity
  # fallback it uses when fluorophore identity is unavailable, run here too
  # so the warning shows up at load time rather than only later at unmixing
  if ( check.collinearity )
    .check.spectra.collinearity(
      spectral.matrix, exclude = af.param, threshold = collinearity.threshold,
      context = paste0( "in '", spectra.file, "'" )
    )

  if ( remove.af )
    spectral.matrix <- spectral.matrix[ !rownames( spectral.matrix ) %in% af.param, ]

  return( spectral.matrix )
}
