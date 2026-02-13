# unmix_ols.r

#' @title Unmix OLS
#'
#' @description
#' Performs spectral unmixing using ordinary least squares
#'
#' @param raw.data Expression data from raw fcs files. Cells in rows and
#' detectors in columns. Columns should be fluorescent data only and must
#' match the columns in spectra.
#' @param spectra Spectral signatures of fluorophores, normalized between 0
#' and 1, with fluorophores in rows and detectors in columns.
#' @param weights Dummy argument to allow dynamic switching between OLS and WLS.
#' Default is `NULL`. Values passed to `weights` will be ignored.
#'
#' @return Unmixed data with cells in rows and fluorophores in columns.
#'
#' @export

unmix.ols <- function( raw.data, spectra, weights = NULL ) {

  # check for data/spectra column matching
  raw.data.cols <- colnames( raw.data )
  spectra.cols <- colnames( spectra )

  if ( !identical( raw.data.cols, spectra.cols ) ) {

    # ensure both actually have the same columns before reordering
    if ( all( spectra.cols %in% raw.data.cols ) &&
         length( spectra.cols ) == length( raw.data.cols ) ) {
      # reorder raw.data to match the order of spectra
      raw.data <- raw.data[, spectra.cols]
      message( "Columns reordered to match spectra." )
    } else {
      stop( "Column names in spectra and raw.data do not match perfectly;
           cannot reorder by name alone." )
    }
  }

  sv <- svd( t( spectra ) )
  # solve unmixing matrix via singular value decomposition
  # more stable than normal equations
  # this is the Moore-Penrose pseudoinverse
  unmixing.matrix <- sv$v %*% ( t( sv$u ) / sv$d )

  unmixed.data <- tcrossprod( raw.data, unmixing.matrix )
  colnames( unmixed.data ) <- rownames( spectra )
  return( unmixed.data )

}
