# unmix_nnls.r

#' @title unmix.nnls
#'
#' @description
#' Performs spectral unmixing using non negative least squares. This adds a
#' positive constraint on the contribution of the markers, i.e. the resulting
#' coefficients must be greater than or equal to zero.
#'
#' @param raw.data Expression data from raw fcs files. Cells in rows and
#' detectors in columns. Columns should be fluorescent data only and must
#' match the columns in spectra.
#' @param spectra Spectral signatures of fluorophores, normalized between 0
#' and 1, with fluorophores in rows and detectors in columns.
#' @param weights Dummy argument to allow dynamic switching between nnls and WLS.
#' Default is `NULL`. Values passed to `weights` will be ignored.
#'
#' @return Unmixed data with cells in rows and fluorophores in columns.
#'
#' @export

unmix.nnls <- function( raw.data, spectra, weights = NULL ) {

  unmixed.data <- RcppML::nnls(
    crossprod( t(spectra) ), tcrossprod( spectra, raw.data ),
    fast_nnls = FALSE, cd_maxit = 19)
  unmixed.data <- t( unmixed.data )

  colnames( unmixed.data ) <- rownames( spectra )

  return( unmixed.data )

}
