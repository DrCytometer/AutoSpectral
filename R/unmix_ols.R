# unmix_ols.r

#' @title unmix.ols
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

  spectra <- t( spectra )

  unmixing.matrix <- solve( crossprod( spectra ) ) %*% t( spectra )
  # unmixing.matrix <- solve( crossprod( spectra ), t( spectra ) )  # alt.

  unmixed.data <- raw.data %*% t( unmixing.matrix )
  # unmixed.data <- tcrossprod(raw.data, unmixing.matrix)  # alt.

  colnames( unmixed.data ) <- colnames( spectra )

  return( unmixed.data )

}
