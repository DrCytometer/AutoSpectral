# unmix_ols_pos.r

#' @title unmix.ols_pos
#'
#' @description
#' Performs spectral unmixing using ordinary least squares on that are
#' previously rendered positive or zero.
#'
#' @param raw.data Expression data from raw fcs files. Cells in rows and
#' detectors in columns. Columns should be fluorescent data only and must
#' match the columns in spectra.
#' @param spectra Spectral signatures of fluorophores, normalized between 0
#' and 1, with fluorophores in rows and detectors in columns.
#' @param weights Dummy argument to allow dynamic switching between ols_pos and WLS.
#' Default is `NULL`. Values passed to `weights` will be ignored.
#'
#' @return Unmixed data with cells in rows and fluorophores in columns.
#'
#' @export

unmix.ols_pos <- function( raw.data, spectra, weights = NULL ) {

  raw.data <- raw.data
  raw.data[ raw.data < 0 ] <- 0

  unmixed.data <- solve( crossprod( t( spectra ) ), tcrossprod( spectra, raw.data))
  unmixed.data <- t( unmixed.data )

  return( unmixed.data )

}
