# unmix_ols_fast.r

#' @title Unmix OLS Fast
#'
#' @description
#' Faster solver for per-cell optimization workflow. Performs spectral unmixing
#' using ordinary least squares.
#'
#' @param raw.data Expression data from raw FCS files. Cells in rows and
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

unmix.ols.fast <- function( raw.data, spectra, weights = NULL ) {

  XtX <- tcrossprod( spectra )

  unmixing.matrix <- solve.default( XtX, spectra )

  raw.data %*% t( unmixing.matrix )

}
