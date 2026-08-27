# unmix_wls_fast.r

#' @title Unmix WLS Fast
#'
#' @description
#' Faster solver for per-cell optimization workflow. Performs spectral unmixing
#' using weighted least squares.
#'
#' @param raw.data Expression data from raw FCS files. Cells in rows and
#' detectors in columns. Columns should be fluorescent data only and must
#' match the columns in spectra.
#' @param spectra Spectral signatures of fluorophores, normalized between 0
#' and 1, with fluorophores in rows and detectors in columns.
#' @param weights The weighting values for the detectors. No checks are
#' performed on supplied weights; `unmix.wls` should be used for most cases.
#' Default is `NULL`, in which case weighting is computed from channel means
#' of `raw.data`, floored at `noise.floor`.
#' @param noise.floor Numeric, default `125`. Lower clamp on mean channel
#' signal before inversion when `weights = NULL`. Signal units, same
#' convention as `noise.floor` elsewhere in the package. Ignored if `weights`
#' is supplied.
#'
#' @return A matrix containing unmixed data with cells in rows and fluorophores in columns.
#'
#' @export

unmix.wls.fast <- function( raw.data, spectra, weights = NULL, noise.floor = 125 ) {

  if ( is.null( weights ) ) {
    weights <- 1 / pmax( abs( colMeans( raw.data ) ), noise.floor )
  }

  Sw  <- sweep( spectra, 2, weights, "*" )
  XtX <- Sw %*% t( spectra )

  # (S W Sᵀ)⁻¹ S W
  unmixing.matrix <- solve( XtX, Sw )

  raw.data %*% t( unmixing.matrix )

}
