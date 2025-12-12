# unmix_wls.r

#' @title Unmix Using Weighted Least Squares
#'
#' @description
#' This function performs unmixing of raw data using weighted least squares,
#' AKA WLS, based on the provided spectra. Weighting is by channel power.
#'
#' @param raw.data Expression data from raw fcs files. Cells in rows and
#' detectors in columns. Columns should be fluorescent data only and must
#' match the columns in spectra.
#' @param spectra Spectral signatures of fluorophores, normalized between 0
#' and 1, with fluorophores in rows and detectors in columns.
#' @param weights Optional numeric vector of weights, one per fluorescent
#' detector. Default is `NULL`, in which case weighting will be done by
#' channel means.
#'
#' @return A matrix containing unnmixed data with cells in rows and
#' fluorophores in columns.
#'
#' @export

unmix.wls <- function( raw.data, spectra, weights = NULL ) {

  spectra <- t( spectra )

  if ( is.null( weights ) ) {
    # weights are inverse of channel variances (mean if Poisson)
    weights <- pmax( abs( colMeans( raw.data ) ), 1e-6 )
    weights <- 1 / weights
    W <- diag( weights )

  } else {
    if ( !is.numeric( weights ) )
      stop( "Weights must be a numeric vector." )

    if ( length( weights ) != nrow( spectra ) )
      stop( "Mismatch between supplied weights and detectors in spectra" )

    if ( length( weights ) != ncol( raw.data ) )
      stop( "Mismatch between supplied weights and detectors in raw.data" )

    W <- diag( as.numeric( weights ) )
  }

  # Weighted LS solution: (M^T W M)^{-1} M^T W
  unmixing.matrix <- solve( t( spectra ) %*% W %*% spectra ) %*%
    ( t( spectra ) %*% W )
  # unmixing.matrix <- solve( t( spectra ) %*% W %*% spectra,
  #   t( spectra ) %*% W )  # alt.

  unmixed.data <- raw.data %*% t( unmixing.matrix )
  # unmixed.data <- tcrossprod(raw.data, unmixing.matrix)  # alt.

  colnames( unmixed.data ) <- colnames( spectra )

  return( unmixed.data )

}
