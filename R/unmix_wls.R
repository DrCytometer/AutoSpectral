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

  # set up weights correctly
  if ( is.null( weights ) ) {
    # weights are inverse of channel variances (mean if Poisson)
    weights <- pmax( abs( colMeans( raw.data ) ), 1e-6 )
    weights <- 1 / weights
    W.half <- diag( sqrt( weights ) )

  } else {
    if ( !is.numeric( weights ) )
      stop( "Weights must be a numeric vector." )

    if ( length( weights ) != ncol( spectra ) )
      stop( "Mismatch between supplied weights and detectors in spectra" )

    if ( length( weights ) != ncol( raw.data ) )
      stop( "Mismatch between supplied weights and detectors in raw.data" )

    W.half <- diag( as.numeric( sqrt( weights ) ) )
  }

  # transpose spectra to be detectors x fluorophores
  A <- t( spectra )

  # Weighted design--solving via singular value decomposition for stability
  A.star <- W.half %*% A

  # SVD of weighted design
  sv <- svd( A.star )

  # Weighted pseudoinverse A+_W
  A.pinv.w <- sv$v %*% ( t( sv$u ) / sv$d )

  # Multiply by W.half: (A.star)^+ W.half = (A^T W A)^-1 A^T W
  unmixing.matrix <- A.pinv.w %*% W.half

  unmixed.data <- tcrossprod( raw.data, unmixing.matrix )
  colnames( unmixed.data ) <- rownames( spectra )
  return( unmixed.data )

}
