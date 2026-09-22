# compute_unmix_rmse.R

#' @title Compute Unmixing Reconstruction RMSE
#'
#' @description
#' Computes the root-mean-square error between raw detector data and its
#' reconstruction from an unmixing result, `fitted = unmixed %*% spectra`
#' (plus any fitted autofluorescence contribution, if applicable). This is a
#' direct measure of how well a given unmixing solution explains the raw
#' signal, independent of any FCS header keyword.
#'
#' @param raw.data Numeric matrix, cells in rows and detectors in columns.
#' @param fitted Numeric matrix of the same dimensions as `raw.data`,
#' giving the reconstructed signal for each cell and detector.
#'
#' @return A single numeric value, the RMSE across every cell and detector.
#'
#' @export

compute.unmix.rmse <- function( raw.data, fitted ) {

  raw.data <- as.matrix( raw.data )
  fitted <- as.matrix( fitted )

  if ( !identical( dim( raw.data ), dim( fitted ) ) ) {
    stop( "`raw.data` and `fitted` must have the same dimensions.", call. = FALSE )
  }

  sqrt( mean( ( raw.data - fitted ) ^ 2 ) )
}
