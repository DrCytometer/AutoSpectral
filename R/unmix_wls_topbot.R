# unmix_wls_qtopbot.r

#' @title Unmix Using Weighted Least Squares
#'
#' @description
#' This function performs unmixing of raw data using weighted least squares,
#' AKA WLS, based on the provided spectra. Weighting is by channel. The top
#' value of each channel is considered here as a measure of importance of the
#' channel. Channels with low top values have less importance, but more than in
#' Ordinary Least Squares.
#'
#' @param raw.data Expression data from raw fcs files. Cells in rows and
#' detectors in columns. Columns should be fluorescent data only and must
#' match the columns in spectra.
#' @param spectra Spectral signatures of fluorophores, normalized between 0
#' and 1, with fluorophores in rows and detectors in columns.
#' @param weights Optional numeric vector of weights, one per fluorescent
#' detector. Default is `NULL`, in which case weighting will be done by
#' channel means.
#' @param qtop quantile to identify the top value of each channel. Maximum is
#' an unstable measure.
#' @param qbot quantile to identify the constant to add to the denominator of
#' the weight formula.
#'
#' @return A matrix containing unnmixed data with cells in rows and
#' fluorophores in columns.
#'
#' @export

unmix.wls_qtopbot <- function(
    raw.data, spectra, weights = NULL, qtop = 0.99, qbot = 0.25
) {

  if ( is.null( weights ) ) {

    # get the top value of each channel; maximum is unstable
    tops <- apply( raw.data, 2, quantile, probs = qtop )
    # converts tops to weights: perfer linear to square, purely arbitrary
    # weights <- (weights)**2 / (weights**2 + quantile(weights**2, probs = qbot))
    weights <- (tops) / (tops + quantile(tops, probs = qbot))

  }

  return( unmix.wls( raw.data, spectra, weights ) )

}
