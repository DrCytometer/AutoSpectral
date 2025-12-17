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
#' detector. Default is `NULL`, in which case weighting aka scaling is
#' based on channel top values and a bottom value among the top values.
#' @param qtop quantile to identify the top value of each channel. Maximum is
#' an unstable measure, so the 99th percentile is used by default.
#' @param qbot quantile to identify the constant to add to moderate the
#' weighting. For the channel with the lowest top, 0 results in a weight of
#' about 2 whereas 1 results in a large weight. The channel with the highest
#' top value has always a weight of 1. WLS results in much larger weights than
#' those obtained with qbot = 1.
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
    weights <- ( tops + quantile( tops, probs = qbot ) ) / tops
    # set minimum to 1, i.e. no change
    # the minimum is the channel with the highest top
    weights <- weights / min(weights)

  }

  return( unmix.wls( raw.data, spectra, weights ) )

}
