# run_downsample.r

#' @title Run Downsample
#'
#' @description
#' This function runs the downsampling process on a list of samples, using the
#' specified peak channels and parameters.
#'
#' @param clean.expr.data List containing cleaned expression data.
#' @param downsample.sample Vector of sample names to be downsampled.
#' @param peak.channels Named vector mapping samples to their peak channels.
#' @param negative.n Number of negative events to select, default `500`.
#' @param positive.n Number of positive events to select, default `1000`.
#' @param verbose Logical, default is `TRUE`. Set to `FALSE` to suppress messages.
#'
#' @return A list containing the downsampled expression data for each sample.

run.downsample <- function( clean.expr.data, downsample.sample, peak.channels,
                            negative.n = 500, positive.n = 1000,
                            verbose = TRUE ) {

  downsample.expr <- lapply( downsample.sample, function( sample.name ){

    downsample.control(
      clean.expr.data,
      sample.name,
      peak.channels,
      negative.n,
      positive.n,
      verbose
    )
  } )

  names( downsample.expr ) <- downsample.sample

  return( downsample.expr )
}
