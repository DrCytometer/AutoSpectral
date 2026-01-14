# run_universal_negative.r

#' @title Run Universal Negative
#' @description
#' This function processes universal negative samples to generate expression
#' data based on specified parameters.
#'
#' @param clean.expr A matrix containing cleaned expression data.
#' @param univ.sample A character vector specifying the names of universal
#' negative samples.
#' @param universal.negatives A list containing universal negative control
#' parameters.
#' @param scatter.param A character vector specifying the scatter parameters.
#' @param peak.channels A character vector specifying the peak channels.
#' @param downsample A numeric value indicating the downsampling factor.
#' @param negative.n A numeric value indicating the number of negative events.
#' @param positive.n A numeric value indicating the number of positive events.
#' @param spectral.channel A character vector specifying the spectral channels.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#' @param control.type A character string specifying the
#' type of control: `beads` or `cells`
#' @param scatter.match A logical value indicating whether scatter matching
#' is performed.
#' @param intermediate.figures Logical, if `TRUE` returns additional figures to
#' show the inner workings of the cleaning, including definition of low-AF cell
#' gates on the PCA-unmixed unstained and spectral ribbon plots of the AF
#' exclusion from the unstained. Default is `FALSE` to speed up processing.
#' @param main.figures Logical, if `TRUE` creates the main figures to show the
#' impact of intrusive autofluorescent event removal and scatter-matching for
#' the negatives.
#' @param verbose Logical, default is `TRUE`. Set to `FALSE` to suppress messages.
#'
#' @return A list containing the processed expression data for each universal
#' negative sample.

run.universal.negative <- function( clean.expr, univ.sample,
                                    universal.negatives,
                                    scatter.param,
                                    peak.channels, downsample,
                                    negative.n, positive.n,
                                    spectral.channel, asp,
                                    control.type,
                                    scatter.match,
                                    intermediate.figures = FALSE,
                                    main.figures = TRUE,
                                    verbose = TRUE ) {

  univ.expr <- lapply( univ.sample, function( sample.name ) {

    get.universal.negative(
      clean.expr,
      sample.name,
      universal.negatives,
      scatter.param,
      peak.channels,
      downsample,
      negative.n,
      positive.n,
      spectral.channel,
      asp,
      control.type,
      scatter.match,
      intermediate.figures,
      main.figures,
      verbose
    )

  } )

  names( univ.expr ) <- univ.sample

  rm( clean.expr )

  return( univ.expr )
}
