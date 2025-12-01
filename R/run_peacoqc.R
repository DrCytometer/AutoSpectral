# run_peacoqc.r

#' @title Run PeacoQC
#'
#' @description
#' This function runs PeacoQC to remove flow fluctuation errors from expression
#' data using parallel processing if specified.
#'
#' @importFrom future plan multisession
#' @importFrom future.apply future_lapply
#' @importFrom flowWorkspace flowjo_biexp
#'
#' @param expr.data A list containing the expression data for each sample.
#' @param spectral.channel A character vector specifying the spectral channels.
#' @param all.channels A character vector specifying all channels.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#' @param figures Logical, if `TRUE` creates the main figures to show the
#' impact of intrusive autofluorescent event removal and scatter-matching for
#' the negatives.
#' @param parallel Logical, default is `FALSE`, in which case parallel processing
#' will not be used. Parallel processing will likely be faster when many small
#' files are read in. If the data is larger, parallel processing may not
#' accelerate the process much.
#' @param verbose Logical, default is `TRUE`. Set to `FALSE` to suppress messages.
#'
#' @return A list containing the cleaned expression data for each sample.

run.peacoQC <- function( expr.data, spectral.channel, all.channels, asp,
                         figures = TRUE, parallel = FALSE, verbose = TRUE ) {

  # set up parallel processing
  if ( parallel ){
    plan( multisession, workers = asp$worker.process.n )
    options( future.globals.maxSize = asp$max.memory.n )
    lapply.function <- future_lapply
  } else {
    lapply.function <- lapply.sequential
  }

  # define parameters for peacoQC
  biexp.transform <- flowWorkspace::flowjo_biexp(
    channelRange = asp$default.transformation.param$length,
    maxValue = asp$default.transformation.param$max.range,
    pos = asp$default.transformation.param$pos,
    neg = asp$default.transformation.param$neg,
    widthBasis = asp$default.transformation.param$width,
    inverse = FALSE
  )

  transform.inv <- flowWorkspace::flowjo_biexp(
    channelRange = asp$default.transformation.param$length,
    maxValue = asp$default.transformation.param$max.range,
    pos = asp$default.transformation.param$pos,
    neg = asp$default.transformation.param$neg,
    widthBasis = asp$default.transformation.param$width,
    inverse = TRUE
  )

  time.param <- asp$default.time.parameter

  # run peacoQC to remove flow fluctuation errors
  clean.expr <- lapply.function( names( expr.data ), function( sample.name ) {
    do.peacoQC( expr.data[[ sample.name ]], sample.name,
                spectral.channel, biexp.transform, transform.inv,
                asp$figure.peacoqc.dir, time.param, all.channels,
                asp$peacoqc.method, figures, verbose )

  }, future.seed = asp$gate.downsample.seed )

  # note that PeacoQC will only run on MAD for low n samples like controls
  # do.peacoQC is therefore set to run with MAD only
  names( clean.expr ) <- names( expr.data )

  rm( expr.data )

  # issue warning if fewer than 500 events for any sample
  clean.expr.n <- sapply( clean.expr, nrow )

  low.sample.n <- which( clean.expr.n < asp$min.cell.warning.n )

  if ( any( clean.expr.n < asp$min.cell.warning.n ) ){
    warning( paste( "\033[31m", "Warning! Fewer than", asp$min.cell.warning.n,
                "gated events in", names( low.sample.n ), "\033[0m", "\n"  ) )
  }

  return( clean.expr )

}
