# run_peacoqc.r

#' @title Run PeacoQC
#'
#' @description
#' This function runs PeacoQC to remove flow fluctuation errors from expression
#' data using parallel processing if specified.
#'
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
#' will not be used. Set to `TRUE` to run in parallel.
#' @param threads Number of cores to use for parallel processing, default is `1`.
#' @param verbose Logical, default is `TRUE`. Set to `FALSE` to suppress messages.
#'
#' @return A list containing the cleaned expression data for each sample.

run.peacoQC <- function( expr.data,
                         spectral.channel,
                         all.channels,
                         asp,
                         figures = TRUE,
                         parallel = FALSE,
                         threads = 1,
                         verbose = TRUE ) {

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

  # construct arguments list
  args.list <- list(
    spectral.channel = spectral.channel,
    biexp.transform = biexp.transform,
    transform.inv = transform.inv,
    output.dir = asp$figure.peacoqc.dir,
    time.param = asp$default.time.parameter,
    all.channels = all.channels,
    method = asp$peacoqc.method,
    figures = figures,
    verbose = verbose
  )

  # set up parallel processing
  if ( parallel ) {
    internal.functions <- c( "do.peacoQC" )
    exports <- c( "args.list", "expr.data", internal.functions )
    result <- create.parallel.lapply(
      asp,
      exports,
      parallel = parallel,
      threads = threads,
      export.env = environment(),
      allow.mclapply.mac = TRUE
    )
    lapply.function <- result$lapply
  } else {
    lapply.function <- lapply
    result <- list( cleanup = NULL )
  }

  # run peacoQC to remove flow fluctuation errors
  clean.expr <- tryCatch( {
    lapply.function( names( expr.data ), function( f ) {
      do.call( do.peacoQC, c( list( expr.data[[ f ]], f ), args.list ) )
    } )
  }, finally = {
    # clean up cluster when done
    if ( !is.null( result$cleanup ) ) result$cleanup()
  } )

  # note that PeacoQC will only run on MAD for low n samples like controls
  # do.peacoQC is therefore set to run with MAD only
  names( clean.expr ) <- names( expr.data )

  rm( expr.data )

  # issue warning if fewer than 500 events for any sample
  clean.expr.n <- sapply( clean.expr, nrow )

  low.sample.n <- which( clean.expr.n < asp$min.cell.warning.n )

  if ( any( clean.expr.n < asp$min.cell.warning.n ) ) {
    warning(
      paste0(
        "\033[31m",
        "Warning! Fewer than ",
        asp$min.cell.warning.n,
        " gated events in ",
        names( low.sample.n ),
        "\033[0m",
        "\n"
      )
    )
  }

  return( clean.expr )
}
