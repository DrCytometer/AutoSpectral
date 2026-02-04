# run_af_removal.r

#' @title Run Autofluorescence Removal
#'
#' @description
#' This function runs the autofluorescence removal process on a list of samples,
#' using the specified parameters and settings.
#'
#' @param clean.expr List containing cleaned expression data.
#' @param af.removal.sample Vector of sample names for which autofluorescence
#' removal is to be performed.
#' @param spectral.channel Vector of spectral channel names.
#' @param peak.channel Vector of peak detection channels for fluorophores.
#' @param universal.negative Name of the universal negative control.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#' @param scatter.param Vector of scatter parameters.
#' @param negative.n Integer. Number of events to include in the downsampled
#' negative population. Default is `500`.
#' @param positive.n Integer. Number of events to include in the downsampled
#' positive population. Default is `1000`.
#' @param scatter.match Logical, default is `TRUE`. Whether to select negative
#' events based on scatter profiles matching the positive events. Defines a
#' region of FSC and SSC based on the distribution of selected positive events.
#' @param intermediate.figures Logical, if `TRUE` returns additional figures to
#' show the inner workings of the cleaning, including definition of low-AF cell
#' gates on the PCA-unmixed unstained and spectral ribbon plots of the AF
#' exclusion from the unstained.
#' @param main.figures Logical, if `TRUE` creates the main figures to show the
#' impact of intrusive autofluorescent event removal and scatter-matching for
#' the negatives.
#' @param parallel Logical, default is `FALSE`, in which case parallel processing
#' will not be used. Set to `TRUE` to run in parallel.
#' @param threads Number of cores to use for parallel processing, default is `1`.
#' @param verbose Logical, default is `TRUE`. Set to `FALSE` to suppress messages.
#'
#' @return A list containing the expression data with autofluorescent events
#' removed for each sample.

run.af.removal <- function(
    clean.expr,
    af.removal.sample,
    spectral.channel,
    peak.channel,
    universal.negative,
    asp,
    scatter.param,
    negative.n = 500,
    positive.n = 1000,
    scatter.match = TRUE,
    intermediate.figures = FALSE,
    main.figures = TRUE,
    parallel = FALSE,
    threads = 1,
    verbose = TRUE
  ) {

  # construct arguments list
  args.list <- list(
    clean.expr = clean.expr,
    spectral.channel = spectral.channel,
    peak.channel = peak.channel,
    universal.negative = universal.negative,
    asp = asp,
    scatter.param = scatter.param,
    negative.n = negative.n,
    positive.n = positive.n,
    scatter.match = scatter.match,
    main.figures = main.figures,
    intermediate.figures = intermediate.figures,
    verbose = verbose
  )

  # set up parallel processing
  if ( parallel ) {
    if ( verbose ) {
      message( "\033[34mIdentifying and removing autofluorescence contamination \033[0m" )
      message( "\033[34mCheck `figure_clean_controls` for plots \033[0m" )
    }

    internal.functions <- c( "remove.af" )
    exports <- c( "args.list", "af.removal.sample", internal.functions )
    result <- create.parallel.lapply(
      asp,
      exports,
      parallel = parallel,
      threads = threads,
      export.env = environment()
    )
    lapply.function <- result$lapply
  } else {
    lapply.function <- lapply
    result <- list( cleanup = NULL )
  }

  # main loop
  af.remove.expr <- tryCatch( {
    lapply.function( af.removal.sample, function( s ) {
      do.call( remove.af, c( list( s ), args.list ) )
    } )
  }, finally = {
    # clean up cluster when done
    if ( !is.null( result$cleanup ) ) result$cleanup()
  } )

  names( af.remove.expr ) <- af.removal.sample

  return( af.remove.expr )
}
