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
#' events based on scatter profiles matching the positive events.
#' @param k.neighbors Numeric, number of scatter-matched unstained events to
#' pair with every positive event for background determination. Default is `3`.
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
#' @param diagnostics.env Optional environment, default `NULL`. If supplied,
#' `remove.af()` populates it (keyed by sample name) with the objects used to
#' identify and exclude intrusive autofluorescence for each cell-based
#' AF-removal sample: `af.peak.channel`, `fluor.peak`, `af.boundaries`,
#' `expr.data.pos`/`expr.data.neg` (spectral channels only), and the
#' resulting gate indices. Intended for diagnostic/manuscript figures (see
#' `plot.spectra.legacy.steps()`); has no effect on the cleaning result.
#' Capture is unreliable when `parallel = TRUE`.
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
    k.neighbors = 3L,
    intermediate.figures = FALSE,
    main.figures = TRUE,
    parallel = FALSE,
    threads = 1,
    verbose = TRUE,
    diagnostics.env = NULL
) {

  if ( parallel && !is.null( diagnostics.env ) )
    warning(
      "diagnostics.env capture is unreliable under parallel = TRUE (forked ",
      "worker processes do not propagate environment mutations back to the ",
      "caller); diagnostics may be incomplete. Use parallel = FALSE for ",
      "guaranteed-complete diagnostics.", call. = FALSE
    )

  # only ship the expression data being cleaned
  needed.samples <- unique( c( af.removal.sample, universal.negative[ af.removal.sample ] ) )

  # construct arguments list
  args.list <- list(
    clean.expr = clean.expr[ needed.samples ],
    spectral.channel = spectral.channel,
    peak.channel = peak.channel,
    universal.negative = universal.negative,
    asp = asp,
    scatter.param = scatter.param,
    negative.n = negative.n,
    positive.n = positive.n,
    scatter.match = scatter.match,
    k.neighbors = k.neighbors,
    main.figures = main.figures,
    intermediate.figures = intermediate.figures,
    verbose = verbose,
    diagnostics.env = diagnostics.env
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

  # guard against a parallel backend silently returning the wrong number of
  # results, or a per-task failure masquerading as valid expression data
  if ( length( af.remove.expr ) != length( af.removal.sample ) )
    stop(
      "run.af.removal(): parallel backend returned ", length( af.remove.expr ),
      " results for ", length( af.removal.sample ), " requested samples.",
      call. = FALSE
    )

  names( af.remove.expr ) <- af.removal.sample

  is.bad.result <- function( x ) {
    if ( is.null( x ) || inherits( x, "try-error" ) ) return( TRUE )
    if ( !is.matrix( x ) && !is.data.frame( x ) ) return( TRUE )
    nrow( x ) == 0
  }
  failed <- vapply( af.remove.expr, is.bad.result, logical( 1 ) )

  if ( any( failed ) )
    stop(
      "run.af.removal(): worker failed for sample(s): ",
      paste( af.removal.sample[ failed ], collapse = ", " ),
      call. = FALSE
    )

  return( af.remove.expr )
}
