# unmix_poisson.r

#' @title Unmix Using Poisson Regression
#'
#' @description
#' This function performs unmixing of raw data using Poisson regression, with
#' iterative reweighted least squares (IRLS) and fallback methods for cells
#' that fail to converge.
#'
#' @param raw.data Matrix containing raw data to be unmixed.
#' @param spectra Matrix containing spectra information.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#' @param initial.weights Optional numeric vector of weights, one per fluorescent
#' detector. Default is `NULL`, in which case weighting will be done by
#' channel means.
#' @param parallel Logical, default is `TRUE`, which enables parallel processing
#' for per-cell unmixing.
#' @param threads Numeric, default is `NULL`, in which case `asp$worker.process.n`
#' will be used if `parallel=TRUE`.
#'
#' @return A matrix containing the unmixed data.
#'
#' @export

unmix.poisson <- function(
    raw.data,
    spectra,
    asp,
    initial.weights = NULL,
    parallel = TRUE,
    threads = NULL
) {

  if ( parallel & is.null( threads ) ) threads <- asp$worker.process.n

  # initialize with WLS unmixing (approximates Poisson weighting)
  unmixed.data <- unmix.wls( raw.data, spectra, initial.weights )

  # handle zeros and negative values
  spectra.t <- t( abs( spectra ) )
  raw.data[ raw.data <= 0 ] <- 1e-6
  unmixed.data[ unmixed.data <= 0 ] <- 1e-6

  # set up parallel processing
  exports <- c( "raw.data", "unmixed.data", "spectra.t", "asp" )

  if ( parallel ) {
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

  # unmix
  unmixed.data <- tryCatch( {
    lapply.function( 1:nrow( raw.data ), function( cell ) {
      raw.data.cell <- raw.data[ cell, ]
      unmixed.data.cell <- unmixed.data[ cell, ]

      # fit glm model with Poisson distribution using identity link function
      fit <- tryCatch(
        suppressWarnings(
          stats::glm.fit(
            spectra.t, raw.data.cell,
            start = unmixed.data.cell,
            family = stats::poisson( link = "identity" ),
            control = stats::glm.control( maxit = asp$rlm.iter.max ),
            intercept = FALSE
          )
        )
      )

      if ( !inherits( fit, "try-error" ) && fit$converged ) {
        return( stats::coef( fit ) )
      } else {
        return( unmixed.data.cell )
      }
    } )
  }, finally = {
    # clean up cluster when done
    if ( !is.null( result$cleanup ) ) result$cleanup()
  } )

  unmixed.data <- do.call( rbind, unmixed.data )

  colnames( unmixed.data ) <- rownames( spectra )

  return( unmixed.data )
}
