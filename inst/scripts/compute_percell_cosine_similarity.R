# compute_percell_cosine_similarity.R

#' @title Per-Cell Cosine Similarity to a Reference Spectrum
#'
#' @description
#' Computes the cosine similarity between each event's raw detector vector
#' and a single fixed reference spectrum, with no unmixing involved. This is
#' the "OLS" comparison arm when assessing how closely a single-stained
#' control's raw signal follows its library spectrum: OLS unmixing assumes
#' every event shares that one spectrum, so this is exactly the assumption
#' being checked.
#'
#' Mathematically, for event \eqn{i}, the reported value is
#' \eqn{\cos(\theta_i) = \dfrac{r_i \cdot s}{\|r_i\| \|s\|}}, where \eqn{r_i}
#' is the raw detector vector for event \eqn{i} and \eqn{s} is the reference
#' spectrum. In lay terms: ignoring how bright the cell is, does the shape
#' of what it lit up across every detector match the textbook shape for
#' this fluorophore?
#'
#' @importFrom lsa cosine
#'
#' @param raw.data Numeric matrix, cells in rows and detectors in columns.
#' @param reference.spectrum Numeric vector, the reference spectrum, in the
#' same detector order as the columns of `raw.data`.
#'
#' @return A numeric vector, one cosine similarity value per row of
#' `raw.data`.
#'
#' @export

per.cell.cosine.similarity <- function( raw.data, reference.spectrum ) {

  reference.spectrum <- as.numeric( reference.spectrum )

  vapply(
    seq_len( nrow( raw.data ) ),
    function( i ) lsa::cosine( reference.spectrum, as.numeric( raw.data[ i, ] ) ),
    numeric( 1 )
  )
}


#' @title Per-Cell Cosine Similarity to the Best-Fit Fluorophore Spectrum
#'
#' @description
#' For a single-stained control, assigns each event its own autofluorescence
#' spectrum (`assign.af.fluorophores()`, then `unmix.af.fwl()`) and, where
#' the target fluorophore's precomputed spectral variants
#' (`get.spectral.variants()`) reduce the residual fit against that event's
#' AF-subtracted signal, its own variant of the target fluorophore's
#' spectrum. Returns, for every event, the cosine similarity between its
#' AF-subtracted detector vector and whichever target-fluorophore spectrum -
#' reference or variant - fit it best.
#'
#' This is the "per-cell fluorophore optimization" counterpart to
#' `per.cell.cosine.similarity()`'s single fixed reference spectrum: instead
#' of asking how well one textbook spectrum describes every event, it asks
#' how well the *best available* spectrum for that event - chosen from a
#' library of spectral variants already measured for this fluorophore -
#' describes it, once that event's own autofluorescence has been removed.
#'
#' For event \eqn{i}, restricted to its own positive fluorophores
#' \eqn{P_i} (from `spectra.variants$thresholds`), each candidate spectrum
#' \eqn{v} for the target fluorophore is scored by the L1 residual
#' \eqn{\| r_i - \hat{r}_i(v) \|_1} of re-unmixing event \eqn{i}'s
#' AF-subtracted signal \eqn{r_i} against \eqn{P_i} with the target row set
#' to \eqn{v}; the variant minimizing this residual (the reference spectrum
#' is always variant 1, so it is never worse off than plain OLS) is kept,
#' and its cosine similarity to \eqn{r_i} is returned. In lay terms: rather
#' than judging every cell against one average fluorophore signature, each
#' cell is allowed to pick whichever measured version of that signature it
#' actually looks most like, and the reported number is how good that best
#' match is.
#'
#' @importFrom lsa cosine
#'
#' @param raw.data Numeric matrix, cells in rows and detectors in columns,
#' for a single-stained control of `target.fluor`.
#' @param spectra Numeric matrix, fluorophore reference spectra
#' (fluorophores in rows, detectors in columns). An `AF` row, if present, is
#' dropped.
#' @param af.spectra Numeric matrix of autofluorescence spectra, as prepared
#' by `get.af.spectra()`.
#' @param spectra.variants The list returned by `get.spectral.variants()`,
#' carrying `$variants` (named list of per-fluorophore variant-spectra
#' matrices) and `$thresholds` (named vector of unmixed positivity
#' thresholds).
#' @param target.fluor Character, the name of the fluorophore in `spectra`
#' (and in `spectra.variants$variants`) that `raw.data` is a single-stained
#' control for.
#' @param asp The AutoSpectral parameter list, used to set up parallel
#' workers via `create.parallel.lapply()`.
#' @param parallel Logical, default `TRUE`.
#' @param threads Numeric or `NULL`. Default `NULL` uses
#' `asp$worker.process.n`.
#'
#' @seealso
#' * [per.cell.cosine.similarity()]
#' * [assign.af.fluorophores()]
#' * [unmix.af.fwl()]
#' * [get.spectral.variants()]
#'
#' @return A numeric vector, one cosine similarity value per row of
#' `raw.data`.
#'
#' @export

compute.percell.variant.cosine.similarity <- function(
    raw.data,
    spectra,
    af.spectra,
    spectra.variants,
    target.fluor,
    asp,
    parallel = TRUE,
    threads = NULL
) {

  if ( "AF" %in% rownames( spectra ) )
    spectra <- spectra[ rownames( spectra ) != "AF", , drop = FALSE ]

  fluorophores <- rownames( spectra )

  if ( !( target.fluor %in% fluorophores ) )
    stop( "`target.fluor` is not a row of `spectra`.", call. = FALSE )

  target.variants <- spectra.variants$variants[[ target.fluor ]]
  if ( is.null( target.variants ) )
    stop(
      paste0( "No spectral variants found for `target.fluor` = \"", target.fluor,
              "\" in `spectra.variants$variants`." ),
      call. = FALSE
    )

  pos.thresholds <- spectra.variants$thresholds

  # per-cell autofluorescence assignment and subtraction
  af.index <- assign.af.fluorophores( raw.data, spectra, af.spectra )
  af.fit <- unmix.af.fwl(
    raw.data = raw.data, spectra = spectra, af.spectra = af.spectra,
    af.index = af.index, return.fitted.af = TRUE
  )

  remaining.raw <- raw.data - af.fit$fitted.af
  baseline.error <- rowSums( abs(
    raw.data - ( af.fit$fluorophores %*% spectra ) - af.fit$fitted.af
  ) )

  threads <- if ( isTRUE( parallel ) ) {
    if ( is.null( threads ) ) asp$worker.process.n else as.integer( threads )
  } else {
    1L
  }

  worker.exports <- c(
    "spectra", "fluorophores", "target.fluor", "target.variants",
    "pos.thresholds", "remaining.raw", "baseline.error", "af.fit"
  )

  result.setup <- create.parallel.lapply(
    asp,
    exports = worker.exports,
    parallel = parallel,
    threads = threads,
    export.env = environment()
  )
  lapply.function <- result.setup$lapply

  process.cell <- function( cell ) {

    cell.raw <- remaining.raw[ cell, ]
    cell.abundance <- af.fit$fluorophores[ cell, fluorophores ]
    error.final <- baseline.error[ cell ]

    pos.fluors <- fluorophores[ cell.abundance >= pos.thresholds[ fluorophores ] ]

    if ( !( target.fluor %in% pos.fluors ) )
      return( lsa::cosine( spectra[ target.fluor, ], cell.raw ) )

    cell.spectra.final <- spectra
    cell.spectra.curr <- spectra[ pos.fluors, , drop = FALSE ]

    for ( v in seq_len( nrow( target.variants ) ) ) {

      cell.spectra.curr[ target.fluor, ] <- target.variants[ v, ]

      unmixed.curr <- unmix.ols( matrix( cell.raw, nrow = 1 ), cell.spectra.curr )
      error.curr <- sum( abs(
        cell.raw - as.numeric( unmixed.curr %*% cell.spectra.curr )
      ) )

      if ( error.curr < error.final ) {
        error.final <- error.curr
        cell.spectra.final[ target.fluor, ] <- target.variants[ v, ]
      }
    }

    lsa::cosine( cell.spectra.final[ target.fluor, ], cell.raw )
  }

  similarity <- tryCatch(
    expr = unlist( lapply.function( seq_len( nrow( raw.data ) ), process.cell ) ),
    finally = {
      if ( !is.null( result.setup$cleanup ) ) result.setup$cleanup()
    }
  )

  similarity
}
