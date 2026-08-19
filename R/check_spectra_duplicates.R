# check_spectra_duplicates.r

#' @title Check Spectra for Duplicate Fluorophores
#'
#' @description
#' Verifies that a spectral reference matrix intended for unmixing contains
#' exactly one row per fluorophore. Multiple controls per fluorophore are
#' permitted upstream (during control file validation, spectrum extraction,
#' and QC comparison) via `sample`-based disambiguation, but two or more rows
#' representing the *same* fluorophore must never reach the unmixing solvers:
#' near-identical reference columns are highly collinear and will badly
#' degrade unmixing conditioning. This function is the single enforcement
#' point that catches that before any unmixing math runs, and is called
#' automatically from `unmix.fcs()`, `unmix.folder()`, `unmix.autospectral()`,
#' `unmix.autospectral.joint()`, `unmix.autospectral.rcpp()`, `unmix.gls()`
#' and `unmix.af.gls()`.
#'
#' If `spectra` carries a `"fluorophore"` attribute (set automatically by
#' `get.fluorophore.spectra()` and `get.spectra.automated()`), true
#' fluorophore identity per row is checked exactly, and a duplicate stops
#' execution. If that attribute is absent -- notably, after `spectra` has
#' been round-tripped through a CSV via `read.spectra()`, which cannot carry
#' the attribute, and where a `sample`-disambiguated rowname is unique by
#' construction and so cannot reveal a duplicate either -- this falls back to
#' a pairwise cosine-similarity screen as a cheap, imperfect surrogate. That
#' screen cannot distinguish "two controls for the same fluorophore" from
#' "two genuinely similar fluorophores", so it warns rather than stopping.
#'
#' @param spectra A matrix of spectral reference data (fluorophores x
#' detectors), as used by `unmix.fcs()` and related functions.
#' @param collinearity.threshold Numeric, default `0.95`. Cosine similarity
#' above which two rows trigger the fallback warning, used only when
#' `spectra` has no `"fluorophore"` attribute. Matches the package default
#' for `asp$similarity.warning.n`.
#'
#' @return Invisibly returns `TRUE` if no duplicates (or, in fallback mode, no
#' high-similarity pairs) are found. Stops with an informative error when an
#' exact duplicate is found via the attribute.
#'
#' @export

check.spectra.duplicates <- function( spectra, collinearity.threshold = 0.95 ) {

  if ( is.null( spectra ) ) return( invisible( TRUE ) )

  fluorophore <- attr( spectra, "fluorophore" )

  if ( !is.null( fluorophore ) ) {

    dup.fluor <- unique( fluorophore[ duplicated( fluorophore ) ] )

    if ( length( dup.fluor ) > 0 ) {
      offending.rows <- rownames( spectra )[ fluorophore %in% dup.fluor ]
      stop(
        paste0(
          "Spectral reference matrix contains more than one row for the ",
          "following fluorophore(s): ", paste( dup.fluor, collapse = ", " ), " ",
          "(rows: ", paste( offending.rows, collapse = ", " ), "). ",
          "Multiple controls per fluorophore are only permitted upstream, for ",
          "QC/comparison purposes (`allow.duplicate.controls = TRUE` in ",
          "`get.fluorophore.spectra()` / `get.spectra.automated()`). Reduce ",
          "the reference library to a single row per fluorophore (e.g. by ",
          "choosing the better-conditioned control, or averaging) before ",
          "unmixing."
        ),
        call. = FALSE
      )
    }

  } else {

    # No fluorophore identity available. Fall back to a cosine-similarity
    # screen; see the class-level description for why this warns rather
    # than stopping.
    .check.spectra.collinearity(
      spectra, exclude = "AF", threshold = collinearity.threshold,
      context = "in the spectral reference matrix"
    )
  }

  invisible( TRUE )
}

#' @title Screen a Spectral Matrix for Suspiciously Similar Rows
#' @description
#' Internal helper shared by `check.spectra.duplicates()` (as its fallback
#' when fluorophore identity is unavailable) and `read.spectra()` (as an
#' early, load-time check). Not exported.
#' @keywords internal
.check.spectra.collinearity <- function( spectra, exclude = "AF",
                                         threshold = 0.95, context = "" ) {

  check.rows <- rownames( spectra )[ !rownames( spectra ) %in% exclude ]
  if ( length( check.rows ) < 2 ) return( invisible( TRUE ) )

  sim.mat  <- cosine.similarity( spectra[ check.rows, , drop = FALSE ] )
  uniq.sim <- sim.mat * lower.tri( sim.mat )
  sim.idx  <- which( uniq.sim > threshold, arr.ind = TRUE )

  if ( nrow( sim.idx ) > 0 ) {
    f1   <- rownames( sim.mat )[ sim.idx[ , 1 ] ]
    f2   <- colnames( sim.mat )[ sim.idx[ , 2 ] ]
    vals <- sim.mat[ sim.idx ]
    warning(
      paste0(
        "High cosine similarity (> ", threshold, ") detected between rows",
        if ( nzchar( context ) ) paste0( " ", context ) else "", ": ",
        paste( sprintf( "%s <-> %s (%.4f)", f1, f2, vals ), collapse = "; " ), ". ",
        "This may mean two controls for the same fluorophore are both ",
        "present (only one should be used for unmixing), or that two ",
        "distinct fluorophores simply have similar spectra. ",
        "Inspect before unmixing."
      ),
      call. = FALSE
    )
    return( invisible( FALSE ) )
  }

  invisible( TRUE )
}
