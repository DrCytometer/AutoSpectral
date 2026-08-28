#' @title Spectral-Location Alignment Between Variability and Mismatch
#'
#' @description
#' For each fluorophore common to both \code{variability.mad} and
#' \code{mismatch.dist}, computes the cosine similarity between that
#' fluorophore's per-detector variability profile (e.g. the denoised MAD
#' profile from \code{assess.variability.mad()}) and the per-detector
#' magnitude of its bead-vs-cell mismatch (\code{abs()} of, e.g.,
#' \code{bead.cell.dist()}). Both profiles are non-negative by
#' construction, so the resulting cosine similarity is bounded in
#' \code{[0, 1]} and reflects purely where in detector space each
#' quantity is concentrated -- a value near 1 means the fluorophore's
#' variant-to-variant variability and its bead-vs-cell mismatch peak at
#' the same detectors, regardless of either quantity's overall magnitude;
#' a value near 0 means they are concentrated in different, non-
#' overlapping parts of the spectrum.
#'
#' \code{variability.mad} and \code{mismatch.dist} are expected to come
#' from independent computations and are not guaranteed to share the same
#' detector set or column order, so the comparison is aligned explicitly
#' by detector name rather than position. A warning is issued if the two
#' detector sets differ.
#'
#' @param variability.mad Numeric matrix, fluorophores in rows and
#'   detectors in columns, as returned by \code{assess.variability.mad()}.
#' @param mismatch.dist Numeric matrix, fluorophores in rows and detectors
#'   in columns, as returned by \code{bead.cell.dist()}. Values are used
#'   as \code{abs(mismatch.dist)} so that alignment reflects spectral
#'   location rather than the sign of the mismatch.
#'
#' @return A one-column numeric matrix of cosine similarity values
#'   (column name \code{"Alignment"}), one row per fluorophore common to
#'   both inputs.
#'
#' @export

assess.variability.alignment <- function(
    variability.mad,
    mismatch.dist
) {

  fluorophores <- intersect( rownames( variability.mad ), rownames( mismatch.dist ) )

  if ( length( fluorophores ) == 0 ) {
    stop(
      "No fluorophore names in common between `variability.mad` and ",
      "`mismatch.dist`.",
      call. = FALSE
    )
  }

  var.channels  <- colnames( variability.mad )
  dist.channels <- colnames( mismatch.dist )
  common.channels <- intersect( var.channels, dist.channels )

  if ( length( common.channels ) == 0 ) {
    stop(
      "No shared detector channels between `variability.mad` and ",
      "`mismatch.dist`.",
      call. = FALSE
    )
  }

  if ( !setequal( var.channels, dist.channels ) ) {
    warning(
      sprintf(
        paste(
          "Variability and mismatch detector sets differ (%d shared of %d",
          "variability / %d mismatch channels). Comparing on the shared",
          "channel(s) only. Missing from mismatch: %s. Missing from variability: %s."
        ),
        length( common.channels ), length( var.channels ), length( dist.channels ),
        paste( setdiff( var.channels, dist.channels ), collapse = ", " ),
        paste( setdiff( dist.channels, var.channels ), collapse = ", " )
      ),
      call. = FALSE
    )
  }

  alignment <- lapply(
    fluorophores,
    function( fluor ) {
      tryCatch(
        {
          var.profile  <- variability.mad[ fluor, common.channels ]
          dist.profile <- abs( mismatch.dist[ fluor, common.channels ] )

          sim <- cosine.similarity( rbind( var.profile, dist.profile ) )
          sim[ lower.tri( sim ) ]
        },
        error = function( e ) {
          message(
            sprintf(
              "[ERROR] Fluorophore %s failed: %s",
              fluor,
              conditionMessage( e )
            )
          )
          return( NULL )
        }
      )
    }
  )

  names( alignment ) <- fluorophores

  result <- do.call( rbind, alignment )
  colnames( result ) <- "Alignment"

  result
}
