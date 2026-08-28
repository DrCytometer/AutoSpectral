#' @title Assess Spectral Angle Between Reference and Test Spectra
#'
#' @description
#' As \code{assess.mismatch()}, but reports the spectral angle (in degrees)
#' between each fluorophore's reference and test spectra rather than their
#' cosine similarity. The angle is \code{acos()} of the same cosine value
#' \code{assess.mismatch()} returns, so the two metrics are a strictly
#' monotonic transform of one another -- but unlike cosine similarity,
#' spectral angle increases with divergence, so it can be plotted directly
#' against other divergence-increasing metrics (e.g. mismatch distance)
#' without an inverted axis.
#'
#' @param reference.variants Named list of spectral variant matrices (one
#'   per fluorophore), as returned in the \code{variants} element of
#'   \code{get.spectral.variants()}. Row 1 of each matrix is treated as the
#'   reference spectrum.
#' @param test.variants Named list of spectral variant matrices in the same
#'   format as \code{reference.variants}, to be compared against it.
#'
#' @return A one-column numeric matrix of spectral angle values in degrees
#'   (column name \code{"Angle"}), one row per fluorophore common to both
#'   lists.
#'
#' @export

assess.mismatch.angle <- function(
    reference.variants,
    test.variants
) {

  cosine.mat <- assess.mismatch( reference.variants, test.variants )

  # clamp for floating-point drift just outside [-1, 1] before acos()
  angle.mat <- acos( pmin( pmax( cosine.mat, -1 ), 1 ) ) * 180 / pi
  colnames( angle.mat ) <- "Angle"

  angle.mat
}
