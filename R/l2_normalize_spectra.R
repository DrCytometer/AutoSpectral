# l2_normalize_spectra.R

#' @title Renormalize Spectra to Unit L2 (Euclidean) Norm
#'
#' @description
#' AutoSpectral's stored reference spectra and spectral variants are
#' L-infinity normalized (each spectrum's peak detector fixed to exactly
#' 1.0), the convention used throughout unmixing. That convention
#' implicitly asserts zero measurement variance at whichever detector
#' happens to be the peak, pushing all of a spectrum's real variability
#' into its off-peak channels and biasing any distance- or variability-
#' based comparison built on top of it. Because renormalization is just a
#' positive per-row rescaling, this function can be applied directly to
#' already L-infinity-normalized spectra (or the raw unnormalized
#' measurements) to obtain the L2-normalized equivalent, without needing
#' access to the original unnormalized data.
#'
#' @param spectra Numeric matrix with spectra in rows (fluorophores or
#'   spectral variants) and detectors in columns.
#'
#' @return A numeric matrix of the same dimensions and dimnames as
#'   `spectra`, with every row divided by its own Euclidean norm. Rows
#'   with a (near-)zero norm are returned unchanged rather than divided
#'   by ~0. Attributes other than `dim`/`dimnames` on the input (e.g.
#'   `noise.floor`, `spillover.spread`) are not preserved.
#'
#' @export

l2.normalize.spectra <- function( spectra ) {
  
  mat <- as.matrix( spectra )
  row.norm  <- sqrt( rowSums( mat^2 ) )
  safe.norm <- ifelse( row.norm > 1e-9, row.norm, 1 )
  
  out <- mat / safe.norm
  dimnames( out ) <- dimnames( mat )
  out
}