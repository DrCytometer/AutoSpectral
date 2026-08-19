# correct_spectra_mean_delta.R

#' @title Correct Reference Spectra for Systematic Variant Offset
#'
#' @description
#' Shifts each fluorophore's reference spectrum by the mean delta captured
#' in `variant.basis` (see `build.variant.basis()`), then re-normalises to
#' L-infinity = 1. A non-zero mean delta means the control-derived
#' reference is systematically off-centre relative to the variant
#' population -- a bias in `spectra`, not a source of variance -- so
#' folding it into the covariance would inflate uncertainty while leaving
#' the bias itself uncorrected. This applies the correction once, up
#' front, so `unmix.gls()` sees a centred reference and `variant.basis`
#' describes only the residual spread around it.
#'
#' The `mean.delta` entries in the returned `variant.basis` are zeroed for
#' every fluorophore corrected, since the mean is now baked into `spectra`;
#' calling this function again on its own output is therefore a no-op.
#'
#' @param spectra Numeric matrix (fluorophores x detectors).
#' @param variant.basis List from `build.variant.basis()`.
#' @param fluorophores Optional character vector restricting the correction
#'   to named rows of `spectra`. Default `NULL` (all rows present in both
#'   `spectra` and `variant.basis`).
#'
#' @return A list with `spectra` (corrected, L-infinity re-normalised) and
#'   `variant.basis` (mean.delta zeroed for corrected fluorophores).
#'
#' @export

correct.spectra.mean.delta <- function( spectra, variant.basis, fluorophores = NULL ) {
  
  target <- if ( is.null( fluorophores ) ) rownames( spectra ) else fluorophores
  target <- intersect( target, rownames( spectra ) )
  target <- intersect( target, names( variant.basis ) )
  
  for ( fl in target ) {
    
    md <- variant.basis[[ fl ]]$mean.delta
    if ( is.null( md ) || !any( md != 0 ) ) next
    
    shifted <- spectra[ fl, ] + md[ colnames( spectra ) ]
    peak    <- max( abs( shifted ) )
    if ( peak > 0 ) shifted <- shifted / peak
    
    spectra[ fl, ] <- shifted
    variant.basis[[ fl ]]$mean.delta[] <- 0
  }
  
  list( spectra = spectra, variant.basis = variant.basis )
}