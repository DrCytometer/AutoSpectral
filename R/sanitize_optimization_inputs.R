# optimize_unmix.r

#' @title Sanitize Optimization Inputs
#'
#' @description
#' Checks the delta list and delta norms for empty, NULL or redundant content,
#' removes fluorophores from the list of ones to be optimized if they fail checks.
#'
#' @param spectra Numeric matrix (fluors x detectors)
#' @param optimize.fluors Character vector of fluorophores present in variants
#' @param variants List of variant matrices per fluorophore
#' @param delta.norms List of delta norms per fluorophore
#'
#' @return Updated list of fluorophores to be optimized
#'
#' @export

sanitize.optimization.inputs <- function(
    spectra,
    optimize.fluors,
    variants,
    delta.norms
) {

  # identify which fluors actually have meaningful variants
  valid.optimization.fluors <- c()

  for ( fl in optimize.fluors ) {
    # does the fluorophore exist in the variants list?
    if ( !fl %in% names( variants ) ) next

    # is the delta.norm valid (not NULL, not zero)?
    if ( is.null( delta.norms[[ fl ]] ) || all( delta.norms[[ fl ]] < 1e-12 ) ) {
      message(
        paste( "Skipping optimization for", fl, "- No spectral variation detected." )
      )
      next
    }

    valid.optimization.fluors <- c( valid.optimization.fluors, fl )
  }

  return( valid.optimization.fluors )
}
