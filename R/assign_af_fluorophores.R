# assign_af_fluorophores.r

#' @title Assign AF Spectrum By Fluorophore Projection
#'
#' @description
#' Projects the autofluorescence spectral variants into fluorophore unmixed
#' space to determine which best fits each cell (event). A fast approximation for
#' brute force sequential unmixing method in early versions of AutoSpectral.
#' Provides essentially identical results to minimization of fluorophore signal
#' (dist0 method). Substantially faster. Use L1 (absolute value) minimization,
#' which works better than the standard L2 (squared error).
#'
#' @param raw.data Expression data from raw FCS files. Cells in rows and
#' detectors in columns. Columns should be fluorescent data only and must
#' match the columns in spectra.
#' @param spectra Spectral signatures of fluorophores, normalized between 0
#' and 1, with fluorophores in rows and detectors in columns.
#' @param af.spectra Spectral signatures of autofluorescences, normalized
#' between 0 and 1, with fluorophores in rows and detectors in columns.
#'
#' @return Row indices for best-fitting AF spectra (from `af.spectra`)
#'
#' @export

assign.af.fluorophores <- function(
    raw.data,
    spectra,
    af.spectra
) {

  # how many AF spectra do we have?
  af.n <- nrow( af.spectra )

  # calculate pseudoinverse
  S <- t( spectra )
  XtX <- tcrossprod( spectra )
  unmixing.matrix <- solve.default( XtX, spectra )

  # how much each AF variant looks like each fluorophore
  v.library <- unmixing.matrix %*% t( af.spectra )

  # calculate the 'residual AF' (the part of AF fluorophores can't explain)
  r.library <- t( af.spectra ) - ( S %*% v.library )

  # predicted AF intensity
  numerator <- raw.data %*% r.library

  # denominator (vector of length af.n)
  denominator <- colSums( r.library^2 )

  # k_matrix (cell.n x af.n): estimated AF intensity per cell/variant
  k.matrix <- sweep( numerator, 2, denominator, "/" )

  # initial unmix (no AF)
  unmixed <- raw.data %*% t( unmixing.matrix )

  # Initialize a matrix to store the 'error' (abs sum of fluors) for each AF choice
  error.matrix <- matrix(
    0,
    nrow = nrow( raw.data ),
    ncol = af.n
  )

  for ( i in seq_len( af.n ) ) {
    # predicted intensity of AF variant i
    k_i <- k.matrix[ , i, drop = FALSE ]

    # dye-leakage signature of AF variant i
    v_i <- t( v.library[, i, drop = FALSE ] )

    # adjusted fluorophore values for all cells using variant i
    adjusted.fluors <- unmixed - ( k_i %*% v_i )

    # worst-case' scenario: sum of absolute fluorophore signals
    error.matrix[ , i ] <- rowSums( abs( adjusted.fluors ) )
  }

  # maximum score corresponds to the minimum L1 error (max of negative is min)
  return( max.col( -error.matrix ) )
}
