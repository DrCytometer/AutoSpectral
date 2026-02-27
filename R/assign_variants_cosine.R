# assign_variants.r

#' @title Assign Variant Spectrum By Fluorophore Projection
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
#' @param fluorophore Name of the fluorophore to be tested for variation.
#' @param variant.spectra Spectral signatures of fluorophore variants, normalized
#' between 0 and 1, with fluorophores in rows and detectors in columns.
#'
#' @return Row indices for best-fitting variant spectra (from `variant.spectra`)
#'
#' @export

assign.variants.cosine <- function(
    raw.data,
    spectra,
    fluorophore,
    variant.spectra
) {

  # how many variant spectra do we have?
  var.n <- nrow( variant.spectra )

  # drop this fluorophore from the spectral mixing matrix
  spectra <- spectra[ !rownames( spectra ) == fluorophore, , drop = FALSE ]

  # calculate pseudoinverse
  S <- t( spectra )
  XtX <- tcrossprod( spectra )
  unmixing.matrix <- solve.default( XtX, spectra )

  # how much each variant looks like every other fluorophore
  v.library <- unmixing.matrix %*% t( variant.spectra )

  # calculate the 'residual fluor' (the part of this fluor other fluorophores can't explain)
  r.library <- t( variant.spectra ) - ( S %*% v.library )

  # predicted fluor intensity
  numerator <- raw.data %*% r.library

  # denominator (vector of length var.n)
  denominator <- colSums( r.library^2 )

  # k_matrix (cell.n x var.n): estimated fluor intensity per cell/variant
  k.matrix <- sweep( numerator, 2, denominator, "/" )

  # similarity weights
  sim.weights <- cosine.similarity( rbind( variant.spectra[ 1, ], spectra ) )
  sim.weights <- sim.weights[ 1, ]
  sim.weights <- sim.weights[ -1 ]

  # initial unmix (no fluor)
  unmixed <- raw.data %*% t( unmixing.matrix )

  # Initialize a matrix to store the 'error' (abs sum of other fluors) for each variant choice
  error.matrix <- matrix(
    0,
    nrow = nrow( raw.data ),
    ncol = var.n
  )

  for ( i in seq_len( var.n ) ) {
    # predicted intensity of fluor variant i
    k_i <- k.matrix[ , i, drop = FALSE ]

    # dye-leakage signature of fluor variant i
    v_i <- t( v.library[, i, drop = FALSE ] )

    # weight by similarity so we downplay minor fluctuations
    v_i <- v_i * sim.weights

    # adjusted fluorophore values for all cells using variant i
    adjusted.fluors <- unmixed - ( k_i %*% v_i )

    # worst-case' scenario: sum of absolute fluorophore signals
    error.matrix[ , i ] <- rowSums( abs( adjusted.fluors ) )
  }

  # maximum score corresponds to the minimum L1 error (max of negative is min)
  return( max.col( -error.matrix ) )
}
