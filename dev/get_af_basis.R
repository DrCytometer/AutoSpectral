# get_af_basis.r

#' @title Build A Continuous Autofluorescence Basis
#'
#' @description
#' Constructs a low-dimensional basis for autofluorescence that replaces the
#' discrete AF library in the unmixing solve. The basis is derived from the
#' part of the AF library that the fluorophore panel cannot explain, because
#' only that component is identifiable and only that component reaches the
#' fluorophore abundances. Ranking by out-of-span energy rather than total
#' energy means the leading directions are the ones that actually matter, and
#' it makes the resulting normal equations diagonal, so the per-cell solve is a
#' single projection with no matrix inverse.
#'
#' @param af.spectra Autofluorescence spectra, variants in rows and detectors
#' in columns. Prepare using `get.af.spectra`.
#' @param spectra Fluorophore spectral signatures, fluorophores in rows and
#' detectors in columns.
#' @param n.components Integer, number of basis directions to retain. Default
#' `NULL`, in which case the count is chosen from `var.explained`.
#' @param var.explained Numeric in (0, 1], default `0.99`. Fraction of the
#' out-of-span variance of the AF library the retained basis must capture.
#' Ignored when `n.components` is supplied.
#'
#' @return A list with `basis` (detectors x components), `sigma` (singular
#' values of the projected library), `directions` (the orthonormal out-of-span
#' directions, detectors x components), `var.explained` (cumulative fraction
#' captured) and `n.components`.
#'
#' @export

get.af.basis.test <- function(
    af.spectra,
    spectra,
    n.components  = NULL,
    var.explained = 0.99
) {

  af.spectra <- as.matrix( af.spectra )
  spectra    <- as.matrix( spectra )

  if ( "AF" %in% rownames( spectra ) )
    spectra <- spectra[ rownames( spectra ) != "AF", , drop = FALSE ]

  unmixing.matrix <- solve.default( tcrossprod( spectra ), spectra )

  # part of each AF spectrum no combination of fluorophores can explain
  r.library <- t( af.spectra ) -
    t( spectra ) %*% ( unmixing.matrix %*% t( af.spectra ) )

  sv    <- svd( r.library )
  power <- sv$d^2
  cumul <- cumsum( power ) / sum( power )

  q <- if ( !is.null( n.components ) )
    min( as.integer( n.components ), length( sv$d ) ) else
      which( cumul >= var.explained )[ 1 ]
  if ( is.na( q ) ) q <- length( sv$d )

  keep  <- seq_len( q )
  basis <- t( af.spectra ) %*% sv$v[ , keep, drop = FALSE ]   # D x q
  rownames( basis ) <- colnames( af.spectra )
  colnames( basis ) <- paste0( "AF", keep )

  list(
    basis         = basis,
    sigma         = sv$d[ keep ],
    directions    = sv$u[ , keep, drop = FALSE ],
    var.explained = cumul[ q ],
    n.components  = q
  )
}
