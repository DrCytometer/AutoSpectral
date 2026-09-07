# unmix_af_basis.r

#' @title Unmix With A Continuous Autofluorescence Basis
#'
#' @description
#' Solves `y = t(spectra) %*% f + basis %*% k` for every cell, fitting
#' autofluorescence as a free combination of basis directions rather than
#' selecting one row from a discrete library. Because the basis is built from
#' the panel-oblique part of the AF library, the normal equations for `k` are
#' diagonal and the per-cell solve is a single projection.
#'
#' Fits are unconstrained, so a cell whose reconstructed AF spectrum goes
#' meaningfully negative is flagged in `negative` and is a candidate for
#' falling back to the discrete library.
#'
#' @param raw.data Expression data from raw FCS files. Cells in rows and
#' detectors in columns.
#' @param spectra Fluorophore spectral signatures, fluorophores in rows and
#' detectors in columns.
#' @param af.basis An AF basis from `get.af.basis`.
#' @param af.spectra Optional AF library used only to report a nearest-library
#' index for each cell, for compatibility with the discrete workflow. Default
#' `NULL`.
#' @param return.fitted.af Logical, default `FALSE`. Whether to return the
#' fitted autofluorescence in detector space.
#' @param negative.tol Numeric, default `0.05`. A cell is flagged when the most
#' negative detector of its reconstructed AF spectrum falls below
#' `-negative.tol` times its largest.
#'
#' @return A list with `fluorophores` (cells x fluorophores), `k` (cells x
#' components), `af` (total fitted AF scale per cell), `negative` (logical
#' vector) and, when requested, `fitted.af` and `af.index`.
#'
#' @export

unmix.af.basis <- function(
    raw.data,
    spectra,
    af.basis,
    af.spectra       = NULL,
    return.fitted.af = FALSE,
    negative.tol     = 0.05
) {

  if ( "AF" %in% rownames( spectra ) )
    spectra <- spectra[ rownames( spectra ) != "AF", , drop = FALSE ]

  unmixing.matrix <- solve.default( tcrossprod( spectra ), spectra )

  B     <- af.basis$basis          # D x q
  sigma <- af.basis$sigma          # q
  W     <- af.basis$directions     # D x q, orthonormal and out-of-span

  # k = (B' Pperp B)^-1 B' Pperp y, and B' Pperp B is diag(sigma^2), so the
  # solve is a projection onto W followed by a scalar division.
  k <- sweep( raw.data %*% W, 2, sigma, "/" )      # cells x q

  fluorophores <- ( raw.data %*% t( unmixing.matrix ) ) -
    k %*% t( unmixing.matrix %*% B )
  colnames( fluorophores ) <- rownames( spectra )

  fitted.af <- k %*% t( B )                        # cells x detectors

  af.max <- apply( fitted.af, 1, max )
  af.min <- apply( fitted.af, 1, min )
  negative <- af.min < -negative.tol * pmax( af.max, 1e-12 )

  out <- list(
    fluorophores = fluorophores,
    k            = k,
    af           = af.max,
    negative     = negative
  )

  if ( return.fitted.af ) out$fitted.af <- fitted.af

  if ( !is.null( af.spectra ) ) {
    # descriptive nearest-library row by cosine similarity, so downstream code
    # expecting an AF index keeps working
    af.norm  <- af.spectra / sqrt( rowSums( af.spectra^2 ) )
    fit.norm <- fitted.af / pmax( sqrt( rowSums( fitted.af^2 ) ), 1e-12 )
    out$af.index <- max.col( fit.norm %*% t( af.norm ) )
  }

  out
}
