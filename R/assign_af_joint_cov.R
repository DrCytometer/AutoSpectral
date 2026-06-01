# assign_af_joint_cov.R

#' @title Assign AF Spectrum By Joint Covariance-Weighted Error
#'
#' @description
#' Assigns each cell to the best-fitting autofluorescence spectral variant using
#' a joint scoring criterion that multiplies two proportional error terms: a
#' covariance-weighted fluorophore error and a raw-space residual error. The
#' covariance of the AF spectra library is propagated into fluorophore space via
#' the unmixing matrix to derive per-channel error weights, giving channels where
#' AF variation matters most a proportionally larger influence on the assignment
#' decision. Multiplying the two terms rewards variants that achieve large
#' improvements on either axis, without requiring an explicit mixing parameter.
#'
#' More principled than \code{assign.af.fluorophores} (plain L1) or
#' \code{assign.af.residuals} (simple residual dot-product) when the AF library
#' contains spectrally diverse variants, because the covariance weights naturally
#' downweight channels where all AF variants look similar and upweight channels
#' where they diverge.
#'
#' @param raw.data Expression data from raw FCS files. Cells in rows and
#' detectors in columns. Columns should be fluorescent data only and must
#' match the columns in \code{spectra}.
#' @param spectra Spectral signatures of fluorophores, normalized between 0
#' and 1, with fluorophores in rows and detectors in columns.
#' @param af.spectra Spectral signatures of autofluorescences, normalized
#' between 0 and 1, with AF variants in rows and detectors in columns. Prepare
#' using \code{get.af.spectra}.
#'
#' @return Integer vector of length \code{nrow(raw.data)} giving the row index
#' (into \code{af.spectra}) of the best-fitting AF variant for each cell.
#'
#' @export

assign.af.joint.cov <- function(
    raw.data,
    spectra,
    af.spectra
) {

  # drop AF row from spectra if present (mirrors assign.af.fluorophores pattern)
  if ( "AF" %in% rownames( spectra ) )
    spectra <- spectra[ rownames( spectra ) != "AF", , drop = FALSE ]

  af.n   <- nrow( af.spectra )
  cell.n <- nrow( raw.data )
  S      <- t( spectra )   # detectors × fluorophores

  # ---- Unmixing matrix (OLS pseudo-inverse) ----
  XtX             <- tcrossprod( spectra )
  unmixing.matrix <- solve.default( XtX, spectra )

  # ---- Covariance-based fluorophore error weights ----
  # Propagate AF spectral covariance into fluorophore space:
  #   Σ_f = U %*% Σ_AF %*% t(U)   [fluorophore.n × fluorophore.n]
  # The diagonal gives per-channel variance; take sqrt for SD-scale weights.
  af.cov           <- stats::cov( af.spectra )
  fluor.cov        <- unmixing.matrix %*% af.cov %*% t( unmixing.matrix )
  af.error.weights <- as.vector( sqrt( abs( diag( fluor.cov ) ) ) )

  # ---- AF projection library ----
  # v.library: how each AF variant projects into fluorophore space [fluorophore.n × af.n]
  v.library <- unmixing.matrix %*% t( af.spectra )

  # r.library: residual AF signal (detector space, orthogonal to fluorophores)
  #            [detector.n × af.n]
  r.library <- t( af.spectra ) - ( S %*% v.library )

  # ---- Estimated AF intensity per cell per variant (k matrix) ----
  numerator                <- raw.data %*% r.library
  denominator              <- colSums( r.library^2 )
  denominator[ denominator == 0 ] <- 1e-10
  k.matrix                 <- sweep( numerator, 2, denominator, "/" )
  k.matrix[ k.matrix < 0 ] <- 0

  # ---- Initial unmix (no AF) and baseline errors ----
  unmixed        <- raw.data %*% t( unmixing.matrix )
  unmixed.nonneg <- unmixed
  unmixed.nonneg[ unmixed.nonneg < 0 ] <- 0
  resids.initial <- raw.data - ( unmixed.nonneg %*% spectra )

  # baseline covariance-weighted L1 fluorophore error (+ small constant for stability)
  base.e.fluor <- rowSums( sweep( abs( unmixed ), 2, af.error.weights, "*" ) ) + 1e-6

  # baseline L2 residual error
  base.e.resid <- sqrt( rowSums( resids.initial^2 ) ) + 1e-6

  # ---- Per-variant multiplicative proportional scoring ----
  # score = p_fluor * p_resid  (both proportional reductions vs baseline)
  # Lower score = better fit. Multiplication rewards large improvements on
  # either axis without requiring an explicit mixing parameter.
  prop.score.matrix <- matrix( 0, nrow = cell.n, ncol = af.n )

  for ( j in seq_len( af.n ) ) {
    k_j <- k.matrix[ , j ]
    v_j <- v.library[ , j ]

    # adjusted fluorophores after removing AF variant j
    adj.fluors <- unmixed - tcrossprod( k_j, v_j )

    # covariance-weighted L1 fluorophore error for variant j
    e.fluor.j <- rowSums( sweep( abs( adj.fluors ), 2, af.error.weights, "*" ) )

    # adjusted raw-space residual after removing AF variant j
    r_j        <- r.library[ , j ]
    adj.resids <- resids.initial - tcrossprod( k_j, r_j )
    e.resid.j  <- sqrt( rowSums( adj.resids^2 ) )

    # proportional scores relative to baseline; product is the joint score
    p.fluor <- e.fluor.j / base.e.fluor
    p.resid <- e.resid.j / base.e.resid

    prop.score.matrix[ , j ] <- p.fluor * p.resid
  }

  # variant that minimises the joint proportional score
  return( max.col( -prop.score.matrix, ties.method = "first" ) )
}
