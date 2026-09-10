# assign_af_joint_cov_l2.R

#' @title Assign AF Spectrum By Joint Covariance-Weighted Squared Error
#'
#' @description
#' Assigns each cell to the best-fitting autofluorescence spectral variant using
#' a joint scoring criterion that multiplies two proportional error terms: a
#' covariance-weighted fluorophore error and a raw-space residual error, both
#' measured as squared (L2) deviations. The covariance of the AF spectra
#' library is propagated into fluorophore space via the unmixing matrix to
#' derive per-channel error weights, giving channels where AF variation
#' matters most a proportionally larger influence on the assignment decision.
#' Multiplying the two terms rewards variants that achieve large improvements
#' on either axis, without requiring an explicit mixing parameter.
#'
#' Because both error terms are quadratic in the per-cell, per-variant AF
#' abundance, each can be expanded into a baseline term, a cross term, and a
#' curvature term. All three are computed as single matrix products across
#' every cell and every variant simultaneously, so no per-variant loop is
#' required. This makes the function substantially faster than the L1
#' (\code{abs}-based) formulation in \code{assign.af.joint.cov}, at the cost
#' of being somewhat less robust to outlier channels, since squared error
#' weights large deviations more heavily than L1.
#'
#' @param raw.data Expression data from raw FCS files. Cells in rows and
#' detectors in columns. Columns should be fluorescent data only and must
#' match the columns in \code{spectra}.
#' @param spectra Spectral signatures of fluorophores, normalized between 0
#' and 1, with fluorophores in rows and detectors in columns.
#' @param af.spectra Spectral signatures of autofluorescences, normalized
#' between 0 and 1, with AF variants in rows and detectors in columns. Prepare
#' using \code{get.af.spectra}.
#' @param return.scores Logical, default \code{FALSE}. If `\code{TRUE}`, also
#' returns the unmixed data and scores for each AF variant per cell.
#'
#' @return Integer vector of length \code{nrow(raw.data)} giving the row index
#' (into \code{af.spectra}) of the best-fitting AF variant for each cell.
#'
#' @export

assign.af.joint.cov.l2 <- function(
    raw.data,
    spectra,
    af.spectra,
    return.scores = FALSE
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
  numerator      <- raw.data %*% r.library
  raw.denominator <- colSums( r.library^2 )
  
  # Identifiability guard: an AF variant lying almost inside the fluorophore
  # span has a vanishing out-of-span residual direction, so its abundance k
  # is not identifiable from the residual and the raw ratio explodes. Floor
  # each self-dot at a fraction of the largest, capping the relative
  # amplification of near-in-span variants. This floored value is used only
  # for computing k; the error expansion below uses the true, unfloored
  # ‖r_j‖² so the quadratic matches the actual squared residual.
  denominator <- pmax( raw.denominator, 0.01 * max( raw.denominator, 1e-10 ) )
  k.matrix                 <- sweep( numerator, 2, denominator, "/" )
  k.matrix[ k.matrix < 0 ] <- 0
  
  # ---- Initial unmix (no AF) and baseline errors ----
  unmixed        <- raw.data %*% t( unmixing.matrix )
  unmixed.nonneg <- unmixed
  unmixed.nonneg[ unmixed.nonneg < 0 ] <- 0
  resids.initial <- raw.data - ( unmixed.nonneg %*% spectra )
  
  # baseline covariance-weighted squared (L2) fluorophore error
  base.e.fluor <- rowSums( sweep( unmixed^2, 2, af.error.weights, "*" ) ) + 1e-6
  
  # baseline squared (L2) residual error
  base.e.resid <- rowSums( resids.initial^2 ) + 1e-6
  
  # ---- Vectorised per-variant squared-error expansion ----
  # For variant j, the adjusted fluorophores are (unmixed - k_j * v_j).
  # Expanding the weighted squared norm gives a quadratic in k_j:
  #   e.fluor_j = base.e.fluor - 2 k_j (Uw %*% v_j) + k_j^2 c_j
  # where Uw is `unmixed` weighted column-wise by af.error.weights and
  # c_j = sum_f w_f v_jf^2. All variants are handled in one pair of matrix
  # products rather than a loop. The residual term follows the same pattern
  # with r.library replacing v.library and an unweighted dot product.
  Uw          <- sweep( unmixed, 2, af.error.weights, "*" )
  cross.fluor <- Uw %*% v.library                                          # cell.n × af.n
  c.fluor     <- colSums( sweep( v.library^2, 1, af.error.weights, "*" ) ) # af.n
  
  cross.resid <- resids.initial %*% r.library   # cell.n × af.n
  c.resid     <- raw.denominator                # af.n, true ‖r_j‖²
  
  e.fluor <- base.e.fluor - 2 * k.matrix * cross.fluor + sweep( k.matrix^2, 2, c.fluor, "*" )
  e.resid <- base.e.resid - 2 * k.matrix * cross.resid + sweep( k.matrix^2, 2, c.resid, "*" )
  
  # numerical floor: expansion is exact in theory but floating point can push
  # near-zero values slightly negative
  e.fluor[ e.fluor < 0 ] <- 0
  e.resid[ e.resid < 0 ] <- 0
  
  # proportional scores relative to baseline; product is the joint score
  p.fluor <- e.fluor / base.e.fluor
  p.resid <- e.resid / base.e.resid
  
  prop.score.matrix <- p.fluor * p.resid
  
  # variant that minimises the joint proportional score
  af.index <- max.col( -prop.score.matrix, ties.method = "first" )
  
  if ( !return.scores ) return( af.index )
  
  # k.matrix and v.library let a caller reconstruct the implied AF abundance
  # and fluorophore correction for any candidate, not only the winner
  list(
    af.index  = af.index,
    scores    = prop.score.matrix,
    k.matrix  = k.matrix,
    v.library = v.library,
    unmixed   = unmixed
  )
}