# calculate_optimize_necessity.R

#' @title Calculate Optimization Necessity Scores for Spectral Variants
#'
#' @description
#' Determines, for each fluorophore with computed spectral variants, whether
#' per-cell spectral optimisation is likely to produce a meaningful improvement
#' in unmixing quality. Fluorophores whose spectral variation does not overlap
#' with the detector channels occupied by other panel members cannot improve
#' unmixing of those channels regardless of how many variants are tested, so
#' optimising them wastes computation without benefit.
#'
#' The score is computed in the same geometric space used by the C++ joint
#' pipeline (`unmix_autospectral_joint_cpp`):
#'
#' 1. For each fluorophore `fl`, the pseudoinverse `U_nof` of the remaining
#'    fluorophores' spectra is computed. This is the same matrix used inside
#'    `FluorPrecomp` to build `pc.w_leakage`.
#' 2. The empirical covariance of `fl`'s delta spectra (variant spread) is
#'    propagated through `U_nof`: `leakage_cov = U_nof %*% delta_cov %*% t(U_nof)`.
#'    The diagonal of `leakage_cov` gives the variance of the unmixing error
#'    induced in each other fluorophore channel by `fl`'s spectral uncertainty.
#' 3. The score for `fl` is the sum of the standard deviations of that leakage:
#'    `score(fl) = sum( sqrt( abs( diag( leakage_cov ) ) ) )`.
#'    Scores are normalised to `[0, 1]` relative to the highest-scoring
#'    fluorophore in the panel.
#'
#' Optionally, if a representative stained sample is provided (as a matrix of
#' per-fluorophore MFIs, typically the median positive signal), scores are
#' additionally weighted by `mu[fl]`, the fluorophore's own brightness. This
#' down-weights fluorophores that are geometrically capable of cross-channel
#' leakage but are always dim in the actual experiment being unmixed. The MFI
#' weighting is applied *after* normalisation so that the geometric score
#' remains interpretable on its own.
#'
#' @param spectra Numeric matrix of spectral signatures (fluorophores x
#'   detectors, values normalised 0-1). Must not include an `"AF"` row.
#' @param delta.list Named list of delta matrices, one per fluorophore. Each
#'   matrix is (n_variants x n_detectors) and contains the element-wise
#'   difference between each variant spectrum and the base spectrum for that
#'   fluorophore. Produced by `get.spectral.variants()` as
#'   `spectra.variants$delta.list`.
#' @param mu Named numeric vector of per-fluorophore MFI values from a
#'   representative stained sample (one value per fluorophore in `spectra`).
#'   Values should be positive; negative or zero values are clamped to zero.
#'   Supply `NULL` (the default) to use the geometric score alone without
#'   brightness weighting.
#' @param threshold Numeric scalar in `[0, 1]`, default `0.01`. Normalised score
#'   below which optimisation is considered unnecessary and the fluorophore is
#'   excluded from the `optimize.recommended` output. A value of `0.01` means
#'   fluorophores scoring below 1% of the top-scoring fluorophore are skipped.
#' @param ridge Numeric scalar, default `1e-4`. Ridge regularisation added to
#'   the diagonal of `delta_cov` before propagation, matching the treatment
#'   inside the C++ precomputation. Prevents degenerate covariance estimates
#'   when a fluorophore has very few variants.
#' @param verbose Logical, default `TRUE`. When `TRUE`, prints a table of
#'   normalised scores and the resulting recommendation for each fluorophore.
#'
#' @return A named list with three elements:
#' \describe{
#'   \item{`scores.raw`}{Named numeric vector of raw (unnormalised) leakage
#'     propagation scores, one per fluorophore in `delta.list`.}
#'   \item{`scores.norm`}{Named numeric vector of scores normalised to `[0, 1]`.}
#'   \item{`optimize.recommended`}{Named logical vector: `TRUE` if the
#'     fluorophore's normalised score meets or exceeds `threshold`, `FALSE`
#'     otherwise. Fluorophores absent from `delta.list` are not included.}
#' }
#'
#' @seealso [get.spectral.variants()], [sanitize.optimization.inputs()]
#'
#' @export

calculate.optimize.necessity <- function(
    spectra,
    delta.list,
    mu        = NULL,
    threshold = 0.01,
    ridge     = 1e-4,
    verbose   = TRUE
) {

  # -------------------------------------------------------------------------
  # Input validation
  # -------------------------------------------------------------------------
  if ( !is.matrix( spectra ) )
    stop( "`spectra` must be a numeric matrix (fluorophores x detectors).",
          call. = FALSE )

  if ( "AF" %in% rownames( spectra ) )
    stop( "`spectra` must not contain an 'AF' row. Remove it before calling.",
          call. = FALSE )

  if ( !is.list( delta.list ) || is.null( names( delta.list ) ) )
    stop( "`delta.list` must be a named list of delta matrices.",
          call. = FALSE )

  if ( !is.numeric( threshold ) || length( threshold ) != 1 ||
       threshold < 0 || threshold > 1 )
    stop( "`threshold` must be a single numeric value in [0, 1].",
          call. = FALSE )

  if ( !is.null( mu ) ) {
    if ( !is.numeric( mu ) || is.null( names( mu ) ) )
      stop( "`mu` must be a named numeric vector of per-fluorophore MFI values.",
            call. = FALSE )
    mu <- pmax( mu, 0 )
  }

  fluorophores <- rownames( spectra )
  D            <- ncol( spectra )

  # only score fluorophores that are both in spectra and have a delta matrix
  score.fluors <- intersect( names( delta.list ), fluorophores )

  if ( length( score.fluors ) == 0 ) {
    warning(
      "No fluorophores found in both `spectra` and `delta.list`. ",
      "Returning empty scores.",
      call. = FALSE
    )
    return( list(
      scores.raw           = numeric(0),
      scores.norm          = numeric(0),
      optimize.recommended = logical(0)
    ) )
  }

  # -------------------------------------------------------------------------
  # Core scoring loop
  # -------------------------------------------------------------------------
  scores.raw <- stats::setNames( rep( 0, length( score.fluors ) ), score.fluors )

  for ( fl in score.fluors ) {

    delta <- delta.list[[ fl ]]   # n_variants x D

    # guard: skip if delta is empty, NULL or all-zero
    if ( is.null( delta ) || !is.matrix( delta ) ||
         nrow( delta ) == 0 || all( abs( delta ) < 1e-12 ) ) {
      next
    }

    other.fl <- fluorophores[ fluorophores != fl ]
    if ( length( other.fl ) == 0 ) next

    # pseudoinverse of the other fluorophores' spectra, U_nof = (S_nof S_nof^T)^-1 S_nof
    # shape: (F-1) x D  — matches pc.U_nof in C++ FluorPrecomp
    S_nof <- spectra[ other.fl, , drop = FALSE ]
    U_nof <- tryCatch(
      solve( S_nof %*% t( S_nof ), S_nof ),
      error = function( e ) {
        # fallback: Moore-Penrose via SVD if S_nof S_nof^T is singular
        sv <- svd( S_nof )
        tol <- max( dim( S_nof ) ) * .Machine$double.eps * sv$d[ 1 ]
        d.inv <- ifelse( sv$d > tol, 1 / sv$d, 0 )
        t( sv$v %*% diag( d.inv, length( d.inv ) ) %*% t( sv$u ) )
      }
    )

    # empirical covariance of delta spectra with ridge regularisation
    # matches `delta_cov.diag() += 1e-4` in C++ Section 2
    if ( nrow( delta ) > 1 ) {
      delta.cov <- stats::cov( delta )
    } else {
      # single variant: use outer product as a rank-1 cov approximation
      delta.cov <- t( delta ) %*% delta / max( nrow( delta ) - 1, 1 )
    }
    diag( delta.cov ) <- diag( delta.cov ) + ridge

    # propagate covariance through U_nof:
    # leakage_cov = U_nof %*% delta_cov %*% t(U_nof)  [(F-1) x (F-1)]
    # matches `leakage_cov = mat(U_nof * delta_cov * U_nof.t())` in C++
    leakage.cov <- U_nof %*% delta.cov %*% t( U_nof )

    # score = sum of leakage standard deviations across other-fluor channels
    # = sum( sqrt( abs( diag( leakage_cov ) ) ) )
    scores.raw[ fl ] <- sum( sqrt( abs( diag( leakage.cov ) ) ) )
  }

  # -------------------------------------------------------------------------
  # Optional MFI brightness weighting
  # -------------------------------------------------------------------------
  # Applied before normalisation so the ranking reflects the combination of
  # geometric leakage potential AND actual signal in this experiment.
  # Only applied to fluorophores where mu is available.
  if ( !is.null( mu ) ) {
    matched <- intersect( score.fluors, names( mu ) )
    if ( length( matched ) < length( score.fluors ) ) {
      missing.fl <- setdiff( score.fluors, matched )
      warning(
        paste(
          "The following fluorophores are missing from `mu` and will use",
          "geometric scores only:", paste( missing.fl, collapse = ", " )
        ),
        call. = FALSE
      )
    }
    for ( fl in matched ) {
      scores.raw[ fl ] <- scores.raw[ fl ] * mu[ fl ]
    }
  }

  # -------------------------------------------------------------------------
  # Normalise to [0, 1]
  # -------------------------------------------------------------------------
  max.score <- max( scores.raw, na.rm = TRUE )

  if ( max.score > 0 ) {
    scores.norm <- scores.raw / max.score
  } else {
    # all scores are zero — every fluorophore is geometrically independent
    # default to FALSE (no optimisation needed) for all
    scores.norm <- scores.raw
  }

  scores.norm[ is.na( scores.norm ) ] <- 0

  # -------------------------------------------------------------------------
  # Threshold to produce recommendation
  # -------------------------------------------------------------------------
  optimize.recommended <- scores.norm >= threshold

  # -------------------------------------------------------------------------
  # Verbose reporting
  # -------------------------------------------------------------------------
  if ( verbose ) {
    score.df <- data.frame(
      Fluorophore = names( scores.norm ),
      Score       = round( scores.norm, 4 ),
      Optimize    = optimize.recommended,
      stringsAsFactors = FALSE
    )
    score.df <- score.df[ order( score.df$Score, decreasing = TRUE ), ]

    message( paste0( "\033[34m", "Variant optimization necessity scores:", "\033[0m" ) )
    message( paste0(
      sprintf( "  %-20s  %6s  %s", "Fluorophore", "Score", "Optimize" ),
      collapse = ""
    ) )
    for ( k in seq_len( nrow( score.df ) ) ) {
      flag <- if ( score.df$Optimize[ k ] ) "\033[32mYES\033[0m" else "\033[31mNO\033[0m"
      message( sprintf(
        "  %-20s  %6.4f  %s",
        score.df$Fluorophore[ k ],
        score.df$Score[ k ],
        flag
      ) )
    }
  }

  return( list(
    scores.raw           = scores.raw,
    scores.norm          = scores.norm,
    optimize.recommended = optimize.recommended
  ) )
}
