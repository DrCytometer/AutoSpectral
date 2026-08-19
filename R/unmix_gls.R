# unmix_gls.R

#' @title Build Variant Covariance Basis
#'
#' @description
#' Converts the per-fluorophore variant deltas from `get.spectral.variants()`
#' into a low-rank covariance basis for `unmix.gls()`. Each fluorophore's
#' contribution to the detector-space covariance is
#' \eqn{x_f^2 \sum_j \lambda_{fj} b_{fj} b_{fj}^\top}.
#'
#' The *mean* delta is returned separately rather than folded into the
#' covariance: a systematic offset from the reference spectrum is a
#' correction to `spectra`, not a source of variance, and treating it as
#' variance both inflates the noise model and leaves the bias uncorrected.
#'
#' @param spectra.variants Either the full list returned by
#'   `get.spectral.variants()` (containing `$delta.list`), or a named list of
#'   delta matrices (variants x detectors).
#' @param rank Integer, maximum number of components retained per
#'   fluorophore. Default `2`.
#' @param var.explained Numeric in (0, 1]. Components are dropped once this
#'   fraction of the delta variance is captured. Default `0.99`.
#' @param min.lambda Numeric, eigenvalues below this fraction of the leading
#'   eigenvalue are dropped. Default `1e-4`.
#'
#' @return A named list, one entry per fluorophore with at least two
#'   variants, each containing `basis` (components x detectors), `lambda`
#'   (numeric, one per component), and `mean.delta` (numeric, length
#'   detectors). Fluorophores with insufficient variants are omitted.
#'
#' @export

build.variant.basis <- function(
    spectra.variants,
    rank            = 6L,
    var.explained   = 0.99,
    min.lambda      = 1e-4,
    pooled.fallback = TRUE,
    verbose         = TRUE
) {

  delta.list <- if ( !is.null( spectra.variants$delta.list ) )
    spectra.variants$delta.list else spectra.variants

  if ( !is.list( delta.list ) || is.null( names( delta.list ) ) )
    stop( "`spectra.variants` must be a named list of delta matrices or the ",
          "output of `get.spectral.variants()`.", call. = FALSE )

  fit.one <- function( delta ) {

    delta <- as.matrix( delta )
    if ( nrow( delta ) < 2L ) return( NULL )

    mean.delta <- colMeans( delta )
    centred    <- sweep( delta, 2L, mean.delta, "-" )

    sv  <- svd( centred )
    lam <- sv$d^2 / ( nrow( centred ) - 1L )

    if ( lam[ 1L ] <= 0 ) return( NULL )

    keep <- which( lam >= min.lambda * lam[ 1L ] )
    keep <- keep[ keep <= min( rank, length( lam ) ) ]
    if ( length( keep ) == 0L ) return( NULL )

    cum <- cumsum( lam[ keep ] ) / sum( lam )
    if ( any( cum >= var.explained ) )
      keep <- keep[ seq_len( which( cum >= var.explained )[ 1L ] ) ]

    basis <- t( sv$v[ , keep, drop = FALSE ] )
    colnames( basis ) <- colnames( delta )

    list( basis = basis, lambda = lam[ keep ], mean.delta = mean.delta,
          n.variants = nrow( delta ), source = "fitted" )
  }

  out <- lapply( delta.list, fit.one )
  names( out ) <- names( delta.list )

  missing.fl <- names( out )[ vapply( out, is.null, logical( 1 ) ) ]

  if ( pooled.fallback && length( missing.fl ) > 0L ) {

    fitted.fl <- setdiff( names( out ), missing.fl )

    if ( length( fitted.fl ) == 0L ) {
      warning( "No fluorophore had enough variants to fit a basis; cannot ",
               "build a pooled fallback.", call. = FALSE )
    } else {

      # Pool *shape*, not scale: normalise each fluorophore's centred delta
      # rows to unit length before stacking, so a fluorophore with large raw
      # deltas does not dominate the pooled directions. The result describes
      # the typical cross-detector pattern of variant noise, rescaled onto
      # the target fluorophore's own variance level below.
      pooled.rows <- do.call( rbind, lapply( fitted.fl, function( fl ) {
        delta      <- as.matrix( delta.list[[ fl ]] )
        mean.delta <- colMeans( delta )
        centred    <- sweep( delta, 2L, mean.delta, "-" )
        row.norm   <- sqrt( rowSums( centred^2 ) )
        centred[ row.norm > 0, , drop = FALSE ] / row.norm[ row.norm > 0 ]
      } ) )

      pooled.sv    <- svd( pooled.rows )
      pooled.lam   <- pooled.sv$d^2 / ( nrow( pooled.rows ) - 1L )
      pooled.keep  <- seq_len( min( rank, length( pooled.lam ) ) )
      pooled.basis <- t( pooled.sv$v[ , pooled.keep, drop = FALSE ] )
      colnames( pooled.basis ) <- colnames( pooled.rows )

      # Scale: median leading eigenvalue among fitted fluorophores. This
      # assumes a fluorophore with no usable variant data wobbles about as
      # much as a typical panel member -- untested, and deliberately
      # conservative rather than zero.
      lam1.median <- stats::median( vapply(
        out[ fitted.fl ], function( x ) x$lambda[ 1L ], numeric( 1 ) ) )
      pooled.lam.scaled <- pooled.lam[ pooled.keep ] *
        ( lam1.median / pooled.lam[ pooled.keep[ 1L ] ] )

      for ( fl in missing.fl ) {
        out[[ fl ]] <- list(
          basis      = pooled.basis,
          lambda     = pooled.lam.scaled,
          mean.delta = stats::setNames( rep( 0, ncol( pooled.basis ) ),
                                        colnames( pooled.basis ) ),
          n.variants = if ( is.null( delta.list[[ fl ]] ) ) 0L else
            nrow( as.matrix( delta.list[[ fl ]] ) ),
          source     = "pooled.fallback"
        )
      }

      if ( verbose )
        message( sprintf(
          "%d fluorophore(s) had insufficient variants for their own basis; ",
          length( missing.fl ) ), "assigned a pooled fallback: ",
          paste( missing.fl, collapse = ", " ) )
    }
  }

  out <- out[ !vapply( out, is.null, logical( 1 ) ) ]

  if ( verbose && length( out ) > 0L ) {
    ranks <- vapply( out, function( x ) length( x$lambda ), integer( 1 ) )
    message( sprintf(
      "Variant basis: %d fluorophore(s), rank %d-%d (median %.0f), total k = %d",
      length( ranks ), min( ranks ), max( ranks ), stats::median( ranks ),
      sum( ranks ) ) )
  }

  out
}


# -----------------------------------------------------------------------------
# Internal: apply Sigma^{-1} to Z, where Sigma = diag(d) + U diag(c) U'
# -----------------------------------------------------------------------------

.sigma.solve <- function( Z, d.vec, U, c.vec, ridge = 1e-8, diagnostics = FALSE ) {

  Dinv.Z <- Z / d.vec

  if ( is.null( U ) || length( c.vec ) == 0L ) {
    if ( !diagnostics ) return( Dinv.Z )
    return( list( W = Dinv.Z, logdet = sum( log( d.vec ) ) ) )
  }

  keep <- which( abs( c.vec ) > 0 )
  if ( length( keep ) == 0L ) {
    if ( !diagnostics ) return( Dinv.Z )
    return( list( W = Dinv.Z, logdet = sum( log( d.vec ) ) ) )
  }

  # Absorb |c| into the columns: Sigma = D + V S V', S = diag( sign( c ) ).
  # The inner matrix then contains S rather than diag( 1 / c ). The 1 / c
  # form loses precision whenever the low-rank part dominates the diagonal,
  # which is the normal regime for bright cells once the x^2 variant term
  # is active.
  V <- sweep( U[ , keep, drop = FALSE ], 2L,
              sqrt( abs( c.vec[ keep ] ) ), "*" )
  s.sign <- sign( c.vec[ keep ] )

  DiV <- V / d.vec
  M   <- diag( s.sign, nrow = length( s.sign ) ) + crossprod( V, DiV )
  diag( M ) <- diag( M ) + ridge * max( abs( diag( M ) ), 1 )

  W <- Dinv.Z - DiV %*% solve( M, crossprod( DiV, Z ) )

  if ( !diagnostics ) return( W )

  # Matrix determinant lemma: det(D + V S V') = det(D) det(S) det(M), and
  # det(S) = +-1 since S is diagonal +-1, so it drops out of the log.
  logdet <- sum( log( d.vec ) ) +
    as.numeric( determinant( M, logarithm = TRUE )$modulus )

  list( W = W, logdet = logdet )
}


#' @title Unmix by Generalised Least Squares
#'
#' @description
#' Per-cell GLS unmixing under a composite noise model:
#' \deqn{\Sigma_i = diag(read.var) + diag(\mu_i)/\kappa
#'   + \frac{1}{\bar\kappa}\sum_f x_{if} m_f [diag(p_f) - p_f p_f^\top]
#'   + \sum_f x_{if}^2 \Sigma_f^{var} + \alpha_i^2 \Sigma_{k_i}^{AF}}
#'
#' where \eqn{m_f = \sum_d S_{fd}} is the spectral mass of fluorophore f.
#' The spillover (multinomial), spectral-variant and AF terms are low rank,
#' so \eqn{\Sigma_i^{-1}} is applied via the Woodbury identity without ever
#' forming a D x D inverse. Abundances are clamped at zero **only** when
#' constructing \eqn{\Sigma_i} (a variance cannot be negative); the returned
#' estimate is unconstrained, preserving continuity through zero.
#'
#' The spillover term and shot noise (`diag(mu_i)/kappa`) both derive from
#' photon partitioning; if `kappa` was estimated on cell or bead data by
#' `estimate.noise.model()` (the normal case), the fitted slope has already
#' absorbed the spillover diagonal in full, and `include.spillover = TRUE`
#' double-counts it. Leave `include.spillover = FALSE` unless `kappa` came
#' from a source measured free of spillover (e.g. `BDSPECTRAL QSPE`).
#'
#' @param raw.data Numeric matrix (events x detectors).
#' @param spectra Numeric matrix (fluorophores x detectors), no `"AF"` row.
#' @param noise.model List from `estimate.noise.model()`.
#' @param variant.basis Optional list from `build.variant.basis()`.
#' @param af.spectra Optional AF dictionary. When supplied with `af.index`,
#'   each cell's assigned AF spectrum is appended to `spectra` for that cell.
#' @param af.index Optional integer vector, length `nrow(raw.data)`, giving
#'   the row of `af.spectra` assigned to each cell.
#' @param af.basis Optional list, one entry per AF dictionary node, each
#'   containing `basis` (components x detectors) and `lambda` (numeric).
#'   Supplies the within-node AF covariance term
#'   \eqn{\alpha_i^2 \Sigma_{k_i}^{AF}} (section 2.3 of
#'   `CONTEXT_GLS_unmixing.md`). This does not need to be built separately:
#'   pass `attr(af.spectra, "af.model")$nodes` directly when `af.spectra`
#'   was produced by `get.af.spectra(..., return.model = TRUE)`, which
#'   already computes exactly this structure, keyed by
#'   `rownames(af.spectra)`. Ignored unless `af.spectra` and `af.index` are
#'   also supplied. Default `NULL`, which omits the term entirely -- AF
#'   mismatch then appears as unmodelled residual rather than covariance.
#' @param n.iter Integer, number of GLS refinement iterations. Default `2`.
#'   The covariance depends on the abundances, so one iteration is a plug-in
#'   estimate from OLS and two is usually enough.
#' @param method `"woodbury"` (default) or `"dense"`. `"dense"` forms
#'   \eqn{\Sigma_i} explicitly and is intended as a correctness reference.
#' @param include.spillover Logical, whether to include the multinomial
#'   spillover term. Default `FALSE` (see description).
#' @param spillover.kappa Optional scalar \eqn{\bar\kappa} for the
#'   multinomial term. Defaults to `noise.model$kappa.pooled`. A scalar is
#'   used deliberately: photon partitioning across detectors is a property
#'   of the emission spectrum, not of per-detector gain.
#' @param gain.cv Numeric, coefficient of variation of multiplicative gain
#'   fluctuation (illumination/transit-time variation), applied as a rank-1
#'   term \eqn{\code{gain.cv}^2\, \mu_i \mu_i^\top}. Distinct from shot
#'   noise and spillover, which are photon-counting effects; `gain.cv`
#'   captures gain variation at fixed photon count. If `curvature.coef`
#'   from `estimate.noise.model()` comes back positive and similar in
#'   magnitude across detectors on a bead run, that is this term's expected
#'   signature, and `sqrt(curvature.coef)` at a representative detector is
#'   a reasonable starting estimate. Default `0` (off): unvalidated on real
#'   data as of this writing.
#' @param active.threshold Numeric. Fluorophores below this fraction of the
#'   cell's largest abundance are excluded from the low-rank terms, keeping
#'   the Woodbury rank well below the detector count. Default `1e-3`.
#' @param ridge Numeric ridge added to the `XtX` solve. Default `1e-8`.
#' @param inner.ridge Numeric ridge added inside `.sigma.solve()`'s Woodbury
#'   inner matrix `M`. Kept separate from `ridge` because the two solves see
#'   different conditioning: `M` is `k x k` in the (typically small) active
#'   rank, `XtX` is `F x F` in the number of fluorophores. Default `1e-10`.
#' @param return.variance Logical. When `TRUE`, also returns per-cell
#'   estimated standard errors, i.e. \eqn{\sqrt{diag((S\Sigma^{-1}S^\top)^{-1})}}.
#'   Default `FALSE`.
#' @param return.fit.stats Logical. When `TRUE`, also returns per-cell
#'   goodness-of-fit: `chisq` (\eqn{r^\top\Sigma^{-1}r}, using the plug-in
#'   \eqn{\Sigma} from the final iteration, evaluated at the converged
#'   \eqn{\hat x} -- the same plug-in convention `se` already uses), `df`
#'   (detectors minus active fluorophores, including AF), `logdet`
#'   (\eqn{\log|\Sigma_i|}, via the matrix determinant lemma at no extra
#'   solve), and `p` (upper-tail chi-square p-value). A small `p` flags a
#'   cell whose residual is inconsistent with the fitted noise model -- the
#'   direct test proposed in `CONTEXT_GLS_unmixing.md` section 5 as a
#'   replacement for `calculate.optimize.necessity()`. Default `FALSE`.
#' @param verbose Logical. Default `TRUE`.
#'
#' @return A matrix (cells x fluorophores), or when `return.variance` or
#'   `return.fit.stats` is `TRUE`, a list with `unmixed` and whichever of
#'   `se`, `chisq`, `df`, `logdet`, `p` were requested.
#'
#' @export

unmix.gls <- function(
    raw.data,
    spectra,
    noise.model,
    variant.basis     = NULL,
    af.spectra        = NULL,
    af.index          = NULL,
    af.basis          = NULL,
    n.iter            = 2L,
    method            = c( "woodbury", "dense" ),
    include.spillover = FALSE,
    spillover.kappa   = NULL,
    gain.cv           = 0,
    active.threshold  = 1e-3,
    ridge             = 1e-8,
    inner.ridge       = 1e-10,
    return.variance   = FALSE,
    return.fit.stats  = FALSE,
    verbose           = TRUE
) {

  method    <- match.arg( method )
  raw.data  <- as.matrix( raw.data )
  det.names <- colnames( spectra )

  # spectral reference matrix must contain exactly one row per fluorophore
  # before any unmixing method runs
  check.spectra.duplicates( spectra )

  if ( !all( det.names %in% colnames( raw.data ) ) )
    stop( "`raw.data` does not contain all detectors present in `spectra`.",
          call. = FALSE )
  raw.data <- raw.data[ , det.names, drop = FALSE ]

  if ( "AF" %in% rownames( spectra ) )
    stop( "`spectra` must not contain an 'AF' row; pass AF via `af.spectra`.",
          call. = FALSE )

  read.var <- as.numeric( noise.model$read.var[ det.names ] )
  kappa    <- as.numeric( noise.model$counts.per.unit[ det.names ] )
  kappa.p  <- if ( is.null( spillover.kappa ) )
    noise.model$kappa.pooled else spillover.kappa

  if ( anyNA( read.var ) || anyNA( kappa ) )
    stop( "`noise.model` does not cover all detectors in `spectra`.",
          call. = FALSE )

  cell.n  <- nrow( raw.data )
  fluor.n <- nrow( spectra )
  use.af  <- !is.null( af.spectra ) && !is.null( af.index )

  if ( use.af ) {
    af.spectra <- as.matrix( af.spectra )[ , det.names, drop = FALSE ]
    if ( length( af.index ) != cell.n )
      stop( "`af.index` must have one entry per event.", call. = FALSE )
  }

  # emission probability vectors (rows sum to 1), for the multinomial term
  .prob.rows <- function( m ) {
    m <- pmax( m, 0 )
    m / pmax( rowSums( m ), 1e-12 )
  }
  P.fluor <- .prob.rows( spectra )
  P.af    <- if ( use.af ) .prob.rows( af.spectra ) else NULL

  # ---------------------------------------------------------------------------
  # Starting values: OLS, grouped by AF assignment
  # ---------------------------------------------------------------------------

  group <- if ( use.af ) af.index else rep( 1L, cell.n )
  x.mat <- matrix( 0, nrow = cell.n, ncol = fluor.n + as.integer( use.af ),
                   dimnames = list( NULL, c( rownames( spectra ),
                                             if ( use.af ) "AF" ) ) )

  for ( g in unique( group ) ) {
    idx <- which( group == g )
    S.g <- if ( use.af )
      rbind( spectra, AF = af.spectra[ g, ] ) else spectra
    x.mat[ idx, ] <- unmix.ols.fast( raw.data[ idx, , drop = FALSE ], S.g )
  }

  se.mat <- if ( return.variance )
    matrix( NA_real_, nrow = cell.n, ncol = ncol( x.mat ),
            dimnames = dimnames( x.mat ) ) else NULL

  chisq.vec  <- if ( return.fit.stats ) rep( NA_real_, cell.n ) else NULL
  df.vec     <- if ( return.fit.stats ) rep( NA_integer_, cell.n ) else NULL
  logdet.vec <- if ( return.fit.stats ) rep( NA_real_, cell.n ) else NULL

  # index of variant bases aligned to the rows of spectra
  vb.idx <- if ( is.null( variant.basis ) ) NULL else
    match( rownames( spectra ), names( variant.basis ) )

  if ( !is.null( vb.idx ) && anyNA( vb.idx ) )
    warning( "No variant basis for: ",
             paste( rownames( spectra )[ is.na( vb.idx ) ], collapse = ", " ),
             ". These fluorophores are modelled with zero spectral ",
             "uncertainty and will absorb unexplained residual as a sink. ",
             "Rebuild `variant.basis` with `build.variant.basis(..., ",
             "pooled.fallback = TRUE)` to avoid this.", call. = FALSE )

  # index of the AF covariance basis aligned to rows of af.spectra
  ab.idx <- if ( is.null( af.basis ) || !use.af ) NULL else {
    if ( !is.null( names( af.basis ) ) && !is.null( rownames( af.spectra ) ) )
      match( rownames( af.spectra ), names( af.basis ) )
    else
      seq_len( nrow( af.spectra ) )
  }

  if ( verbose )
    message( sprintf( "GLS unmixing %d events, %d iterations (%s)",
                      cell.n, n.iter, method ) )

  # ---------------------------------------------------------------------------
  # Per-cell refinement
  # ---------------------------------------------------------------------------

  for ( i in seq_len( cell.n ) ) {

    y   <- raw.data[ i, ]
    g   <- group[ i ]
    S.i <- if ( use.af ) rbind( spectra, AF = af.spectra[ g, ] ) else spectra
    P.i <- if ( use.af ) rbind( P.fluor, AF = P.af[ g, ] ) else P.fluor
    m.i <- rowSums( pmax( S.i, 0 ) )
    k.i <- nrow( S.i )

    x <- x.mat[ i, ]

    for ( it in seq_len( n.iter ) ) {

      xc <- pmax( x, 0 )
      mu <- as.vector( xc %*% S.i )

      d.vec <- read.var + pmax( mu, 0 ) / kappa

      u.list <- list()
      c.vals <- numeric( 0 )

      x.max <- max( xc )
      act   <- if ( x.max > 0 ) which( xc >= active.threshold * x.max ) else integer( 0 )

      for ( f in act ) {

        if ( include.spillover ) {
          p.f   <- P.i[ f, ]
          n.f   <- xc[ f ] * m.i[ f ] / kappa.p
          d.vec <- d.vec + n.f * p.f
          u.list[[ length( u.list ) + 1L ]] <- p.f
          c.vals <- c( c.vals, -n.f )
        }

        if ( !is.null( vb.idx ) && f <= fluor.n && !is.na( vb.idx[ f ] ) ) {
          vb <- variant.basis[[ vb.idx[ f ] ]]
          for ( j in seq_along( vb$lambda ) ) {
            u.list[[ length( u.list ) + 1L ]] <- vb$basis[ j, ]
            c.vals <- c( c.vals, xc[ f ]^2 * vb$lambda[ j ] )
          }
        }
      }

      if ( use.af && !is.null( ab.idx ) && !is.na( ab.idx[ g ] ) && xc[ k.i ] > 0 ) {
        ab <- af.basis[[ ab.idx[ g ] ]]
        for ( j in seq_along( ab$lambda ) ) {
          u.list[[ length( u.list ) + 1L ]] <- ab$basis[ j, ]
          c.vals <- c( c.vals, xc[ k.i ]^2 * ab$lambda[ j ] )
        }
      }

      if ( gain.cv > 0 ) {
        u.list[[ length( u.list ) + 1L ]] <- mu
        c.vals <- c( c.vals, gain.cv^2 )
      }

      U <- if ( length( u.list ) > 0L ) do.call( cbind, u.list ) else NULL

      Z <- cbind( t( S.i ), y )

      sigma.out <- if ( method == "woodbury" ) {
        .sigma.solve( Z, d.vec, U, c.vals, ridge = inner.ridge,
                      diagnostics = return.fit.stats )
      } else {
        Sigma <- diag( d.vec, nrow = length( d.vec ) )
        if ( !is.null( U ) ) Sigma <- Sigma + U %*% ( c.vals * t( U ) )
        solve( Sigma, Z )
      }

      W <- if ( return.fit.stats && method == "woodbury" ) sigma.out$W else sigma.out

      XtX <- S.i %*% W[ , seq_len( k.i ), drop = FALSE ]
      Xty <- S.i %*% W[ , k.i + 1L, drop = FALSE ]

      diag( XtX ) <- diag( XtX ) + ridge * mean( abs( diag( XtX ) ) )

      x <- as.vector( solve( XtX, Xty ) )
    }

    x.mat[ i, ] <- x

    if ( return.variance )
      se.mat[ i, ] <- sqrt( pmax( diag( solve( XtX ) ), 0 ) )

    if ( return.fit.stats ) {
      Sigma.inv.r     <- W[ , k.i + 1L ] - as.vector( W[ , seq_len( k.i ), drop = FALSE ] %*% x )
      r               <- y - as.vector( x %*% S.i )
      chisq.vec[ i ]  <- sum( r * Sigma.inv.r )
      df.vec[ i ]     <- length( d.vec ) - k.i
      logdet.vec[ i ] <- if ( method == "woodbury" ) sigma.out$logdet else
        as.numeric( determinant(
          diag( d.vec, nrow = length( d.vec ) ) +
            if ( is.null( U ) ) 0 else U %*% ( c.vals * t( U ) ),
          logarithm = TRUE )$modulus )
    }

    if ( verbose && cell.n >= 20000L && i %% 10000L == 0L )
      message( sprintf( "  %d / %d", i, cell.n ) )
  }

  if ( return.variance || return.fit.stats ) {
    out <- list( unmixed = x.mat )
    if ( return.variance )  out$se     <- se.mat
    if ( return.fit.stats ) {
      out$chisq  <- chisq.vec
      out$df     <- df.vec
      out$logdet <- logdet.vec
      out$p      <- stats::pchisq( chisq.vec, df.vec, lower.tail = FALSE )
    }
    out
  } else {
    x.mat
  }
}
