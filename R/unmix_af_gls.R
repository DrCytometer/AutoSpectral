# unmix_af_gls.R

#' @title Unmix Autofluorescence by GLS with Per-Node Covariance
#'
#' @description
#' Selects an autofluorescence dictionary entry and abundance per cell by
#' maximising a Gaussian log-likelihood in which the AF spectrum itself is
#' uncertain:
#' \deqn{\Sigma_i = diag(read.var) + diag(\mu_i)/\kappa + \alpha_i^2 \Sigma_k^{AF}
#'   + \sum_f x_{if}^2 \Sigma_f^{var}}
#'
#' The AF covariance term softens the discreteness of the dictionary without
#' adding per-cell degrees of freedom: a cell whose true AF falls between two
#' nodes is no longer forced to attribute the mismatch to photon noise. This
#' is deliberately distinct from a continuous AF model, which is known to
#' inflate variance past the Cramer-Rao bound.
#'
#' Applies the same low-rank Woodbury solve as `unmix.gls()` (`.sigma.solve()`)
#' rather than a dense Cholesky factorisation of the full detector-space
#' covariance. `af.model$nodes[[k]]` is already in the same `basis`/`lambda`
#' shape `build.variant.basis()` produces, so the AF term for the candidate
#' under test is folded into the same low-rank accumulator as fluorophore
#' variant uncertainty -- no D x D matrix is ever formed. Cost per candidate
#' per iteration is `O(D k^2 + k^3)` with `k` the active low-rank total
#' (typically well under 20), rather than `O(D^3)`.
#'
#' @param raw.data Numeric matrix (events x detectors).
#' @param af.spectra AF dictionary (nodes x detectors), from
#'   `get.af.spectra(return.model = TRUE)`.
#' @param af.model Per-node model. Defaults to the `"af.model"` attribute of
#'   `af.spectra`. When supplied explicitly, must have one entry per row of
#'   `af.spectra`, in the same order.
#' @param spectra Optional fluorophore spectra (fluorophores x detectors). When
#'   supplied, fluorophores are fitted jointly with AF.
#' @param variant.basis Optional list from `build.variant.basis()`, applied to
#'   the fluorophore terms exactly as in `unmix.gls()`. Ignored unless
#'   `spectra` is supplied. Default `NULL`.
#' @param noise.model Optional list with `read.var` and `counts.per.unit`. When
#'   `NULL`, a flat diagonal is used and only the AF covariance shapes the fit.
#' @param use.af.covariance Logical, include the `alpha^2 Sigma_k^AF` term.
#'   Set `FALSE` for the A/B comparison. Default `TRUE`.
#' @param use.node.prior Logical, add `log(prior_k)` to the score. Default `TRUE`.
#' @param use.abundance.prior Logical, add the per-node lognormal prior on
#'   abundance. Default `FALSE`.
#' @param include.spillover Logical, include the multinomial spillover term
#'   for active fluorophores and the AF row itself. Default `FALSE` -- see
#'   `unmix.gls()` for why this double-counts shot noise under the usual
#'   `kappa` estimation route.
#' @param spillover.kappa Optional scalar spillover `kappa`. Defaults to
#'   `noise.model$kappa.pooled`.
#' @param gain.cv Numeric, multiplicative-gain coefficient of variation; see
#'   `unmix.gls()`. Default `0`.
#' @param active.threshold Numeric, fluorophores below this fraction of the
#'   cell's largest abundance are excluded from the low-rank terms. Default
#'   `1e-3`.
#' @param n.candidates Integer, number of dictionary entries scored per cell,
#'   shortlisted by the joint covariance-weighted score. `Inf` scores all.
#'   Default `5`.
#' @param n.iter Integer, GLS refinement iterations per candidate. Default `2`.
#' @param ridge Numeric ridge added to the `XtX` solve. Default `1e-8`.
#' @param inner.ridge Numeric ridge added inside the Woodbury inner matrix.
#'   Default `1e-10`.
#' @param bend.max.p Numeric in `(0, 1]`. The per-cell AF refinement is
#'   applied only to cells whose chi-square fit p-value falls below this,
#'   i.e. cells whose residual is inconsistent with the fitted noise model
#'   and therefore carry autofluorescence the dictionary cannot represent.
#'   Cells that fit are left on their discrete dictionary entry. The
#'   refinement can absorb genuine dye signal, so restricting it to misfit
#'   cells confines that exposure to the small population that gains from
#'   it. Requires a fitted `noise.model` to be meaningful, since the
#'   p-value is only calibrated against a real noise model. `1` bends every
#'   cell, which is not recommended: an ungated refinement
#'   absorbs dim fluorophore signal across the whole sample, while gating it
#'   to misfit cells retains the autofluorescence gain and removes the cost.
#'   Default `1e-6`.
#' @param return.af.spectra Logical, return the per-cell AF spectrum refined
#'   within the winning node's covariance basis (the best linear unbiased
#'   predictor of the shape deviation), rather than only the dictionary entry.
#'   Adds no free parameters: the deviation is determined by `Sigma` and the
#'   node's `lambda`, which the likelihood already uses. Costs an
#'   events x detectors matrix. Default `FALSE`.
#' @param verbose Logical. Default `TRUE`.
#'
#' @return A list with `af.index`, `af.spectrum.name`, `af.abundance`,
#'   `loglik`, `chisq`, `df`, `logdet`, `p`, `unmixed` (when `spectra`
#'   supplied), `se`, and `af.spectrum` (when `return.af.spectra = TRUE`).
#'   `af.spectrum` rows carry the same scale convention as the dictionary
#'   rows they refine, so the reconstructed AF contribution for a cell is
#'   `af.abundance * af.spectrum`; the refined rows are not re-normalised,
#'   since re-normalising would require rescaling `af.abundance` to match.
#'
#' @export

unmix.af.gls <- function(
    raw.data,
    af.spectra,
    af.model            = attr( af.spectra, "af.model" ),
    spectra             = NULL,
    variant.basis       = NULL,
    noise.model         = NULL,
    use.af.covariance   = TRUE,
    use.node.prior      = TRUE,
    use.abundance.prior = FALSE,
    include.spillover   = FALSE,
    spillover.kappa     = NULL,
    gain.cv             = 0,
    active.threshold    = 1e-3,
    n.candidates        = 5L,
    n.iter              = 2L,
    ridge               = 1e-8,
    inner.ridge         = 1e-10,
    return.af.spectra   = FALSE,
    bend.max.p          = 1e-6,
    verbose             = TRUE
) {

  if ( is.null( af.model ) &&
       ( use.af.covariance || use.node.prior || use.abundance.prior ) )
    stop( "No AF model supplied. Re-run get.af.spectra() with ",
          "`return.model = TRUE`, or set `use.af.covariance = use.node.prior ",
          "= use.abundance.prior = FALSE`.", call. = FALSE )

  if ( !is.null( af.model ) &&
       ( length( af.model$nodes ) != nrow( af.spectra ) ||
         !identical( names( af.model$nodes ), rownames( af.spectra ) ) ) )
    stop( "`af.model` does not align with `af.spectra`: node count or row ",
          "names differ. `af.model$nodes[[k]]` is indexed positionally ",
          "against `af.spectra` row `k`; a mismatch here silently scores ",
          "cells against the wrong node's covariance. Re-attach ",
          "`attr(af.spectra, \"af.model\")` from the same `get.af.spectra()` ",
          "call that produced this `af.spectra`, without reordering or ",
          "subsetting either one afterward.", call. = FALSE )

  # spectral reference matrix (if supplied) must contain exactly one row per
  # fluorophore before any unmixing method runs
  if ( !is.null( spectra ) ) check.spectra.duplicates( spectra )

  det.names  <- colnames( af.spectra )
  raw.data   <- as.matrix( raw.data )[ , det.names, drop = FALSE ]
  af.spectra <- as.matrix( af.spectra )

  cell.n <- nrow( raw.data )
  node.n <- nrow( af.spectra )
  det.n  <- length( det.names )

  if ( is.null( noise.model ) ) {

    # A flat unit diagonal is not neutral: alpha^2 Sigma_k^AF is then 5-7
    # orders of magnitude larger, so the AF covariance directions are given
    # effectively infinite variance and discarded rather than down-weighted.
    # Fall back to a per-detector scale taken from the lower half of the data,
    # which is commensurate with the AF term.

    read.var <- apply( raw.data, 2, function( v ) {
      v <- v[ is.finite( v ) ]
      q <- stats::quantile( v, c( 0.1587, 0.5 ), names = FALSE )
      max( ( q[ 2 ] - q[ 1 ] )^2, .Machine$double.eps )
    } )
    kappa        <- rep( Inf, det.n )
    kappa.pooled <- Inf

    warning( "`noise.model` is NULL; using a per-detector scale derived from ",
             "the data. This is a diagnostic fallback only -- supply a fitted ",
             "noise model for meaningful results.", call. = FALSE )

  } else {

    read.var     <- as.numeric( noise.model$read.var[ det.names ] )
    kappa        <- as.numeric( noise.model$counts.per.unit[ det.names ] )
    kappa.pooled <- noise.model$kappa.pooled

    if ( anyNA( read.var ) || anyNA( kappa ) )
      stop( "`noise.model` does not cover all detectors in `af.spectra`. ",
            "Missing: ", paste( det.names[ is.na( read.var ) | is.na( kappa ) ],
                                collapse = ", " ), call. = FALSE )

    if ( any( read.var <= 0 ) || any( kappa <= 0 ) )
      stop( "`noise.model` contains non-positive read.var or ",
            "counts.per.unit at one or more detectors; Sigma would not be ",
            "positive definite there. Check that detector names match ",
            "between `noise.model` and `af.spectra`.", call. = FALSE )
  }

  kappa.p <- if ( is.null( spillover.kappa ) ) kappa.pooled else spillover.kappa

  fluor.n <- if ( is.null( spectra ) ) 0L else nrow( spectra )
  S.base  <- if ( is.null( spectra ) ) NULL else
    as.matrix( spectra )[ , det.names, drop = FALSE ]

  .prob.rows <- function( m ) {
    m <- pmax( m, 0 )
    m / pmax( rowSums( m ), 1e-12 )
  }
  P.fluor <- if ( !is.null( S.base ) ) .prob.rows( S.base ) else NULL

  # index of variant bases aligned to the rows of spectra
  vb.idx <- if ( is.null( variant.basis ) || is.null( S.base ) ) NULL else
    match( rownames( S.base ), names( variant.basis ) )

  if ( !is.null( vb.idx ) && anyNA( vb.idx ) )
    warning( "No variant basis for: ",
             paste( rownames( S.base )[ is.na( vb.idx ) ], collapse = ", " ),
             ". These fluorophores are modelled with zero spectral ",
             "uncertainty.", call. = FALSE )

  # ---------------------------------------------------------------------------
  # Candidate shortlist and starting values
  # ---------------------------------------------------------------------------
  # Shortlisting on a raw least-squares fit to each AF spectrum is valid only
  # for unstained data; on a stained sample that fit is dominated by
  # fluorophore signal. Use the established joint scoring instead:
  # covariance-weighted fluorophore error times residual error, with AF
  # intensity estimated in the space orthogonal to the fluorophores.

  # Scalar least-squares score, computed unconditionally: minimising it is
  # equivalent to maximising cosine similarity to the raw vector, so it is
  # an exact shortlister for the closest-matching AF shape. On stained data
  # it is polluted by fluorophore signal, so it supplements rather than
  # replaces the joint score below.
  aa        <- rowSums( af.spectra^2 )
  ya        <- raw.data %*% t( af.spectra )
  score.rss <- rowSums( raw.data^2 ) - sweep( ya^2, 2L, aa, "/" )

  if ( is.null( S.base ) ) {
    score <- score.rss
    joint <- NULL
  } else {
    joint <- assign.af.joint.cov(
      raw.data      = raw.data,
      spectra       = S.base,
      af.spectra    = af.spectra,
      return.scores = TRUE
    )
    score <- joint$scores
  }

  n.cand <- min( n.candidates, node.n )

  log.prior <- if ( use.node.prior && !is.null( af.model ) )
    log( af.model$prior ) else rep( 0, node.n )

  best.k      <- integer( cell.n )
  best.alpha  <- numeric( cell.n )
  best.ll     <- rep( -Inf, cell.n )
  best.chisq  <- rep( NA_real_, cell.n )
  best.df     <- rep( NA_integer_, cell.n )
  best.logdet <- rep( NA_real_, cell.n )
  best.x      <- matrix( NA_real_, cell.n, fluor.n + 1L )
  best.se     <- matrix( NA_real_, cell.n, fluor.n + 1L )
  best.af.spectrum <- if ( return.af.spectra )
    matrix( NA_real_, cell.n, det.n, dimnames = list( NULL, det.names ) ) else NULL

  if ( verbose )
    message( sprintf( "AF GLS: %d events, %d nodes, %d candidates each",
                      cell.n, node.n, n.cand ) )

  for ( i in seq_len( cell.n ) ) {

    y    <- raw.data[ i, ]
    # Union shortlist: the joint score protects fluorophore space but does
    # not always rank the closest-matching AF shape highly, while the
    # scalar least-squares score is exactly the cosine criterion. Taking
    # candidates from both guarantees the closest node is always scored by
    # the likelihood, which remains the arbiter. On heavily stained cells
    # the scalar half spends its slots poorly, but the joint half is
    # unaffected.
    cand <- unique( c(
      order( score[ i, ] )[ seq_len( ceiling( n.cand / 2 ) ) ],
      order( score.rss[ i, ] )[ seq_len( ceiling( n.cand / 2 ) ) ]
    ) )

    for ( k in cand ) {

      S.i <- if ( is.null( S.base ) )
        af.spectra[ k, , drop = FALSE ] else
          rbind( S.base, AF = af.spectra[ k, ] )
      m.i <- rowSums( pmax( S.i, 0 ) )
      p.n <- nrow( S.i )

      P.i <- if ( include.spillover ) {
        p.af.k <- pmax( af.spectra[ k, ], 0 )
        p.af.k <- p.af.k / max( sum( p.af.k ), 1e-12 )
        if ( is.null( P.fluor ) ) rbind( AF = p.af.k ) else
          rbind( P.fluor, AF = p.af.k )
      } else NULL

      nd <- if ( use.af.covariance ) af.model$nodes[[ k ]] else NULL

      # start from the joint-scoring solution rather than a fresh OLS: the
      # fluorophore estimate with AF variant k already removed, plus the AF
      # abundance implied by the orthogonal residual
      x <- if ( is.null( joint ) ) {
        tryCatch( as.vector( solve( tcrossprod( S.i ), S.i %*% y ) ),
                  error = function( e ) rep( 0, p.n ) )
      } else {
        c( joint$unmixed[ i, ] - joint$k.matrix[ i, k ] * joint$v.library[ , k ],
           joint$k.matrix[ i, k ] )
      }

      W <- NULL; XtX <- NULL; sigma.out <- NULL; iter.ok <- TRUE

      # always build and solve Sigma at least once, so that n.iter = 0
      # scores the starting values rather than skipping the cell entirely
      for ( it in seq_len( max( n.iter, 1L ) ) ) {

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

          if ( f <= fluor.n && !is.null( vb.idx ) && !is.na( vb.idx[ f ] ) ) {
            vb <- variant.basis[[ vb.idx[ f ] ]]
            for ( j in seq_along( vb$lambda ) ) {
              u.list[[ length( u.list ) + 1L ]] <- vb$basis[ j, ]
              c.vals <- c( c.vals, xc[ f ]^2 * vb$lambda[ j ] )
            }
          }
        }

        if ( !is.null( nd ) && xc[ p.n ] > 0 ) {
          for ( j in seq_along( nd$lambda ) ) {
            u.list[[ length( u.list ) + 1L ]] <- nd$basis[ j, ]
            c.vals <- c( c.vals, xc[ p.n ]^2 * nd$lambda[ j ] )
          }
        }

        if ( gain.cv > 0 ) {
          u.list[[ length( u.list ) + 1L ]] <- mu
          c.vals <- c( c.vals, gain.cv^2 )
        }

        U <- if ( length( u.list ) > 0L ) do.call( cbind, u.list ) else NULL

        Z <- cbind( t( S.i ), y )
        sigma.out <- tryCatch(
          .sigma.solve( Z, d.vec, U, c.vals, ridge = inner.ridge, diagnostics = TRUE ),
          error = function( e ) NULL )

        if ( is.null( sigma.out ) ) { iter.ok <- FALSE; break }

        W   <- sigma.out$W
        XtX <- S.i %*% W[ , seq_len( p.n ), drop = FALSE ]
        diag( XtX ) <- diag( XtX ) + ridge * mean( abs( diag( XtX ) ) )

        if ( n.iter > 0L )
          x <- tryCatch(
            as.vector( solve( XtX, S.i %*% W[ , p.n + 1L, drop = FALSE ] ) ),
            error = function( e ) x )
      }

      if ( !iter.ok ) next

      Sigma.inv.r <- W[ , p.n + 1L ] - as.vector( W[ , seq_len( p.n ), drop = FALSE ] %*% x )
      r           <- y - as.vector( x %*% S.i )
      chisq.k     <- sum( r * Sigma.inv.r )
      df.k        <- det.n - p.n

      # Best linear unbiased predictor of this cell's AF shape deviation
      # within the candidate node's basis. The dictionary entry is a
      # centroid; cells whose true AF falls between nodes, or beyond the
      # brightest node the SOM could place, are otherwise forced to charge
      # the mismatch to photon noise. The deviation is shrunk by the node's
      # own variance against this cell's noise, so it vanishes for dim cells
      # and is only substantial where the residual is large relative to
      # Sigma.
      af.shape.k <- NULL
      if ( return.af.spectra ) {

        af.shape.k <- af.spectra[ k, ]

        fit.p <- if ( is.finite( chisq.k ) && df.k > 0L )
          stats::pchisq( chisq.k, df.k, lower.tail = FALSE ) else 0

        if ( !is.null( nd ) && x[ p.n ] > 0 && fit.p < bend.max.p ) {
          t.hat      <- x[ p.n ] * nd$lambda * as.vector( nd$basis %*% Sigma.inv.r )
          af.shape.k <- af.shape.k + as.vector( crossprod( nd$basis, t.hat ) )
        }
      }

      ll <- -0.5 * ( sigma.out$logdet + chisq.k ) + log.prior[ k ]

      if ( use.abundance.prior && !is.null( nd ) ) {
        if ( is.finite( nd$log.alpha.sd ) && nd$log.alpha.sd > 0 && x[ p.n ] > 0 )
          ll <- ll + stats::dnorm( log( x[ p.n ] ), nd$log.alpha.mu,
                                   nd$log.alpha.sd, log = TRUE )
      }

      if ( ll > best.ll[ i ] ) {
        best.ll[ i ]     <- ll
        best.k[ i ]      <- k
        best.alpha[ i ]  <- x[ p.n ]
        best.chisq[ i ]  <- chisq.k
        best.df[ i ]     <- df.k
        best.logdet[ i ] <- sigma.out$logdet
        best.x[ i, ]     <- x
        best.se[ i, ]    <- tryCatch(
          sqrt( pmax( diag( solve( XtX ) ), 0 ) ),
          error = function( e ) rep( NA_real_, p.n ) )
        if ( return.af.spectra ) best.af.spectrum[ i, ] <- af.shape.k
      }
    }

    if ( verbose && cell.n >= 20000L && i %% 10000L == 0L )
      message( sprintf( "  %d / %d", i, cell.n ) )
  }

  n.failed <- sum( best.k == 0L )
  if ( n.failed > 0L )
    warning( n.failed, " of ", cell.n, " cell(s) had no candidate that ",
             "produced a finite score; their `af.index` is 0 and ",
             "`unmixed`/`se` are NA. This should be rare -- check for ",
             "non-finite values in `raw.data`.", call. = FALSE )

  colnames( best.x ) <- colnames( best.se ) <-
    c( if ( fluor.n > 0 ) rownames( S.base ), "AF" )

  list(
    af.index         = best.k,
    af.spectrum.name = ifelse( best.k == 0L, NA_character_,
                               rownames( af.spectra )[ pmax( best.k, 1L ) ] ),
    af.abundance     = best.alpha,
    loglik           = best.ll,
    chisq            = best.chisq,
    df               = best.df,
    logdet           = best.logdet,
    p                = stats::pchisq( best.chisq, best.df, lower.tail = FALSE ),
    unmixed          = if ( fluor.n > 0 ) best.x else NULL,
    se               = best.se,
    af.spectrum      = best.af.spectrum
  )
}
