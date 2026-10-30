# noise_model_fit.R

# Detector noise model: Var(y_d) = sigma0_d^2 + mu_d / kappa_d (+ c_d mu_d^2),
# fitted jointly over groups of events (single-stain control positives, the
# unstained cells, bead negatives), each projected by OLS onto its own design.
# `estimate.noise.model()` builds the groups; this file holds the estimator.
#
# A group is a list with
#   name    character
#   kind    "control", "unstained" or "beads"
#   raw     events x detectors matrix
#   design  rows x detectors matrix, or NULL for a mean-only group
#   offset  length-D vector subtracted before projecting and added back to
#           the fitted signal (bead background), or NULL
#   bin.row design row whose coefficient is the binning variable (1)
#   cls     particle class, "cells" or "beads"
#
# With heteroscedastic noise the OLS residual of detector d mixes in the
# variance of every detector: E[r_d^2] = sum_e M_de V_e(mu_e), with
# M = (I - H) o (I - H) and H the detector-space hat matrix of the group's
# design. Events are split into quantile bins of the group's own abundance,
# and each bin gives one equation per detector: its robust residual variance
# against the bin-mean fitted signal. All groups, bins and detectors form one
# non-negative least-squares problem for s = sigma0^2 and g = 1 / kappa (and
# optionally c, and one c_g per control), weighted by n / v^2 (the variance
# of a variance), with IRLS passes on the fitted variance. The first solve
# uses the squared normalised MAD of each bin; later rounds standardise the
# residuals by the current model before taking the MAD, which removes the
# scale-mixture bias of a MAD over events of different brightness.
#
# Every row of a (group, bin) block is a column-scaled copy of M, so the
# normal equations are accumulated block by block from M' W M; the stacked
# design matrix is never formed.
#
# Afterwards s_d is clamped above at the lowest bin variance any group shows
# at d (divided by M_dd, widened by the bin variance's sampling SD) and below
# at floor.frac of the median s, so that a dark detector on bead events
# cannot take an unbounded weight; g is refitted with s fixed when a clamp
# applied. A detector with no shot-noise information (g = 0) takes the median
# g and is flagged.
#
# Options:
#   per.control.quadratic  one non-negative c_g mu^2 per control group in the
#     joint fit, absorbing control-specific excess variance so that it does
#     not bias kappa; stored for calibration and QC, not used for unmixing.
#   sigma0.by.class  one sigma0 per particle class, kappa shared.
#   sigma0.fallback = "laser"  a clamped or floored detector takes the median
#     sigma0 of the unflagged detectors on its laser, per class (with two
#     classes, the other class's value scaled by the per-laser ratio when the
#     other class identifies it); fewer than 3 such detectors: the class's
#     global median. Never above the detector's upper bound.

.NOISE.FLAG.SIGMA0.HIGH     <- 1L
.NOISE.FLAG.SIGMA0.LOW      <- 2L
.NOISE.FLAG.KAPPA.POOLED    <- 4L
.NOISE.FLAG.SIGMA0.FALLBACK <- 8L
.NOISE.CLASSES   <- c( "cells", "beads" )
.NOISE.TINY      <- 1e-12
.NOISE.MAD.VAR.REL <- 2.72

#' @title Default Noise Model Estimator Settings
#'
#' @description
#' Settings of the joint noise-model fit used by `estimate.noise.model()`.
#'
#' @param n.bins Integer, quantile bins per group. Default `20`.
#' @param min.bin.events Integer, fewest events in a bin. Default `50`.
#' @param fit.quadratic Logical, fit a per-detector `c mu^2` term. Default
#'   `FALSE`.
#' @param n.irls Integer, IRLS passes per solve. Default `2`.
#' @param n.restandardise Integer, rounds that recompute the bin variances
#'   from residuals standardised by the current model. Default `2`.
#' @param floor.frac Numeric, lower clamp on `sigma0^2` as a fraction of its
#'   median. Default `0.01`.
#' @param bound.sds Numeric, sampling SDs added to the upper clamp on
#'   `sigma0^2`. Default `1`.
#' @param per.control.quadratic Logical. Default `FALSE`.
#' @param sigma0.by.class Logical. Default `FALSE`.
#' @param sigma0.fallback `"floor"` (default) or `"laser"`.
#'
#' @return A named list.
#'
#' @noRd

.noise.settings <- function(
    n.bins                = 20L,
    min.bin.events        = 50L,
    fit.quadratic         = FALSE,
    n.irls                = 2L,
    n.restandardise       = 2L,
    floor.frac            = 0.01,
    bound.sds             = 1,
    per.control.quadratic = FALSE,
    sigma0.by.class       = FALSE,
    sigma0.fallback       = c( "floor", "laser" )
) {
  list(
    n.bins = as.integer( n.bins ), min.bin.events = as.integer( min.bin.events ),
    fit.quadratic = isTRUE( fit.quadratic ), n.irls = as.integer( n.irls ),
    n.restandardise = as.integer( n.restandardise ), floor.frac = floor.frac,
    bound.sds = bound.sds, per.control.quadratic = isTRUE( per.control.quadratic ),
    sigma0.by.class = isTRUE( sigma0.by.class ),
    sigma0.fallback = match.arg( sigma0.fallback )
  )
}

# Median of each column of a finite matrix: one sort of every column at once.
.col.medians <- function( x ) {
  m <- nrow( x )
  if ( m == 0L ) return( rep( NA_real_, ncol( x ) ) )
  s <- matrix( x[ order( col( x ), x ) ], m )
  if ( m %% 2L == 1L ) s[ ( m + 1L ) %/% 2L, ] else ( s[ m %/% 2L, ] + s[ m %/% 2L + 1L, ] ) / 2
}

# Squared normalised MAD of each column.
.robust.var <- function( x ) {
  x <- as.matrix( x )
  ( 1.4826 * .col.medians( abs( sweep( x, 2, .col.medians( x ), "-" ) ) ) )^2
}

.noise.flag.text <- function( flags ) {
  labels <- c( "sigma0 clamped high", "sigma0 raised to floor", "kappa pooled",
               "sigma0 from same laser" )
  bits <- c( .NOISE.FLAG.SIGMA0.HIGH, .NOISE.FLAG.SIGMA0.LOW,
             .NOISE.FLAG.KAPPA.POOLED, .NOISE.FLAG.SIGMA0.FALLBACK )
  vapply( flags, function( f )
    paste( labels[ bitwAnd( as.integer( f ), bits ) > 0 ], collapse = ", " ), character( 1 ) )
}

# ---------------------------------------------------------------------------
# Projection and binning
# ---------------------------------------------------------------------------

# Orthonormal rows (k x D) spanning `design`, dropping directions below rtol
# of the largest singular value.
.orthonormal.rows <- function( design, rtol = 1e-8 ) {
  sv   <- svd( design )
  keep <- sv$d > rtol * max( sv$d[ 1 ], .NOISE.TINY )
  t( sv$v[ , keep, drop = FALSE ] )
}

# (I - H) o (I - H); for a mean-only group (1 - 1/n) I.
.mixing.matrix <- function( q, n.detectors, n.events = 0L ) {
  if ( is.null( q ) )
    return( diag( n.detectors ) * ( 1 - 1 / max( n.events, 2 ) ) )
  resid.maker <- diag( n.detectors ) - crossprod( q )
  resid.maker * resid.maker
}

.pinv <- function( a ) {
  sv   <- svd( a )
  keep <- sv$d > 1e-15 * max( sv$d )
  sv$v[ , keep, drop = FALSE ] %*% ( t( sv$u[ , keep, drop = FALSE ] ) / sv$d[ keep ] )
}

.noise.project <- function( group, keep.centred = TRUE ) {
  raw <- as.matrix( group$raw )
  storage.mode( raw ) <- "double"
  n <- nrow( raw )
  d <- ncol( raw )
  offset <- if ( is.null( group$offset ) ) numeric( d ) else as.numeric( group$offset )
  centred <- sweep( raw, 2, offset, "-" )

  if ( is.null( group$design ) ) {
    mean.c <- colMeans( centred )
    fit    <- matrix( mean.c + offset, n, d, byrow = TRUE )
    return( list( coef.bin = numeric( n ), fit = fit,
                  resid = sweep( centred, 2, mean.c, "-" ),
                  centred = if ( keep.centred ) centred else NULL, q = NULL,
                  M = .mixing.matrix( NULL, d, n ) ) )
  }

  design <- as.matrix( group$design )
  q      <- .orthonormal.rows( design )
  within <- ( centred %*% t( q ) ) %*% q
  coef   <- centred %*% .pinv( design )
  bin.row <- if ( is.null( group$bin.row ) ) 1L else group$bin.row
  list( coef.bin = coef[ , bin.row ], fit = sweep( within, 2, offset, "+" ),
        resid = centred - within, centred = if ( keep.centred ) centred else NULL,
        q = q, M = .mixing.matrix( q, d ) )
}

# Event indices of each quantile bin of the binning variable; one bin for a
# mean-only group; bins with fewer than min.bin.events dropped.
.bin.members <- function( proj, n.bins, min.bin.events, mean.only ) {
  n <- nrow( proj$fit )
  if ( mean.only || n < 2L * min.bin.events )
    return( if ( n >= min.bin.events ) list( seq_len( n ) ) else list() )
  nb    <- as.integer( max( 1L, min( n.bins, n %/% min.bin.events ) ) )
  edges <- stats::quantile( proj$coef.bin, seq( 0, 1, length.out = nb + 1L ), names = FALSE )
  which.bin <- pmin( pmax( findInterval( proj$coef.bin, edges ) - 1L, 0L ), nb - 1L )
  members <- lapply( seq_len( nb ) - 1L, function( b ) which( which.bin == b ) )
  members[ lengths( members ) >= min.bin.events ]
}

# Per bin and detector: the residual variance, the mean fitted signal
# (clamped at 0) and its mean square. With `pred` (the current model's
# residual variance per event), the variance is the bin mean of pred times
# the squared normalised MAD of resid / sqrt(pred).
.bin.stats <- function( proj, members, pred = NULL ) {
  mu <- pmax( proj$fit, 0 )
  d  <- ncol( mu )
  if ( length( members ) == 0L ) {
    empty <- matrix( 0, 0, d )
    return( list( n = integer( 0 ), var = empty, mu = empty, mu2 = empty, centre = numeric( 0 ) ) )
  }
  v <- if ( is.null( pred ) ) {
    lapply( members, function( m ) .robust.var( proj$resid[ m, , drop = FALSE ] ) )
  } else {
    lapply( members, function( m )
      colMeans( pred[ m, , drop = FALSE ] ) *
        .robust.var( proj$resid[ m, , drop = FALSE ] / sqrt( pred[ m, , drop = FALSE ] ) ) )
  }
  list(
    n      = lengths( members ),
    var    = do.call( rbind, v ),
    mu     = do.call( rbind, lapply( members, function( m ) colMeans( mu[ m, , drop = FALSE ] ) ) ),
    mu2    = do.call( rbind, lapply( members, function( m ) colMeans( mu[ m, , drop = FALSE ]^2 ) ) ),
    centre = vapply( members, function( m ) stats::median( proj$coef.bin[ m ] ), numeric( 1 ) )
  )
}

# ---------------------------------------------------------------------------
# Rows, layout and normal equations
# ---------------------------------------------------------------------------

# One row per group, bin and detector (the stacked problem's targets and
# bookkeeping), and one block per (group, bin).
.noise.rows <- function( projections, bins ) {
  b <- n <- grp <- bn <- det <- mu.own <- m.diag <- list()
  blocks <- list()
  pos <- 0L
  for ( gi in seq_along( projections ) ) {
    proj <- projections[[ gi ]]
    bi   <- bins[[ gi ]]
    d    <- nrow( proj$M )
    for ( k in seq_along( bi$n ) ) {
      idx <- pos + seq_len( d )
      pos <- pos + d
      blocks[[ length( blocks ) + 1L ]] <- list(
        group = gi, rows = idx, mu = bi$mu[ k, ], mu2 = bi$mu2[ k, ],
        v = as.vector( proj$M %*% bi$mu2[ k, ] ) )
      b[[ length( b ) + 1L ]]           <- bi$var[ k, ]
      n[[ length( n ) + 1L ]]           <- rep( as.numeric( bi$n[ k ] ), d )
      grp[[ length( grp ) + 1L ]]       <- rep( gi, d )
      bn[[ length( bn ) + 1L ]]         <- rep( k, d )
      det[[ length( det ) + 1L ]]       <- seq_len( d )
      mu.own[[ length( mu.own ) + 1L ]] <- bi$mu[ k, ]
      m.diag[[ length( m.diag ) + 1L ]] <- diag( proj$M )
    }
  }
  if ( length( b ) == 0L )
    stop( "No group has enough events for a noise model bin.", call. = FALSE )
  list( b = unlist( b ), n = unlist( n ), group = unlist( grp ), bin = unlist( bn ),
        det = unlist( det ), mu.own = unlist( mu.own ), m.diag = unlist( m.diag ),
        blocks = blocks )
}

# Column layout: K sigma0 blocks, g, optional c, then one c_g per control.
.noise.layout <- function( d, k, cls.of.group, quadratic, controls ) {
  s.cols <- lapply( seq_len( k ), function( j ) ( j - 1L ) * d + seq_len( d ) )
  g.cols <- k * d + seq_len( d )
  c.cols <- if ( quadratic ) ( k + 1L ) * d + seq_len( d ) else integer( 0 )
  cg.start <- ( k + 1L + as.integer( quadratic ) ) * d
  cg.col <- integer( length( cls.of.group ) )
  cg.col[ controls ] <- cg.start + seq_along( controls )
  list( d = d, k = k, cls.of.group = cls.of.group, quadratic = quadratic, controls = controls,
        s.cols = s.cols, g.cols = g.cols, c.cols = c.cols, cg.col = cg.col,
        n.cols = cg.start + length( controls ) )
}

# Normal equations G = A' W A and h = A' W b of the stacked problem, with W
# the row weights (zero for rows left out).
.noise.normal.equations <- function( rows, projections, lay, w ) {
  G <- matrix( 0, lay$n.cols, lay$n.cols )
  h <- numeric( lay$n.cols )
  b <- ifelse( w > 0, rows$b, 0 )
  d <- lay$d
  gc <- lay$g.cols
  cc <- lay$c.cols
  # only the upper triangle of blocks is accumulated; the lower is mirrored at the end
  for ( blk in rows$blocks ) {
    wb <- w[ blk$rows ]
    if ( !any( wb > 0 ) ) next
    M  <- projections[[ blk$group ]]$M
    P  <- crossprod( M, wb * M )
    u  <- as.vector( crossprod( M, wb * b[ blk$rows ] ) )
    sc <- lay$s.cols[[ lay$cls.of.group[ blk$group ] ]]
    mu <- blk$mu
    G[ sc, sc ] <- G[ sc, sc ] + P
    G[ sc, gc ] <- G[ sc, gc ] + P * rep( mu, each = d )
    G[ gc, gc ] <- G[ gc, gc ] + P * tcrossprod( mu )
    h[ sc ] <- h[ sc ] + u
    h[ gc ] <- h[ gc ] + mu * u
    if ( lay$quadratic ) {
      mu2 <- blk$mu2
      G[ sc, cc ] <- G[ sc, cc ] + P * rep( mu2, each = d )
      G[ gc, cc ] <- G[ gc, cc ] + P * tcrossprod( mu, mu2 )
      G[ cc, cc ] <- G[ cc, cc ] + P * tcrossprod( mu2 )
      h[ cc ] <- h[ cc ] + mu2 * u
    }
    j <- lay$cg.col[ blk$group ]
    if ( j > 0L ) {
      r <- as.vector( crossprod( M, wb * blk$v ) )
      G[ sc, j ] <- G[ sc, j ] + r
      G[ gc, j ] <- G[ gc, j ] + mu * r
      if ( lay$quadratic ) G[ cc, j ] <- G[ cc, j ] + blk$mu2 * r
      G[ j, j ] <- G[ j, j ] + sum( wb * blk$v^2 )
      h[ j ]    <- h[ j ] + sum( wb * blk$v * b[ blk$rows ] )
    }
  }
  # mirror: every block above sits at rows < columns in the layout order
  # (sigma0 blocks, g, c, c_g), except the symmetric diagonal blocks
  upper.only <- upper.tri( G )
  diag.blocks <- matrix( FALSE, lay$n.cols, lay$n.cols )
  for ( sc in lay$s.cols ) diag.blocks[ sc, sc ] <- TRUE
  diag.blocks[ gc, gc ] <- TRUE
  if ( lay$quadratic ) diag.blocks[ cc, cc ] <- TRUE
  lower <- t( G )
  take  <- !diag.blocks & !upper.only
  G[ take ] <- lower[ take ]
  list( G = G, h = h )
}

# Predicted residual variance of every row for a solution vector.
.noise.row.predict <- function( rows, projections, lay, sol ) {
  pred <- numeric( length( rows$b ) )
  g  <- sol[ lay$g.cols ]
  cq <- if ( lay$quadratic ) sol[ lay$c.cols ] else 0
  for ( blk in rows$blocks ) {
    M <- projections[[ blk$group ]]$M
    s <- sol[ lay$s.cols[[ lay$cls.of.group[ blk$group ] ]] ]
    val <- as.vector( M %*% ( s + g * blk$mu + cq * blk$mu2 ) )
    j <- lay$cg.col[ blk$group ]
    if ( j > 0L ) val <- val + sol[ j ] * blk$v
    pred[ blk$rows ] <- val
  }
  pred
}

.noise.solve.rows <- function( rows, projections, lay, ok, n.irls ) {
  w <- ifelse( ok, rows$n / pmax( rows$b, .NOISE.TINY )^2, 0 )
  ne  <- .noise.normal.equations( rows, projections, lay, w )
  sol <- .nnls.normal( ne$G, ne$h )
  for ( pass in seq_len( max( 0L, n.irls ) ) ) {
    pred <- .noise.row.predict( rows, projections, lay, sol )
    w    <- ifelse( ok, rows$n / pmax( pred, .NOISE.TINY )^2, 0 )
    ne   <- .noise.normal.equations( rows, projections, lay, w )
    sol  <- .nnls.normal( ne$G, ne$h )
  }
  sol
}

# Split a solution: s (K x D), g, c, c_g (one per group).
.noise.split <- function( sol, lay, n.groups ) {
  s  <- matrix( unlist( lapply( lay$s.cols, function( cc ) sol[ cc ] ) ), nrow = lay$k, byrow = TRUE )
  cq <- if ( lay$quadratic ) sol[ lay$c.cols ] else numeric( lay$d )
  cg <- numeric( n.groups )
  cg[ lay$cg.col > 0L ] <- sol[ lay$cg.col[ lay$cg.col > 0L ] ]
  list( s = s, g = sol[ lay$g.cols ], c = cq, cg = cg )
}

.noise.join <- function( s, g, cq, cg, lay ) {
  sol <- numeric( lay$n.cols )
  for ( j in seq_len( lay$k ) ) sol[ lay$s.cols[[ j ]] ] <- s[ j, ]
  sol[ lay$g.cols ] <- g
  if ( lay$quadratic ) sol[ lay$c.cols ] <- cq
  sol[ lay$cg.col[ lay$cg.col > 0L ] ] <- cg[ lay$cg.col > 0L ]
  sol
}

# ---------------------------------------------------------------------------
# Clamps and fallback
# ---------------------------------------------------------------------------

# Upper bound on sigma0_d^2: V_d <= v_bd / M_dd in every bin. Each bin's
# variance has a relative sampling SD of about sqrt(2.72 / n), and the
# smallest of N such estimates lies about sqrt(2 ln N) SDs low, so each
# bound is divided by 1 - z * rel.sd with z = sqrt(2 ln N_d) + extra.sds
# before the minimum is taken; bins with z * rel.sd >= 0.9 are skipped.
.sigma0.upper.bound <- function( rows, ok, d, extra.sds ) {
  rel    <- sqrt( .NOISE.MAD.VAR.REL / pmax( rows$n, 1 ) )
  n.rows <- tabulate( rows$det[ ok ], nbins = d )
  z      <- sqrt( 2 * log( pmax( n.rows, 2 ) ) ) + extra.sds
  shrink <- 1 - z[ rows$det ] * rel
  usable <- ok & shrink > 0.1
  bound  <- ifelse( usable, rows$b / pmax( rows$m.diag, .NOISE.TINY ) / ifelse( usable, shrink, 1 ), Inf )
  upper  <- rep( Inf, d )
  agg    <- tapply( bound, rows$det, min )
  upper[ as.integer( names( agg ) ) ] <- agg
  upper
}

# sigma0^2 (K x D) with each flagged detector of each class replaced from
# the same laser, capped at its upper bound; returns list(s, applied).
.laser.fallback <- function( s, flagged, upper, lasers ) {
  k.n <- nrow( s )
  d   <- ncol( s )
  las <- if ( is.null( lasers ) ) rep( "", d ) else as.character( lasers )
  sd  <- sqrt( pmax( s, 0 ) )
  out <- sd
  applied <- matrix( FALSE, k.n, d )
  for ( k in seq_len( k.n ) ) {
    ok.k <- !flagged[ k, ] & sd[ k, ] > 0
    if ( !any( ok.k ) ) next
    global.med <- stats::median( sd[ k, ok.k ] )
    for ( i in which( flagged[ k, ] ) ) {
      same  <- las == las[ i ]
      value <- NULL
      for ( j in seq_len( k.n ) ) {
        if ( j == k || flagged[ j, i ] || sd[ j, i ] <= 0 ) next
        both <- !flagged[ k, ] & !flagged[ j, ] & sd[ k, ] > 0 & sd[ j, ] > 0
        ref  <- if ( sum( both & same ) >= 3L ) both & same else both
        if ( any( ref ) ) {
          value <- sd[ j, i ] * stats::median( sd[ k, ref ] / sd[ j, ref ] )
          break
        }
      }
      if ( is.null( value ) ) {
        peers <- ok.k & same
        peers[ i ] <- FALSE
        value <- if ( sum( peers ) >= 3L ) stats::median( sd[ k, peers ] ) else global.med
      }
      out[ k, i ] <- if ( is.finite( upper[ k, i ] ) ) min( value, sqrt( upper[ k, i ] ) ) else value
      applied[ k, i ] <- TRUE
    }
  }
  list( s = out^2, applied = applied )
}

# sum_e M_de V_e(mu_e) per event of one group, with its sigma0 and its c_g.
.group.resid.var <- function( proj, sigma0.sq, g, cq, c.group = 0 ) {
  mu  <- pmax( proj$fit, 0 )
  var <- sweep( sweep( mu, 2, g, "*" ), 2, sigma0.sq, "+" )
  if ( any( cq != 0 ) || c.group != 0 )
    var <- var + sweep( mu * mu, 2, cq + c.group, "*" )
  pmax( var, .NOISE.TINY ) %*% t( proj$M )
}

# ---------------------------------------------------------------------------
# Fit
# ---------------------------------------------------------------------------

#' @title Fit the Detector Noise Model to Groups of Events
#'
#' @description
#' Joint mixing NNLS over every group (see the header of
#' `noise_model_fit.R`).
#'
#' @param groups List of groups.
#' @param detectors Character vector of detector names.
#' @param settings List from `.noise.settings()`.
#' @param lasers Optional character vector, one laser per detector, for the
#'   same-laser fallback.
#'
#' @return A noise model list: `detectors`, `sigma0.sq`, `g` (= 1 / kappa),
#'   `c`, `flags`, `classes`, `sigma0.sq.by.class` and `flags.by.class`
#'   (with a sigma0 per class), `group.c`, `fit.table`, `groups`, `settings`.
#'
#' @noRd

.fit.noise.model <- function( groups, detectors, settings = .noise.settings(), lasers = NULL ) {

  s <- settings
  d <- length( detectors )
  n.groups <- length( groups )
  group.cls <- vapply( groups, function( gr ) if ( is.null( gr$cls ) ) "cells" else gr$cls, character( 1 ) )

  classes <- if ( s$sigma0.by.class ) {
    cl <- .NOISE.CLASSES[ .NOISE.CLASSES %in% group.cls ]
    if ( length( cl ) == 0L ) .NOISE.CLASSES[ 1 ] else cl
  } else {
    if ( n.groups > 0L ) group.cls[ 1 ] else .NOISE.CLASSES[ 1 ]
  }
  cls.of.group <- ifelse( group.cls %in% classes, match( group.cls, classes ), 1L )
  k.n <- length( classes )
  controls <- if ( s$per.control.quadratic )
    which( vapply( groups, function( gr ) gr$kind == "control", logical( 1 ) ) ) else integer( 0 )
  lay <- .noise.layout( d, k.n, cls.of.group, s$fit.quadratic, controls )

  projections <- lapply( groups, .noise.project, keep.centred = FALSE )
  members <- mapply( function( p, gr ) .bin.members( p, s$n.bins, s$min.bin.events, is.null( gr$design ) ),
                     projections, groups, SIMPLIFY = FALSE )

  solve.bins <- function( bins ) {
    rows <- .noise.rows( projections, bins )
    ok   <- is.finite( rows$b ) & rows$b > 0
    if ( !any( ok ) )
      stop( "Every binned residual variance is zero or not finite.", call. = FALSE )
    sol <- .noise.solve.rows( rows, projections, lay, ok, s$n.irls )
    c( list( rows = rows, ok = ok ), .noise.split( sol, lay, n.groups ) )
  }

  floored <- function( s.mat ) {
    out <- s.mat
    for ( k in seq_len( k.n ) ) {
      positive <- s.mat[ k, s.mat[ k, ] > 0 ]
      out[ k, ] <- pmax( s.mat[ k, ], if ( length( positive ) ) s$floor.frac * stats::median( positive ) else .NOISE.TINY )
    }
    out
  }

  bins <- mapply( .bin.stats, projections, members, SIMPLIFY = FALSE )
  fit  <- solve.bins( bins )
  for ( round in seq_len( max( 0L, s$n.restandardise ) ) ) {
    s.cur <- floored( fit$s )
    bins <- lapply( seq_len( n.groups ), function( i )
      .bin.stats( projections[[ i ]], members[[ i ]],
                  .group.resid.var( projections[[ i ]], s.cur[ cls.of.group[ i ], ], fit$g, fit$c, fit$cg[ i ] ) ) )
    fit <- solve.bins( bins )
  }
  rows  <- fit$rows
  ok    <- fit$ok
  s.mat <- fit$s
  g     <- fit$g
  cq    <- fit$c
  cg    <- fit$cg

  flags <- matrix( 0L, k.n, d )
  upper <- matrix( Inf, k.n, d )
  for ( k in seq_len( k.n ) ) {
    ok.k <- ok & cls.of.group[ rows$group ] == k
    if ( !any( ok.k ) ) next
    upper[ k, ] <- .sigma0.upper.bound( rows, ok.k, d, s$bound.sds )
    high <- s.mat[ k, ] > upper[ k, ]
    s.mat[ k, high ] <- upper[ k, high ]
    flags[ k, high ] <- bitwOr( flags[ k, high ], .NOISE.FLAG.SIGMA0.HIGH )
    positive <- s.mat[ k, s.mat[ k, ] > 0 ]
    floor.k <- s$floor.frac * stats::median( if ( length( positive ) ) positive else rows$b[ ok.k ] )
    low <- s.mat[ k, ] < floor.k
    s.mat[ k, ] <- pmax( s.mat[ k, ], floor.k )
    flags[ k, low ] <- bitwOr( flags[ k, low ], .NOISE.FLAG.SIGMA0.LOW )
  }

  clamped <- matrix( bitwAnd( flags, .NOISE.FLAG.SIGMA0.HIGH + .NOISE.FLAG.SIGMA0.LOW ) > 0, k.n, d )
  if ( s$sigma0.fallback == "laser" && any( clamped ) ) {
    fb <- .laser.fallback( s.mat, clamped, upper, lasers )
    s.mat <- fb$s
    flags[ fb$applied ] <- bitwOr( flags[ fb$applied ], .NOISE.FLAG.SIGMA0.FALLBACK )
  }

  # refit g (and c, c_g) with sigma0 fixed: the rest of the problem with the
  # fixed part moved to the target, h_rest - G_rest,s s
  if ( any( clamped ) ) {
    sol.cur <- .noise.join( s.mat, g, cq, cg, lay )
    current <- .noise.row.predict( rows, projections, lay, sol.cur )
    w   <- ifelse( ok, rows$n / pmax( current, .NOISE.TINY )^2, 0 )
    ne  <- .noise.normal.equations( rows, projections, lay, w )
    s.idx <- unlist( lay$s.cols )
    rest  <- setdiff( seq_len( lay$n.cols ), s.idx )
    h.rest <- ne$h[ rest ] - as.vector( ne$G[ rest, s.idx, drop = FALSE ] %*% sol.cur[ s.idx ] )
    sub <- numeric( lay$n.cols )
    sub[ rest ] <- .nnls.normal( ne$G[ rest, rest, drop = FALSE ], h.rest )
    parts <- .noise.split( sub, lay, n.groups )
    g  <- parts$g
    cq <- parts$c
    cg <- parts$cg
  }

  g.pos <- g[ g > 0 ]
  if ( length( g.pos ) == 0L )
    stop( "No detector shows shot noise: the controls do not span enough intensity.", call. = FALSE )
  pooled <- g <= 1e-6 * stats::median( g.pos )
  g[ pooled ] <- stats::median( g.pos )
  flags[ , pooled ] <- bitwOr( flags[ , pooled ], .NOISE.FLAG.KAPPA.POOLED )

  cq <- pmax( cq, 0 )
  group.names <- vapply( groups, function( gr ) gr$name, character( 1 ) )
  group.c <- stats::setNames( cg[ controls ], group.names[ controls ] )

  sol.final <- .noise.join( s.mat, g, cq, cg, lay )
  fit.table <- data.frame(
    group     = group.names[ rows$group ],
    bin       = rows$bin,
    detector  = detectors[ rows$det ],
    mu        = rows$mu.own,
    observed  = rows$b,
    predicted = .noise.row.predict( rows, projections, lay, sol.final ),
    n         = as.integer( rows$n ),
    centre    = unlist( lapply( bins, function( bi ) rep( bi$centre, each = d ) ) ),
    stringsAsFactors = FALSE
  )

  flags.all <- apply( flags, 2, function( f ) Reduce( bitwOr, f, 0L ) )

  model <- list(
    detectors = detectors,
    sigma0.sq = stats::setNames( s.mat[ 1, ], detectors ),
    g         = stats::setNames( g, detectors ),
    c         = stats::setNames( cq, detectors ),
    flags     = stats::setNames( as.integer( flags.all ), detectors ),
    classes   = classes,
    sigma0.sq.by.class = NULL,
    flags.by.class     = NULL,
    group.c   = group.c,
    fit.table = fit.table,
    groups    = data.frame(
      name   = group.names,
      kind   = vapply( groups, function( gr ) gr$kind, character( 1 ) ),
      class  = group.cls,
      n.fit  = vapply( groups, function( gr ) nrow( gr$raw ), integer( 1 ) ),
      n.bins = vapply( bins, function( bi ) length( bi$n ), integer( 1 ) ),
      c.group = cg,
      stringsAsFactors = FALSE ),
    settings  = s
  )
  if ( k.n > 1L ) {
    model$sigma0.sq.by.class <- s.mat
    model$flags.by.class     <- flags
    dimnames( model$sigma0.sq.by.class ) <- list( classes, detectors )
    dimnames( model$flags.by.class )     <- list( classes, detectors )
  }
  model
}

# ---------------------------------------------------------------------------
# Model use
# ---------------------------------------------------------------------------

# sigma0^2 of a particle class; the first class when cls is NULL or not fitted.
.noise.sigma0.sq.for <- function( model, cls = NULL ) {
  if ( is.null( cls ) || is.null( model$sigma0.sq.by.class ) || !cls %in% model$classes )
    return( model$sigma0.sq )
  model$sigma0.sq.by.class[ match( cls, model$classes ), ]
}

# sigma0^2 m^2 + mu+ / kappa + c mu+^2 per event and detector.
.noise.variance <- function( model, mu, sigma0.mult = 1, cls = NULL ) {
  mu.pos <- pmax( mu, 0 )
  var <- sweep( mu.pos, 2, model$g, "*" )
  var <- sweep( var, 2, sigma0.mult^2 * .noise.sigma0.sq.for( model, cls ), "+" )
  if ( any( model$c != 0 ) ) var <- var + sweep( mu.pos * mu.pos, 2, model$c, "*" )
  pmax( var, .NOISE.TINY )
}

# A model with the same variance phi at every detector and signal.
.flat.noise.model <- function( detectors, phi ) {
  d <- length( detectors )
  list( detectors = detectors, sigma0.sq = stats::setNames( rep( phi, d ), detectors ),
        g = stats::setNames( numeric( d ), detectors ), c = stats::setNames( numeric( d ), detectors ),
        flags = stats::setNames( integer( d ), detectors ), classes = .NOISE.CLASSES[ 1 ],
        source = "flat" )
}

# Poisson weights with a floor: V_d = scale g_d (floor + mu+).
.poisson.noise.model <- function( detectors, g, floor, scale = 1 ) {
  d <- length( detectors )
  g <- as.numeric( g ) * scale
  list( detectors = detectors, sigma0.sq = stats::setNames( g * floor, detectors ),
        g = stats::setNames( g, detectors ), c = stats::setNames( numeric( d ), detectors ),
        flags = stats::setNames( integer( d ), detectors ), classes = .NOISE.CLASSES[ 1 ],
        source = "poisson" )
}

# Whether the negatives are clamped at zero: no detector has min.share or
# more of its values at or below zero (not clamped below min.events events).
.negatives.clamped <- function( negatives, min.share = 0.01, min.events = 200L ) {
  negatives <- Filter( function( x ) !is.null( x ) && nrow( x ) > 0L, negatives )
  if ( length( negatives ) == 0L ) return( list( clamped = FALSE, share = NULL ) )
  stacked <- do.call( rbind, negatives )
  share   <- colMeans( stacked <= 0 )
  if ( nrow( stacked ) < min.events ) return( list( clamped = FALSE, share = share ) )
  list( clamped = all( share < min.share ), share = share )
}
