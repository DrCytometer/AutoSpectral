# noise_model_calibrate.R

# Predictive check of a fitted noise model on held-out events, each set
# projected with the same design as its fitting group: per detector and
# intensity tranche, the squared normalised MAD of r_d / sqrt(sum_e M_de V_e)
# (1 when the model is right), and per set the share of events whose GLS
# chi-square has p below p.cut. Each set uses its class's sigma0 and, for a
# control, its own c_g.

# V(mu) of one set: the model's variance plus the set's own c_g mu+^2.
.noise.set.variance <- function( model, mu, sigma0.mult = 1, cls = NULL, c.group = 0 ) {
  var <- .noise.variance( model, mu, sigma0.mult, cls )
  if ( c.group != 0 ) var <- var + c.group * pmax( mu, 0 )^2
  var
}

# sum_e M_de V_e(mu_e) per event and detector.
.noise.predicted.resid.var <- function( model, proj, sigma0.mult = 1, cls = NULL, c.group = 0 ) {
  .noise.set.variance( model, proj$fit, sigma0.mult, cls, c.group ) %*% t( proj$M )
}

# Per-event chi-square of the GLS fit on the set's design, with the variance
# at the OLS fitted signal: list(chi2, df). Mean-only sets use the residual
# from the set mean.
.noise.gls.chi2 <- function( model, proj, sigma0.mult = 1, cls = NULL, c.group = 0 ) {
  v <- .noise.set.variance( model, proj$fit, sigma0.mult, cls, c.group )
  d <- ncol( v )
  if ( is.null( proj$q ) )
    return( list( chi2 = rowSums( proj$resid^2 / v ), df = d ) )
  q <- proj$q
  k <- nrow( q )
  w <- 1 / v
  y <- proj$centred
  chi <- numeric( nrow( v ) )
  for ( i in seq_len( nrow( v ) ) ) {
    qw    <- sweep( q, 2, w[ i, ], "*" )
    theta <- solve( tcrossprod( qw, q ), as.vector( qw %*% y[ i, ] ) )
    r     <- y[ i, ] - as.vector( theta %*% q )
    chi[ i ] <- sum( r * r * w[ i, ] )
  }
  list( chi2 = chi, df = d - k )
}

# (D x T) squared normalised MAD of z in quantile tranches of each detector's
# mu, with counts and median levels.
.noise.tranche.ratios <- function( z, mu, n.tranches, min.events ) {
  d <- ncol( z )
  ratio <- level <- matrix( NA_real_, d, n.tranches )
  count <- matrix( 0L, d, n.tranches )
  for ( j in seq_len( d ) ) {
    edges <- stats::quantile( mu[ , j ], seq( 0, 1, length.out = n.tranches + 1L ), names = FALSE )
    which.t <- pmin( pmax( findInterval( mu[ , j ], edges ) - 1L, 0L ), n.tranches - 1L )
    for ( t.i in seq_len( n.tranches ) ) {
      sel <- which.t == t.i - 1L
      count[ j, t.i ] <- sum( sel )
      if ( count[ j, t.i ] >= min.events ) {
        ratio[ j, t.i ] <- ( 1.4826 * stats::median( abs( z[ sel, j ] ) ) )^2
        level[ j, t.i ] <- stats::median( mu[ sel, j ] )
      }
    }
  }
  list( ratio = ratio, count = count, level = level )
}

#' @title Calibrate a Noise Model on Held-Out Events
#'
#' @description
#' Observed / predicted residual variance on held-out sets (see the header of
#' `noise_model_calibrate.R`). `ratio` (D x T) pools every set and
#' `ratio.controls` the single-stain controls only, each over quantile
#' tranches of the detector's own fitted signal; z-scores are centred per set
#' and detector first. `in.band` and `in.band.controls` are the shares of
#' filled cells inside `band`.
#'
#' @param model Noise model list from `.fit.noise.model()`.
#' @param groups List of held-out groups (same designs as the fit).
#' @param n.tranches Integer. Default `10`.
#' @param band Numeric length 2. Default `c(0.8, 1.25)`.
#' @param min.events Integer, fewest events in a tranche. Default `30`.
#' @param sigma0.mult Numeric multiplier on sigma0. Default `1`.
#' @param p.cut Numeric chi-square p-value cut. Default `1e-3`.
#'
#' @return A list: `groups` (data frame, one row per set: events, median
#'   per-detector ratio, share of detectors in band, share of events with
#'   chi-square p below `p.cut`, median chi-square / df), `per.detector`
#'   (sets x detectors matrix of ratios), `ratio`, `count`, `level`,
#'   `in.band` and the same with `.controls`, `band`, `n.tranches`,
#'   `sigma0.mult`, `per.control.c`.
#'
#' @noRd

.calibrate.noise.model <- function( model, groups, n.tranches = 10L, band = c( 0.8, 1.25 ),
                                    min.events = 30L, sigma0.mult = 1, p.cut = 1e-3 ) {

  zs <- mus <- kinds <- list()
  rows <- list()
  per.det <- list()
  for ( gr in groups ) {
    proj <- .noise.project( gr )
    cls  <- if ( is.null( gr$cls ) ) "cells" else gr$cls
    c.group <- if ( gr$kind == "control" && gr$name %in% names( model$group.c ) )
      model$group.c[[ gr$name ]] else 0
    p <- .noise.predicted.resid.var( model, proj, sigma0.mult, cls, c.group )
    z <- proj$resid / sqrt( p )
    z <- sweep( z, 2, .col.medians( z ), "-" )
    pd <- ( 1.4826 * .col.medians( abs( z ) ) )^2
    chi <- .noise.gls.chi2( model, proj, sigma0.mult, cls, c.group )
    df  <- max( as.integer( chi$df ), 1L )
    zs[[ length( zs ) + 1L ]]   <- z
    mus[[ length( mus ) + 1L ]] <- pmax( proj$fit, 0 )
    kinds[[ length( kinds ) + 1L ]] <- rep( gr$kind, nrow( z ) )
    per.det[[ gr$name ]] <- pd
    rows[[ length( rows ) + 1L ]] <- data.frame(
      name = gr$name, kind = gr$kind, n = nrow( z ),
      median.ratio = stats::median( pd ),
      in.band = mean( pd >= band[ 1 ] & pd <= band[ 2 ] ),
      chi2.p.small = mean( stats::pchisq( chi$chi2, df, lower.tail = FALSE ) < p.cut ),
      chi2.per.df = stats::median( chi$chi2 ) / df,
      stringsAsFactors = FALSE )
  }

  out <- list( groups = do.call( rbind, rows ),
               per.detector = do.call( rbind, per.det ),
               band = band, n.tranches = n.tranches, sigma0.mult = sigma0.mult,
               per.control.c = length( model$group.c ) > 0L )
  if ( length( zs ) == 0L ) return( out )
  if ( !is.null( out$per.detector ) ) colnames( out$per.detector ) <- model$detectors

  z.all    <- do.call( rbind, zs )
  mu.all   <- do.call( rbind, mus )
  kind.all <- unlist( kinds )
  for ( suffix in c( "", ".controls" ) ) {
    sel <- if ( suffix == "" ) rep( TRUE, nrow( z.all ) ) else kind.all == "control"
    if ( sum( sel ) < n.tranches * min.events ) {
      out[[ paste0( "ratio", suffix ) ]]   <- NULL
      out[[ paste0( "in.band", suffix ) ]] <- NA_real_
      next
    }
    tr <- .noise.tranche.ratios( z.all[ sel, , drop = FALSE ], mu.all[ sel, , drop = FALSE ],
                                 n.tranches, min.events )
    rownames( tr$ratio ) <- rownames( tr$count ) <- rownames( tr$level ) <- model$detectors
    filled <- is.finite( tr$ratio )
    out[[ paste0( "ratio", suffix ) ]] <- tr$ratio
    out[[ paste0( "count", suffix ) ]] <- tr$count
    out[[ paste0( "level", suffix ) ]] <- tr$level
    out[[ paste0( "in.band", suffix ) ]] <- if ( any( filled ) )
      mean( tr$ratio[ filled ] >= band[ 1 ] & tr$ratio[ filled ] <= band[ 2 ] ) else NA_real_
  }
  out
}

# Observed / predicted residual variance of a projected group, median over detectors.
.noise.dispersion.scale <- function( model, proj ) {
  observed  <- .robust.var( proj$resid )
  predicted <- .col.medians( .noise.predicted.resid.var( model, proj ) )
  stats::median( observed / pmax( predicted, .NOISE.TINY ) )
}
