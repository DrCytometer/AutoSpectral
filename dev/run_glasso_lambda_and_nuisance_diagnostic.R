# run_glasso_lambda_and_nuisance_diagnostic.R
#
# Two questions about correct.spectra.glasso() on concatenated single stained
# controls, answered without re-running the function.
#
# Stage 1 asks where on the lambda grid the two-way holdout rule actually
# lands for each target, and how much of the target-negative population's
# variance the best fit on that grid explains. A row that lands at grid
# point 1 has had every coefficient set to zero, which is not a statement
# that no spillover exists; it is the selection rule refusing to spend any
# error budget at all.
#
# Stage 2 takes one fluorophore's own population and refits its signature
# three ways - with the active set the lasso selected, with the whole panel,
# and with the target alone - so the effect of the nuisance set on the
# recovered row can be read off directly against ground truth.
#
# Requires devtools::load_all() so the non-exported .glasso.* helpers are
# callable. Expects the objects a correct.spectra.glasso() run already
# builds: `g.fit` (its return value), `spectra.start` (what it was given),
# `spectra.truth` (the reference the corrected rows should move toward, or
# NULL), plus `flow.control`, `asp`, and the two file paths.

library( AutoSpectral )

# ---------------------------------------------------------------------------
# Inputs - edit these
# ---------------------------------------------------------------------------

g.fit         <- glasso.fit
spectra.start <- cell.spectra
spectra.truth <- bd.spectra          # or NULL
fs.path       <- fully.stained.path
un.path       <- unstained.path
af.name       <- "AF"
bg.mode       <- "af.deconv"
diag.dir      <- file.path( asp$fix.unmixing.dir, "glasso_diagnostic" )

if ( ! dir.exists( diag.dir ) ) dir.create( diag.dir, recursive = TRUE )

fluorophores <- setdiff( rownames( spectra.start ), af.name )

# ---------------------------------------------------------------------------
# Stage 0: what the run already tells us, no recomputation
# ---------------------------------------------------------------------------

cat( "\n=== Stage 0: selection summary from the completed run ===\n" )

cat( sprintf( "coefficients selected: %d of %d pairs\n",
              sum( g.fit$coefficient.log$selected ),
              nrow( g.fit$coefficient.log ) ) )

cat( sprintf( "targets with an all-zero row: %d of %d\n",
              sum( g.fit$lambda.log$n.selected == 0 ),
              nrow( g.fit$lambda.log ) ) )

cat( sprintf( "active sets of size one (target only): %d of %d\n",
              sum( lengths( g.fit$active.set ) == 1 ),
              length( g.fit$active.set ) ) )

cat( sprintf( "reversed donors found anywhere: %d\n",
              sum( lengths( g.fit$reversed.donors ) ) ) )

cat( sprintf( "off-diagonal spillover, max |value|: %.5f\n",
              max( abs( g.fit$spillover[ row( g.fit$spillover ) !=
                                           col( g.fit$spillover ) ] ) ) ) )

print( g.fit$convergence.log )

if ( ! is.null( g.fit$signature.log ) )
  print( table( g.fit$signature.log$reason ) )

# ---------------------------------------------------------------------------
# Stage 1: where the lambda rule lands, per target
# ---------------------------------------------------------------------------

cat( "\n=== Stage 1: lambda grid position and explained variance ===\n" )

read.gated <- function( file.name, gate.polygon = NULL, label ) {

  expr.data <- readFCS( file.name,
                        columns = flow.control$scatter.and.channel.spectral )
  gate.data <- expr.data[ , flow.control$scatter.parameter ]

  if ( is.null( gate.polygon ) )
    gate.polygon <- do.gate(
      gate.data, viability.gate = FALSE, large.gate = TRUE, samp = label,
      scatter.and.channel.label = flow.control$scatter.and.channel.label,
      control.type = "cells", asp )

  keep <- which( sp::point.in.polygon( gate.data[ , 1 ], gate.data[ , 2 ],
                                       gate.polygon$x, gate.polygon$y ) != 0 )

  list( data = expr.data[ keep, flow.control$spectral.channel, drop = FALSE ],
        gate = gate.polygon )
}

un.in <- read.gated( un.path, NULL, "unstained raw" )
fs.in <- read.gated( fs.path, un.in$gate, "fully stained raw" )

unstained.raw <- un.in$data
stained.raw   <- fs.in$data

spectra.fluor <- spectra.start[ fluorophores, , drop = FALSE ]

background.basis <- if ( identical( bg.mode, "af.deconv" ) )
  as.matrix( get.af.basis( unstained.raw, n.pc = "auto", verbose = FALSE ) ) else
    if ( identical( bg.mode, "af.row" ) )
      spectra.start[ af.name, , drop = FALSE ] else NULL

design <- if ( is.null( background.basis ) ) spectra.fluor else
  rbind( background.basis, spectra.fluor )

project <- function( y ) {
  coefs <- unmix.ols.fast( y, design )
  colnames( coefs ) <- rownames( design )
  coefs[ , fluorophores, drop = FALSE ]
}

stained.abundance   <- project( stained.raw )
unstained.abundance <- project( unstained.raw )

thresholds <- 1.3 * apply( unstained.abundance, 2, stats::quantile,
                           probs = 0.99, names = FALSE )
names( thresholds ) <- fluorophores

lambda.rows <- list()

for ( target in fluorophores ) {

  sources <- setdiff( fluorophores, target )
  neg.idx <- which( stained.abundance[ , target ] <= thresholds[ target ] )

  if ( length( neg.idx ) < 200 ) next
  if ( length( neg.idx ) > 20000 ) neg.idx <- sample( neg.idx, 20000 )

  x.mat <- stained.abundance[ neg.idx, sources, drop = FALSE ]
  y.vec <- stained.abundance[ neg.idx, target ]

  sel <- .glasso.select.lambda( x.mat, y.vec, n.lambda = 40L,
                                lambda.min.ratio = 1e-3, margin.frac = 0.05 )

  total <- sel$path$total
  best  <- which.min( total )

  lambda.rows[[ length( lambda.rows ) + 1L ]] <- data.frame(
    target            = target,
    n.negative        = length( neg.idx ),
    grid.index.chosen = which( sel$path$lambda == sel$lambda )[ 1 ],
    grid.index.best   = best,
    r.squared.best    = 1 - total[ best ] / total[ 1 ],
    n.selected        = sum( sel$beta != 0 ),
    max.abs.beta      = max( abs( sel$beta ) ),
    n.selected.at.min = sum( abs( sel$beta ) > 0 ),
    row.names         = NULL, stringsAsFactors = FALSE )
}

lambda.table <- do.call( rbind, lambda.rows )
print( lambda.table )

utils::write.csv( lambda.table,
                  file.path( diag.dir, "glasso_lambda_grid_position.csv" ),
                  row.names = FALSE )

cat( sprintf(
  "\ntargets landing at grid index 1 (all coefficients zeroed): %d of %d\n",
  sum( lambda.table$grid.index.chosen == 1 ), nrow( lambda.table ) ) )
cat( sprintf( "median R-squared at the grid minimum: %.4f\n",
              stats::median( lambda.table$r.squared.best ) ) )

# ---------------------------------------------------------------------------
# Stage 2: the nuisance set's effect on one recovered signature
# ---------------------------------------------------------------------------

cat( "\n=== Stage 2: same population, three nuisance sets ===\n" )

compensation <- solve( g.fit$spillover )
unmixed.comp <- stained.abundance %*% compensation

background <- if ( is.null( background.basis ) ) 0 else {
  coefs <- unmix.ols.fast( stained.raw, design )
  colnames( coefs ) <- rownames( design )
  b <- coefs[ , rownames( background.basis ), drop = FALSE ]
  b[ b < 0 ] <- 0
  b %*% background.basis
}

residual <- stained.raw - background

thr.comp <- 1.3 * apply( unstained.abundance %*% compensation, 2,
                         stats::quantile, probs = 0.99, names = FALSE )
names( thr.comp ) <- fluorophores

dyn.range <- pmax( apply( unmixed.comp, 2, stats::quantile, probs = 0.999,
                          names = FALSE ) - thr.comp, .Machine$double.eps )

score    <- sweep( pmax( sweep( unmixed.comp, 2, thr.comp, "-" ), 0 ), 2,
                   dyn.range, "/" )
dominant <- max.col( score, ties.method = "first" )
dominant[ score[ cbind( seq_len( nrow( score ) ), dominant ) ] <= 0 ] <- 0L

background.idx <- which( dominant == 0L )

angle.to <- function( a, b ) {
  cs <- sum( a * b ) / max( sqrt( sum( a^2 ) ) * sqrt( sum( b^2 ) ),
                            .Machine$double.eps )
  180 / pi * acos( pmin( 1, pmax( -1, cs ) ) )
}

nuisance.rows <- list()

for ( f in seq_along( fluorophores ) ) {

  j   <- fluorophores[ f ]
  idx <- which( dominant == f )

  if ( length( idx ) < 200 ) next

  sets <- list(
    selected    = g.fit$active.set[[ j ]],
    panel       = fluorophores,
    target.only = j )

  for ( set.name in names( sets ) ) {

    candidate <- tryCatch( extract.raw.signature(
      raw.data       = residual[ idx, , drop = FALSE ],
      spectra        = spectra.start[ fluorophores, , drop = FALSE ],
      abundance      = unmixed.comp[ idx, , drop = FALSE ],
      target         = j,
      active         = sets[[ set.name ]],
      intercept      = TRUE,
      multivariate   = TRUE,
      ridge          = 1e-6,
      n.levels       = 60L,
      min.bin.events = 50L,
      min.events     = 200L,
      background.raw = if ( length( background.idx ) > 0 )
        residual[ background.idx, , drop = FALSE ] else NULL ),
      error = function( e ) NULL )

    if ( is.null( candidate ) ) next

    st <- candidate$stats

    nuisance.rows[[ length( nuisance.rows ) + 1L ]] <- data.frame(
      fluorophore  = j,
      nuisance.set = set.name,
      n.events     = length( idx ),
      n.active     = st$n.active,
      joint        = st$joint,
      n.bins       = st$n.bins,
      resid.rel    = st$resid.rel,
      explained.total = st$explained.total,
      vif.target   = st$vif.target,
      deg.change   = st$deg.change,
      deg.to.truth = if ( is.null( spectra.truth ) ||
                          ! j %in% rownames( spectra.truth ) ) NA_real_ else
                            angle.to( candidate$signature,
                                      spectra.truth[ j, colnames( spectra.start ) ] ),
      deg.start    = if ( is.null( spectra.truth ) ||
                          ! j %in% rownames( spectra.truth ) ) NA_real_ else
                            angle.to( spectra.start[ j, ],
                                      spectra.truth[ j, colnames( spectra.start ) ] ),
      row.names    = NULL, stringsAsFactors = FALSE )
  }
}

nuisance.table <- do.call( rbind, nuisance.rows )
print( nuisance.table )

utils::write.csv( nuisance.table,
                  file.path( diag.dir, "glasso_nuisance_set_comparison.csv" ),
                  row.names = FALSE )

if ( ! all( is.na( nuisance.table$deg.to.truth ) ) )
  print( stats::aggregate( cbind( deg.to.truth, deg.start ) ~ nuisance.set,
                           data = nuisance.table, FUN = median ) )
