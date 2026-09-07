# run_signature_correction_diagnostics.R
#
# Staged diagnostics for the signature-error-correction workflow, run against
# the bead / cell single-stained control pair.
#
# Stages A-D are standalone: they define their own helpers and do not depend on
# any pending change to the package, so they can be run before editing anything.
# Stage E exercises the package loop and should only be run once the changes to
# correct.unmixing.signatures(), cluster.unmixed.events(),
# restrict.cluster.unmixing() and fit.signature.error.model() are in.
#
# Run top to bottom. Every stage prints a table; those tables are the output to
# report back.

# ---------------------------------------------------------------------------
# 0. Setup - EDIT THIS SECTION
# ---------------------------------------------------------------------------

asp.dir <- "/Users/oliverburton/Bioinformatics/AutoSpectral"
asp.rcpp.dir <- "/Users/oliverburton/Bioinformatics/AutoSpectralRcpp"
devtools::load_all(asp.rcpp.dir)
devtools::load_all(asp.dir)
asp <- get.autospectral.param( cytometer = "aurora" )
cell.dir <- "./PBMC direct stain"
cell.file <-"./PBMC direct stain/cell_control_file_corrections.csv"
bd.dir <- "./BD beads"
bd.file <- "./BD beads/bead_control_file_corrections.csv"
cell.spectra <- read.spectra("Cells_spectra_corrections.csv")
bd.spectra <- read.spectra("BD_spectra_corrections.csv")
concat.cells <- readFCS("./concatenated_fcs/Concatenated_cells.fcs")
concat.beads <- readFCS("./concatenated_fcs/Concatenated_BD_beads.fcs")


# Independently-measured spectra per particle type, fluorophores x detectors,
# INCLUDING the autofluorescence row(s) if present.
diag.spectra <- list(
  Cells = cell.spectra,
  Beads = bd.spectra
)

# Concatenated single-stained raw data per particle type.
diag.raw <- list(
  Cells = concat.cells,
  Beads = concat.beads
)

# Unstained control per particle type, raw. Used for positivity thresholds and
# for the negative-population spread reference.
diag.unstained <- list(
  Cells = readFCS( file.path( cell.dir, "A10 unstained_010_Cells.fcs" ) ),
  Beads = readFCS( file.path( bd.dir, "A10 unstained_010_Beads_1.fcs" ) )
)

# Name of the autofluorescence row in `diag.spectra`, or NULL if the tables are
# fluorophore-only. It is kept in the design at all times rather than being
# thresholded out.
af.name <- "AF"

# Positivity percentile on the unstained control, and how many events to use.
diag.threshold.probs <- 0.995

# Variant lists from get.spectral.variants(), one per particle type. Only
# $spillover.spread is used.
diag.variants <- list(
  Cells = readRDS( "./figure_spectral_variants/Spectral_variants_cells.rds" ),
  Beads = readRDS( "./figure_spectral_variants/Spectral_variants_beads.rds" )
)
diag.max.events      <- 200000L
diag.cofactor        <- 500
diag.seed            <- 42L

write.csv(diag.variants$Cells[["spillover.spread"]],
          file = file.path(cell.dir, "cell_spillover_spread.csv"))
write.csv(diag.variants$Beads[["spillover.spread"]],
          file = file.path(bd.dir, "bead_spillover_spread.csv"))
# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

#' Row-wise cosine similarity of two matched spectra matrices.
#' @noRd
.diag.row.cosine <- function( a, b ) {
  stopifnot( identical( dim( a ), dim( b ) ) )
  num <- rowSums( a * b )
  den <- sqrt( rowSums( a^2 ) ) * sqrt( rowSums( b^2 ) )
  num / den
}

#' Scale-matched residual angle in degrees, and the fraction of the initial
#' error removed. Cosine similarity saturates above 0.99 on real spectra and
#' cannot resolve differences at the scale that matters here; the angle is the
#' same quantity on a readable scale, and recovered.fraction is linear in it
#' (0 = no change, 1 = exact recovery, negative = moved away).
#' @noRd
.diag.angle <- function( a, b ) {
  cs <- rowSums( a * b ) / ( sqrt( rowSums( a^2 ) ) * sqrt( rowSums( b^2 ) ) )
  180 / pi * acos( pmin( 1, pmax( -1, cs ) ) )
}

#' @noRd
.diag.recovered.fraction <- function( start, corrected, truth ) {
  ang.start <- .diag.angle( start, truth )
  ang.end   <- .diag.angle( corrected, truth )
  ( ang.start - ang.end ) / ang.start
}

#' L-infinity renormalisation, matching the package convention.
#' @noRd
.diag.renorm <- function( s ) {
  row.max <- apply( s, 1, max )
  row.max[ row.max <= 0 ] <- 1
  s / row.max
}

#' Split a spectral error into components parallel and perpendicular to the
#' row space of a reference spectra matrix.
#' @noRd
.diag.error.split <- function( error, reference ) {
  ref  <- as.matrix( reference )
  proj <- t( ref ) %*% MASS::ginv( ref %*% t( ref ) ) %*% ref
  par  <- error %*% proj
  list( parallel = par, perpendicular = error - par )
}

#' Per-cluster centroids in linear units, given a fixed cluster assignment.
#' @noRd
.diag.centroids <- function( cluster.id, unmixed, raw.data, min.size = 5 ) {

  tab  <- table( cluster.id )
  keep <- as.integer( names( tab )[ tab >= min.size ] )

  x.clust <- t( vapply( keep, function( cl )
    colMeans( unmixed[ which( cluster.id == cl ), , drop = FALSE ] ),
    numeric( ncol( unmixed ) ) ) )

  y.clust <- t( vapply( keep, function( cl )
    colMeans( raw.data[ which( cluster.id == cl ), , drop = FALSE ] ),
    numeric( ncol( raw.data ) ) ) )

  colnames( x.clust ) <- colnames( unmixed )
  colnames( y.clust ) <- colnames( raw.data )
  rownames( x.clust ) <- rownames( y.clust ) <- keep

  list( x.clust = x.clust, y.clust = y.clust,
        size = as.integer( tab[ as.character( keep ) ] ) )
}

#' Restricted unmix of a matrix of centroids against a chosen fluorophore set.
#' @noRd
.diag.restricted.unmix <- function( y, spectra, active ) {
  sub <- spectra[ active, , drop = FALSE ]
  y %*% t( solve.default( tcrossprod( sub ), sub ) )
}

#' Autofluorescence principal components from an unstained control, matching
#' the basis get.spectral.variants() builds: top right singular vectors of the
#' raw (uncentred) unstained matrix, so the leading component is the mean
#' background direction rather than a deviation from it.
#' @noRd
.diag.af.pcs <- function( unstained, n.pc = 4L, max.events = 20000L ) {
  
  dat <- as.matrix( unstained )
  if ( nrow( dat ) > max.events )
    dat <- dat[ sample.int( nrow( dat ), max.events ), , drop = FALSE ]
  
  sv  <- svd( dat, nu = 0, nv = n.pc )
  pcs <- t( sv$v )
  dimnames( pcs ) <- list( paste0( "AFPC", seq_len( nrow( pcs ) ) ), colnames( dat ) )
  pcs
}

#' Scatter-matched per-event background subtraction, mirroring the kNN stage of
#' get.fluor.variants(). Each event has the mean spectrum of its k nearest
#' unstained neighbours in scatter space subtracted. Matching on FSC/SSC rather
#' than on the spectral profile is the point: the amount subtracted is set by an
#' independent measurement channel, so it cannot be absorbed back into the
#' fluorophore abundance the way a spectrally-derived background estimate can.
#' @noRd
.diag.bg.subtract <- function( raw.data, scatter, unstained, unstained.scatter,
                               k.neighbors = 20L, max.reference = 50000L ) {
  
  ref.idx <- seq_len( nrow( unstained ) )
  if ( length( ref.idx ) > max.reference )
    ref.idx <- sample( ref.idx, max.reference )
  
  knn.idx <- FNN::knnx.index(
    data  = as.matrix( unstained.scatter[ ref.idx, , drop = FALSE ] ),
    query = as.matrix( scatter ),
    k     = k.neighbors
  )
  
  bg <- matrix( 0, nrow( raw.data ), ncol( raw.data ),
                dimnames = dimnames( raw.data ) )
  
  for ( ki in seq_len( k.neighbors ) )
    bg <- bg + as.matrix( unstained[ ref.idx, , drop = FALSE ] )[ knn.idx[ , ki ], , drop = FALSE ]
  
  raw.data - bg / k.neighbors
}

#' Hotspot scale (sqrt of the pseudoinverse of the cosine-similarity matrix) for
#' a candidate restricted basis, reported for one target row. Used to decide
#' per dye whether AF components can safely be added to that dye's projection
#' basis: a dye whose spectrum is close to the AF span will have real signal
#' removed by the projection rather than background.
#' @noRd
.diag.basis.hotspot <- function( basis, target ) {
  h <- sqrt( abs( MASS::ginv( cosine.similarity( basis ) ) ) )
  dimnames( h ) <- list( rownames( basis ), rownames( basis ) )
  max( h[ target, ] )
}

#' Main-population scatter gate: keep events whose 2D scatter density exceeds
#' a fraction of the modal density. On bead files this removes the noise and
#' debris fraction, whose distinct autofluorescence profile otherwise enters
#' the dominance populations, the global-mean background estimate, and the
#' evaluation medians; on cells it removes debris and most aggregates.
#' Density-relative rather than polygonal, so it needs no per-run landmarks.
#' @noRd
.diag.main.gate <- function( scatter, gate.level = 0.1, grid.n = 128L,
                             max.events = 100000L ) {
  
  sc <- as.matrix( scatter[ , 1:2, drop = FALSE ] )
  
  fit.idx <- seq_len( nrow( sc ) )
  if ( length( fit.idx ) > max.events )
    fit.idx <- sample( fit.idx, max.events )
  
  bw <- vapply( 1:2, function( k ) {
    b <- tryCatch( MASS::bandwidth.nrd( sc[ fit.idx, k ] ),
                   error = function( e ) 0 )
    if ( !is.finite( b ) || b <= 0 )
      b <- 4 * 1.06 * stats::sd( sc[ fit.idx, k ] ) *
        length( fit.idx )^( -0.2 )
    max( b, .Machine$double.eps )
  }, numeric( 1 ) )
  
  kde <- MASS::kde2d( sc[ fit.idx, 1 ], sc[ fit.idx, 2 ],
                      h = bw, n = grid.n )
  
  ix   <- findInterval( sc[ , 1 ], kde$x, all.inside = TRUE )
  iy   <- findInterval( sc[ , 2 ], kde$y, all.inside = TRUE )
  dens <- kde$z[ cbind( ix, iy ) ]
  
  dens >= gate.level * max( kde$z )
}

# ---------------------------------------------------------------------------
# STAGE A - algebra sanity check on noise-free synthetic data
# ---------------------------------------------------------------------------
# Confirms that regressing the cluster residual on cluster abundance estimates
# the spectral error itself, so that spectra + beta is the correct update. No
# clustering, no noise, no autofluorescence: if this stage does not behave as
# stated, the problem is in the linear algebra rather than in the data.

cat( "\n===================== STAGE A: algebra sanity =====================\n" )

set.seed( diag.seed )

a.f <- 6L
a.d <- 20L
a.c <- 400L

a.true <- matrix( abs( stats::rnorm( a.f * a.d ) ), nrow = a.f )
a.true <- .diag.renorm( a.true )
dimnames( a.true ) <- list( paste0( "F", seq_len( a.f ) ), paste0( "D", seq_len( a.d ) ) )

a.wrong <- pmax( a.true - 0.03 * matrix( stats::rnorm( a.f * a.d ), nrow = a.f ), 0 )
dimnames( a.wrong ) <- dimnames( a.true )
a.error <- a.true - a.wrong

a.x <- matrix( abs( stats::rnorm( a.c * a.f ) ) * 1000, nrow = a.c )
colnames( a.x ) <- rownames( a.true )
a.y <- a.x %*% a.true

a.xhat  <- unmix.ols.fast( a.y, a.wrong )
a.resid <- a.y - a.xhat %*% a.wrong

cat( sprintf( "  max |r %%*%% t(spectra)| (expected ~0): %.3e\n",
              max( abs( a.resid %*% t( a.wrong ) ) ) ) )

a.fit  <- stats::lm.fit( x = cbind( 1, a.xhat ), y = a.resid )
a.beta <- stats::coef( a.fit )[ -1, , drop = FALSE ]
dimnames( a.beta ) <- dimnames( a.wrong )

a.split <- .diag.error.split( a.error, a.wrong )

a.table <- data.frame(
  fluorophore   = rownames( a.true ),
  cos.start     = .diag.row.cosine( a.wrong, a.true ),
  cos.increase  = .diag.row.cosine( .diag.renorm( pmax( a.wrong + a.beta, 0 ) ), a.true ),
  cos.decrease  = .diag.row.cosine( .diag.renorm( pmax( a.wrong - a.beta, 0 ) ), a.true ),
  cos.ceiling   = .diag.row.cosine(
    .diag.renorm( pmax( a.wrong + a.split$perpendicular, 0 ) ), a.true ),
  row.names     = NULL
)

cat( sprintf( "  ||E|| = %.4f   ||E.parallel|| = %.4f   ||E.perpendicular|| = %.4f\n",
              norm( a.error, "F" ), norm( a.split$parallel, "F" ),
              norm( a.split$perpendicular, "F" ) ) )
cat( sprintf( "  cor( beta, E.perpendicular ) = %.4f\n",
              stats::cor( as.vector( a.beta ), as.vector( a.split$perpendicular ) ) ) )
cat( "\n" )
print( a.table, digits = 6 )
cat( "\n  EXPECTED: cos.increase > cos.start on every row, cos.decrease < cos.start,\n" )
cat( "  and cos.increase close to cos.ceiling.\n" )

# ---------------------------------------------------------------------------
# STAGE B - identifiability ceiling on the real bead/cell pair
# ---------------------------------------------------------------------------
# How much of the real bead-versus-cell spectral difference lies inside the row
# space of the starting spectra, and is therefore invisible to any residual-based
# correction that unmixes against the full panel. cos.ceiling is the best cosine
# similarity a single full-panel correction step could reach.

cat( "\n============ STAGE B: identifiability ceiling, real data ============\n" )

b.shared.fluor <- intersect( rownames( diag.spectra$Cells ), rownames( diag.spectra$Beads ) )
b.shared.det   <- intersect( colnames( diag.spectra$Cells ), colnames( diag.spectra$Beads ) )

cat( sprintf( "  Shared rows: %d (%s)\n", length( b.shared.fluor ),
              paste( b.shared.fluor, collapse = ", " ) ) )
cat( sprintf( "  Shared detectors: %d\n\n", length( b.shared.det ) ) )

b.spectra <- lapply( diag.spectra, function( s )
  as.matrix( s[ b.shared.fluor, b.shared.det, drop = FALSE ] ) )

b.report <- lapply( c( "Cells", "Beads" ), function( target ) {

  wrong <- if ( target == "Cells" ) "Beads" else "Cells"

  s.start <- b.spectra[[ wrong ]]
  s.true  <- b.spectra[[ target ]]
  e.mat   <- s.true - s.start
  split   <- .diag.error.split( e.mat, s.start )

  out <- data.frame(
    target        = target,
    fluorophore   = b.shared.fluor,
    cos.start     = .diag.row.cosine( s.start, s.true ),
    err.norm      = sqrt( rowSums( e.mat^2 ) ),
    frac.parallel = sqrt( rowSums( split$parallel^2 ) ) / sqrt( rowSums( e.mat^2 ) ),
    cos.ceiling   = .diag.row.cosine(
      .diag.renorm( pmax( s.start + split$perpendicular, 0 ) ), s.true ),
    row.names     = NULL
  )
  out$headroom <- out$cos.ceiling - out$cos.start

  cat( sprintf( "-- %s raw data, starting from %s spectra --\n", target, wrong ) )
  print( out, digits = 5 )
  cat( sprintf( "  panel totals: ||E|| = %.4f, parallel fraction = %.3f\n\n",
                norm( e.mat, "F" ),
                norm( split$parallel, "F" ) / norm( e.mat, "F" ) ) )
  out
} )
names( b.report ) <- c( "Cells", "Beads" )

cat( "  INTERPRETATION: frac.parallel is the share of the error that a full-panel\n" )
cat( "  residual correction cannot see. headroom is all the full-panel path can win.\n" )
cat( "  The restricted / single-positive path is not bound by this ceiling.\n" )

# ---------------------------------------------------------------------------
# STAGE C - cluster structure and threshold audit
# ---------------------------------------------------------------------------
# Two questions. First, does the clustering give every fluorophore a usable
# number of populated abundance levels? Density-driven clustering (SOM) spends
# its nodes where the events are, which in pooled single-stained controls is the
# shared negative population, so a rare positive population can end up with one
# node or none. Stratification assigns each event to the fluorophore it is most
# strongly positive for relative to that fluorophore's own threshold, then bins
# by abundance, guaranteeing equal representation. Second, does a flat
# positivity threshold survive spillover spread, or does it misclassify bright
# events as double-positive?

cat( "\n============ STAGE C: cluster and threshold audit ============\n" )

c.results <- list()

for ( target in c( "Cells", "Beads" ) ) {

  wrong   <- if ( target == "Cells" ) "Beads" else "Cells"
  s.start <- as.matrix( diag.spectra[[ wrong ]][ , b.shared.det, drop = FALSE ] )
  fluors  <- rownames( s.start )
  panel   <- setdiff( fluors, af.name )

  scatter.det <- grep( "^(FSC|SSC)", colnames( diag.raw[[ target ]] ), value = TRUE )
  
  set.seed( diag.seed )
  keep.idx <- if ( nrow( diag.raw[[ target ]] ) > diag.max.events )
    sample.int( nrow( diag.raw[[ target ]] ), diag.max.events ) else
      seq_len( nrow( diag.raw[[ target ]] ) )
  
  raw.all     <- diag.raw[[ target ]][ keep.idx, b.shared.det, drop = FALSE ]
  scatter.all <- diag.raw[[ target ]][ keep.idx, scatter.det, drop = FALSE ]
  
  gate.keep   <- .diag.main.gate( scatter.all )
  cat( sprintf( "  main-population gate, %s: kept %d of %d events (%.1f%%)\n",
                target, sum( gate.keep ), length( gate.keep ),
                100 * mean( gate.keep ) ) )
  raw.all     <- raw.all[ gate.keep, , drop = FALSE ]
  scatter.all <- scatter.all[ gate.keep, , drop = FALSE ]
  
  unst         <- diag.unstained[[ target ]][ , b.shared.det, drop = FALSE ]
  unst.scatter <- diag.unstained[[ target ]][ , scatter.det, drop = FALSE ]
  
  unst.keep    <- .diag.main.gate( unst.scatter )
  cat( sprintf( "  main-population gate, %s unstained: kept %d of %d events (%.1f%%)\n",
                target, sum( unst.keep ), length( unst.keep ),
                100 * mean( unst.keep ) ) )
  unst         <- unst[ unst.keep, , drop = FALSE ]
  unst.scatter <- unst.scatter[ unst.keep, , drop = FALSE ]
  unst.scatter <- diag.unstained[[ target ]][ , scatter.det, drop = FALSE ]
  unst.um      <- unmix.ols.fast( unst, s.start )
  thr          <- apply( unst.um, 2, stats::quantile, probs = diag.threshold.probs )
  
  unmixed <- unmix.ols.fast( raw.all, s.start )
  
  ss.mat  <- diag.variants[[ target ]]$spillover.spread
  af.only <- if ( !is.null( af.name ) && af.name %in% fluors ) af.name else NULL
  
  # ---- threshold audit, before any clustering ----
  # A fluorophore with pct.above near zero either has no positive events in the
  # concatenated file at all, or a threshold sitting above its own positive
  # population. The two need different fixes, so separate them here rather than
  # discovering it downstream as an empty cluster group.
  
  thr.table <- data.frame(
    fluorophore = panel,
    threshold   = thr[ panel ],
    q99.9       = apply( unmixed[ , panel, drop = FALSE ], 2, stats::quantile,
                         probs = 0.999 ),
    max.unmixed = apply( unmixed[ , panel, drop = FALSE ], 2, max ),
    pct.above   = 100 * colMeans(
      sweep( unmixed[ , panel, drop = FALSE ], 2, thr[ panel ], ">" ) ),
    row.names   = NULL
  )
  
  # ---- clustering: density-driven SOM versus abundance stratification ----
  
  cl.som <- cluster.unmixed.events(
    unmixed = unmixed, raw.data = raw.all, asp = asp,
    method = "som", som.dim = 20, verbose = FALSE
  )$cluster.id
  
  cl.strat <- cluster.unmixed.events(
    unmixed            = unmixed,
    raw.data           = raw.all,
    asp                = asp,
    method             = "stratify",
    unmixed.thresholds = thr,
    n.levels           = 10L,
    stratify.exclude   = af.only,
    verbose            = FALSE
  )$cluster.id
  
  # How many of ten abundance bands between a fluorophore's threshold and its
  # brightest cluster centroid actually contain a centroid. Bands are asinh
  # spaced so a decade at the dim end counts the same as a decade at the bright
  # end, and the threshold is floored at zero because a threshold taken against
  # the wrong spectra can come out negative.
  level.count <- function( cluster.id ) {
    cen <- .diag.centroids( cluster.id, unmixed, raw.all )
    vapply( panel, function( j ) {
      v <- cen$x.clust[ , j ]
      v <- v[ v > thr[ j ] ]
      if ( length( v ) < 2 ) return( 0L )
      lo <- asinh( max( thr[ j ], 0 ) / diag.cofactor )
      hi <- asinh( max( v ) / diag.cofactor )
      if ( !is.finite( lo ) || !is.finite( hi ) || hi <= lo ) return( 0L )
      brk <- seq( lo, hi, length.out = 11 )
      length( unique( cut( asinh( v / diag.cofactor ), breaks = brk,
                           include.lowest = TRUE ) ) )
    }, integer( 1 ) )
  }
  
  # Events assigned to each fluorophore's stratification group, reproducing the
  # dominance rule used inside cluster.unmixed.events().
  strat.score    <- sweep( unmixed[ , panel, drop = FALSE ], 2,
                           pmax( thr[ panel ], .Machine$double.eps ), "/" )
  strat.dominant <- max.col( strat.score, ties.method = "first" )
  strat.top      <- strat.score[ cbind( seq_len( nrow( strat.score ) ), strat.dominant ) ]
  strat.dominant[ strat.top <= 1 ] <- 0L
  
  strat.events <- vapply( seq_along( panel ), function( f )
    sum( strat.dominant == f ), integer( 1 ) )
  
  # ---- active-set purity, on the stratified centroids ----
  
  cen <- .diag.centroids( cl.strat, unmixed, raw.all )
  
  mask.flat <- sweep( cen$x.clust, 2, thr, ">" )
  
  thr.var  <- get.spread.thresholds(
    unmixed = cen$x.clust, thresholds = thr,
    spillover.spread = ss.mat, spread.kappa = 2, verbose = FALSE
  )
  mask.var <- cen$x.clust > thr.var
  
  if ( !is.null( af.only ) ) {
    mask.flat[ , af.only ] <- TRUE
    mask.var[ , af.only ]  <- TRUE
  }
  
  n.active.flat <- rowSums( mask.flat[ , panel, drop = FALSE ] )
  n.active.var  <- rowSums( mask.var[ , panel, drop = FALSE ] )
  
  c.table <- data.frame(
    target          = target,
    fluorophore     = panel,
    levels.som      = level.count( cl.som ),
    levels.strat    = level.count( cl.strat ),
    strat.events    = strat.events,
    single.pos.flat = vapply( panel, function( j )
      sum( n.active.flat == 1 & mask.flat[ , j ] ), integer( 1 ) ),
    single.pos.var  = vapply( panel, function( j )
      sum( n.active.var == 1 & mask.var[ , j ] ), integer( 1 ) ),
    row.names       = NULL
  )
  
  cat( sprintf( "-- %s raw data, starting from %s spectra --\n", target, wrong ) )
  print( thr.table, digits = 4 )
  cat( sprintf( "\n  events used: %d;  SOM clusters: %d;  stratified clusters: %d\n",
                nrow( raw.all ),
                length( unique( stats::na.omit( cl.som ) ) ),
                nrow( cen$x.clust ) ) )
  cat( sprintf( "  stratified cluster size: median %.0f, min %d, max %d\n",
                stats::median( cen$size ), min( cen$size ), max( cen$size ) ) )
  cat( sprintf( "  events not dominant for any fluorophore (background): %d (%.1f%%)\n",
                sum( strat.dominant == 0 ),
                100 * mean( strat.dominant == 0 ) ) )
  cat( sprintf( "  clusters by active-set size, flat threshold: %s\n",
                paste( utils::capture.output( table( n.active.flat ) ), collapse = " | " ) ) )
  cat( sprintf( "  clusters by active-set size, spread-scaled:  %s\n",
                paste( utils::capture.output( table( n.active.var ) ), collapse = " | " ) ) )
  print( c.table )
  cat( "\n" )
  
  c.results[[ target ]] <- list(
    table = c.table, spectra = s.start, raw = raw.all, unmixed = unmixed,
    scatter = scatter.all, unstained = unst, unstained.scatter = unst.scatter,
    thresholds = thr, ss.mat = ss.mat, panel = panel, fluors = fluors,
    cluster.id = cl.strat, strat.dominant = strat.dominant
  )
}

cat( "  INTERPRETATION: levels.som vs levels.strat says how much abundance\n" )
cat( "  resolution density-driven clustering was giving away. strat.events near\n" )
cat( "  zero with pct.above also near zero means the control is missing or empty;\n" )
cat( "  strat.events near zero with pct.above healthy means the dominance rule is\n" )
cat( "  losing that dye to a brighter neighbour. single.pos.flat vs single.pos.var\n" )
cat( "  says how many clean single-positive clusters a flat threshold is losing.\n" )

# ---------------------------------------------------------------------------
# STAGE D - one-shot single-positive correction
# ---------------------------------------------------------------------------
# For each fluorophore, take only the events that are positive for it and
# nothing else, bin them into abundance levels, unmix each bin against a
# restricted design (that fluorophore plus autofluorescence), and regress the
# residual on abundance. The slope is the correction to that fluorophore's
# spectrum. This is not bound by the Stage B ceiling: with a rank-1 or rank-2
# design the only unidentifiable direction is the row's own scale, which
# renormalisation removes.

cat( "\n============ STAGE D: single-positive correction ============\n" )

#' Per-fluorophore spectral correction from dominance-assigned populations.
#'
#' Each event is assigned to the fluorophore it is most strongly positive for
#' relative to that fluorophore's own dynamic range, so a dye that is never
#' strictly single-positive (because a spectrally collinear partner is always
#' above threshold too) still gets a population. Within a group, any co-active
#' dye is carried in the restricted design as a nuisance term; only the
#' dominant dye's slope is used.
#'
#' Background is removed from the raw data before fitting rather than modelled
#' as a row in the design, so the dye's abundance and the background amount are
#' no longer estimated from the same spectral vector.
#'
#' @noRd
.diag.dominance.fit <- function(
    raw.data, spectra, thresholds, panel, af.name,
    spillover.spread  = NULL,
    spread.kappa      = 2,
    n.levels          = 10L,
    n.iter            = 6L,
    min.events        = 200L,
    min.span          = 5,
    min.explained     = 0.5,
    min.gain          = 0.002,
    step.grid         = c( 0, 0.03125, 0.0625, 0.125, 0.25, 0.5, 1 ),
    n.split.trials    = 1L,
    min.split.frac    = 0.6,
    max.step          = 0.15,
    max.span.drift    = 1.10,
    max.bg.alignment  = -0.9,
    nuisance.frac     = 0.5,
    background.n      = 5000L,
    bg.mode           = c( "scatter.knn", "global.mean", "none" ),
    scatter           = NULL,
    unstained         = NULL,
    unstained.scatter = NULL,
    k.neighbors       = 20L
) {
  
  bg.mode <- match.arg( bg.mode )
  
  unmixed <- unmix.ols.fast( raw.data, spectra )
  
  if ( is.null( spillover.spread ) ) {
    above <- sweep( unmixed[ , panel, drop = FALSE ], 2, thresholds[ panel ], ">" )
  } else {
    threshold.matrix <- get.spread.thresholds(
      unmixed          = unmixed,
      thresholds       = thresholds,
      spillover.spread = spillover.spread,
      spread.kappa     = spread.kappa,
      verbose          = FALSE
    )
    above <- unmixed[ , panel, drop = FALSE ] >
      threshold.matrix[ , panel, drop = FALSE ]
  }
  
  # Dominance is scored as a fraction of each dye's own dynamic range above its
  # threshold, not as a multiple of the threshold. Thresholds vary by hundreds
  # of fold across a panel and can be negative when taken against the wrong
  # spectra, so dividing by them hands every dim or ambiguous event to whichever
  # dye happens to have the smallest threshold.
  excess <- pmax( sweep( unmixed[ , panel, drop = FALSE ],
                         2, thresholds[ panel ], "-" ), 0 )
  
  dyn.range <- apply( unmixed[ , panel, drop = FALSE ], 2, stats::quantile,
                      probs = 0.999 ) - thresholds[ panel ]
  dyn.range <- pmax( dyn.range, .Machine$double.eps )
  
  score    <- sweep( excess, 2, dyn.range, "/" )
  dominant <- max.col( score, ties.method = "first" )
  top      <- score[ cbind( seq_len( nrow( score ) ), dominant ) ]
  dominant[ top <= 0 ] <- 0L
  
  spectra.new    <- spectra
  background.idx <- which( dominant == 0L )
  
  # Scatter-matched subtraction is only meaningful when scatter predicts
  # background, which is true of cells and false of beads: uniform particles
  # have uninformative scatter, so kNN matching draws an essentially arbitrary
  # unstained event and replaces a near-constant background with added variance.
  if ( bg.mode == "scatter.knn" && !is.null( scatter ) &&
       !is.null( unstained ) && !is.null( unstained.scatter ) ) {
    
    raw.data <- .diag.bg.subtract(
      raw.data, scatter, unstained, unstained.scatter, k.neighbors = k.neighbors )
    
  } else if ( bg.mode != "none" && length( background.idx ) >= min.events ) {
    
    raw.data <- sweep( raw.data, 2,
                       colMeans( raw.data[ background.idx, , drop = FALSE ] ), "-" )
    
  }
  
  if ( length( background.idx ) > background.n )
    background.idx <- sample( background.idx, background.n )
  
  fit.log    <- list()
  span.first <- stats::setNames( rep( NA_real_, length( panel ) ), panel )
  
  for ( iter in seq_len( n.iter ) ) {
    
    for ( f in seq_along( panel ) ) {
      
      j   <- panel[ f ]
      idx <- which( dominant == f )
      if ( length( idx ) < min.events ) next
      
      # co-active dyes carried as nuisance predictors, not discarded
      co.frac  <- colMeans( above[ idx, , drop = FALSE ] )
      nuisance <- setdiff( panel[ co.frac > nuisance.frac ], j )
      active   <- c( j, nuisance )
      
      # Slope from a chosen subset of the group's events: bin by abundance,
      # anchor with the background population, regress the restricted residual
      # on the dominant dye's abundance.
      fit.slope <- function( use.idx ) {
        
        if ( length( use.idx ) < min.events %/% 2 ) return( NULL )
        
        y.use <- raw.data[ use.idx, , drop = FALSE ]
        x.use <- .diag.restricted.unmix( y.use, spectra.new, active )[ , 1 ]
        
        brk <- unique( stats::quantile(
          x.use, probs = seq( 0, 1, length.out = n.levels + 1 ) ) )
        if ( length( brk ) < 3 ) return( NULL )
        
        bin   <- as.integer( cut( x.use, breaks = brk, include.lowest = TRUE ) )
        y.bin <- t( vapply( sort( unique( bin ) ), function( b )
          colMeans( y.use[ bin == b, , drop = FALSE ] ), numeric( ncol( y.use ) ) ) )
        
        if ( length( background.idx ) >= min.events )
          y.bin <- rbind(
            colMeans( raw.data[ background.idx, , drop = FALSE ] ), y.bin )
        
        if ( nrow( y.bin ) < length( active ) + 3 ) return( NULL )
        
        x.bin <- .diag.restricted.unmix( y.bin, spectra.new, active )
        r.bin <- y.bin - x.bin %*% spectra.new[ active, , drop = FALSE ]
        
        fit   <- stats::lm.fit( x = cbind( 1, x.bin ), y = r.bin )
        slope <- stats::coef( fit )[ 2, ]
        slope[ !is.finite( slope ) ] <- 0
        
        ss.res <- sum( fit$residuals^2 )
        ss.tot <- sum( sweep( r.bin, 2, colMeans( r.bin ) )^2 )
        
        # At its brightest, how much of the background-subtracted signal does
        # this dye's own term account for? A control where the answer is small
        # is not really a control for this dye: the regression is then fitting
        # whatever else is in the tube, and the dye is too dim to be
        # contributing much error to the panel in the first place.
        top.bin   <- which.max( x.bin[ , 1 ] )
        top.norm  <- sqrt( sum( y.bin[ top.bin, ]^2 ) )
        explained <- if ( top.norm > 0 )
          sqrt( sum( ( x.bin[ top.bin, 1 ] * spectra.new[ j, ] )^2 ) ) / top.norm else 0
        
        list(
          slope     = slope,
          r.sq      = if ( ss.tot > 0 ) 1 - ss.res / ss.tot else 0,
          x.span    = max( x.bin[ , 1 ] ) - min( x.bin[ , 1 ] ),
          explained = explained,
          n.bins    = nrow( r.bin )
        )
      }
      
      full <- fit.slope( idx )
      if ( is.null( full ) ) next
      
      # Background-confound statistic. A genuine spectral error gives the
      # event-level regression an intercept and a slope with independent
      # physical origins. A common-mode background residual whose magnitude
      # tracks brightness is one physical direction split in two by the
      # regression, so the intercept and slope come out anti-collinear.
      # Strong anti-alignment marks a correction that would walk the row
      # into the background direction rather than fix its shape.
      y.evt <- raw.data[ idx, , drop = FALSE ]
      x.evt <- .diag.restricted.unmix( y.evt, spectra.new, active )
      r.evt <- y.evt - x.evt %*% spectra.new[ active, , drop = FALSE ]
      
      fit.evt   <- stats::lm.fit( x = cbind( 1, x.evt[ , 1 ] ), y = r.evt )
      alpha.evt <- stats::coef( fit.evt )[ 1, ]
      beta.evt  <- stats::coef( fit.evt )[ 2, ]
      alpha.evt[ !is.finite( alpha.evt ) ] <- 0
      beta.evt[ !is.finite( beta.evt ) ]   <- 0
      
      alpha.norm <- sqrt( sum( alpha.evt^2 ) )
      beta.norm  <- sqrt( sum( beta.evt^2 ) )
      bg.align   <- if ( alpha.norm > 0 && beta.norm > 0 )
        sum( alpha.evt * beta.evt ) / ( alpha.norm * beta.norm ) else NA_real_
      
      # Held-out step search. Fit on one half, pick the step size that most
      # reduces the restricted residual on the other. This bounds the step but
      # cannot on its own tell a corrected row from a repurposed one, since the
      # residual objective constrains the span rather than the rows - hence the
      # magnitude and drift gates below.
      #
      # A single 50/50 split is a noisy instrument for a fluorophore whose
      # true correction is real but small relative to per-event noise: the
      # default n.split.trials = 1 reproduces the original fixed,
      # un-reseeded rep_len(c(1L, 2L), length(idx)) split exactly, for
      # backward compatibility. Raising n.split.trials instead draws that
      # many independent random 50/50 splits and requires only
      # min.split.frac of them to agree a step helps, taking the median of
      # the agreeing splits' step and gain - trading one noisy verdict for
      # a vote.
      residual.gain <- function( fit.idx, test.idx ) {
        
        fs <- fit.slope( fit.idx )
        if ( is.null( fs ) ) return( NULL )
        
        y.test <- raw.data[ test.idx, , drop = FALSE ]
        base   <- sqrt( sum( y.test^2 ) )
        if ( base <= 0 ) return( NULL )
        
        obj <- vapply( step.grid, function( t ) {
          s.try <- spectra.new
          s.try[ j, ] <- pmax( s.try[ j, ] + t * fs$slope, 0 )
          if ( max( s.try[ j, ] ) <= 0 ) return( Inf )
          s.try <- .diag.renorm( s.try )
          x.try <- .diag.restricted.unmix( y.test, s.try, active )
          sqrt( sum( ( y.test - x.try %*% s.try[ active, , drop = FALSE ] )^2 ) ) / base
        }, numeric( 1 ) )
        
        best <- which.min( obj )
        list( t = step.grid[ best ], gain = obj[ 1 ] - obj[ best ] )
      }
      
      split.trial <- function( half ) {
        
        idx.a <- idx[ half == 1L ]
        idx.b <- idx[ half == 2L ]
        
        gain.ab <- residual.gain( idx.a, idx.b )
        gain.ba <- residual.gain( idx.b, idx.a )
        
        if ( is.null( gain.ab ) || is.null( gain.ba ) )
          return( c( t = 0, gain = NA_real_ ) )
        
        c( t = min( gain.ab$t, gain.ba$t ),
           gain = min( gain.ab$gain, gain.ba$gain ) )
      }
      
      if ( n.split.trials <= 1L ) {
        
        one   <- split.trial( rep_len( c( 1L, 2L ), length( idx ) ) )
        t.hat <- one[ "t" ]
        gain  <- one[ "gain" ]
        
      } else {
        
        trials <- t( vapply( seq_len( n.split.trials ), function( trial )
          split.trial( sample( rep_len( c( 1L, 2L ), length( idx ) ) ) ),
          numeric( 2 ) ) )
        
        passed <- trials[ , "t" ] > 0
        
        if ( any( passed ) && mean( passed ) >= min.split.frac ) {
          t.hat <- stats::median( trials[ passed, "t" ] )
          gain  <- stats::median( trials[ passed, "gain" ] )
        } else {
          t.hat <- 0
          gain  <- NA_real_
        }
      }
      
      slope     <- full$slope
      step.norm <- sqrt( sum( slope^2 ) )
      row.norm  <- sqrt( sum( spectra.new[ j, ]^2 ) )
      rel.step  <- step.norm / row.norm
      
      if ( is.na( span.first[ j ] ) ) span.first[ j ] <- full$x.span
      span.drift <- full$x.span / span.first[ j ]
      
      # A correction larger than the row it corrects is wrong on its face,
      # whatever it does to the residual, so a step that would need scaling down
      # to fit is rejected outright rather than applied in miniature - repeatedly
      # applying a shrunken version of a wrong direction looks like convergence
      # and is worse than not acting. Span drift catches the same pathology from
      # the other side: a row rotating toward whatever else the control contains
      # makes its own apparent abundance grow, feeding the next iteration.
      accepted <- full$x.span > min.span * abs( thresholds[ j ] ) &&
        is.finite( full$explained ) && full$explained > min.explained &&
        t.hat > 0 && is.finite( gain ) && gain > min.gain &&
        is.finite( rel.step ) && rel.step > 0 && rel.step <= max.step &&
        is.finite( span.drift ) && span.drift <= max.span.drift &&
        ( !is.finite( bg.align ) || bg.align > max.bg.alignment )
      
      if ( accepted )
        spectra.new[ j, ] <- pmax( spectra.new[ j, ] + t.hat * slope, 0 )
      
      fit.log[[ length( fit.log ) + 1L ]] <- data.frame(
        iter        = iter,
        fluorophore = j,
        n.events    = length( idx ),
        n.nuisance  = length( nuisance ),
        n.bins      = full$n.bins,
        x.span      = full$x.span,
        explained   = full$explained,
        r.squared   = full$r.sq,
        bg.align    = bg.align,
        t.hat       = t.hat,
        gain        = gain,
        rel.step    = rel.step,
        span.drift  = span.drift,
        accepted    = accepted,
        row.names   = NULL
      )
    }
    
    spectra.new <- .diag.renorm( spectra.new )
  }
  
  list( spectra = spectra.new, fit.log = fit.log, dominant = dominant )
}

for ( target in c( "Cells", "Beads" ) ) {
  
  ctx     <- c.results[[ target ]]
  af.pcs  <- .diag.af.pcs( ctx$unstained )
  
  hot <- vapply( ctx$panel, function( j )
    .diag.basis.hotspot( rbind( af.pcs, ctx$spectra[ j, , drop = FALSE ] ), j ),
    numeric( 1 ) )
  
  cat( sprintf( "-- %s: hotspot scale of [AF PCs; dye] --\n", target ) )
  print( round( sort( hot, decreasing = TRUE ), 2 ) )
  cat( "\n" )
}

d.results <- list()

for ( target in c( "Cells", "Beads" ) ) {

  ctx     <- c.results[[ target ]]
  s.true  <- as.matrix( diag.spectra[[ target ]][ rownames( ctx$spectra ), b.shared.det,
                                                  drop = FALSE ] )

  d.fit <- .diag.dominance.fit(
    raw.data          = ctx$raw,
    spectra           = ctx$spectra,
    thresholds        = ctx$thresholds,
    panel             = ctx$panel,
    af.name           = af.name,
    spillover.spread  = ctx$ss.mat,
    spread.kappa      = 2,
    bg.mode           = if ( target == "Cells" ) "scatter.knn" else "global.mean",
    scatter           = ctx$scatter,
    unstained         = ctx$unstained,
    unstained.scatter = ctx$unstained.scatter
  )

  ceiling.col <- b.report[[ target ]]$cos.ceiling
  names( ceiling.col ) <- b.report[[ target ]]$fluorophore

  d.table <- data.frame(
    target      = target,
    fluorophore = ctx$panel,
    cos.start   = .diag.row.cosine( ctx$spectra[ ctx$panel, , drop = FALSE ],
                                    s.true[ ctx$panel, , drop = FALSE ] ),
    cos.after   = .diag.row.cosine( d.fit$spectra[ ctx$panel, , drop = FALSE ],
                                    s.true[ ctx$panel, , drop = FALSE ] ),
    cos.ceiling = ceiling.col[ ctx$panel ],
    row.names   = NULL
  )
  
  fit.table <- do.call( rbind, d.fit$fit.log )
  
  d.table$n.events <- vapply( ctx$panel, function( j ) {
    rows <- fit.table$n.events[ fit.table$fluorophore == j ]
    if ( length( rows ) == 0 ) 0L else as.integer( rows[ 1 ] )
  }, integer( 1 ) )
  
  d.table$n.steps <- vapply( ctx$panel, function( j )
    sum( fit.table$accepted[ fit.table$fluorophore == j ] ), integer( 1 ) )
  
  d.table$span.drift <- vapply( ctx$panel, function( j ) {
    v <- fit.table$x.span[ fit.table$fluorophore == j ]
    if ( length( v ) < 2 ) 1 else v[ length( v ) ] / v[ 1 ]
  }, numeric( 1 ) )
  d.table$deg.start <- .diag.angle( ctx$spectra[ ctx$panel, , drop = FALSE ],
                                    s.true[ ctx$panel, , drop = FALSE ] )
  d.table$deg.after <- .diag.angle( d.fit$spectra[ ctx$panel, , drop = FALSE ],
                                    s.true[ ctx$panel, , drop = FALSE ] )
  d.table$recovered <- .diag.recovered.fraction(
    ctx$spectra[ ctx$panel, , drop = FALSE ],
    d.fit$spectra[ ctx$panel, , drop = FALSE ],
    s.true[ ctx$panel, , drop = FALSE ]
  )

  cat( sprintf( "-- %s raw data --\n", target ) )
  print( d.table, digits = 6 )
  cat( "\n  per-fit diagnostics (first iteration and every accepted step):\n" )
  print( fit.table[ fit.table$iter == 1 | fit.table$accepted, ], digits = 4 )
  cat( "\n  steps accepted per fluorophore:\n" )
  print( table( fit.table$fluorophore[ fit.table$accepted ] ) )
  cat( "\n" )

  d.results[[ target ]] <- list( table = d.table, spectra = d.fit$spectra,
                                 fit.log = d.fit$fit.log,
                                 dominant = d.fit$dominant )
}

cat( "  INTERPRETATION: delta > 0 on most rows, and cos.after exceeding the\n" )
cat( "  full-panel cos.ceiling, is the result that would justify making the\n" )
cat( "  restricted path the default. n.events near zero means the threshold or the\n" )
cat( "  single-positive definition is too strict for that dye.\n" )

# ---------------------------------------------------------------------------
# STAGE E - downstream effect on unmixed data
# ---------------------------------------------------------------------------
# Cosine similarity to a reference spectrum is a proxy. What a user sees is the
# unmixed abundances, so measure those directly. Leakage needs no ground truth
# and is therefore the metric that survives into production; abundance error
# against the correct spectra is the benchmark available only here.

cat( "\n============ STAGE E: downstream effect on unmixed data ============\n" )

e.results <- list()

for ( target in c( "Cells", "Beads" ) ) {
  
  ctx <- c.results[[ target ]]
  fl  <- rownames( ctx$spectra )
  
  s.start <- ctx$spectra
  s.corr  <- d.results[[ target ]]$spectra[ fl, , drop = FALSE ]
  s.true  <- as.matrix( diag.spectra[[ target ]][ fl, b.shared.det, drop = FALSE ] )
  
  # The AF row is not corrected by this method and its reference value is not
  # trustworthy in either table, so hold it fixed across all three unmixes and
  # let the contrast isolate the fluorophore rows.
  if ( !is.null( af.name ) && af.name %in% fl )
    s.true[ af.name, ] <- s.start[ af.name, ]
  
  um <- list(
    start     = unmix.ols.fast( ctx$raw, s.start ),
    corrected = unmix.ols.fast( ctx$raw, s.corr ),
    true      = unmix.ols.fast( ctx$raw, s.true )
  )
  
  dom <- d.results[[ target ]]$dominant
  
  # Leakage: total off-target abundance within a dye's own population, as a
  # fraction of that dye's own signal. In a single-stained control every other
  # channel should read zero, so this measures unmixing error without needing
  # to know the correct spectra.
  leak <- vapply( seq_along( ctx$panel ), function( f ) {
    
    idx <- which( dom == f )
    if ( length( idx ) < 100 ) return( rep( NA_real_, 3 ) )
    
    j   <- ctx$panel[ f ]
    off <- setdiff( ctx$panel, j )
    
    vapply( um, function( u ) {
      on <- stats::median( u[ idx, j ] )
      if ( !is.finite( on ) || on <= 0 ) return( NA_real_ )
      sum( abs( apply( u[ idx, off, drop = FALSE ], 2, stats::median ) ) ) / on
    }, numeric( 1 ) )
    
  }, numeric( 3 ) )
  
  # Per-event abundance error against the correct spectra, scaled by that
  # channel's own spread, so a channel is judged against the resolution it
  # actually offers rather than in raw units.
  agree <- vapply( ctx$panel, function( j ) {
    sc <- stats::mad( um$true[ , j ] )
    if ( !is.finite( sc ) || sc <= 0 ) return( c( NA_real_, NA_real_ ) )
    c( stats::median( abs( um$start[ , j ]     - um$true[ , j ] ) ) / sc,
       stats::median( abs( um$corrected[ , j ] - um$true[ , j ] ) ) / sc )
  }, numeric( 2 ) )
  
  e.table <- data.frame(
    target      = target,
    fluorophore = ctx$panel,
    leak.start  = leak[ 1, ],
    leak.corr   = leak[ 2, ],
    leak.true   = leak[ 3, ],
    err.start   = agree[ 1, ],
    err.corr    = agree[ 2, ],
    row.names   = NULL
  )
  e.table$leak.gain <- e.table$leak.start - e.table$leak.corr
  
  cat( sprintf( "-- %s raw data --\n", target ) )
  print( e.table, digits = 3 )
  cat( sprintf( "\n  median leakage: start %.3f -> corrected %.3f (correct spectra %.3f)\n",
                stats::median( e.table$leak.start, na.rm = TRUE ),
                stats::median( e.table$leak.corr,  na.rm = TRUE ),
                stats::median( e.table$leak.true,  na.rm = TRUE ) ) )
  cat( sprintf( "  median abundance error vs correct spectra: %.2f -> %.2f MADs\n\n",
                stats::median( e.table$err.start, na.rm = TRUE ),
                stats::median( e.table$err.corr,  na.rm = TRUE ) ) )
  
  e.results[[ target ]] <- list( table = e.table, unmixed = um )
}

cat( "  INTERPRETATION: leak.corr should sit between leak.start and leak.true;\n" )
cat( "  leak.true is the floor set by the spectra themselves. err.corr below\n" )
cat( "  err.start means the abundances a user reads off are closer to the ones\n" )
cat( "  the correct spectra would have given.\n" )

# Stages G-I of the signature-error-correction diagnostics.
#
# Requires the environment left by run_signature_correction_diagnostics_1.R:
# c.results, diag.spectra, diag.variants, b.shared.det, af.name, and the
# .diag.* helper functions. Run top to bottom; the printed tables are the
# output to report back.
#
# Stage G  Descent audit. For each dominance population, decomposes the
#          first derivative of the held-out residual objective at step zero
#          into its slope and intercept terms. A dye with a strong
#          abundance-correlated residual but a non-descending objective
#          (Jprime0 >= 0) is being vetoed by the intercept term, not by a
#          wrong slope: the line-search gate is rejecting a correct
#          correction because it does not also absorb the population's
#          constant background misfit.
#
# Stage H  Estimator comparison and curvature. Runs alternating least
#          squares (ALS) on each dominance population to convergence, with
#          the abundance, the background vector and the dominant row all
#          re-fit jointly, and compares one-step slope, converged ALS, and
#          ground truth. Also measures whether the residual direction
#          drifts across abundance bins, and whether that drift lies in the
#          span of the dye's own spectral-variant deltas (the tandem
#          variant-mixture hypothesis).
#
# Stage I  Restricted identifiability ceiling. With ground truth available,
#          measures per dye how much of the true spectral error lies inside
#          the span of the restricted design's nuisance rows, and is
#          therefore invisible to any residual-based correction using that
#          design. Also reports the nuisance sets under flat and
#          spread-scaled co-activity masks.

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

diag2.min.events    <- 200L
diag2.nuisance.frac <- 0.5
diag2.spread.kappa  <- 2
diag2.n.bins        <- 8L
diag2.background.n  <- 3000L
diag2.als.max.iter  <- 30L
diag2.als.tol       <- 1e-6

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

#' Dominance assignment, background handling, and co-activity masks for one
#' particle type, replicating the Stage D configuration so every stage here
#' sees the same populations Stage D fitted.
#' @noRd
.diag2.population <- function( target ) {
  
  ctx        <- c.results[[ target ]]
  panel      <- ctx$panel
  thresholds <- ctx$thresholds
  
  unmixed <- unmix.ols.fast( ctx$raw, ctx$spectra )
  
  excess <- pmax( sweep( unmixed[ , panel, drop = FALSE ],
                         2, thresholds[ panel ], "-" ), 0 )
  
  dyn.range <- apply( unmixed[ , panel, drop = FALSE ], 2, stats::quantile,
                      probs = 0.999 ) - thresholds[ panel ]
  dyn.range <- pmax( dyn.range, .Machine$double.eps )
  
  score    <- sweep( excess, 2, dyn.range, "/" )
  dominant <- max.col( score, ties.method = "first" )
  top      <- score[ cbind( seq_len( nrow( score ) ), dominant ) ]
  dominant[ top <= 0 ] <- 0L
  
  background.idx <- which( dominant == 0L )
  
  raw.use <- ctx$raw
  if ( target == "Cells" ) {
    raw.use <- .diag.bg.subtract(
      ctx$raw, ctx$scatter, ctx$unstained, ctx$unstained.scatter )
  } else if ( length( background.idx ) >= diag2.min.events ) {
    raw.use <- sweep( ctx$raw, 2,
                      colMeans( ctx$raw[ background.idx, , drop = FALSE ] ), "-" )
  }
  
  above.flat <- sweep( unmixed[ , panel, drop = FALSE ],
                       2, thresholds[ panel ], ">" )
  
  above.spread <- NULL
  if ( !is.null( ctx$ss.mat ) ) {
    threshold.matrix <- get.spread.thresholds(
      unmixed          = unmixed,
      thresholds       = thresholds,
      spillover.spread = ctx$ss.mat,
      spread.kappa     = diag2.spread.kappa,
      verbose          = FALSE
    )
    above.spread <- unmixed[ , panel, drop = FALSE ] >
      threshold.matrix[ , panel, drop = FALSE ]
  }
  
  list( ctx = ctx, panel = panel, thresholds = thresholds,
        dominant = dominant, background.idx = background.idx,
        raw.use = raw.use, above.flat = above.flat,
        above.spread = above.spread )
}

#' Nuisance set for one dominance population under a given co-activity mask.
#' @noRd
.diag2.nuisance <- function( pop, idx, j, mask ) {
  co.frac <- colMeans( mask[ idx, , drop = FALSE ] )
  setdiff( pop$panel[ co.frac > diag2.nuisance.frac ], j )
}

#' Alternating least squares refit of one spectra row on a dominance
#' population. The abundance vector, the constant background vector, and the
#' dominant dye's row are re-fit jointly; nuisance rows stay fixed. This is
#' the estimator whose first linearised step equals the intercept-model
#' regression slope: iterating it to convergence removes the need for a step
#' size, a direction search, or a learning rate.
#' @noRd
.diag2.als <- function( y, spectra, j, nuisance,
                        max.iter = diag2.als.max.iter,
                        tol      = diag2.als.tol ) {
  
  active <- c( j, nuisance )
  s      <- spectra
  alpha  <- rep( 0, ncol( y ) )
  
  for ( it in seq_len( max.iter ) ) {
    
    y.adj <- sweep( y, 2, alpha )
    x     <- .diag.restricted.unmix( y.adj, s, active )
    r     <- y.adj - x %*% s[ active, , drop = FALSE ]
    alpha <- alpha + colMeans( r )
    
    xj    <- x[ , 1 ]
    denom <- sum( xj^2 )
    if ( denom <= 0 ) break
    
    nuisance.part <- if ( length( nuisance ) > 0 )
      x[ , -1, drop = FALSE ] %*% s[ nuisance, , drop = FALSE ] else 0
    
    row.new <- colSums( xj * ( sweep( y, 2, alpha ) - nuisance.part ) ) / denom
    row.new <- pmax( row.new, 0 )
    if ( max( row.new ) <= 0 ) break
    row.new <- row.new / max( row.new )
    
    step <- sqrt( sum( ( row.new - s[ j, ] )^2 ) )
    s[ j, ] <- row.new
    if ( step < tol ) break
  }
  
  list( spectra = s, alpha = alpha, n.iter = it )
}

#' Cosine between two vectors, guarded.
#' @noRd
.diag2.vcos <- function( a, b ) {
  na <- sqrt( sum( a^2 ) )
  nb <- sqrt( sum( b^2 ) )
  if ( na <= 0 || nb <= 0 ) return( NA_real_ )
  sum( a * b ) / ( na * nb )
}

#' Fraction of a vector's norm lying in the row span of a matrix.
#' @noRd
.diag2.span.fraction <- function( v, basis ) {
  if ( is.null( basis ) || nrow( basis ) == 0 ) return( 0 )
  proj <- t( basis ) %*% MASS::ginv( basis %*% t( basis ) ) %*% basis
  nv   <- sqrt( sum( v^2 ) )
  if ( nv <= 0 ) return( NA_real_ )
  sqrt( sum( ( v %*% proj )^2 ) ) / nv
}

# ---------------------------------------------------------------------------
# STAGE G - descent audit of the residual line search
# ---------------------------------------------------------------------------
# With the restricted design S_A and a perturbation of row j along beta, the
# envelope theorem gives the derivative of the re-fit residual objective at
# step zero as
#
#   J'(0) = -2 [ (sum x_j) (alpha . beta) + (sum x_j^2) ||beta||^2 + eps ]
#
# where alpha and beta are the intercept and slope of the per-detector
# regression of the restricted residual on the dominant abundance. The slope
# term always favours stepping; the intercept term vetoes it whenever the
# constant background misfit is anti-correlated with the correction. A dye
# with high r.squared and t.hat = 0 in Stage D should show term.intercept
# negative and larger in magnitude than term.slope here.

cat( "\n============ STAGE G: descent audit ============\n" )

g.results <- list()

for ( target in c( "Cells", "Beads" ) ) {
  
  pop    <- .diag2.population( target )
  s.true <- as.matrix( diag.spectra[[ target ]][ rownames( pop$ctx$spectra ),
                                                 b.shared.det, drop = FALSE ] )
  
  g.rows  <- list()
  g.store <- list()
  
  for ( f in seq_along( pop$panel ) ) {
    
    j   <- pop$panel[ f ]
    idx <- which( pop$dominant == f )
    if ( length( idx ) < diag2.min.events ) next
    
    nuisance <- .diag2.nuisance( pop, idx, j, pop$above.flat )
    active   <- c( j, nuisance )
    
    y  <- pop$raw.use[ idx, , drop = FALSE ]
    x  <- .diag.restricted.unmix( y, pop$ctx$spectra, active )
    xj <- x[ , 1 ]
    r  <- y - x %*% pop$ctx$spectra[ active, , drop = FALSE ]
    
    fit   <- stats::lm.fit( x = cbind( 1, xj ), y = r )
    alpha <- stats::coef( fit )[ 1, ]
    beta  <- stats::coef( fit )[ 2, ]
    alpha[ !is.finite( alpha ) ] <- 0
    beta[ !is.finite( beta ) ]   <- 0
    
    s1 <- sum( xj )
    s2 <- sum( xj^2 )
    
    term.slope     <- s2 * sum( beta^2 )
    term.intercept <- s1 * sum( alpha * beta )
    jprime.emp     <- -2 * sum( xj * ( r %*% beta ) )
    
    beta.noint <- colSums( xj * r ) / s2
    e.true     <- s.true[ j, ] - pop$ctx$spectra[ j, ]
    
    g.rows[[ j ]] <- data.frame(
      target         = target,
      fluorophore    = j,
      n.events       = length( idx ),
      n.nuisance     = length( nuisance ),
      cos.alpha.beta = .diag2.vcos( alpha, beta ),
      term.slope     = term.slope,
      term.intercept = term.intercept,
      jprime0        = -2 * ( term.intercept + term.slope ),
      jprime0.emp    = jprime.emp,
      descends       = jprime.emp < 0,
      cos.beta.E     = .diag2.vcos( beta, e.true ),
      cos.betaNI.E   = .diag2.vcos( beta.noint, e.true ),
      row.names      = NULL
    )
    
    g.store[[ j ]] <- list( idx = idx, nuisance = nuisance,
                            beta = beta, alpha = alpha )
  }
  
  g.table <- do.call( rbind, g.rows )
  cat( sprintf( "-- %s --\n", target ) )
  print( g.table, digits = 4 )
  cat( "\n" )
  
  g.results[[ target ]] <- list( table = g.table, store = g.store, pop = pop,
                                 s.true = s.true )
}

cat( "  INTERPRETATION: descends = FALSE with cos.beta.E clearly positive is\n" )
cat( "  the smoking gun: the slope points at the true error, and the residual\n" )
cat( "  line search rejects it because term.intercept outweighs term.slope.\n" )
cat( "  cos.betaNI.E below cos.beta.E confirms the no-intercept (descent)\n" )
cat( "  direction is the background-contaminated one.\n" )

# ---------------------------------------------------------------------------
# STAGE H - ALS estimator, concordance, and residual curvature
# ---------------------------------------------------------------------------
# Three corrections per dye: the one-step intercept-model slope (what Stage D
# applies when accepted), the converged ALS refit, and ground truth. The
# concordance statistic cos(step, als) is a candidate production gate: two
# estimators with different failure modes agreeing on a direction. The
# curvature block bins the population by abundance and asks whether the
# residual direction drifts across bins, and whether that drift lies in the
# span of the dye's own spectral-variant deltas.

cat( "\n============ STAGE H: ALS, concordance, curvature ============\n" )

h.results <- list()

for ( target in c( "Cells", "Beads" ) ) {
  
  gres   <- g.results[[ target ]]
  pop    <- gres$pop
  s.true <- gres$s.true
  
  delta.list <- diag.variants[[ target ]]$delta.list
  
  h.rows <- list()
  h.spectra.als <- pop$ctx$spectra
  
  bg.idx <- pop$background.idx
  if ( length( bg.idx ) > diag2.background.n )
    bg.idx <- sample( bg.idx, diag2.background.n )
  
  for ( j in names( gres$store ) ) {
    
    st       <- gres$store[[ j ]]
    idx      <- st$idx
    nuisance <- st$nuisance
    active   <- c( j, nuisance )
    
    s.start.row <- pop$ctx$spectra[ j, ]
    e.true      <- s.true[ j, ] - s.start.row
    
    # one-step estimator: current Stage D update at unit step
    row.step <- pmax( s.start.row + st$beta, 0 )
    row.step <- row.step / max( row.step )
    
    # converged ALS, with background events appended to anchor alpha
    y.fit <- rbind( pop$raw.use[ idx, , drop = FALSE ],
                    pop$raw.use[ bg.idx, , drop = FALSE ] )
    als   <- .diag2.als( y.fit, pop$ctx$spectra, j, nuisance )
    row.als <- als$spectra[ j, ]
    h.spectra.als[ j, ] <- row.als
    
    u.step <- row.step - s.start.row
    u.als  <- row.als  - s.start.row
    
    # curvature: residual direction per abundance bin
    y  <- pop$raw.use[ idx, , drop = FALSE ]
    x  <- .diag.restricted.unmix( y, pop$ctx$spectra, active )
    xj <- x[ , 1 ]
    r  <- y - x %*% pop$ctx$spectra[ active, , drop = FALSE ]
    
    brk <- unique( stats::quantile(
      xj, probs = seq( 0, 1, length.out = diag2.n.bins + 1 ) ) )
    
    cos.dir.min    <- NA_real_
    drift.in.delta <- NA_real_
    
    if ( length( brk ) >= 4 ) {
      
      bin   <- as.integer( cut( xj, breaks = brk, include.lowest = TRUE ) )
      bins  <- sort( unique( bin ) )
      r.bin <- t( vapply( bins, function( b )
        colMeans( r[ bin == b, , drop = FALSE ] ), numeric( ncol( r ) ) ) )
      
      dirs <- sweep( r.bin, 2, st$alpha )
      dir.norm <- sqrt( rowSums( dirs^2 ) )
      keep <- which( dir.norm > 0 )
      
      if ( length( keep ) >= 3 ) {
        
        dirs <- dirs[ keep, , drop = FALSE ] / dir.norm[ keep ]
        ref  <- dirs[ nrow( dirs ), ]
        
        cos.bin     <- as.vector( dirs %*% ref )
        cos.dir.min <- min( cos.bin )
        
        d.mat <- delta.list[[ j ]]
        if ( !is.null( d.mat ) && nrow( d.mat ) > 0 ) {
          shared <- intersect( colnames( d.mat ), colnames( r ) )
          if ( length( shared ) == ncol( r ) ) {
            d.mat <- as.matrix( d.mat[ , colnames( r ), drop = FALSE ] )
            drift <- sweep( dirs, 2, ref )
            drift.frac <- apply( drift[ -nrow( drift ), , drop = FALSE ], 1,
                                 .diag2.span.fraction, basis = d.mat )
            drift.in.delta <- stats::median( drift.frac, na.rm = TRUE )
          }
        }
      }
    }
    
    ang <- function( v ) .diag.angle( matrix( v, nrow = 1 ),
                                      matrix( s.true[ j, ], nrow = 1 ) )
    deg.start <- ang( s.start.row )
    
    h.rows[[ j ]] <- data.frame(
      target         = target,
      fluorophore    = j,
      deg.start      = deg.start,
      deg.step       = ang( row.step ),
      deg.als        = ang( row.als ),
      rec.step       = ( deg.start - ang( row.step ) ) / deg.start,
      rec.als        = ( deg.start - ang( row.als ) ) / deg.start,
      cos.step.als   = .diag2.vcos( u.step, u.als ),
      cos.als.E      = .diag2.vcos( u.als, e.true ),
      als.iter       = als$n.iter,
      cos.dir.min    = cos.dir.min,
      drift.in.delta = drift.in.delta,
      row.names      = NULL
    )
  }
  
  h.table <- do.call( rbind, h.rows )
  cat( sprintf( "-- %s --\n", target ) )
  print( h.table, digits = 4 )
  cat( "\n" )
  
  h.results[[ target ]] <- list( table = h.table, spectra.als = h.spectra.als )
}

cat( "  INTERPRETATION: rec.als above rec.step on the dyes the line search\n" )
cat( "  vetoed (t.hat = 0 in Stage D) confirms the gate, not the estimator,\n" )
cat( "  was the blocker. cos.step.als high on recovering dyes and low on\n" )
cat( "  diverging ones would validate concordance as a production gate.\n" )
cat( "  cos.dir.min well below 1 with drift.in.delta high confirms the tandem\n" )
cat( "  variant-mixture hypothesis: the residual direction drifts with\n" )
cat( "  brightness, and the drift lives in the dye's own variant span.\n" )

# ---------------------------------------------------------------------------
# STAGE I - restricted identifiability ceiling per dye
# ---------------------------------------------------------------------------
# The restricted residual is orthogonal to every active row, so the component
# of a dye's true error lying in the span of its nuisance rows is invisible
# to the restricted path, exactly as the full-panel path was blind to the
# full row space. For tandem pairs the nuisance span sits on top of the
# donor/acceptor rebalancing direction, so this ceiling is expected to bind
# hardest on precisely the dyes that under-recover.

cat( "\n============ STAGE I: restricted identifiability ceiling ============\n" )

for ( target in c( "Cells", "Beads" ) ) {
  
  gres   <- g.results[[ target ]]
  pop    <- gres$pop
  s.true <- gres$s.true
  
  i.rows <- list()
  
  for ( j in names( gres$store ) ) {
    
    st       <- gres$store[[ j ]]
    idx      <- st$idx
    e.true   <- s.true[ j, , drop = FALSE ] - pop$ctx$spectra[ j, , drop = FALSE ]
    
    nuis.flat <- st$nuisance
    nuis.spread <- if ( !is.null( pop$above.spread ) )
      .diag2.nuisance( pop, idx, j, pop$above.spread ) else NULL
    
    basis.flat <- if ( length( nuis.flat ) > 0 )
      pop$ctx$spectra[ nuis.flat, , drop = FALSE ] else NULL
    basis.spread <- if ( !is.null( nuis.spread ) && length( nuis.spread ) > 0 )
      pop$ctx$spectra[ nuis.spread, , drop = FALSE ] else NULL
    
    i.rows[[ j ]] <- data.frame(
      target           = target,
      fluorophore      = j,
      err.norm         = sqrt( sum( e.true^2 ) ),
      frac.in.own.row  = .diag2.span.fraction(
        e.true, pop$ctx$spectra[ j, , drop = FALSE ] ),
      frac.in.nuis     = .diag2.span.fraction( e.true, basis.flat ),
      frac.in.nuis.spr = if ( is.null( pop$above.spread ) ) NA_real_ else
        .diag2.span.fraction( e.true, basis.spread ),
      n.nuis.flat      = length( nuis.flat ),
      n.nuis.spread    = if ( is.null( nuis.spread ) ) NA_integer_ else
        length( nuis.spread ),
      nuisance.flat    = paste( nuis.flat, collapse = ", " ),
      row.names        = NULL
    )
  }
  
  i.table <- do.call( rbind, i.rows )
  cat( sprintf( "-- %s --\n", target ) )
  print( i.table, digits = 3 )
  cat( "\n" )
}

cat( "  INTERPRETATION: frac.in.nuis is the share of the true error the\n" )
cat( "  restricted design structurally cannot recover for that dye; expect it\n" )
cat( "  high for the tandems that stall. frac.in.nuis.spr below frac.in.nuis\n" )
cat( "  quantifies how much recoverability the spread-scaled co-activity mask\n" )
cat( "  buys by removing phantom nuisances. frac.in.own.row is scale error\n" )
cat( "  that renormalisation absorbs and can be ignored.\n" )

# ---------------------------------------------------------------------------
# STAGE J - max.step ceiling sweep
# ---------------------------------------------------------------------------
# CONTEXT_A10_A11_session_summary.md section 7 flagged, and explicitly
# deferred, that correct.unmixing.signatures()'s max.step = 0.08 ceiling
# "almost certainly has the same structural problem as fix.my.unmix()'s
# max.coefficient" - a single fixed cap on a quantity whose plausible
# magnitude scales with how oblique/collinear a dye's error is. That cap was
# raised 0.2 -> 0.5 on directly analogous fix.my.unmix() evidence and is now
# the settled default everywhere it is used. This has never been tested on
# the max.step side. The restricted-design rescue that inherits this same
# cap is net-positive but noisier (71% would.help, n=17) than the primary
# arm - consistent with, not proof of, some capped-but-real corrections
# being turned away. This stage tests it directly, on the same benchmark
# Stage D already scored against ground truth.
#
# Unlike fix.my.unmix()'s angle gate, this loop has no partial-step
# (clamped) fallback: a step is either accepted whole or rejected outright.
# Raising max.step therefore lets a large step through unmodified, not
# scaled down - watch for regressions, not just gains.

cat( "\n============ STAGE J: max.step ceiling sweep ============\n" )

j.step.grid <- c( 0.08, 0.15, 0.25, 0.5 )

j.results <- list()

for ( target in c( "Cells", "Beads" ) ) {
  
  ctx    <- c.results[[ target ]]
  s.true <- as.matrix( diag.spectra[[ target ]][ rownames( ctx$spectra ),
                                                 b.shared.det, drop = FALSE ] )
  
  j.rows <- list()
  
  for ( ms in j.step.grid ) {
    
    fit <- .diag.dominance.fit(
      raw.data          = ctx$raw,
      spectra           = ctx$spectra,
      thresholds        = ctx$thresholds,
      panel             = ctx$panel,
      af.name           = af.name,
      spillover.spread  = ctx$ss.mat,
      spread.kappa      = 2,
      max.step          = ms,
      bg.mode           = if ( target == "Cells" ) "scatter.knn" else "global.mean",
      scatter           = ctx$scatter,
      unstained         = ctx$unstained,
      unstained.scatter = ctx$unstained.scatter
    )
    
    fit.table <- do.call( rbind, fit$fit.log )
    
    j.rows[[ length( j.rows ) + 1L ]] <- data.frame(
      target      = target,
      max.step    = ms,
      fluorophore = ctx$panel,
      deg.start   = .diag.angle( ctx$spectra[ ctx$panel, , drop = FALSE ],
                                 s.true[ ctx$panel, , drop = FALSE ] ),
      deg.after   = .diag.angle( fit$spectra[ ctx$panel, , drop = FALSE ],
                                 s.true[ ctx$panel, , drop = FALSE ] ),
      recovered   = .diag.recovered.fraction(
        ctx$spectra[ ctx$panel, , drop = FALSE ],
        fit$spectra[ ctx$panel, , drop = FALSE ],
        s.true[ ctx$panel, , drop = FALSE ] ),
      n.steps     = vapply( ctx$panel, function( j )
        sum( fit.table$accepted[ fit.table$fluorophore == j ] ), integer( 1 ) ),
      max.rel.step.seen = vapply( ctx$panel, function( j ) {
        v <- fit.table$rel.step[ fit.table$fluorophore == j ]
        if ( length( v ) == 0 ) NA_real_ else max( v, na.rm = TRUE )
      }, numeric( 1 ) ),
      row.names   = NULL
    )
  }
  
  j.table <- do.call( rbind, j.rows )
  cat( sprintf( "-- %s --\n", target ) )
  print( j.table, digits = 3 )
  cat( "\n" )
  
  agg <- stats::aggregate( recovered ~ max.step, j.table, mean, na.rm = TRUE )
  cat( "  mean recovered fraction by max.step:\n" )
  print( agg )
  cat( "\n" )
  
  j.results[[ target ]] <- j.table
}

cat( "  INTERPRETATION: compare each row's recovered at max.step > 0.08\n" )
cat( "  against Stage D's own max.step = 0.08 result (d.results). A dye whose\n" )
cat( "  recovered fraction improves as max.step rises AND whose\n" )
cat( "  max.rel.step.seen at 0.08 sat at or near 0.08 (the row was being\n" )
cat( "  capped, not genuinely small) is the direct analogue of the\n" )
cat( "  max.coefficient finding - a real, identifiable correction the\n" )
cat( "  ceiling was blocking outright. A dye whose recovered fraction\n" )
cat( "  WORSENS as max.step rises means the fixed step is not the primary\n" )
cat( "  constraint here, and the priority shifts to porting a partial-step\n" )
cat( "  fallback (like fix.my.unmix()'s angle clamp) before raising the\n" )
cat( "  ceiling blind.\n" )

# ---------------------------------------------------------------------------
# STAGE K (corrected) - held-out split stability
# ---------------------------------------------------------------------------
# Two bugs in the first version of this stage, found by comparing its
# output against Stage D/J's actual per-fit log rather than trusting the
# pass-rate summary on its own:
#
#   1. `above` (used to decide which co-active fluorophores are carried as
#      nuisance terms) was computed from FLAT thresholds. Stage D/J's real
#      calls all pass `spillover.spread = ctx$ss.mat`, so production always
#      uses get.spread.thresholds()'s WIDENED thresholds for this decision.
#      For Beads specifically, three fluorophores (BV711, RB705,
#      PerCP-eFluor 710) show n.nuisance = 1 in Stage D's own iteration-1
#      log; the first version of this stage could easily have assembled a
#      different restricted design for them.
#   2. For Cells, .diag.bg.subtract() subsamples its unstained reference
#      pool down to max.reference = 50000 whenever the pool exceeds that
#      (it does here: ~102,000 events after gating), via an un-reseeded
#      sample() call. The first version of this stage called
#      .diag.bg.subtract() once, fresh, right after its own
#      set.seed(diag.seed) - a different reference draw than whatever the
#      original Stage D run happened to use, applied identically across
#      all 40 "trials" for every Cells fluorophore. That has nothing to do
#      with the held-out split and would explain a whole-target, largely
#      unanimous shift far more readily than resampling noise would -
#      exactly the pattern seen on PerCP-eFluor 710 (t.hat = 1.0 in
#      Stage D, 0/40 here).
#
# Neither bug explains Beads' widespread flip on its own - that one is
# still open. This version fixes both known issues and reports n.nuisance
# and gain quantiles, not just a pass count, so a real small-N instability
# can be told apart from a further hidden difference from Stage D's setup.

cat( "\n============ STAGE K (corrected): held-out split stability ============\n" )

k.targets <- list(
  Cells = c( "Spark UV 387", "BUV395", "BUV805", "eFluor 450", "FITC",
             "PerCP-eFluor 710", "BV711" ),
  Beads = c( "Spark UV 387", "BUV395", "eFluor 450", "BV711", "RB705",
             "BUV805", "Spark Violet 538", "PE" )
)
k.n.splits <- 40L

set.seed( diag.seed )
k.results <- list()

for ( target in names( k.targets ) ) {
  
  ctx <- c.results[[ target ]]
  dom <- d.results[[ target ]]$dominant
  
  above <- if ( is.null( ctx$ss.mat ) ) {
    sweep( ctx$unmixed[ , ctx$panel, drop = FALSE ], 2,
           ctx$thresholds[ ctx$panel ], ">" )
  } else {
    threshold.matrix <- get.spread.thresholds(
      unmixed = ctx$unmixed, thresholds = ctx$thresholds,
      spillover.spread = ctx$ss.mat, spread.kappa = 2, verbose = FALSE )
    ctx$unmixed[ , ctx$panel, drop = FALSE ] > threshold.matrix[ , ctx$panel, drop = FALSE ]
  }
  
  background.idx.full <- which( dom == 0L )
  
  raw.bg <- if ( target == "Cells" )
    .diag.bg.subtract( ctx$raw, ctx$scatter, ctx$unstained, ctx$unstained.scatter,
                       max.reference = 200000L ) else
                         sweep( ctx$raw, 2,
                                colMeans( ctx$raw[ background.idx.full, , drop = FALSE ] ), "-" )
  
  background.idx <- background.idx.full
  if ( length( background.idx ) > 5000L )
    background.idx <- sample( background.idx, 5000L )
  
  k.rows <- list()
  
  for ( j in k.targets[[ target ]] ) {
    
    f   <- match( j, ctx$panel )
    idx <- which( dom == f )
    if ( length( idx ) < 200L ) next
    
    co.frac  <- colMeans( above[ idx, , drop = FALSE ] )
    nuisance <- setdiff( ctx$panel[ co.frac > 0.5 ], j )
    active   <- c( j, nuisance )
    
    trials <- t( vapply( seq_len( k.n.splits ), function( k ) {
      half <- sample( rep_len( c( 1L, 2L ), length( idx ) ) )
      .diag.split.trial( half, idx, raw.bg, ctx$spectra, active, j,
                         background.idx )
    }, numeric( 2 ) ) )
    colnames( trials ) <- c( "t", "gain" )
    
    k.rows[[ j ]] <- data.frame(
      target        = target,
      fluorophore   = j,
      n.nuisance    = length( nuisance ),
      n.splits.pass = sum( trials[ , "t" ] > 0, na.rm = TRUE ),
      n.splits      = k.n.splits,
      median.t      = stats::median( trials[ , "t" ] ),
      gain.q25      = unname( stats::quantile( trials[ , "gain" ], 0.25, na.rm = TRUE ) ),
      gain.median   = stats::median( trials[ , "gain" ], na.rm = TRUE ),
      gain.q75      = unname( stats::quantile( trials[ , "gain" ], 0.75, na.rm = TRUE ) ),
      row.names     = NULL
    )
  }
  
  k.table <- do.call( rbind, k.rows )
  cat( sprintf( "-- %s --\n", target ) )
  print( k.table, digits = 4 )
  cat( "\n" )
  
  k.results[[ target ]] <- k.table
}

cat( "  INTERPRETATION: check n.nuisance here against Stage D's iter=1\n" )
cat( "  fit.log for the same fluorophore - a mismatch means the restricted\n" )
cat( "  design itself differs from what Stage D evaluated, which explains a\n" )
cat( "  result on its own before any conclusion about split-noise or data\n" )
cat( "  volume is drawn. Where n.nuisance matches and the pass rate is still\n" )
cat( "  near-unanimous (0/40 or 40/40) rather than a mixed rate like\n" )
cat( "  BUV805's 19/40, that argues for a further, still-unidentified\n" )
cat( "  difference between this stage and Stage D's actual run rather than\n" )
cat( "  genuine split-to-split noise, which should produce intermediate\n" )
cat( "  pass rates, not unanimous ones, for a fluorophore near the boundary.\n" )
# ---------------------------------------------------------------------------
# STAGE L - repeated-split estimator vs. ground truth
# ---------------------------------------------------------------------------
cat( "\n============ STAGE L: repeated-split vs. ground truth ============\n" )

for ( target in c( "Cells", "Beads" ) ) {
  
  ctx    <- c.results[[ target ]]
  s.true <- as.matrix( diag.spectra[[ target ]][ rownames( ctx$spectra ),
                                                 b.shared.det, drop = FALSE ] )
  
  fit.new <- .diag.dominance.fit(
    raw.data          = ctx$raw,
    spectra           = ctx$spectra,
    thresholds        = ctx$thresholds,
    panel             = ctx$panel,
    af.name           = af.name,
    spillover.spread  = ctx$ss.mat,
    spread.kappa      = 2,
    max.step          = 0.15,
    n.split.trials    = 15L,
    min.split.frac    = 0.6,
    bg.mode           = if ( target == "Cells" ) "scatter.knn" else "global.mean",
    scatter           = ctx$scatter,
    unstained         = ctx$unstained,
    unstained.scatter = ctx$unstained.scatter
  )
  
  l.table <- data.frame(
    target      = target,
    fluorophore = ctx$panel,
    deg.start   = .diag.angle( ctx$spectra[ ctx$panel, , drop = FALSE ],
                               s.true[ ctx$panel, , drop = FALSE ] ),
    deg.new     = .diag.angle( fit.new$spectra[ ctx$panel, , drop = FALSE ],
                               s.true[ ctx$panel, , drop = FALSE ] ),
    recovered   = .diag.recovered.fraction(
      ctx$spectra[ ctx$panel, , drop = FALSE ],
      fit.new$spectra[ ctx$panel, , drop = FALSE ],
      s.true[ ctx$panel, , drop = FALSE ] ),
    row.names   = NULL
  )
  
  cat( sprintf( "-- %s (n.split.trials = 15, min.split.frac = 0.6) --\n", target ) )
  print( l.table, digits = 4 )
  cat( "\n" )
}

cat( "  INTERPRETATION: compare deg.new / recovered here against Stage D's\n" )
cat( "  deg.after / recovered (n.split.trials = 1 implicitly). Improvement on\n" )
cat( "  the five stuck Cells fluorophores without regression elsewhere is the\n" )
cat( "  result that justifies porting n.split.trials/min.split.frac into\n" )
cat( "  correct_unmixing_signatures.R (section 5c above). Watch NovaFluor\n" )
cat( "  Blue 610-30S and PerCP specifically for regression, since they're the\n" )
cat( "  two fluorophores already known to be fragile at this benchmark.\n" )

# ---------------------------------------------------------------------------
# STAGE M - direction check on Beads' reproducibly-accepted-but-wrong
#   corrections
# ---------------------------------------------------------------------------
# Stage K/L together rule out split noise for BUV805 and Spark Violet 538
# (Beads): the held-out step search accepts a step on essentially every one
# of 40 random splits, with a tight, reproducible, above-threshold gain
# (~0.0101, ~0.0094) - yet the accepted correction moves the row AWAY from
# ground truth (recovered -0.10, -1.57), unchanged by n.split.trials = 15.
# The held-out objective only checks whether restricted-design residual
# norm goes down on held-out data; it says nothing about whether the
# fitted slope points toward this dye's true error or toward something
# else that happens to reduce residual norm on this specific restricted
# population. This checks what the fitted slope is actually aligned with.

cat( "\n============ STAGE M: fitted-direction check ============\n" )

m.targets <- list( Beads = c( "BUV805", "Spark Violet 538", "BUV395" ) )
# BUV395 is a contrast case: also Beads, also low-gain, but REJECTED
# (gain never clears min.gain) rather than accepted-wrong - useful to see
# whether its slope direction looks any different from the two problem
# dyes.

for ( target in names( m.targets ) ) {
  
  ctx    <- c.results[[ target ]]
  dom    <- d.results[[ target ]]$dominant
  s.true <- as.matrix( diag.spectra[[ target ]][ rownames( ctx$spectra ),
                                                 colnames( ctx$spectra ),
                                                 drop = FALSE ] )
  
  above <- if ( is.null( ctx$ss.mat ) ) {
    sweep( ctx$unmixed[ , ctx$panel, drop = FALSE ], 2,
           ctx$thresholds[ ctx$panel ], ">" )
  } else {
    threshold.matrix <- get.spread.thresholds(
      unmixed = ctx$unmixed, thresholds = ctx$thresholds,
      spillover.spread = ctx$ss.mat, spread.kappa = 2, verbose = FALSE )
    ctx$unmixed[ , ctx$panel, drop = FALSE ] > threshold.matrix[ , ctx$panel, drop = FALSE ]
  }
  
  background.idx <- which( dom == 0L )
  raw.bg <- sweep( ctx$raw, 2,
                   colMeans( ctx$raw[ background.idx, , drop = FALSE ] ), "-" )
  if ( length( background.idx ) > 5000L )
    background.idx <- sample( background.idx, 5000L )
  
  m.rows <- list()
  
  for ( j in m.targets[[ target ]] ) {
    
    f   <- match( j, ctx$panel )
    idx <- which( dom == f )
    
    co.frac  <- colMeans( above[ idx, , drop = FALSE ] )
    nuisance <- setdiff( ctx$panel[ co.frac > 0.5 ], j )
    active   <- c( j, nuisance )
    
    slope.hat <- local( {
      y.use <- raw.bg[ idx, , drop = FALSE ]
      x.use <- .diag.restricted.unmix( y.use, ctx$spectra, active )[ , 1 ]
      brk <- unique( stats::quantile( x.use, probs = seq( 0, 1, length.out = 11 ) ) )
      bin <- as.integer( cut( x.use, breaks = brk, include.lowest = TRUE ) )
      idx.by.bin <- split( seq_along( bin ), bin )
      y.bin <- t( vapply( idx.by.bin, function( ii )
        colMeans( y.use[ ii, , drop = FALSE ] ), numeric( ncol( y.use ) ) ) )
      y.bin <- rbind( colMeans( raw.bg[ background.idx, , drop = FALSE ] ), y.bin )
      x.bin <- .diag.restricted.unmix( y.bin, ctx$spectra, active )
      r.bin <- y.bin - x.bin %*% ctx$spectra[ active, , drop = FALSE ]
      fit   <- stats::lm.fit( x = cbind( 1, x.bin ), y = r.bin )
      slope <- stats::coef( fit )[ 2, ]
      slope[ !is.finite( slope ) ] <- 0
      slope
    } )
    
    slope.unit  <- slope.hat / sqrt( sum( slope.hat^2 ) )
    e.true      <- s.true[ j, ] - ctx$spectra[ j, ]
    e.true.unit <- e.true / sqrt( sum( e.true^2 ) )
    mean.bright <- colMeans( ctx$spectra[ ctx$panel, , drop = FALSE ] )
    mean.unit   <- mean.bright / sqrt( sum( mean.bright^2 ) )
    
    cos.to.other <- vapply( setdiff( ctx$panel, j ), function( k ) {
      row.k <- ctx$spectra[ k, ]
      sum( slope.unit * row.k ) / sqrt( sum( row.k^2 ) )
    }, numeric( 1 ) )
    best.other <- names( which.max( abs( cos.to.other ) ) )
    
    m.rows[[ j ]] <- data.frame(
      target                = target,
      fluorophore           = j,
      cos.slope.true        = sum( slope.unit * e.true.unit ),
      cos.slope.meanbright  = sum( slope.unit * mean.unit ),
      best.other.match      = best.other,
      cos.slope.best.other  = unname( cos.to.other[ best.other ] ),
      row.names             = NULL
    )
  }
  
  m.table <- do.call( rbind, m.rows )
  cat( sprintf( "-- %s --\n", target ) )
  print( m.table, digits = 3 )
  cat( "\n" )
}

cat( "  INTERPRETATION: cos.slope.true near +1 means the fitted correction\n" )
cat( "  really does point toward this dye's true error. A low or negative\n" )
cat( "  value on BUV805/Spark Violet 538, alongside their reproducible\n" )
cat( "  held-out gain, would show the objective is being satisfied by\n" )
cat( "  something other than the real spectral error. cos.slope.meanbright\n" )
cat( "  near +-1 flags a generic scale/brightness confound; a high\n" )
cat( "  cos.slope.best.other flags an unrecognised co-varying partner that\n" )
cat( "  nuisance.frac = 0.5 didn't catch as a nuisance term. BUV395 is the\n" )
cat( "  contrast case - if its cos.slope.true also comes out low, that's\n" )
cat( "  evidence the direction problem isn't specific to the two dyes that\n" )
cat( "  happened to clear min.gain, and cuts across the whole restricted-\n" )
cat( "  design method on this substrate.\n" )

# ---------------------------------------------------------------------------
# STAGE N - background-anchor and high-abundance leverage check
# ---------------------------------------------------------------------------
# Two candidate explanations for BUV805/Spark Violet 538's reproducible,
# wrong-direction slope on Beads, neither yet tested: (1) the single
# background-anchor bin has outsized leverage in a ~10-11-bin regression,
# and (2) a small number of high-abundance events (plausible aggregate/
# doublet contamination, which beads are more prone to than cells) are
# driving the trend. Refits the same slope three ways - baseline (as
# Stage M), without the background anchor, and after trimming the
# brightest 2% of events - and reports cos.slope.true for each.

cat( "\n============ STAGE N: anchor and leverage check ============\n" )

n.targets <- list( Beads = c( "BUV805", "Spark Violet 538", "BUV395" ) )

.diag.slope.variant <- function( idx, raw.data, spectra, active, j,
                                 background.idx, use.anchor = TRUE,
                                 trim.top = 0 ) {
  
  y.use <- raw.data[ idx, , drop = FALSE ]
  x.use <- .diag.restricted.unmix( y.use, spectra, active )[ , 1 ]
  
  if ( trim.top > 0 ) {
    keep  <- x.use <= stats::quantile( x.use, 1 - trim.top )
    y.use <- y.use[ keep, , drop = FALSE ]
    x.use <- x.use[ keep ]
  }
  
  brk <- unique( stats::quantile( x.use, probs = seq( 0, 1, length.out = 11 ) ) )
  bin <- as.integer( cut( x.use, breaks = brk, include.lowest = TRUE ) )
  idx.by.bin <- split( seq_along( bin ), bin )
  y.bin <- t( vapply( idx.by.bin, function( ii )
    colMeans( y.use[ ii, , drop = FALSE ] ), numeric( ncol( y.use ) ) ) )
  
  if ( use.anchor )
    y.bin <- rbind( colMeans( raw.data[ background.idx, , drop = FALSE ] ), y.bin )
  
  x.bin <- .diag.restricted.unmix( y.bin, spectra, active )
  r.bin <- y.bin - x.bin %*% spectra[ active, , drop = FALSE ]
  fit   <- stats::lm.fit( x = cbind( 1, x.bin ), y = r.bin )
  slope <- stats::coef( fit )[ 2, ]
  slope[ !is.finite( slope ) ] <- 0
  slope
}

for ( target in names( n.targets ) ) {
  
  ctx    <- c.results[[ target ]]
  dom    <- d.results[[ target ]]$dominant
  s.true <- as.matrix( diag.spectra[[ target ]][ rownames( ctx$spectra ),
                                                 colnames( ctx$spectra ),
                                                 drop = FALSE ] )
  
  above <- if ( is.null( ctx$ss.mat ) ) {
    sweep( ctx$unmixed[ , ctx$panel, drop = FALSE ], 2,
           ctx$thresholds[ ctx$panel ], ">" )
  } else {
    threshold.matrix <- get.spread.thresholds(
      unmixed = ctx$unmixed, thresholds = ctx$thresholds,
      spillover.spread = ctx$ss.mat, spread.kappa = 2, verbose = FALSE )
    ctx$unmixed[ , ctx$panel, drop = FALSE ] > threshold.matrix[ , ctx$panel, drop = FALSE ]
  }
  
  background.idx <- which( dom == 0L )
  raw.bg <- sweep( ctx$raw, 2,
                   colMeans( ctx$raw[ background.idx, , drop = FALSE ] ), "-" )
  if ( length( background.idx ) > 5000L )
    background.idx <- sample( background.idx, 5000L )
  
  n.rows <- list()
  
  for ( j in n.targets[[ target ]] ) {
    
    f   <- match( j, ctx$panel )
    idx <- which( dom == f )
    
    co.frac  <- colMeans( above[ idx, , drop = FALSE ] )
    nuisance <- setdiff( ctx$panel[ co.frac > 0.5 ], j )
    active   <- c( j, nuisance )
    
    e.true.unit <- {
      e <- s.true[ j, ] - ctx$spectra[ j, ]
      e / sqrt( sum( e^2 ) )
    }
    
    cos.of <- function( slope )
      sum( ( slope / sqrt( sum( slope^2 ) ) ) * e.true.unit )
    
    s.base    <- .diag.slope.variant( idx, raw.bg, ctx$spectra, active, j,
                                      background.idx, use.anchor = TRUE,  trim.top = 0 )
    s.noanchor <- .diag.slope.variant( idx, raw.bg, ctx$spectra, active, j,
                                       background.idx, use.anchor = FALSE, trim.top = 0 )
    s.trim    <- .diag.slope.variant( idx, raw.bg, ctx$spectra, active, j,
                                      background.idx, use.anchor = TRUE,  trim.top = 0.02 )
    
    n.rows[[ j ]] <- data.frame(
      target             = target,
      fluorophore        = j,
      cos.baseline       = cos.of( s.base ),
      cos.no.anchor      = cos.of( s.noanchor ),
      cos.trim.top2pct   = cos.of( s.trim ),
      row.names          = NULL
    )
  }
  
  n.table <- do.call( rbind, n.rows )
  cat( sprintf( "-- %s --\n", target ) )
  print( n.table, digits = 3 )
  cat( "\n" )
}

cat( "  INTERPRETATION: cos.no.anchor far from cos.baseline implicates the\n" )
cat( "  background bin's leverage. cos.trim.top2pct far from cos.baseline\n" )
cat( "  implicates a handful of high-abundance (possible aggregate) events.\n" )
cat( "  If neither variant moves cos.baseline meaningfully for Spark Violet\n" )
cat( "  538, the anti-alignment is a property of the bulk of the population,\n" )
cat( "  not a leverage artifact, and the next step is a direct look at that\n" )
cat( "  dye's bead vs. cell reference spectra rather than at this fitting\n" )
cat( "  procedure.\n" )

# ---------------------------------------------------------------------------
# STAGE O - substrate-specific AF, in-memory test
# ---------------------------------------------------------------------------
# bd.spectra["AF",] and cell.spectra["AF",] are identical (confirmed: their
# difference is exactly zero across every detector), meaning at least one
# substrate's AF row is not actually that substrate's own autofluorescence.
# AF never enters the per-dye restricted regression directly (`active`
# excludes af.name by construction) - but it DOES enter the full-panel
# unmix used to decide which events are dominance-assigned to each dye and
# which are background. Spark Violet 538 is flagged as strongly collinear
# with AF; if the AF row used for Beads is wrong, that collinear pair's
# full-panel unmix is exactly where population assignment goes unstable -
# a population-wide contamination, consistent with Stage N finding no
# leverage- or tail-driven explanation. This derives Beads' own AF
# direction from Beads' own unstained control (top singular vector of the
# raw, uncentred matrix - the mean background direction, same technique
# already used for the Stage D hotspot check) and reruns dominance
# assignment and the Stage M direction check with only that one row
# changed, to see whether the corruption in Spark Violet 538's assigned
# population - and its slope direction - clears up.

cat( "\n============ STAGE O: substrate-specific AF, in-memory test ============\n" )

o.targets <- c( "Spark Violet 538", "BUV805", "BUV395" )

for ( target in "Beads" ) {
  
  ctx    <- c.results[[ target ]]
  s.true <- as.matrix( diag.spectra[[ target ]][ rownames( ctx$spectra ),
                                                 colnames( ctx$spectra ),
                                                 drop = FALSE ] )
  
  af.own <- .diag.af.pcs( ctx$unstained, n.pc = 1L )[ 1, ]
  af.own <- af.own / max( af.own )   # match the row's own L-infinity convention
  
  cos.shared.af <- vapply( o.targets, function( j )
    sum( ctx$spectra[ "AF", ] * ctx$spectra[ j, ] ) /
      ( sqrt( sum( ctx$spectra[ "AF", ]^2 ) ) * sqrt( sum( ctx$spectra[ j, ]^2 ) ) ),
    numeric( 1 ) )
  
  cos.own.af <- vapply( o.targets, function( j )
    sum( af.own * ctx$spectra[ j, ] ) /
      ( sqrt( sum( af.own^2 ) ) * sqrt( sum( ctx$spectra[ j, ]^2 ) ) ),
    numeric( 1 ) )
  
  cat( "-- AF/dye collinearity: shared (wrong) AF vs. Beads' own unstained-derived AF --\n" )
  print( data.frame( fluorophore = o.targets, cos.shared.af = cos.shared.af,
                     cos.own.af = cos.own.af, row.names = NULL ), digits = 3 )
  cat( "\n" )
  
  spectra.fix <- ctx$spectra
  spectra.fix[ "AF", ] <- af.own
  
  unst.um.fix <- unmix.ols.fast( ctx$unstained, spectra.fix )
  thr.fix     <- apply( unst.um.fix, 2, stats::quantile, probs = diag.threshold.probs )
  
  unmixed.fix <- unmix.ols.fast( ctx$raw, spectra.fix )
  
  above.fix <- if ( is.null( ctx$ss.mat ) ) {
    sweep( unmixed.fix[ , ctx$panel, drop = FALSE ], 2, thr.fix[ ctx$panel ], ">" )
  } else {
    threshold.matrix <- get.spread.thresholds(
      unmixed = unmixed.fix, thresholds = thr.fix,
      spillover.spread = ctx$ss.mat, spread.kappa = 2, verbose = FALSE )
    unmixed.fix[ , ctx$panel, drop = FALSE ] > threshold.matrix[ , ctx$panel, drop = FALSE ]
  }
  
  excess    <- pmax( sweep( unmixed.fix[ , ctx$panel, drop = FALSE ], 2,
                            thr.fix[ ctx$panel ], "-" ), 0 )
  dyn.range <- apply( unmixed.fix[ , ctx$panel, drop = FALSE ], 2, stats::quantile,
                      probs = 0.999 ) - thr.fix[ ctx$panel ]
  dyn.range <- pmax( dyn.range, .Machine$double.eps )
  score     <- sweep( excess, 2, dyn.range, "/" )
  dom.fix   <- max.col( score, ties.method = "first" )
  top       <- score[ cbind( seq_len( nrow( score ) ), dom.fix ) ]
  dom.fix[ top <= 0 ] <- 0L
  
  n.compare <- data.frame(
    fluorophore  = o.targets,
    n.old        = vapply( o.targets, function( j )
      sum( d.results[[ target ]]$dominant == match( j, ctx$panel ) ), integer( 1 ) ),
    n.new        = vapply( o.targets, function( j )
      sum( dom.fix == match( j, ctx$panel ) ), integer( 1 ) ),
    row.names    = NULL
  )
  cat( "-- dominance population size, old (shared AF) vs. new (Beads' own AF) --\n" )
  print( n.compare )
  cat( "\n" )
  
  background.idx <- which( dom.fix == 0L )
  raw.bg <- sweep( ctx$raw, 2,
                   colMeans( ctx$raw[ background.idx, , drop = FALSE ] ), "-" )
  if ( length( background.idx ) > 5000L )
    background.idx <- sample( background.idx, 5000L )
  
  o.rows <- list()
  
  for ( j in o.targets ) {
    
    f   <- match( j, ctx$panel )
    idx <- which( dom.fix == f )
    if ( length( idx ) < 200L ) { o.rows[[ j ]] <- NULL; next }
    
    co.frac  <- colMeans( above.fix[ idx, , drop = FALSE ] )
    nuisance <- setdiff( ctx$panel[ co.frac > 0.5 ], j )
    active   <- c( j, nuisance )
    
    slope.hat <- local( {
      y.use <- raw.bg[ idx, , drop = FALSE ]
      x.use <- .diag.restricted.unmix( y.use, spectra.fix, active )[ , 1 ]
      brk <- unique( stats::quantile( x.use, probs = seq( 0, 1, length.out = 11 ) ) )
      bin <- as.integer( cut( x.use, breaks = brk, include.lowest = TRUE ) )
      idx.by.bin <- split( seq_along( bin ), bin )
      y.bin <- t( vapply( idx.by.bin, function( ii )
        colMeans( y.use[ ii, , drop = FALSE ] ), numeric( ncol( y.use ) ) ) )
      y.bin <- rbind( colMeans( raw.bg[ background.idx, , drop = FALSE ] ), y.bin )
      x.bin <- .diag.restricted.unmix( y.bin, spectra.fix, active )
      r.bin <- y.bin - x.bin %*% spectra.fix[ active, , drop = FALSE ]
      fit   <- stats::lm.fit( x = cbind( 1, x.bin ), y = r.bin )
      slope <- stats::coef( fit )[ 2, ]
      slope[ !is.finite( slope ) ] <- 0
      slope
    } )
    
    e.true      <- s.true[ j, ] - ctx$spectra[ j, ]
    e.true.unit <- e.true / sqrt( sum( e.true^2 ) )
    slope.unit  <- slope.hat / sqrt( sum( slope.hat^2 ) )
    
    o.rows[[ j ]] <- data.frame(
      fluorophore     = j,
      n.nuisance.new  = length( nuisance ),
      cos.slope.old   = c( "Spark Violet 538" = -0.824, "BUV805" = 0.350,
                           "BUV395" = 0.352 )[ j ],
      cos.slope.new   = sum( slope.unit * e.true.unit ),
      row.names       = NULL
    )
  }
  
  o.table <- do.call( rbind, o.rows )
  cat( "-- fitted-direction check, with Beads' own AF in place of the shared row --\n" )
  print( o.table, digits = 3 )
  cat( "\n" )
}

cat( "  INTERPRETATION: cos.shared.af high for Spark Violet 538 confirms the\n" )
cat( "  collinearity claim quantitatively. n.new far from n.old means the\n" )
cat( "  dominance population itself was contaminated under the shared AF, not\n" )
cat( "  just the regression on top of it - matching Stage N's population-wide\n" )
cat( "  (not leverage-driven) result. cos.slope.new moving from cos.slope.old\n" )
cat( "  toward +1 for Spark Violet 538 confirms the wrong AF, not the fitting\n" )
cat( "  procedure, was the cause. If cos.slope.new is still far from +1 even\n" )
cat( "  with Beads' own AF, that points to a real, physical spectral\n" )
cat( "  closeness between this dye and bead autofluorescence rather than a\n" )
cat( "  data problem - worth checking against the hotspot-of-[AF;dye] metric\n" )
cat( "  the Stage D preamble already computes for Cells, extended to Beads.\n" )
