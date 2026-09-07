# run_pair_restricted_residual_validation.R
#
# Scores a fix.my.unmix() run using detector-space residuals from a design
# that excludes BOTH members of each candidate collinear pair.
#
# The full-design OLS residual cannot validate fix.my.unmix(): if a row's
# error lies in the row space of the spectra, the projector onto that row
# space is unchanged and the residual is identically unchanged. Removing one
# row is also insufficient, because a contaminating donor left in the design
# reabsorbs the error. Removing both members of the pair makes the error
# out-of-span and therefore visible.
#
# For each pair (j, k) this script builds the restricted design
#   A = panel \ {j, k}  (+ background basis, if supplied)
# regresses the restricted residual jointly on the two abundance proxies, and
# compares the fitted responses against what each candidate spectra matrix
# predicts. It also reports whether the pair is separable at all under that
# design, so an unidentifiable pair is labelled rather than scored.
#
# PREREQUISITES, in the session:
#   prv.raw.data        Numeric matrix, events x detectors. The gated raw
#                       fully stained sample fix.my.unmix() was run on.
#   prv.spectra.before  Numeric matrix, fluorophores x detectors. Input spectra.
#   prv.spectra.after   Numeric matrix, same rows/columns. fix$spectra.
#   prv.af.name         Character or NULL. Row name to exclude from the panel.
#
# OPTIONAL:
#   prv.background.basis  Numeric matrix, components x detectors. The AF basis
#                         used by the run; always kept in the restricted design.
#   prv.thresholds        Named numeric, per-fluorophore positivity thresholds.
#                         Defaults to a per-column quantile of the unmixed data.
#
# Run top to bottom.
prv.raw.data <- lung.stained.data[,colnames(spectra)]
prv.spectra.before <- spectra
prv.spectra.after <- fixed.sig.lung.nods$spectra
prv.af.name <- "AF"
prv.background.basis <- fixed.sig.lung.nods$af.basis
prv.thresholds <- fixed.sig.lung.nods$thresholds.final

# ---------------------------------------------------------------------------
# 0. Setup - EDIT THIS SECTION
# ---------------------------------------------------------------------------

prv.output.dir <- "./pair_restricted_residual_validation_before"

# How many pairs to score, taken in descending order of cosine similarity
# under the starting spectra. The whole panel is F(F-1)/2 pairs; the ones
# that matter are the collinear ones.
prv.n.pairs <- 40L

# Minimum events required in each of the two single-positive populations.
prv.min.events <- 300L

# Where the non-pair rows come from when scoring. "after" isolates the pair,
# so a before/after difference reflects only what happened to s_j and s_k.
# "before" scores the panel as a whole and will mix in every other accepted
# row change. Run both; they answer different questions.
prv.a.from <- "before"

# Number of independent 50/50 event splits used to put an error bar on each
# angle. The spread across splits is the noise floor a real improvement has
# to clear.
prv.n.splits <- 8L

prv.seed <- 42L

if ( !dir.exists( prv.output.dir ) )
  dir.create( prv.output.dir, recursive = TRUE )

set.seed( prv.seed )

# ---------------------------------------------------------------------------
# 1. Input checks
# ---------------------------------------------------------------------------

for ( obj in c( "prv.raw.data", "prv.spectra.before", "prv.spectra.after" ) )
  if ( !exists( obj ) )
    stop( "Missing required object: ", obj, call. = FALSE )

if ( !exists( "prv.af.name" ) )           prv.af.name <- "AF"
if ( !exists( "prv.background.basis" ) )  prv.background.basis <- NULL
if ( !exists( "prv.thresholds" ) )        prv.thresholds <- NULL

prv.detectors <- colnames( prv.spectra.before )

if ( !identical( colnames( prv.spectra.after ), prv.detectors ) )
  stop( "prv.spectra.before and prv.spectra.after have different detectors.",
        call. = FALSE )

if ( !all( prv.detectors %in% colnames( prv.raw.data ) ) )
  stop( "prv.raw.data is missing detectors present in the spectra: ",
        paste( setdiff( prv.detectors, colnames( prv.raw.data ) ),
               collapse = ", " ), call. = FALSE )

prv.raw <- as.matrix( prv.raw.data[ , prv.detectors, drop = FALSE ] )

prv.panel <- setdiff( rownames( prv.spectra.before ), prv.af.name )

prv.shared <- intersect( prv.panel, rownames( prv.spectra.after ) )
if ( length( prv.shared ) < length( prv.panel ) )
  stop( "Fluorophores present before but not after: ",
        paste( setdiff( prv.panel, prv.shared ), collapse = ", " ),
        call. = FALSE )

if ( !is.null( prv.background.basis ) )
  prv.background.basis <-
    as.matrix( prv.background.basis[ , prv.detectors, drop = FALSE ] )

# ---------------------------------------------------------------------------
# 2. Helpers
# ---------------------------------------------------------------------------

# Residual of y against a design, using the same solver the package uses.
# Falls back to a pseudoinverse if the restricted design is singular.
.prv.residual <- function( y, design ) {

  gram <- tcrossprod( design )

  unmixing <- tryCatch( solve.default( gram, design ),
                        error = function( e ) MASS::ginv( gram ) %*% design )

  y - ( y %*% t( unmixing ) ) %*% design
}

# The component of a single spectrum the restricted design cannot explain.
# This is what the restricted residual's response to that row's abundance
# should equal if the spectrum is correct.
.prv.orthogonalise <- function( s, design ) {
  as.vector( .prv.residual( matrix( s, nrow = 1L ), design ) )
}

.prv.angle <- function( a, b ) {

  na <- sqrt( sum( a^2 ) )
  nb <- sqrt( sum( b^2 ) )

  if ( !is.finite( na ) || !is.finite( nb ) || na <= 0 || nb <= 0 )
    return( NA_real_ )

  180 / pi * acos( pmin( 1, pmax( -1, sum( a * b ) / ( na * nb ) ) ) )
}

# Build the restricted design for a pair: every panel row except j and k,
# plus the background basis, which is always retained so that background
# structure is projected out rather than mistaken for pair signal.
.prv.design <- function( spectra, panel, j, k, background ) {

  keep <- setdiff( panel, c( j, k ) )
  out  <- spectra[ keep, , drop = FALSE ]

  if ( !is.null( background ) ) out <- rbind( background, out )

  out
}

# ---------------------------------------------------------------------------
# 3. Abundance proxies and populations, fixed once from the starting spectra
# ---------------------------------------------------------------------------
# Both the regressors and the event sets are taken from the BEFORE run and
# then held fixed. Recomputing them per candidate would change the question
# being asked between the two scorings.

prv.design.full <- if ( is.null( prv.background.basis ) )
  prv.spectra.before[ prv.panel, , drop = FALSE ] else
    rbind( prv.background.basis,
           prv.spectra.before[ prv.panel, , drop = FALSE ] )

prv.unmixed <- unmix.ols.fast( prv.raw, prv.design.full )
colnames( prv.unmixed ) <- rownames( prv.design.full )
prv.unmixed <- prv.unmixed[ , prv.panel, drop = FALSE ]

if ( is.null( prv.thresholds ) ) {
  prv.thresholds <- apply( prv.unmixed, 2, stats::quantile,
                           probs = 0.75, names = FALSE )
  names( prv.thresholds ) <- prv.panel
}

prv.high <- apply( prv.unmixed, 2, stats::quantile, probs = 0.999,
                   names = FALSE )
names( prv.high ) <- prv.panel

prv.above <- sweep( prv.unmixed, 2, prv.thresholds[ prv.panel ], ">" )

# ---------------------------------------------------------------------------
# 4. Candidate pairs, ranked by collinearity under the starting spectra
# ---------------------------------------------------------------------------

prv.cos <- cosine.similarity( prv.spectra.before[ prv.panel, , drop = FALSE ] )

prv.pair.grid <- expand.grid( j = prv.panel, k = prv.panel,
                              stringsAsFactors = FALSE )
prv.pair.grid <- prv.pair.grid[
  match( prv.pair.grid$j, prv.panel ) < match( prv.pair.grid$k, prv.panel ), ]

prv.pair.grid$cosine <- prv.cos[ cbind( prv.pair.grid$j, prv.pair.grid$k ) ]
prv.pair.grid <- prv.pair.grid[ order( -prv.pair.grid$cosine ), ]

prv.pair.grid <- utils::head( prv.pair.grid, prv.n.pairs )

# ---------------------------------------------------------------------------
# 5. Per-pair scoring
# ---------------------------------------------------------------------------
# For each pair, the two single-positive populations plus the events negative
# for both are pooled. The joint regression of the restricted residual on
# (x_j, x_k) with an intercept recovers each row's out-of-span response
# without either absorbing the other; the intercept takes the common-mode
# background the basis did not.

prv.spectra.a <- if ( identical( prv.a.from, "after" ) )
  prv.spectra.after else prv.spectra.before

prv.rows <- list()

for ( p in seq_len( nrow( prv.pair.grid ) ) ) {

  j <- prv.pair.grid$j[ p ]
  k <- prv.pair.grid$k[ p ]

  idx.j  <- which(  prv.above[ , j ] & !prv.above[ , k ] )
  idx.k  <- which( !prv.above[ , j ] &  prv.above[ , k ] )
  idx.bg <- which( !prv.above[ , j ] & !prv.above[ , k ] )

  if ( length( idx.j ) < prv.min.events || length( idx.k ) < prv.min.events ) {

    prv.rows[[ length( prv.rows ) + 1L ]] <- data.frame(
      j = j, k = k, cosine = prv.pair.grid$cosine[ p ],
      n.j = length( idx.j ), n.k = length( idx.k ),
      resid.dim = NA_integer_, separability = NA_real_,
      regressor.corr = NA_real_,
      deg.j.before = NA_real_, deg.j.after = NA_real_,
      deg.k.before = NA_real_, deg.k.after = NA_real_,
      deg.j.sd = NA_real_, deg.k.sd = NA_real_,
      status = "too.few.events",
      row.names = NULL, stringsAsFactors = FALSE )

    next
  }

  # Cap the double-negative pool so it cannot swamp the two populations that
  # carry the pair's signal.
  if ( length( idx.bg ) > 2L * ( length( idx.j ) + length( idx.k ) ) )
    idx.bg <- sample( idx.bg, 2L * ( length( idx.j ) + length( idx.k ) ) )

  idx <- c( idx.j, idx.k, idx.bg )

  design.a <- .prv.design( prv.spectra.a, prv.panel, j, k,
                           prv.background.basis )

  y     <- prv.raw[ idx, , drop = FALSE ]
  r.a   <- .prv.residual( y, design.a )
  x.reg <- cbind( 1, prv.unmixed[ idx, j ], prv.unmixed[ idx, k ] )

  # Model-predicted out-of-span responses, before and after.
  m.j.before <- .prv.orthogonalise( prv.spectra.before[ j, ], design.a )
  m.k.before <- .prv.orthogonalise( prv.spectra.before[ k, ], design.a )
  m.j.after  <- .prv.orthogonalise( prv.spectra.after[  j, ], design.a )
  m.k.after  <- .prv.orthogonalise( prv.spectra.after[  k, ], design.a )

  # Separability: how distinguishable the two rows' responses are inside the
  # restricted residual space. Near 1 means the design cannot tell which of
  # the two carries an error, and the fitted direction is free to flip sign.
  separability <- abs( cos( pi / 180 * .prv.angle( m.j.after, m.k.after ) ) )

  regressor.corr <- abs( stats::cor( prv.unmixed[ idx, j ],
                                     prv.unmixed[ idx, k ] ) )

  resid.dim <- ncol( prv.raw ) - nrow( design.a )

  # Split-half estimates of the observed responses, giving both a point
  # estimate (pooled) and a spread (across splits).
  fit.observed <- function( use ) {

    fit <- stats::lm.fit( x = x.reg[ use, , drop = FALSE ],
                          y = r.a[ use, , drop = FALSE ] )
    cf  <- stats::coef( fit )
    cf[ !is.finite( cf ) ] <- 0

    list( j = cf[ 2, ], k = cf[ 3, ] )
  }

  obs.all <- fit.observed( seq_along( idx ) )

  splits <- vapply( seq_len( prv.n.splits ), function( s ) {

    half <- sample( rep_len( c( TRUE, FALSE ), length( idx ) ) )
    ob   <- fit.observed( which( half ) )

    c( .prv.angle( ob$j, m.j.after ), .prv.angle( ob$k, m.k.after ) )

  }, numeric( 2 ) )

  prv.rows[[ length( prv.rows ) + 1L ]] <- data.frame(
    j              = j,
    k              = k,
    cosine         = prv.pair.grid$cosine[ p ],
    n.j            = length( idx.j ),
    n.k            = length( idx.k ),
    resid.dim      = resid.dim,
    separability   = separability,
    regressor.corr = regressor.corr,
    deg.j.before   = .prv.angle( obs.all$j, m.j.before ),
    deg.j.after    = .prv.angle( obs.all$j, m.j.after ),
    deg.k.before   = .prv.angle( obs.all$k, m.k.before ),
    deg.k.after    = .prv.angle( obs.all$k, m.k.after ),
    deg.j.sd       = stats::sd( splits[ 1, ], na.rm = TRUE ),
    deg.k.sd       = stats::sd( splits[ 2, ], na.rm = TRUE ),
    status         = "scored",
    row.names      = NULL, stringsAsFactors = FALSE )
}

prv.result <- do.call( rbind, prv.rows )

# ---------------------------------------------------------------------------
# 6. Verdicts
# ---------------------------------------------------------------------------
# A pair is only informative if the two rows are separable inside the
# restricted residual space and the two abundances are not themselves
# collinear across the scored events. Where either fails, the sign of the
# fitted response is not determined by the data and the pair is reported as
# unidentifiable rather than as a pass or a fail.

prv.result$identifiable <- with( prv.result,
  status == "scored" & is.finite( separability ) & separability < 0.98 &
    is.finite( regressor.corr ) & regressor.corr < 0.9 )

prv.result$delta.j <- prv.result$deg.j.after - prv.result$deg.j.before
prv.result$delta.k <- prv.result$deg.k.after - prv.result$deg.k.before

# Each member's own status against its own split-half noise floor, kept
# separate so a pair where one row improves and the other worsens cannot be
# collapsed into a single "improved" verdict.
.prv.member.status <- function( delta, sd ) {
  ifelse( !is.finite( delta ) | !is.finite( sd ), "unchanged",
          ifelse( delta < -sd, "improved",
                  ifelse( delta > sd, "worsened", "unchanged" ) ) )
}

prv.result$status.j <- ifelse( prv.result$identifiable,
                               .prv.member.status( prv.result$delta.j, prv.result$deg.j.sd ),
                               NA_character_ )
prv.result$status.k <- ifelse( prv.result$identifiable,
                               .prv.member.status( prv.result$delta.k, prv.result$deg.k.sd ),
                               NA_character_ )

prv.result$verdict <- ifelse( !prv.result$identifiable, "unidentifiable",
                              ifelse(
                                ( prv.result$status.j == "improved" & prv.result$status.k == "worsened" ) |
                                  ( prv.result$status.j == "worsened" & prv.result$status.k == "improved" ),
                                "mixed",
                                ifelse(
                                  prv.result$status.j == "improved" | prv.result$status.k == "improved",
                                  "improved",
                                  ifelse(
                                    prv.result$status.j == "worsened" | prv.result$status.k == "worsened",
                                    "worsened", "unchanged" ) ) ) )

# Panel-level scalar, weighted by the abundance each row actually reaches so
# a large angle on a dye that never gets bright does not dominate. A true
# correction lowers this; a rotation that moves error between two rows
# leaves it roughly where it was.
prv.scored <- prv.result[ prv.result$identifiable, , drop = FALSE ]

prv.panel.score <- function( col.j, col.k ) {
  sum( prv.high[ prv.scored$j ] * prv.scored[[ col.j ]]^2 +
         prv.high[ prv.scored$k ] * prv.scored[[ col.k ]]^2, na.rm = TRUE )
}

prv.summary <- data.frame(
  n.pairs.scored   = nrow( prv.scored ),
  n.unidentifiable = sum( prv.result$verdict == "unidentifiable" ),
  n.improved       = sum( prv.result$verdict == "improved" ),
  n.worsened       = sum( prv.result$verdict == "worsened" ),
  n.mixed          = sum( prv.result$verdict == "mixed" ),
  n.unchanged      = sum( prv.result$verdict == "unchanged" ),
  panel.before     = prv.panel.score( "deg.j.before", "deg.k.before" ),
  panel.after      = prv.panel.score( "deg.j.after",  "deg.k.after" ),
  a.from           = prv.a.from,
  row.names        = NULL, stringsAsFactors = FALSE )

utils::write.csv( prv.result,
                  file.path( prv.output.dir,
                             paste0( "pair_restricted_residual_",
                                     prv.a.from, ".csv" ) ),
                  row.names = FALSE )

utils::write.csv( prv.summary,
                  file.path( prv.output.dir,
                             paste0( "pair_restricted_summary_",
                                     prv.a.from, ".csv" ) ),
                  row.names = FALSE )

print( prv.summary )

cat( "\nWorst remaining pairs after correction:\n" )
print( utils::head( prv.scored[ order( -pmax( prv.scored$deg.j.after,
                                              prv.scored$deg.k.after ) ),
                                c( "j", "k", "cosine", "separability",
                                   "deg.j.after", "deg.k.after",
                                   "verdict" ) ], 10 ) )

cat( "\nPairs the restricted design cannot resolve:\n" )
print( prv.result[ prv.result$verdict == "unidentifiable",
                   c( "j", "k", "cosine", "separability",
                      "regressor.corr", "n.j", "n.k", "status" ) ] )
