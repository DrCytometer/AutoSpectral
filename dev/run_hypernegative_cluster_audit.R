# run_hypernegative_cluster_audit.R
#
# Post-hoc audit of a fitted fix.my.unmix() result, using event clusters
# instead of the raw event-level fraction max.hypernegative.delta checks.
#
# Two questions, both answered from the same clustering:
#
#   1. Does an accepted spillover coefficient push any biologically coherent
#      cluster's channel centroid past its own negative boundary - not
#      just some fraction of noisy individual events?
#
#   2. For a fluorophore whose signature update got rejected with reason
#      "offset", how much of the population extract.raw.signature() fit
#      against is actually mixed with real signal from a suspected partner,
#      cluster by cluster?
#
# Requires fix.my.unmix() results that expose unmixed.final, residual.final,
# threshold.matrix.final, neg.threshold.matrix.final, dominant.final and
# spillover - all five are only present after the accompanying
# fix_my_unmix.R patch.
#
# Run top to bottom. The three printed tables are the output to report back.

# ---------------------------------------------------------------------------
# 0. Setup - EDIT THIS SECTION
# ---------------------------------------------------------------------------

# One fix.my.unmix() result per tissue, named for reporting. Always include
# more than one tissue where possible - a coefficient's effect on cluster
# negativity, and a population's contamination by a partner, should not
# depend on which tissue happened to be unmixed, and a single tissue can't
# tell the two apart from one one another (see the BUV805/PE-Cy7 session
# this script came out of).
fits <- list(
  Lung   = fixed.sig.lung,
  Spleen = fixed.sig.spleen
)

# Pairs to check in section 2. NULL auto-detects every off-diagonal
# |spillover coefficient| at or above `min.coefficient`, separately in each
# tissue's own fitted matrix. Set explicitly (e.g.
# data.frame(source = "BUV805", channel = "PE-Cy7")) to track one pair's
# behaviour across every tissue in `fits`, whatever each tissue's own fit
# for that pair turned out to be.
pairs           <- NULL
min.coefficient <- 0.01

# Fluorophore pair for section 3's population-composition check. Only
# meaningful for a fluorophore whose signature.log shows a real rejection to
# explain - set to NULL to skip section 3 entirely.
target  <- "BUV805"
partner <- "PE-Cy7"

# min.cluster.size is set well above cluster.unmixed.events()'s own default
# of 5: elsewhere in the package a cluster only needs to exist to stabilise
# a regression input, but here the centroid itself is the evidence, and a
# 5-event centroid is not a low-noise one.
som.dim          <- 30
min.cluster.size <- 20
cluster.cofactor <- 500

asp <- asp  # must already be in scope, as passed to fix.my.unmix() itself

# ---------------------------------------------------------------------------
# 1. Cluster each tissue once
# ---------------------------------------------------------------------------

clusters <- lapply( names( fits ), function( tissue ) {

  fit <- fits[[ tissue ]]

  message( sprintf( "Clustering %s (%d events)...",
                    tissue, nrow( fit$unmixed.final ) ) )

  cluster.unmixed.events(
    unmixed          = fit$unmixed.final,
    raw.data         = fit$residual.final,
    asp              = asp,
    method           = "som",
    som.dim          = som.dim,
    min.cluster.size = min.cluster.size,
    cluster.cofactor = cluster.cofactor,
    verbose          = FALSE )
} )
names( clusters ) <- names( fits )

# ---------------------------------------------------------------------------
# 2. Per-cluster negativity check: did this pair's coefficient push a
#    cluster that read non-negative before, past the boundary?
# ---------------------------------------------------------------------------

.audit.negativity.one <- function( tissue ) {

  fit <- fits[[ tissue ]]
  cl  <- clusters[[ tissue ]]

  x.clust <- cl$x.clust
  size    <- cl$cluster.size

  pair.set <- pairs
  if ( is.null( pair.set ) ) {

    sp       <- fit$spillover
    off.diag <- which( abs( sp ) >= min.coefficient & row( sp ) != col( sp ),
                       arr.ind = TRUE )

    if ( nrow( off.diag ) == 0 ) return( NULL )

    pair.set <- data.frame(
      source  = rownames( sp )[ off.diag[ , 1 ] ],
      channel = colnames( sp )[ off.diag[ , 2 ] ],
      stringsAsFactors = FALSE )
  }

  # Directly-measured negative boundary from fix.my.unmix() itself, not the
  # mirrored `-threshold.matrix.final` shortcut - the whole reason a
  # negative boundary is measured directly from the unstained population's
  # own tail, rather than derived by negating the positive one, is that AF
  # inflates the positive tail asymmetrically and makes that mirror too
  # loose. Using the mirror here would silently undo that fix for this
  # audit alone.
  neg.threshold <- fit$neg.threshold.matrix.final

  # Per-cluster mean of the per-event boundary, aligned to x.clust by row
  # name rather than assumed sort order.
  neg.threshold.clust <- t( sapply( rownames( x.clust ), function( id ) {
    idx <- which( as.character( cl$cluster.id ) == id )
    colMeans( neg.threshold[ idx, , drop = FALSE ] )
  } ) )

  rows <- lapply( seq_len( nrow( pair.set ) ), function( i ) {

    source  <- pair.set$source[ i ]
    channel <- pair.set$channel[ i ]

    if ( !( source %in% colnames( x.clust ) ) ||
         !( channel %in% colnames( x.clust ) ) ) return( NULL )

    beta <- fit$spillover[ source, channel ]

    before <- x.clust[ , channel ]
    after  <- before - beta * x.clust[ , source ]
    bound  <- neg.threshold.clust[ , channel ]

    newly.crossed <- ( after < bound ) & ( before >= bound )

    if ( !any( newly.crossed ) ) return( NULL )

    data.frame(
      tissue       = tissue,
      source       = source,
      channel      = channel,
      coefficient  = beta,
      cluster      = rownames( x.clust )[ newly.crossed ],
      cluster.size = size[ newly.crossed ],
      before       = before[ newly.crossed ],
      after        = after[ newly.crossed ],
      boundary     = bound[ newly.crossed ],
      stringsAsFactors = FALSE )
  } )

  do.call( rbind, rows )
}

negativity.audit <- do.call( rbind, lapply( names( fits ), .audit.negativity.one ) )

cat( "\n=== Section 2: clusters newly pushed past their negative boundary ===\n" )
if ( is.null( negativity.audit ) || nrow( negativity.audit ) == 0 ) {
  cat( "None found for any tested pair, in any tissue.\n" )
} else {
  print(
    negativity.audit[ order( negativity.audit$tissue, negativity.audit$source,
                             negativity.audit$channel,
                             -negativity.audit$cluster.size ), ],
    row.names = FALSE )
}

# ---------------------------------------------------------------------------
# 3. Population composition: how much of `target`'s own phase-two
#    population sits in a cluster that also reads positive for `partner`?
# ---------------------------------------------------------------------------

if ( !is.null( target ) && !is.null( partner ) ) {

  .audit.composition.one <- function( tissue ) {

    fit <- fits[[ tissue ]]
    cl  <- clusters[[ tissue ]]

    target.idx <- which( rownames( fit$spillover ) == target )

    if ( length( target.idx ) == 0 )
      stop( sprintf( "`%s` not found in %s's spillover matrix.",
                     target, tissue ), call. = FALSE )

    dom.idx <- which( fit$dominant.final == target.idx )

    if ( length( dom.idx ) == 0 ) return( NULL )

    dom.cluster <- as.character( cl$cluster.id[ dom.idx ] )
    dom.cluster <- dom.cluster[ !is.na( dom.cluster ) ]

    tab <- table( dom.cluster )

    partner.threshold <- mean(
      fit$threshold.matrix.final[ dom.idx, partner ] )

    data.frame(
      tissue            = tissue,
      cluster           = names( tab ),
      n.dominant        = as.integer( tab ),
      partner.centroid  = cl$x.clust[ names( tab ), partner ],
      partner.threshold = partner.threshold,
      stringsAsFactors  = FALSE )
  }

  composition <- do.call( rbind, lapply( names( fits ), .audit.composition.one ) )
  composition$co.positive <- composition$partner.centroid >
    composition$partner.threshold

  cat( sprintf(
    "\n=== Section 3: %s's dominant population, by cluster, vs %s ===\n",
    target, partner ) )
  print(
    composition[ order( composition$tissue, -composition$n.dominant ), ],
    row.names = FALSE )

  cat( sprintf(
    "\n%s of %s's dominant events fall in a cluster whose %s centroid reads positive:\n",
    "Fraction", target, partner ) )
  for ( tissue in names( fits ) ) {
    sub  <- composition[ composition$tissue == tissue, ]
    frac <- sum( sub$n.dominant[ sub$co.positive ] ) / sum( sub$n.dominant )
    cat( sprintf( "  %s: %.1f%% (%d of %d dominant events)\n",
                 tissue, 100 * frac,
                 sum( sub$n.dominant[ sub$co.positive ] ),
                 sum( sub$n.dominant ) ) )
  }
}
