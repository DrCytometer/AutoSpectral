# test_hypernegative_trust_weight.R
#
# Standalone check of the hypernegative trust-weighting formula proposed for
# fix.my.unmix(), isolated from the rest of the pipeline so its shape can be
# inspected and min.hypernegative.frac / max.hypernegative.frac /
# min.hypernegative.events tuned before trusting the defaults on real data.
# Mirrors the block pasted into fix_my_unmix.R exactly - if this and that
# block ever disagree, this file is the one to doubt, since it has no
# compiled estimator or acceptance gates feeding it real numbers.
#
# Run top to bottom. No packages beyond base R are required.

# ---------------------------------------------------------------------------
# 1. The formula itself, isolated
# ---------------------------------------------------------------------------

# source.positive.channel, source.positive.source: numeric vectors, the
#   channel and source abundances for events already restricted to source's
#   own positive population (unmixed.comp[source.positive, channel] and
#   unmixed.comp[source.positive, source] in fix.my.unmix() itself).
# neg.threshold: numeric vector (or scalar), channel's negative boundary for
#   those same events (neg.threshold.matrix[source.positive, channel]).
# slope.use: the candidate coefficient under test.
hypernegative.trust.weight <- function(
    source.positive.channel,
    source.positive.source,
    neg.threshold,
    slope.use,
    min.hypernegative.frac = 0.01,
    max.hypernegative.frac = 0.05
) {

  hypernegative.base <- mean( source.positive.channel < neg.threshold )

  corrected <- source.positive.channel - slope.use * source.positive.source
  hypernegative.after <- mean( corrected < neg.threshold )
  hypernegative.new   <- max( hypernegative.after - hypernegative.base, 0 )

  excess.frac <- max( hypernegative.after - min.hypernegative.frac, 0 )
  frac.range  <- max( max.hypernegative.frac - min.hypernegative.frac,
                      .Machine$double.eps )
  multiplier  <- 1 - min( excess.frac / frac.range, 1 )

  list(
    hypernegative.base  = hypernegative.base,
    hypernegative.after = hypernegative.after,
    hypernegative.new   = hypernegative.new,
    multiplier          = multiplier
  )
}

# ---------------------------------------------------------------------------
# 2. Synthetic BUV805-like source-positive population -- EDIT THIS SECTION
# ---------------------------------------------------------------------------
# Replace with real vectors if you have them: pull unmixed.final[, "BUV805"],
# unmixed.final[, "PE-Cy7"] and threshold.matrix.final[, "BUV805"] from a
# fitted result, subset to the positive mask, and neg.threshold from
# neg.threshold.matrix (all four now returned by fix.my.unmix()).

set.seed( 1 )

n.events      <- 2000    # size of BUV805's own positive population
neg.threshold <- -4127   # PE-Cy7 negative boundary for these events (from the audit)

channel.noise.sd  <- 1200
source.abundance  <- rgamma( n.events, shape = 2, scale = 4000 )
channel.true      <- rnorm( n.events, mean = 0, sd = channel.noise.sd )

# ---------------------------------------------------------------------------
# 3. Case A: no real spillover, slope should be ~ 0 and stay trusted
# ---------------------------------------------------------------------------

cat( "=== Case A: negligible true coefficient ===\n" )
for ( slope in c( 0, 0.002, 0.01 ) ) {
  r <- hypernegative.trust.weight(
    channel.true, source.abundance, neg.threshold, slope )
  cat( sprintf(
    "slope %.4f: base %.4f  after %.4f  new %.4f  multiplier %.4f\n",
    slope, r$hypernegative.base, r$hypernegative.after,
    r$hypernegative.new, r$multiplier ) )
}

# ---------------------------------------------------------------------------
# 4. Case B: escalating over-correction -- multiplier should collapse by 0.05
# ---------------------------------------------------------------------------

cat( "\n=== Case B: escalating over-correction ===\n" )
for ( slope in c( 0.02, 0.03, 0.04, 0.05, 0.08 ) ) {
  r <- hypernegative.trust.weight(
    channel.true, source.abundance, neg.threshold, slope )
  cat( sprintf(
    "slope %.4f: base %.4f  after %.4f  new %.4f  multiplier %.4f\n",
    slope, r$hypernegative.base, r$hypernegative.after,
    r$hypernegative.new, r$multiplier ) )
}

# ---------------------------------------------------------------------------
# 5. Case C: pre-existing artefact from a DIFFERENT coefficient
# ---------------------------------------------------------------------------
# Simulates a channel already pushed down by some other pair's accepted
# coefficient, before this pair's own candidate is even considered. Tests
# that a small, well-supported candidate for THIS pair still gets rejected
# if the absolute fraction is already past the cap for a reason this pair
# did not cause.

cat( "\n=== Case C: pre-existing artefact from elsewhere ===\n" )
channel.already.corrected <- channel.true - 0.06 * source.abundance
for ( slope in c( 0, 0.002, 0.01 ) ) {
  r <- hypernegative.trust.weight(
    channel.already.corrected, source.abundance, neg.threshold, slope )
  cat( sprintf(
    "slope %.4f: base %.4f  after %.4f  new %.4f  multiplier %.4f\n",
    slope, r$hypernegative.base, r$hypernegative.after,
    r$hypernegative.new, r$multiplier ) )
}

# ---------------------------------------------------------------------------
# 6. Shape of the penalty across the full range
# ---------------------------------------------------------------------------

frac.grid <- seq( 0, 0.10, by = 0.002 )
multiplier.grid <- sapply( frac.grid, function( f ) {
  excess <- max( f - 0.01, 0 )
  range  <- max( 0.05 - 0.01, .Machine$double.eps )
  1 - min( excess / range, 1 )
} )

png( "hypernegative_trust_weight_shape.png", width = 900, height = 600 )
plot( frac.grid, multiplier.grid, type = "l", lwd = 2,
      xlab = "hypernegative.after (fraction of source-positive population)",
      ylab = "trust multiplier",
      main = "Hypernegative trust-weighting shape" )
abline( v = 0.01, lty = 2, col = "grey50" )
abline( v = 0.05, lty = 2, col = "red" )
text( 0.01, 0.9, "min.hypernegative.frac", pos = 4, col = "grey50" )
text( 0.05, 0.9, "max.hypernegative.frac (hard reject)", pos = 4, col = "red" )
dev.off()

cat( "\nSaved hypernegative_trust_weight_shape.png\n" )
