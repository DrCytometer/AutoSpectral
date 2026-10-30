# bead_angle_stats.R
#
# Statistics for Figure 2: does the bead type affect the bead-vs-cell spectral
# angle, and is the effect the same on every cytometer?
#
# Response: log2( spectral angle, degrees ). On the log scale a model contrast
# is a log2 ratio of geometric means, so 2^estimate is the fold difference in
# angle between two bead types.
#
# Pairing: within one cytometer every bead type is compared against the SAME
# cell spectrum of each fluorophore, and the same fluorophores (reagents) are
# used on every cytometer. Both are handled with random intercepts:
#   (1 | Fluorophore)            fluorophore effect shared by all cytometers
#   (1 | Fluorophore:Cytometer)  fluorophore effect specific to one cytometer
#                                (e.g. an unusual cell spectrum on one instrument)
#
# Models (lme4 / lmerTest; emmeans for contrasts, Satterthwaite df):
#   cell-means  log2A ~ 0 + Group + RE     Group = Cytometer x Bead (all beads)
#               -> per-cytometer contrasts (each bead vs the reference bead,
#                  and all pairs), Holm-adjusted within each cytometer
#   additive    log2A ~ Bead + Cytometer + RE   (beads on >= min.cytometers)
#               -> pooled bead effect (equal weight per cytometer)
#   LRT (ML)    Cytometer  vs  Bead + Cytometer  vs  Group (= Bead x Cytometer)
#               -> is there a bead effect; does it differ between cytometers
#
# Rank-based checks (no distributional assumptions):
#   per cytometer: Friedman test (fluorophores = blocks, complete cases),
#                  paired Wilcoxon signed-rank tests for every bead pair (Holm)
#   across cytometers: Kendall's W of the bead mean ranks, permutation test.

.random.part <- "(1 | Fluorophore) + (1 | FluorCyto)"

#' Prepare the long angle table for modelling
#' @param d Output of load.bead.angles().
#' @param bead.order,cytometer.order Factor level orders.
prepare.angle.model.data <- function( d, bead.order, cytometer.order ) {
  d <- d[ d$Bead %in% bead.order & d$Cytometer %in% cytometer.order, , drop = FALSE ]
  d$log2A     <- log2( d$Angle )
  d$Cytometer <- factor( d$Cytometer, levels = cytometer.order )
  d$Bead      <- factor( d$Bead, levels = bead.order )
  d$FluorCyto <- factor( paste( d$Fluorophore, d$Cytometer, sep = ":" ) )
  grp <- unique( d[ order( d$Cytometer, d$Bead ), c( "Cytometer", "Bead" ) ] )
  d$Group <- factor( paste( d$Cytometer, d$Bead, sep = " | " ),
                     levels = paste( grp$Cytometer, grp$Bead, sep = " | " ) )
  droplevels( d )
}

#' Tidy an emmeans summary (handles finite and asymptotic df)
.tidy.emm <- function( s, est.col ) {
  s  <- as.data.frame( s )
  lo <- grep( "lower\\.CL|asymp\\.LCL", names( s ), value = TRUE )[ 1 ]
  hi <- grep( "upper\\.CL|asymp\\.UCL", names( s ), value = TRUE )[ 1 ]
  st <- grep( "\\.ratio$", names( s ), value = TRUE )
  out <- data.frame(
    Estimate = s[[ est.col ]],
    SE       = s$SE,
    df       = if ( "df" %in% names( s ) ) s$df else Inf,
    Lower    = s[[ lo ]],
    Upper    = s[[ hi ]],
    stringsAsFactors = FALSE
  )
  if ( length( st ) == 1 ) {
    out$Statistic <- s[[ st ]]
    out$P.value   <- s$p.value
  }
  out
}

#' Build emmeans contrast vectors on the Group levels
#' @param pairs data.frame with columns Cytometer, Bead1, Bead2 (Bead1 / Bead2).
.group.contrasts <- function( pairs, group.levels ) {
  con <- lapply( seq_len( nrow( pairs ) ), function( i ) {
    v <- numeric( length( group.levels ) )
    v[ group.levels == paste( pairs$Cytometer[ i ], pairs$Bead1[ i ], sep = " | " ) ] <-  1
    v[ group.levels == paste( pairs$Cytometer[ i ], pairs$Bead2[ i ], sep = " | " ) ] <- -1
    v
  } )
  names( con ) <- paste0( pairs$Cytometer, ": ", pairs$Bead1, " / ", pairs$Bead2 )
  con
}

#' Fit the mixed models and return contrasts and tests
#'
#' @param d Output of prepare.angle.model.data().
#' @param reference.bead Bead type used as the comparator in the forest plot.
#' @param min.cytometers Bead types measured on fewer cytometers are kept in the
#'   per-cytometer analysis but left out of the pooled model and the LRT.
#' @param conf.level Confidence level of the intervals (unadjusted).
#' @return list of data.frames and the fitted models.
fit.bead.angle.models <- function(
    d, reference.bead = "BD", min.cytometers = 2, conf.level = 0.95
) {

  for ( pkg in c( "lme4", "lmerTest", "emmeans" ) )
    if ( !requireNamespace( pkg, quietly = TRUE ) )
      stop( "Package '", pkg, "' is required: install.packages(\"", pkg, "\")", call. = FALSE )

  if ( !reference.bead %in% levels( d$Bead ) )
    stop( "reference.bead '", reference.bead, "' not found in the data", call. = FALSE )

  f.cell <- stats::as.formula( paste( "log2A ~ 0 + Group +", .random.part ) )

  # ---- per-cytometer model (all bead types) ---------------------------------
  fit.cell <- lmerTest::lmer( f.cell, data = d, REML = TRUE )
  emm.g    <- emmeans::emmeans( fit.cell, ~ Group, lmer.df = "satterthwaite" )

  group.means <- .tidy.emm( summary( emm.g, infer = c( TRUE, FALSE ), level = conf.level ), "emmean" )
  grp.split   <- do.call( rbind, strsplit( levels( d$Group ), " | ", fixed = TRUE ) )
  group.means <- data.frame(
    Cytometer = grp.split[ , 1 ], Bead = grp.split[ , 2 ],
    Geometric.mean.deg = 2^group.means$Estimate,
    Lower.deg = 2^group.means$Lower, Upper.deg = 2^group.means$Upper,
    log2.mean = group.means$Estimate, SE = group.means$SE, df = group.means$df,
    stringsAsFactors = FALSE
  )

  present <- split( as.character( d$Bead ), d$Cytometer )
  present <- lapply( present, function( b ) levels( d$Bead )[ levels( d$Bead ) %in% b ] )

  run.contrasts <- function( pairs ) {
    if ( nrow( pairs ) == 0 ) return( NULL )
    con <- .group.contrasts( pairs, levels( d$Group ) )
    ctr <- emmeans::contrast( emm.g, method = con, adjust = "none" )
    out <- cbind( pairs, .tidy.emm( summary( ctr, infer = c( TRUE, TRUE ), level = conf.level ),
                                    "estimate" ) )
    out$P.holm <- stats::ave( out$P.value, out$Cytometer, FUN = function( p ) stats::p.adjust( p, "holm" ) )
    out
  }

  # each bead vs the reference bead, per cytometer
  vs.ref.pairs <- do.call( rbind, lapply( names( present ), function( cy ) {
    b <- setdiff( present[[ cy ]], reference.bead )
    if ( !reference.bead %in% present[[ cy ]] || length( b ) == 0 ) return( NULL )
    data.frame( Cytometer = cy, Bead1 = b, Bead2 = reference.bead, stringsAsFactors = FALSE )
  } ) )
  per.cyto.vs.ref <- run.contrasts( vs.ref.pairs )

  # every pair of beads, per cytometer
  all.pairs <- do.call( rbind, lapply( names( present ), function( cy ) {
    b <- present[[ cy ]]
    if ( length( b ) < 2 ) return( NULL )
    cb <- utils::combn( b, 2 )
    data.frame( Cytometer = cy, Bead1 = cb[ 2, ], Bead2 = cb[ 1, ], stringsAsFactors = FALSE )
  } ) )
  per.cyto.pairwise <- run.contrasts( all.pairs )

  # ---- pooled analysis (beads on >= min.cytometers) ------------------------
  n.cyto <- tapply( as.character( d$Cytometer ), d$Bead, function( x ) length( unique( x ) ) )
  pooled.beads <- names( n.cyto )[ n.cyto >= min.cytometers ]
  if ( !reference.bead %in% pooled.beads )
    stop( "reference.bead '", reference.bead, "' is on fewer than ", min.cytometers,
          " cytometers", call. = FALSE )
  dp <- droplevels( d[ d$Bead %in% pooled.beads, , drop = FALSE ] )

  f.null <- stats::as.formula( paste( "log2A ~ Cytometer +", .random.part ) )
  f.add  <- stats::as.formula( paste( "log2A ~ Bead + Cytometer +", .random.part ) )

  m.null <- lme4::lmer( f.null, data = dp, REML = FALSE )
  m.add  <- lme4::lmer( f.add,  data = dp, REML = FALSE )
  m.int  <- lme4::lmer( f.cell, data = dp, REML = FALSE )
  # anova.merMod orders the models by number of parameters (null < add < int);
  # column names differ slightly between lme4 versions
  lrt <- as.data.frame( stats::anova( m.null, m.add, m.int ) )
  npar.col <- if ( "npar" %in% names( lrt ) ) "npar" else "Df"
  df.col   <- if ( "Chi Df" %in% names( lrt ) ) "Chi Df" else "Df"
  lrt <- data.frame(
    Model = c( "Cytometer", "Bead + Cytometer", "Bead x Cytometer (cell means)" ),
    Test  = c( NA, "Bead effect", "Bead x cytometer interaction" ),
    n.par = lrt[[ npar.col ]], AIC = lrt$AIC, BIC = lrt$BIC, logLik = lrt$logLik,
    Chisq = lrt$Chisq, Df = lrt[[ df.col ]], P.value = lrt[[ "Pr(>Chisq)" ]],
    stringsAsFactors = FALSE
  )

  fit.add <- lmerTest::lmer( f.add, data = dp, REML = TRUE )
  emm.b   <- emmeans::emmeans( fit.add, ~ Bead, lmer.df = "satterthwaite" )

  pooled.means <- .tidy.emm( summary( emm.b, infer = c( TRUE, FALSE ), level = conf.level ), "emmean" )
  pooled.means <- data.frame(
    Bead = levels( dp$Bead ),
    Geometric.mean.deg = 2^pooled.means$Estimate,
    Lower.deg = 2^pooled.means$Lower, Upper.deg = 2^pooled.means$Upper,
    stringsAsFactors = FALSE
  )

  b.lev <- levels( dp$Bead )
  pooled.con <- lapply( setdiff( b.lev, reference.bead ), function( b )
    as.numeric( b.lev == b ) - as.numeric( b.lev == reference.bead ) )
  names( pooled.con ) <- paste0( setdiff( b.lev, reference.bead ), " / ", reference.bead )
  pooled.vs.ref <- .tidy.emm(
    summary( emmeans::contrast( emm.b, method = pooled.con, adjust = "none" ),
             infer = c( TRUE, TRUE ), level = conf.level ), "estimate" )
  pooled.vs.ref <- data.frame( Cytometer = "Pooled", Bead1 = setdiff( b.lev, reference.bead ),
                               Bead2 = reference.bead, pooled.vs.ref, stringsAsFactors = FALSE )
  pooled.vs.ref$P.holm <- stats::p.adjust( pooled.vs.ref$P.value, "holm" )

  pooled.pairwise <- as.data.frame( summary( graphics::pairs( emm.b, adjust = "holm", reverse = TRUE ),
                                             infer = c( TRUE, TRUE ), level = conf.level ) )

  # consistency of each pooled contrast across cytometers
  pooled.vs.ref$N.cytometers <- NA_integer_
  pooled.vs.ref$N.same.direction <- NA_integer_
  pooled.vs.ref$N.significant.same.direction <- NA_integer_
  for ( i in seq_len( nrow( pooled.vs.ref ) ) ) {
    pc <- per.cyto.vs.ref[ per.cyto.vs.ref$Bead1 == pooled.vs.ref$Bead1[ i ], , drop = FALSE ]
    s  <- sign( pooled.vs.ref$Estimate[ i ] )
    pooled.vs.ref$N.cytometers[ i ] <- nrow( pc )
    pooled.vs.ref$N.same.direction[ i ] <- sum( sign( pc$Estimate ) == s )
    pooled.vs.ref$N.significant.same.direction[ i ] <- sum( sign( pc$Estimate ) == s & pc$P.holm < 0.05 )
  }

  # variance components of the additive model
  vc <- as.data.frame( lme4::VarCorr( fit.add ) )
  vc <- data.frame( Component = vc$grp, Variance = vc$vcov, SD = vc$sdcor,
                    stringsAsFactors = FALSE )
  vc$Component[ vc$Component == "FluorCyto" ] <- "Fluorophore:Cytometer"

  list(
    group.means       = group.means,
    per.cyto.vs.ref   = per.cyto.vs.ref,
    per.cyto.pairwise = per.cyto.pairwise,
    pooled.beads      = pooled.beads,
    pooled.means      = pooled.means,
    pooled.vs.ref     = pooled.vs.ref,
    pooled.pairwise   = pooled.pairwise,
    lrt               = lrt,
    variance          = vc,
    n.obs             = c( all = nrow( d ), pooled = nrow( dp ) ),
    models            = list( cell = fit.cell, additive = fit.add,
                              ml.null = m.null, ml.add = m.add, ml.int = m.int )
  )
}

#' Wide fluorophore x bead matrix of angles for one cytometer (complete cases)
angle.wide <- function( d, cytometer, beads ) {
  x <- d[ d$Cytometer == cytometer & d$Bead %in% beads, , drop = FALSE ]
  w <- stats::reshape( x[ , c( "Fluorophore", "Bead", "Angle" ) ],
                       idvar = "Fluorophore", timevar = "Bead", direction = "wide" )
  rownames( w ) <- w$Fluorophore
  w <- as.matrix( w[ , -1, drop = FALSE ] )
  colnames( w ) <- sub( "^Angle\\.", "", colnames( w ) )
  w <- w[ , intersect( beads, colnames( w ) ), drop = FALSE ]
  w[ stats::complete.cases( w ), , drop = FALSE ]
}

#' Friedman and paired Wilcoxon tests per cytometer
#'
#' @param d Output of prepare.angle.model.data().
#' @return list(friedman, wilcoxon, mean.ranks)
bead.angle.nonparametric <- function( d ) {

  cytos <- levels( d$Cytometer )
  beads <- levels( d$Bead )

  res <- lapply( cytos, function( cy ) {
    w <- angle.wide( as.data.frame( lapply( d, function( x ) if ( is.factor( x ) ) as.character( x ) else x ) ),
                     cy, beads )
    if ( ncol( w ) < 2 || nrow( w ) < 3 ) return( NULL )

    fr <- stats::friedman.test( w )
    rk <- t( apply( w, 1, rank ) )
    mean.rank <- colMeans( rk )

    cb <- utils::combn( colnames( w ), 2 )
    wx <- do.call( rbind, lapply( seq_len( ncol( cb ) ), function( j ) {
      a <- cb[ 2, j ]; b <- cb[ 1, j ]
      wt <- suppressWarnings( stats::wilcox.test( w[ , a ], w[ , b ], paired = TRUE ) )
      lr <- log2( w[ , a ] / w[ , b ] )
      data.frame( Cytometer = cy, Bead1 = a, Bead2 = b, N.fluorophores = nrow( w ),
                  Median.ratio = 2^stats::median( lr ),
                  N.Bead1.larger = sum( lr > 0 ),
                  V = unname( wt$statistic ), P.value = wt$p.value,
                  stringsAsFactors = FALSE )
    } ) )
    wx$P.holm <- stats::p.adjust( wx$P.value, "holm" )

    list(
      friedman = data.frame(
        Cytometer = cy, N.fluorophores = nrow( w ), N.beads = ncol( w ),
        Beads = paste( colnames( w ), collapse = ";" ),
        Chisq = unname( fr$statistic ), Df = unname( fr$parameter ), P.value = fr$p.value,
        Kendall.W.within = unname( fr$statistic ) / ( nrow( w ) * ( ncol( w ) - 1 ) ),
        stringsAsFactors = FALSE ),
      wilcoxon = wx,
      mean.ranks = data.frame( Cytometer = cy, Bead = names( mean.rank ),
                               Mean.rank = unname( mean.rank ), N.fluorophores = nrow( w ),
                               stringsAsFactors = FALSE )
    )
  } )
  res <- res[ !vapply( res, is.null, logical( 1 ) ) ]

  list(
    friedman   = do.call( rbind, lapply( res, `[[`, "friedman" ) ),
    wilcoxon   = do.call( rbind, lapply( res, `[[`, "wilcoxon" ) ),
    mean.ranks = do.call( rbind, lapply( res, `[[`, "mean.ranks" ) )
  )
}

#' Agreement of the bead ranking between cytometers (Kendall's W)
#'
#' For each cytometer that has every bead in `beads`, the beads are ranked
#' within each fluorophore (complete cases) and the mean rank per bead is
#' ranked again (1 = smallest angle). W is Kendall's coefficient of
#' concordance of these rankings across cytometers. P from a permutation test
#' that shuffles the ranking of each cytometer independently.
#'
#' @param d Output of prepare.angle.model.data().
#' @param beads Bead types to rank (must be present on >= 2 cytometers).
bead.rank.concordance <- function( d, beads, n.perm = 10000, seed = 42 ) {

  dd <- as.data.frame( lapply( d, function( x ) if ( is.factor( x ) ) as.character( x ) else x ) )
  cytos <- levels( d$Cytometer )
  has.all <- vapply( cytos, function( cy ) all( beads %in% dd$Bead[ dd$Cytometer == cy ] ), logical( 1 ) )
  cytos <- cytos[ has.all ]
  if ( length( cytos ) < 2 ) {
    warning( "Fewer than 2 cytometers have all of: ", paste( beads, collapse = ", " ), call. = FALSE )
    return( NULL )
  }

  mean.ranks <- t( vapply( cytos, function( cy ) {
    w <- angle.wide( dd, cy, beads )
    colMeans( t( apply( w, 1, rank ) ) )[ beads ]
  }, numeric( length( beads ) ) ) )
  R <- t( apply( mean.ranks, 1, rank ) )

  kendall.w <- function( R ) {
    m <- nrow( R ); n <- ncol( R )
    S <- sum( ( colSums( R ) - m * ( n + 1 ) / 2 )^2 )
    12 * S / ( m^2 * ( n^3 - n ) )
  }

  W <- kendall.w( R )
  set.seed( seed )
  W.perm <- replicate( n.perm, kendall.w( t( apply( R, 1, sample ) ) ) )
  p <- ( 1 + sum( W.perm >= W - 1e-12 ) ) / ( 1 + n.perm )

  list(
    W = W, p.value = p, n.perm = n.perm, cytometers = cytos, beads = beads,
    mean.ranks = mean.ranks, ranks = R
  )
}
