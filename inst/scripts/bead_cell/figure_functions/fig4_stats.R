# fig4_stats.R
#
# Mixed models for Figure 4. They extend the Figure 2 models
# (bead_angle_stats.R): the response is on the log2 scale, the cytometer x
# bead-type combination (`Group`) is a fixed nuisance effect, and fluorophores
# are random intercepts:
#   (1 | Fluorophore)            fluorophore effect shared by all cytometers
#   (1 | Fluorophore:Cytometer)  fluorophore effect specific to one cytometer
# A random term is only included if its levels are repeated in the data (e.g.
# cell variability has one value per fluorophore per cytometer, so
# Fluorophore:Cytometer is then the residual and is left out). Without any
# repeated level, an ordinary linear model is fitted.
#
#   fit.class.models()  categorical fluorophore property (dye class): each
#                       class vs a reference class, per cytometer and pooled;
#                       likelihood-ratio tests for the class effect and the
#                       class x cytometer interaction.
#   fit.slope.models()  continuous fluorophore properties (variability,
#                       messiness, excitation-emission separation): slope of
#                       the log2 response on each predictor, adjusted for the
#                       others, per cytometer and pooled.

# ---- helpers ----------------------------------------------------------------

.fig4.random.terms <- function( d ) {
  rep.f  <- any( duplicated( d$Fluorophore ) )
  rep.fc <- "FluorCyto" %in% names( d ) && any( duplicated( d$FluorCyto ) )
  c( if ( rep.f ) "(1 | Fluorophore)", if ( rep.fc ) "(1 | FluorCyto)" )
}

# fixed-effect nuisance term: Group if it has > 1 level, otherwise nothing
.fig4.nuisance <- function( d, nuisance = "Group" ) {
  if ( is.null( nuisance ) || !nuisance %in% names( d ) ) return( character( 0 ) )
  if ( length( unique( d[[ nuisance ]] ) ) > 1 ) nuisance else character( 0 )
}

.fig4.formula <- function( response, fixed, random ) {
  rhs <- c( fixed, random )
  if ( length( rhs ) == 0 ) rhs <- "1"
  stats::as.formula( paste( response, "~", paste( rhs, collapse = " + " ) ) )
}

# lmer (lmerTest) if the formula has random terms, lm otherwise
.fig4.fit <- function( formula, data, REML = TRUE ) {
  if ( length( lme4::findbars( formula ) ) > 0 ) {
    if ( REML ) lmerTest::lmer( formula, data = data, REML = TRUE )
    else lme4::lmer( formula, data = data, REML = FALSE )
  } else {
    stats::lm( formula, data = data )
  }
}

# likelihood-ratio test of two nested fits (ML for lmer; F test for lm)
.fig4.lrt <- function( m0, m1, label ) {
  if ( inherits( m1, "merMod" ) ) {
    a <- as.data.frame( stats::anova( m0, m1 ) )
    df.col <- if ( "Chi Df" %in% names( a ) ) "Chi Df" else "Df"
    data.frame( Test = label, Statistic = "Chisq", Value = a$Chisq[ 2 ],
                Df = a[[ df.col ]][ 2 ], Df.resid = NA_real_, P.value = a[[ "Pr(>Chisq)" ]][ 2 ],
                stringsAsFactors = FALSE )
  } else {
    a <- as.data.frame( stats::anova( m0, m1 ) )
    data.frame( Test = label, Statistic = "F", Value = a$F[ 2 ], Df = a$Df[ 2 ],
                Df.resid = a$Res.Df[ 2 ], P.value = a[[ "Pr(>F)" ]][ 2 ],
                stringsAsFactors = FALSE )
  }
}

# ---- categorical: dye class --------------------------------------------------

#' Class effect on a log2 response
#'
#' @param d data.frame with the response column, Class, Cytometer,
#'   Fluorophore, FluorCyto and the nuisance column (`Group`).
#' @param response Response column (log2 scale).
#' @param reference.class Comparator class.
#' @param min.class.fluors Classes with fewer distinct fluorophores (pooled,
#'   or within a cytometer for the per-cytometer contrasts) are left out of
#'   the statistics.
#' @param class.order Factor order of the classes.
#' @return list(per.cyto, pooled, lrt, class.n, excluded.classes, models)
fit.class.models <- function(
    d, response = "log2A", reference.class = "Small", min.class.fluors = 3,
    class.order = NULL, nuisance = "Group", conf.level = 0.95
) {
  for ( pkg in c( "lme4", "lmerTest", "emmeans" ) )
    if ( !requireNamespace( pkg, quietly = TRUE ) )
      stop( "Package '", pkg, "' is required", call. = FALSE )

  d <- as.data.frame( d )
  d <- d[ is.finite( d[[ response ]] ) & !is.na( d$Class ), , drop = FALSE ]

  n.fl <- tapply( d$Fluorophore, d$Class, function( x ) length( unique( x ) ) )
  class.n <- data.frame( Class = names( n.fl ), N.fluorophores = as.integer( n.fl ),
                         stringsAsFactors = FALSE )
  excluded <- names( n.fl )[ n.fl < min.class.fluors ]
  if ( reference.class %in% excluded || !reference.class %in% names( n.fl ) )
    stop( "Reference class '", reference.class, "' has < ", min.class.fluors,
          " fluorophores or is absent", call. = FALSE )

  if ( is.null( class.order ) ) class.order <- sort( unique( d$Class ) )
  keep.classes <- c( reference.class, setdiff( intersect( class.order, names( n.fl ) ), c( reference.class, excluded ) ) )
  dm <- d[ d$Class %in% keep.classes, , drop = FALSE ]
  dm$Class <- factor( dm$Class, levels = keep.classes )
  dm <- droplevels( dm )

  contrasts.vs.ref <- function( fit ) {
    emm <- if ( inherits( fit, "merMod" ) ) emmeans::emmeans( fit, ~ Class, lmer.df = "satterthwaite" ) else
      emmeans::emmeans( fit, ~ Class )
    lev <- as.character( as.data.frame( summary( emm ) )$Class )
    oth <- setdiff( lev, reference.class )
    if ( length( oth ) == 0 ) return( NULL )
    con <- lapply( oth, function( k ) as.numeric( lev == k ) - as.numeric( lev == reference.class ) )
    names( con ) <- paste0( oth, " / ", reference.class )
    s <- .tidy.emm( summary( emmeans::contrast( emm, method = con, adjust = "none" ),
                             infer = c( TRUE, TRUE ), level = conf.level ), "estimate" )
    data.frame( Class = oth, Reference = reference.class, s, stringsAsFactors = FALSE )
  }

  # ---- per cytometer ----
  per.cyto <- do.call( rbind, lapply( levels( factor( dm$Cytometer ) ), function( cy ) {
    dc <- droplevels( dm[ dm$Cytometer == cy, , drop = FALSE ] )
    nf <- tapply( dc$Fluorophore, dc$Class, function( x ) length( unique( x ) ) )
    ok <- names( nf )[ nf >= min.class.fluors ]
    if ( !reference.class %in% ok || length( ok ) < 2 ) return( NULL )
    dc <- droplevels( dc[ dc$Class %in% ok, , drop = FALSE ] )
    dc$FluorCyto <- NULL
    f <- .fig4.formula( response, c( "Class", .fig4.nuisance( dc, nuisance ) ), .fig4.random.terms( dc ) )
    fit <- .fig4.fit( f, dc )
    out <- contrasts.vs.ref( fit )
    if ( is.null( out ) ) return( NULL )
    out$P.holm <- stats::p.adjust( out$P.value, "holm" )
    data.frame( Cytometer = cy, out, stringsAsFactors = FALSE )
  } ) )

  # ---- pooled ----
  rnd <- .fig4.random.terms( dm )
  nui <- .fig4.nuisance( dm, nuisance )
  f.add <- .fig4.formula( response, c( "Class", nui ), rnd )
  fit.add <- .fig4.fit( f.add, dm )
  pooled <- contrasts.vs.ref( fit.add )
  pooled$P.holm <- stats::p.adjust( pooled$P.value, "holm" )
  pooled <- data.frame( Cytometer = "Pooled", pooled, stringsAsFactors = FALSE )

  # consistency of each pooled contrast across cytometers
  pooled$N.cytometers <- pooled$N.same.direction <- pooled$N.significant.same.direction <- NA_integer_
  if ( !is.null( per.cyto ) ) for ( i in seq_len( nrow( pooled ) ) ) {
    pc <- per.cyto[ per.cyto$Class == pooled$Class[ i ], , drop = FALSE ]
    s <- sign( pooled$Estimate[ i ] )
    pooled$N.cytometers[ i ] <- nrow( pc )
    pooled$N.same.direction[ i ] <- sum( sign( pc$Estimate ) == s )
    pooled$N.significant.same.direction[ i ] <- sum( sign( pc$Estimate ) == s & pc$P.holm < 0.05 )
  }

  # ---- likelihood-ratio tests (ML) ----
  f.null <- .fig4.formula( response, nui, rnd )
  f.int  <- .fig4.formula( response, c( "Class * Cytometer", nui ), rnd )
  m.null <- .fig4.fit( f.null, dm, REML = FALSE )
  m.add  <- .fig4.fit( f.add,  dm, REML = FALSE )
  lrt <- .fig4.lrt( m.null, m.add, "Class effect" )
  if ( length( unique( dm$Cytometer ) ) > 1 ) {
    m.int <- suppressMessages( .fig4.fit( f.int, dm, REML = FALSE ) )
    lrt <- rbind( lrt, .fig4.lrt( m.add, m.int, "Class x cytometer interaction" ) )
  }

  list( per.cyto = per.cyto, pooled = pooled, lrt = lrt, class.n = class.n,
        excluded.classes = excluded, classes.tested = keep.classes,
        n.obs = nrow( dm ), formula = deparse( f.add ),
        models = list( additive = fit.add ) )
}

# ---- continuous predictors ----------------------------------------------------

#' Slopes of a log2 response on fluorophore properties
#'
#' @param d data.frame with the response, the predictor columns, Cytometer,
#'   Fluorophore, FluorCyto and the nuisance column.
#' @param predictors Predictor columns (already on the analysis scale, e.g.
#'   log2(variability) or separation / 100 nm). All are fitted together, so
#'   each slope is adjusted for the others.
#' @return list(per.cyto, pooled, n.obs, formula, models); slopes in log2 units
#'   of the response per unit of the predictor (2^slope = fold change).
fit.slope.models <- function(
    d, response = "log2A", predictors, nuisance = "Group", conf.level = 0.95
) {
  d <- as.data.frame( d )
  ok <- is.finite( d[[ response ]] ) & Reduce( `&`, lapply( predictors, function( p ) is.finite( d[[ p ]] ) ) )
  d <- d[ ok, , drop = FALSE ]

  slope.table <- function( fit, dd ) {
    cf <- if ( inherits( fit, "merMod" ) ) summary( fit )$coefficients else summary( fit )$coefficients
    cf <- cf[ predictors, , drop = FALSE ]
    df <- if ( "df" %in% colnames( cf ) ) cf[ , "df" ] else rep( stats::df.residual( fit ), nrow( cf ) )
    q <- stats::qt( 1 - ( 1 - conf.level ) / 2, df )
    data.frame(
      Predictor = predictors, Estimate = cf[ , "Estimate" ], SE = cf[ , "Std. Error" ], df = df,
      Lower = cf[ , "Estimate" ] - q * cf[ , "Std. Error" ],
      Upper = cf[ , "Estimate" ] + q * cf[ , "Std. Error" ],
      Statistic = cf[ , grep( "t value", colnames( cf ) ) ],
      P.value = cf[ , grep( "^Pr\\(", colnames( cf ) ) ],
      N.obs = nrow( dd ), N.fluorophores = length( unique( dd$Fluorophore ) ),
      stringsAsFactors = FALSE, row.names = NULL
    )
  }

  per.cyto <- do.call( rbind, lapply( levels( factor( d$Cytometer ) ), function( cy ) {
    dc <- droplevels( d[ d$Cytometer == cy, , drop = FALSE ] )
    dc$FluorCyto <- NULL
    f <- .fig4.formula( response, c( predictors, .fig4.nuisance( dc, nuisance ) ), .fig4.random.terms( dc ) )
    data.frame( Cytometer = cy, slope.table( .fig4.fit( f, dc ), dc ), stringsAsFactors = FALSE )
  } ) )
  per.cyto$P.holm <- stats::ave( per.cyto$P.value, per.cyto$Cytometer,
                                 FUN = function( p ) stats::p.adjust( p, "holm" ) )

  f <- .fig4.formula( response, c( predictors, .fig4.nuisance( d, nuisance ) ), .fig4.random.terms( d ) )
  fit <- .fig4.fit( f, d )
  pooled <- data.frame( Cytometer = "Pooled", slope.table( fit, d ), stringsAsFactors = FALSE )
  pooled$P.holm <- stats::p.adjust( pooled$P.value, "holm" )

  # consistency across cytometers
  pooled$N.same.direction <- vapply( seq_len( nrow( pooled ) ), function( i )
    sum( sign( per.cyto$Estimate[ per.cyto$Predictor == pooled$Predictor[ i ] ] ) == sign( pooled$Estimate[ i ] ) ),
    integer( 1 ) )

  list( per.cyto = per.cyto, pooled = pooled, n.obs = nrow( d ), formula = deparse( f ),
        models = list( pooled = fit ) )
}

#' Spearman correlations of a response with each predictor, per subgroup
#' (descriptive; no adjustment)
spearman.table <- function( d, response, predictors, by = c( "Cytometer", "Bead" ) ) {
  by <- intersect( by, names( d ) )
  key <- interaction( d[ , by, drop = FALSE ], drop = TRUE, sep = " | " )
  do.call( rbind, lapply( split( seq_len( nrow( d ) ), key ), function( i ) {
    do.call( rbind, lapply( predictors, function( p ) {
      ok <- is.finite( d[[ response ]][ i ] ) & is.finite( d[[ p ]][ i ] )
      if ( sum( ok ) < 4 ) return( NULL )
      ct <- suppressWarnings( stats::cor.test( d[[ response ]][ i ][ ok ], d[[ p ]][ i ][ ok ],
                                               method = "spearman", exact = FALSE ) )
      data.frame( d[ i[ 1 ], by, drop = FALSE ], Predictor = p, N = sum( ok ),
                  Spearman.rho = unname( ct$estimate ), P.value = ct$p.value,
                  stringsAsFactors = FALSE, row.names = NULL )
    } ) )
  } ) )
}
