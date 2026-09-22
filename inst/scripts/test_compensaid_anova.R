# test_compensaid_anova.R
#
# Blocked ANOVA with Dunnett's many-to-one comparisons for
# compare.compensaid.folders() output. Replicate identity (the same
# physical specimen unmixed by every operator/method) is the block, so
# specimen-to-specimen differences in flag counts are removed before
# methods are compared.
#
# Requires the Suggested package 'multcomp'. Depends on .resolve.df()
# (defined in test_unmix_comparison.R -- despite the filename, this is a
# production dependency, not test-only code).


#' @title Blocked ANOVA and Dunnett Comparisons of CompensAID Flags Against a
#' Reference
#'
#' @description
#' Fits `metric ~ replicate + folder` to the per-replicate CompensAID
#' summary, with replicate as a fixed blocking factor, then reports (1) the
#' omnibus F-test for any difference among folders and (2) Dunnett's
#' many-to-one comparisons of every other folder against `reference.folder`.
#' Dunnett's adjustment controls the family-wise error rate across the
#' folders compared in one call, exploiting the fact that every comparison
#' shares the same reference.
#'
#' Multiplicity across metrics is not handled here: each call tests one
#' `metric`. Choose one primary metric in advance (typically `"n.flagged"` or
#' `"frac.flagged"`) and treat any others as descriptive. The severity-band
#' columns sum to `n.flagged`, so they are not independent endpoints.
#'
#' With exactly two folders the omnibus test, the Dunnett test and a paired
#' t-test on the same metric give identical p-values.
#'
#' The model assumes the folder effect is the same in every replicate
#' (additivity) and that residual variance is similar across replicates.
#' Counts with a large range across replicates usually satisfy this better
#' after `transform = "sqrt"` or `"log1p"`. Residual diagnostics can be
#' drawn from the returned `model`.
#'
#' @param replicate.summary Either the `replicate.summary` data frame
#' returned by `compare.compensaid.folders()` (or read back in from its
#' `replicate.summary.csv`), or a character path to that CSV.
#' @param metric Character, length one. Column of `replicate.summary` to
#' test. Default `"n.flagged"`.
#' @param reference.folder Character. Folder every other folder is compared
#' against. Default `"AutoSpectral"`.
#' @param transform Character, one of `"none"` (default), `"sqrt"` or
#' `"log1p"`. Applied to `metric` before fitting. Estimates and confidence
#' limits are on the transformed scale.
#' @param conf.level Numeric, default `0.95`. Level of the simultaneous
#' confidence intervals for the folder - reference differences.
#' @param omnibus.csv Character. Path for the omnibus F-test table. Default
#' `"compensaid_comparison_anova_omnibus.csv"`.
#' @param contrast.csv Character. Path for the Dunnett contrast table.
#' Default `"compensaid_comparison_anova_contrasts.csv"`.
#' @param verbose Logical, default `TRUE`.
#'
#' @return Invisibly, a named list with `omnibus` (one-row data frame),
#' `contrasts` (one row per non-reference folder: `estimate` = folder minus
#' reference, `std.error`, `statistic`, Dunnett-adjusted `p.adjusted`, and
#' simultaneous `conf.low`/`conf.high`), and `model` (the `lm` fit).
#'
#' @importFrom stats lm anova relevel confint
#'
#' @seealso [compare.compensaid.folders()], [test.compensaid.comparison()]
#'
#' @export

test.compensaid.comparison.anova <- function(
    replicate.summary,
    metric            = "n.flagged",
    reference.folder  = "AutoSpectral",
    transform         = c( "none", "sqrt", "log1p" ),
    conf.level        = 0.95,
    omnibus.csv       = "compensaid_comparison_anova_omnibus.csv",
    contrast.csv      = "compensaid_comparison_anova_contrasts.csv",
    verbose           = TRUE
) {
  
  transform <- match.arg( transform )
  
  if ( length( metric ) != 1 ) {
    stop( "`metric` must be a single column name.", call. = FALSE )
  }
  
  if ( !requireNamespace( "multcomp", quietly = TRUE ) ) {
    stop(
      "Package 'multcomp' is required but not installed. Install it with `install.packages(\"multcomp\")`.",
      call. = FALSE
    )
  }
  
  summary.df <- .resolve.df( replicate.summary, "replicate.summary" )
  
  required.cols <- c( "folder", "replicate", metric )
  missing.cols  <- setdiff( required.cols, colnames( summary.df ) )
  if ( length( missing.cols ) > 0 ) {
    stop(
      paste0(
        "`replicate.summary` is missing required column(s): ",
        paste( missing.cols, collapse = ", " )
      ),
      call. = FALSE
    )
  }
  
  if ( !reference.folder %in% summary.df$folder ) {
    stop(
      paste0( "`reference.folder` ('", reference.folder, "') not found in `replicate.summary`." ),
      call. = FALSE
    )
  }
  
  df <- data.frame(
    folder    = as.character( summary.df$folder ),
    replicate = as.character( summary.df$replicate ),
    value     = summary.df[[ metric ]],
    stringsAsFactors = FALSE
  )
  df <- df[ is.finite( df$value ), , drop = FALSE ]
  
  # a replicate seen in only one folder carries no between-folder information
  folders.per.replicate <- tapply( df$folder, df$replicate, function( x ) length( unique( x ) ) )
  dropped <- names( folders.per.replicate )[ folders.per.replicate < 2 ]
  if ( length( dropped ) > 0 ) {
    warning(
      paste0(
        "Dropping replicate(s) present in only one folder: ",
        paste( dropped, collapse = ", " )
      ),
      call. = FALSE
    )
    df <- df[ !df$replicate %in% dropped, , drop = FALSE ]
  }
  
  df$value <- switch(
    transform,
    none  = df$value,
    sqrt  = sqrt( df$value ),
    log1p = log1p( df$value )
  )
  
  if ( any( !is.finite( df$value ) ) ) {
    stop( paste0( "`transform = \"", transform, "\"` produced non-finite values." ), call. = FALSE )
  }
  
  df$folder    <- stats::relevel( factor( df$folder ), ref = reference.folder )
  df$replicate <- factor( df$replicate )
  
  n.folders    <- nlevels( df$folder )
  n.replicates <- nlevels( df$replicate )
  
  if ( n.folders < 2 ) {
    stop( "Fewer than two folders have usable data.", call. = FALSE )
  }
  if ( nrow( df ) - ( n.replicates + n.folders - 1 ) < 1 ) {
    stop( "Too few observations for residual degrees of freedom.", call. = FALSE )
  }
  
  # replicate first, so the sequential F-test for folder is adjusted for it
  fit     <- stats::lm( value ~ replicate + folder, data = df )
  aov.tab <- stats::anova( fit )
  
  omnibus <- data.frame(
    metric       = metric,
    transform    = transform,
    n.folders    = n.folders,
    n.replicates = n.replicates,
    n.obs        = nrow( df ),
    df.folder    = aov.tab[ "folder", "Df" ],
    df.residual  = aov.tab[ "Residuals", "Df" ],
    F.statistic  = aov.tab[ "folder", "F value" ],
    p.value      = aov.tab[ "folder", "Pr(>F)" ],
    stringsAsFactors = FALSE
  )
  
  dunnett <- summary(
    multcomp::glht( fit, linfct = multcomp::mcp( folder = "Dunnett" ) )
  )
  dunnett.ci <- stats::confint(
    multcomp::glht( fit, linfct = multcomp::mcp( folder = "Dunnett" ) ),
    level = conf.level
  )
  
  contrasts <- data.frame(
    folder     = levels( df$folder )[ -1 ],
    reference  = reference.folder,
    metric     = metric,
    transform  = transform,
    estimate   = unname( dunnett$test$coefficients ),
    std.error  = unname( dunnett$test$sigma ),
    statistic  = unname( dunnett$test$tstat ),
    p.adjusted = unname( dunnett$test$pvalues ),
    conf.low   = unname( dunnett.ci$confint[ , "lwr" ] ),
    conf.high  = unname( dunnett.ci$confint[ , "upr" ] ),
    stringsAsFactors = FALSE
  )
  
  utils::write.csv( omnibus,   omnibus.csv,  row.names = FALSE )
  utils::write.csv( contrasts, contrast.csv, row.names = FALSE )
  
  if ( verbose ) {
    message( sprintf( "\033[32mWrote omnibus test: %s\033[0m", omnibus.csv ) )
    message( sprintf( "\033[32mWrote Dunnett contrasts: %s\033[0m", contrast.csv ) )
  }
  
  return( invisible( list( omnibus = omnibus, contrasts = contrasts, model = fit ) ) )
}