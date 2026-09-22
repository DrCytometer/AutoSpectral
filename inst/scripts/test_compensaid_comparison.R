# test_compensaid_comparison.R
#
# Paired statistical testing for compare.compensaid.folders() output,
# mirroring test_unmix_comparison.R's two-tier design (test.unmix.comparison()
# vs. test.unmix.comparison.channels()) one level down: a conservative
# replicate-level test on aggregate flagged counts, and a more sensitive,
# non-independent pair-level test on individual SSI values.
#
# Depends on .resolve.df(), .run.paired.tests(), and
# .write.paired.test.outputs() (all defined in test_unmix_comparison.R --
# despite the filename, these are production dependencies, not test-only
# code). test_unmix_comparison.R must be loaded into the package namespace
# alongside this file.


#' @title Test Whether Each Unmixing Method's CompensAID Flags Differ From a
#' Reference, Per Replicate
#'
#' @description
#' Runs one paired test per folder against a chosen reference folder
#' (typically `"AutoSpectral"`), pairing observations by replicate identity
#' -- the same physical specimen, unmixed by every operator/method, as set
#' up by `compare.compensaid.folders()`'s `sample.files` argument. This is
#' the conservative, replicate-level test: each replicate contributes one
#' paired observation per metric, so within-replicate correlation across
#' different marker combinations cannot inflate the result the way it can
#' in `test.compensaid.comparison.pairs()`.
#'
#' One test is run per (folder, metric) pair, where the metrics are the
#' total number of flagged marker combinations (`n.flagged`), the fraction
#' of tested combinations flagged (`frac.flagged`), and the count of flagged
#' combinations in each severity band named in `severity.labels`. P-values
#' are adjusted once across every (folder, metric) test performed in the
#' call (`stats::p.adjust()`).
#'
#' @param replicate.summary Either the `replicate.summary` data frame
#' returned by `compare.compensaid.folders()` (or read back in from its
#' `replicate.summary.csv`), or a character path to that CSV.
#' @param severity.labels Character vector of severity band names, matching
#' the columns present in `replicate.summary`. Default
#' `c("Severe", "Moderate", "Mild")`, matching
#' `compare.compensaid.folders()`'s own default.
#' @param reference.folder Character. Name of the folder every other folder
#' is compared against. Default `"AutoSpectral"`.
#' @param test.method Character, one of `"wilcox"` (paired Wilcoxon
#' signed-rank test, the default) or `"t.test"` (paired t-test).
#' @param p.adjust.method Character, passed to `stats::p.adjust()`. Default
#' `"BH"` (Benjamini-Hochberg). Must be one of `stats::p.adjust.methods`.
#' @param min.pairs Numeric, default `4`. Minimum number of replicates
#' present in both a folder and `reference.folder`, for a given metric,
#' before a test is attempted; below this, the comparison is recorded as
#' `NA` (with a message).
#' @param log.transform Logical, default `FALSE`. Log10-transform `value`
#' before pairing and testing -- see `test.unmix.comparison()` for what this
#' changes about the hypothesis being tested. All values must be strictly
#' positive when `TRUE` -- `stop()`s otherwise (this excludes `n.flagged`
#' and the severity-band counts whenever any replicate has zero flags).
#' @param output.csv Character. Path for the wide adjusted-p-value matrix
#' (folders in rows, metrics in columns). Default
#' `"compensaid_comparison_stats.csv"`.
#' @param detail.csv Character. Path for the long-format detail table (one
#' row per folder/metric, with `n.pairs`, the test statistic, raw
#' `p.value`, and `p.adjusted`). Default
#' `"compensaid_comparison_stats_detail.csv"`.
#' @param verbose Logical, default `TRUE`.
#'
#' @return Invisibly, a named list with `wide` (the folder x metric adjusted
#' p-value matrix also written to `output.csv`) and `detail` (the
#' long-format table also written to `detail.csv`).
#'
#' @importFrom stats wilcox.test t.test p.adjust
#'
#' @seealso [compare.compensaid.folders()], [test.compensaid.comparison.pairs()]
#'
#' @export

test.compensaid.comparison <- function(
    replicate.summary,
    severity.labels   = c( "Severe", "Moderate", "Mild" ),
    reference.folder  = "AutoSpectral",
    test.method       = c( "wilcox", "t.test" ),
    p.adjust.method   = "BH",
    min.pairs         = 4,
    log.transform     = FALSE,
    output.csv        = "compensaid_comparison_stats.csv",
    detail.csv        = "compensaid_comparison_stats_detail.csv",
    verbose           = TRUE
) {

  test.method <- match.arg( test.method )

  if ( !p.adjust.method %in% stats::p.adjust.methods ) {
    stop(
      paste0(
        "`p.adjust.method` must be one of: ",
        paste( stats::p.adjust.methods, collapse = ", " )
      ),
      call. = FALSE
    )
  }

  summary.df <- .resolve.df( replicate.summary, "replicate.summary" )

  required.cols <- c( "folder", "replicate", "n.flagged", "frac.flagged" )
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

  missing.severity.cols <- setdiff( severity.labels, colnames( summary.df ) )
  if ( length( missing.severity.cols ) > 0 ) {
    stop(
      paste0(
        "`replicate.summary` is missing severity column(s): ",
        paste( missing.severity.cols, collapse = ", " ),
        ". Pass the `severity.labels` actually used to build it."
      ),
      call. = FALSE
    )
  }

  # --- assemble one long table: folder, pair.id, metric, value -------------
  # `pair.id` is replicate identity here (one paired observation per
  # replicate per metric).

  metric.cols <- c( "n.flagged", "frac.flagged", severity.labels )

  pairing.rows <- list()
  for ( m in metric.cols ) {
    pairing.rows[[ m ]] <- data.frame(
      folder  = summary.df$folder,
      pair.id = summary.df$replicate,
      metric  = m,
      value   = summary.df[[ m ]],
      stringsAsFactors = FALSE
    )
  }
  pairing.df <- do.call( rbind, pairing.rows )

  detail.df <- .run.paired.tests(
    pairing.df       = pairing.df,
    reference.folder = reference.folder,
    test.method      = test.method,
    min.pairs        = min.pairs,
    log.transform    = log.transform,
    verbose          = verbose
  )

  # --- adjust across every (folder, metric) test performed in this call ----

  detail.df$p.adjusted <- stats::p.adjust( detail.df$p.value, method = p.adjust.method )

  test.folders <- setdiff( unique( pairing.df$folder ), reference.folder )
  metrics      <- unique( pairing.df$metric )

  out <- .write.paired.test.outputs(
    detail.df    = detail.df,
    test.folders = test.folders,
    metrics      = metrics,
    output.csv   = output.csv,
    detail.csv   = detail.csv,
    verbose      = verbose
  )

  return( invisible( out ) )
}


#' @title Test Whether Each Unmixing Method's CompensAID Flags Differ From a
#' Reference, Per Marker Combination
#'
#' @description
#' The finer-grained counterpart to `test.compensaid.comparison()`: runs the
#' same paired-test-per-folder procedure, but directly on the full `results`
#' table's continuous SSI values rather than on `replicate.summary`'s
#' per-replicate counts, pairing observations by the combination of
#' replicate identity, primary fluorophore, and secondary fluorophore.
#'
#' \strong{Caveat:} the paired observations at this grain are not
#' independent -- every marker combination belonging to the same replicate
#' is measured on that one sample's single CompensAID run, so
#' within-replicate correlation across marker combinations inflates the
#' effective sample size the test sees. This will tend to produce smaller,
#' more "significant" p-values than a genuinely independent sample of the
#' same size would justify. Treat this as a more sensitive, exploratory
#' companion to `test.compensaid.comparison()` (which pairs by replicate and
#' is the statistically conservative version), not a replacement for it.
#'
#' @param results Either the long-format `results` data frame returned by
#' `compare.compensaid.folders()` (or read back in from its `output.csv`),
#' or a character path to that CSV.
#' @param reference.folder Character. Name of the folder every other folder
#' is compared against. Default `"AutoSpectral"`.
#' @param test.method Character, one of `"wilcox"` (paired Wilcoxon
#' signed-rank test, the default) or `"t.test"` (paired t-test).
#' @param p.adjust.method Character, passed to `stats::p.adjust()`. Default
#' `"BH"` (Benjamini-Hochberg). Must be one of `stats::p.adjust.methods`.
#' @param min.pairs Numeric, default `4`. Minimum number of (replicate,
#' primary, secondary) combinations present in both a folder and
#' `reference.folder` before a test is attempted; below this, the
#' comparison is recorded as `NA` (with a message).
#' @param log.transform Logical, default `FALSE`. Log10-transform `value`
#' before pairing and testing -- see `test.unmix.comparison()` for what
#' this changes about the hypothesis being tested. SSI values are commonly
#' negative, so `TRUE` will usually `stop()` here unless the data are
#' shifted or filtered first.
#' @param output.csv Character. Path for the wide adjusted-p-value matrix
#' (folders in rows, a single `SSI` column). Default
#' `"compensaid_comparison_stats_pairs.csv"`.
#' @param detail.csv Character. Path for the long-format detail table (one
#' row per folder, with `n.pairs`, the test statistic, raw `p.value`, and
#' `p.adjusted`). Default `"compensaid_comparison_stats_pairs_detail.csv"`.
#' @param verbose Logical, default `TRUE`.
#'
#' @return Invisibly, a named list with `wide` (the folder x metric adjusted
#' p-value matrix also written to `output.csv`) and `detail` (the
#' long-format table also written to `detail.csv`).
#'
#' @importFrom stats wilcox.test t.test p.adjust
#'
#' @seealso [compare.compensaid.folders()], [test.compensaid.comparison()]
#'
#' @export

test.compensaid.comparison.pairs <- function(
    results,
    reference.folder  = "AutoSpectral",
    test.method       = c( "wilcox", "t.test" ),
    p.adjust.method   = "BH",
    min.pairs         = 4,
    log.transform     = FALSE,
    output.csv        = "compensaid_comparison_stats_pairs.csv",
    detail.csv        = "compensaid_comparison_stats_pairs_detail.csv",
    verbose           = TRUE
) {

  test.method <- match.arg( test.method )

  if ( !p.adjust.method %in% stats::p.adjust.methods ) {
    stop(
      paste0(
        "`p.adjust.method` must be one of: ",
        paste( stats::p.adjust.methods, collapse = ", " )
      ),
      call. = FALSE
    )
  }

  results.df <- .resolve.df( results, "results" )

  required.cols <- c(
    "folder", "replicate", "primary.fluorophore", "secondary.fluorophore", "ssi"
  )
  missing.cols <- setdiff( required.cols, colnames( results.df ) )
  if ( length( missing.cols ) > 0 ) {
    stop(
      paste0(
        "`results` is missing required column(s): ",
        paste( missing.cols, collapse = ", " )
      ),
      call. = FALSE
    )
  }

  pairing.df <- data.frame(
    folder = results.df$folder,
    pair.id = paste(
      results.df$replicate, results.df$primary.fluorophore,
      results.df$secondary.fluorophore, sep = " | "
    ),
    metric = "SSI",
    value  = results.df$ssi,
    stringsAsFactors = FALSE
  )

  detail.df <- .run.paired.tests(
    pairing.df       = pairing.df,
    reference.folder = reference.folder,
    test.method      = test.method,
    min.pairs        = min.pairs,
    log.transform    = log.transform,
    verbose          = verbose
  )

  detail.df$p.adjusted <- stats::p.adjust( detail.df$p.value, method = p.adjust.method )

  test.folders <- setdiff( unique( pairing.df$folder ), reference.folder )
  metrics      <- unique( pairing.df$metric )

  out <- .write.paired.test.outputs(
    detail.df    = detail.df,
    test.folders = test.folders,
    metrics      = metrics,
    output.csv   = output.csv,
    detail.csv   = detail.csv,
    verbose      = verbose
  )

  return( invisible( out ) )
}
