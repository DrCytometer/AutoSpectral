# test_unmix_comparison.R

## Internal helper. Accepts either a data frame or a path to a CSV and
## returns a data frame either way, so the `test.unmix.comparison*()`
## functions can be called directly on the list returned by
## `compare.unmix.folders()` or on the CSVs it wrote, in a later session.
##
## @keywords internal
.resolve.df <- function( x, arg.name ) {
  
  if ( is.data.frame( x ) ) return( x )
  
  if ( is.character( x ) && length( x ) == 1 ) {
    if ( !file.exists( x ) ) {
      stop( paste0( "`", arg.name, "`: file not found: '", x, "'" ), call. = FALSE )
    }
    return( utils::read.csv( x, stringsAsFactors = FALSE ) )
  }
  
  stop(
    paste0( "`", arg.name, "` must be a data frame or a path to a CSV." ),
    call. = FALSE
  )
}


## Internal helper. Runs one paired test per (folder, metric) against
## `reference.folder`, pairing rows of `pairing.df` by `pair.id` (a
## fluorophore identity for the per-fluorophore comparison, or an
## on-target/off-target channel-pair identity for the channel-level
## comparison). `pairing.df` must have columns `folder`, `pair.id`,
## `metric`, `value`. Returns the long-format detail table (`folder`,
## `metric`, `n.pairs`, `statistic`, `p.value`, `test.method`) with
## `p.adjusted` NOT yet added -- adjustment is applied once by the caller,
## across whatever rows it assembles, so the family of comparisons being
## corrected for is explicit at the call site rather than buried in here.
##
## @keywords internal
.run.paired.tests <- function(
    pairing.df,
    reference.folder,
    test.method,
    min.pairs,
    log.transform = FALSE,
    verbose
) {
  
  if ( log.transform ) {
    if ( any( pairing.df$value <= 0, na.rm = TRUE ) ) {
      stop(
        paste0(
          "`log.transform = TRUE` but some values are <= 0 (log10 is ",
          "undefined for these); the transform is not applied ",
          "automatically to avoid silently dropping observations from a ",
          "paired test. Filter or offset the data first, or leave ",
          "`log.transform = FALSE`."
        ),
        call. = FALSE
      )
    }
    pairing.df$value <- log10( pairing.df$value )
  }
  
  if ( !reference.folder %in% pairing.df$folder ) {
    stop(
      paste0(
        "`reference.folder` ('", reference.folder, "') not found among the ",
        "folders in the supplied data."
      ),
      call. = FALSE
    )
  }
  
  test.folders <- setdiff( unique( pairing.df$folder ), reference.folder )
  metrics      <- unique( pairing.df$metric )
  
  if ( length( test.folders ) == 0 ) {
    stop( "No folders to test other than `reference.folder`.", call. = FALSE )
  }
  
  detail.rows <- list()
  
  for ( m in metrics ) {
    
    ref.rows <- pairing.df[
      pairing.df$folder == reference.folder & pairing.df$metric == m,
    ]
    
    for ( fld in test.folders ) {
      
      fld.rows <- pairing.df[
        pairing.df$folder == fld & pairing.df$metric == m,
      ]
      
      common.id <- intersect( ref.rows$pair.id, fld.rows$pair.id )
      n.pairs   <- length( common.id )
      
      p.raw     <- NA_real_
      statistic <- NA_real_
      
      if ( n.pairs >= min.pairs ) {
        
        x <- fld.rows$value[ match( common.id, fld.rows$pair.id ) ]
        y <- ref.rows$value[  match( common.id, ref.rows$pair.id ) ]
        
        test.result <- tryCatch(
          {
            if ( test.method == "wilcox" ) {
              stats::wilcox.test( x, y, paired = TRUE, exact = FALSE )
            } else {
              stats::t.test( x, y, paired = TRUE )
            }
          },
          error = function( e ) NULL
        )
        
        if ( !is.null( test.result ) ) {
          p.raw     <- test.result$p.value
          statistic <- unname( test.result$statistic )
        } else {
          warning(
            paste0(
              "Folder '", fld, "', metric '", m,
              "': paired test failed (likely zero variance in the ",
              "differences); recorded as NA."
            ),
            call. = FALSE
          )
        }
        
      } else if ( verbose ) {
        message( sprintf(
          "\033[33mFolder '%s', metric '%s': only %d paired observation(s) (< min.pairs = %d); recorded as NA.\033[0m",
          fld, m, n.pairs, min.pairs
        ) )
      }
      
      detail.rows[[ length( detail.rows ) + 1 ]] <- data.frame(
        folder      = fld,
        metric      = m,
        n.pairs     = n.pairs,
        statistic   = statistic,
        p.value     = p.raw,
        test.method = test.method,
        stringsAsFactors = FALSE
      )
    }
  }
  
  detail.df <- do.call( rbind, detail.rows )
  rownames( detail.df ) <- NULL
  
  detail.df
}


## Internal helper. Builds the folder x metric wide matrix of adjusted
## p-values from a `detail.df` (as returned by `.run.paired.tests()`, after
## `p.adjusted` has been added), writes both tables to CSV, and returns them.
##
## @keywords internal
.write.paired.test.outputs <- function(
    detail.df, test.folders, metrics, output.csv, detail.csv, verbose
) {
  
  utils::write.csv( detail.df, detail.csv, row.names = FALSE )
  if ( verbose ) message( sprintf( "\033[32mWrote detail: %s\033[0m", detail.csv ) )
  
  wide.df <- data.frame( folder = test.folders, stringsAsFactors = FALSE )
  for ( m in metrics ) {
    m.rows         <- detail.df[ detail.df$metric == m, ]
    wide.df[[ m ]] <- m.rows$p.adjusted[ match( wide.df$folder, m.rows$folder ) ]
  }
  
  utils::write.csv( wide.df, output.csv, row.names = FALSE )
  if ( verbose ) message( sprintf( "\033[32mWrote stats: %s\033[0m", output.csv ) )
  
  list( wide = wide.df, detail = detail.df )
}


#' @title Test Whether Each Unmixing Method Differs From a Reference
#'
#' @description
#' For every metric produced by `compare.unmix.folders()`, runs one paired
#' test per folder against a chosen reference folder (typically
#' `"AutoSpectral"`), pairing observations by fluorophore. Six metrics are
#' tested: the five already in the per-fluorophore `summary` table (`SSI`,
#' `Delta.MFI`, `Spillover.ratio`, `FPR`, `Mahalanobis`) plus
#' `Unstained.rSD`, pulled in from the long-format `results` table since it
#' is per-fluorophore but was never folded into `summary` (it has no
#' per-channel spread to summarize away).
#'
#' P-values are adjusted once across every (folder, metric) test performed
#' in the call (`stats::p.adjust()`), not separately per metric or per
#' folder, since that is the full family of comparisons being made in one
#' call. This function only writes CSVs; it does not modify or add to the
#' `compare.unmix.folders()` figures.
#'
#' See also `test.unmix.comparison.channels()`, which runs the same
#' procedure directly on `results` (one off-target channel at a time)
#' rather than on the fluorophore-level `summary`.
#'
#' @param results Either the long-format `results` data frame returned by
#' `compare.unmix.folders()` (or read back in from its `output.csv`), or a
#' character path to that CSV.
#' @param summary Either the per-fluorophore `summary` data frame returned
#' by `compare.unmix.folders()` (or read back in from its `summary.csv`), or
#' a character path to that CSV.
#' @param reference.folder Character. Name of the folder every other folder
#' is compared against. Default `"AutoSpectral"`.
#' @param test.method Character, one of `"wilcox"` (paired Wilcoxon
#' signed-rank test, the default) or `"t.test"` (paired t-test). Wilcoxon is
#' the safer default given the typically small, non-normal, non-negative
#' per-fluorophore metrics involved.
#' @param p.adjust.method Character, passed to `stats::p.adjust()`. Default
#' `"BH"` (Benjamini-Hochberg). Must be one of `stats::p.adjust.methods`.
#' @param min.pairs Numeric, default `4`. Minimum number of fluorophores
#' present in both a folder and `reference.folder`, for a given metric,
#' before a test is attempted; below this, the comparison is recorded as
#' `NA` (with a message) rather than run on too few pairs to be meaningful.
#' @param log.transform Logical, default `FALSE`. Log10-transform `value`
#' before pairing and testing. This changes what is being tested: on the
#' raw scale the test asks whether the typical *absolute difference*
#' between a folder and `reference.folder` is zero; on the log scale it
#' asks whether the typical *ratio* is 1, which is often the more
#' appropriate question when a metric spans a wide range across
#' fluorophores. All values must be strictly positive when `TRUE` --
#' `stop()`s otherwise, rather than silently dropping observations from the
#' test.
#' @param output.csv Character. Path for the wide adjusted-p-value matrix
#' (folders in rows, metrics in columns). Default
#' `"unmix_comparison_stats.csv"`.
#' @param detail.csv Character. Path for the long-format detail table (one
#' row per folder/metric, with `n.pairs`, the test statistic, raw
#' `p.value`, and `p.adjusted`) — this is what explains *why* a cell in
#' `output.csv` is `NA` (too few pairs vs. a failed test). Default
#' `"unmix_comparison_stats_detail.csv"`.
#' @param verbose Logical, default `TRUE`.
#'
#' @return Invisibly, a named list with `wide` (the folder x metric adjusted
#' p-value matrix also written to `output.csv`) and `detail` (the
#' long-format table also written to `detail.csv`).
#'
#' @importFrom stats wilcox.test t.test p.adjust
#'
#' @export

test.unmix.comparison <- function(
    results,
    summary,
    reference.folder = "AutoSpectral",
    test.method       = c( "wilcox", "t.test" ),
    p.adjust.method   = "BH",
    min.pairs         = 4,
    log.transform     = FALSE,
    output.csv        = "unmix_comparison_stats.csv",
    detail.csv        = "unmix_comparison_stats_detail.csv",
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
  summary.df <- .resolve.df( summary, "summary" )
  
  required.results.cols <- c( "folder", "off.target.fluorophore", "metric", "value" )
  missing.results.cols  <- setdiff( required.results.cols, colnames( results.df ) )
  if ( length( missing.results.cols ) > 0 ) {
    stop(
      paste0(
        "`results` is missing required column(s): ",
        paste( missing.results.cols, collapse = ", " )
      ),
      call. = FALSE
    )
  }
  
  required.summary.cols <- c( "folder", "fluorophore", "metric", "value" )
  missing.summary.cols  <- setdiff( required.summary.cols, colnames( summary.df ) )
  if ( length( missing.summary.cols ) > 0 ) {
    stop(
      paste0(
        "`summary` is missing required column(s): ",
        paste( missing.summary.cols, collapse = ", " )
      ),
      call. = FALSE
    )
  }
  
  # --- assemble one long table: folder, pair.id, metric, value -------------
  # `pair.id` is fluorophore identity here (one paired observation per
  # fluorophore per metric).
  
  unstained.rows <- results.df[ results.df$metric == "Unstained.rSD", ]
  unstained.long <- data.frame(
    folder  = unstained.rows$folder,
    pair.id = unstained.rows$off.target.fluorophore,
    metric  = "Unstained.rSD",
    value   = unstained.rows$value,
    stringsAsFactors = FALSE
  )
  
  summary.long <- data.frame(
    folder  = summary.df$folder,
    pair.id = summary.df$fluorophore,
    metric  = summary.df$metric,
    value   = summary.df$value,
    stringsAsFactors = FALSE
  )
  
  pairing.df <- rbind( summary.long, unstained.long )
  
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


#' @title Test Whether Each Unmixing Method Differs From a Reference, Per
#' Off-Target Channel
#'
#' @description
#' The channel-level counterpart to `test.unmix.comparison()`: runs the same
#' paired-test-per-folder procedure, but directly on the full `results`
#' table rather than the per-fluorophore `summary` table, pairing
#' observations by the combination of on-target fluorophore and off-target
#' fluorophore (i.e. by spectral channel pair) rather than by fluorophore
#' alone. This uses every individual off-target-channel measurement instead
#' of first collapsing each fluorophore's off-target channels down to a
#' single rSD (MAD) value.
#'
#' Only applies to the four metrics that have both an on-target and an
#' off-target fluorophore in `results` (`SSI`, `Delta.MFI`,
#' `Spillover.ratio`, `FPR`); `Unstained.rSD` has no on-target fluorophore
#' (there is only one measurement per fluorophore channel, with no further
#' pairing dimension available) and `Mahalanobis` has no off-target
#' fluorophore (it is already a single summary value per fluorophore), so
#' both are excluded automatically by requiring both columns to be
#' non-missing.
#'
#' \strong{Caveat:} the paired observations at this grain are not fully
#' independent -- every off-target channel belonging to the same on-target
#' fluorophore is measured on that fluorophore's single set of positive/
#' negative events, so within-fluorophore correlation inflates the
#' effective sample size the test sees. This will tend to produce smaller,
#' more "significant" p-values than a genuinely independent sample of the
#' same size would justify. Treat this as a more sensitive, exploratory
#' companion to `test.unmix.comparison()` (which pairs by fluorophore and
#' is the statistically conservative version), not a replacement for it.
#'
#' @param results Either the long-format `results` data frame returned by
#' `compare.unmix.folders()` (or read back in from its `output.csv`), or a
#' character path to that CSV.
#' @param reference.folder Character. Name of the folder every other folder
#' is compared against. Default `"AutoSpectral"`.
#' @param test.method Character, one of `"wilcox"` (paired Wilcoxon
#' signed-rank test, the default) or `"t.test"` (paired t-test).
#' @param p.adjust.method Character, passed to `stats::p.adjust()`. Default
#' `"BH"` (Benjamini-Hochberg). Must be one of `stats::p.adjust.methods`.
#' @param min.pairs Numeric, default `4`. Minimum number of on-target/
#' off-target channel pairs present in both a folder and
#' `reference.folder`, for a given metric, before a test is attempted;
#' below this, the comparison is recorded as `NA` (with a message).
#' @param log.transform Logical, default `FALSE`. Log10-transform `value`
#' before pairing and testing -- see `test.unmix.comparison()` for what
#' this changes about the hypothesis being tested. All values must be
#' strictly positive when `TRUE` -- `stop()`s otherwise.
#' @param output.csv Character. Path for the wide adjusted-p-value matrix
#' (folders in rows, metrics in columns). Default
#' `"unmix_comparison_stats_channels.csv"`.
#' @param detail.csv Character. Path for the long-format detail table (one
#' row per folder/metric, with `n.pairs`, the test statistic, raw
#' `p.value`, and `p.adjusted`). Default
#' `"unmix_comparison_stats_channels_detail.csv"`.
#' @param verbose Logical, default `TRUE`.
#'
#' @return Invisibly, a named list with `wide` (the folder x metric adjusted
#' p-value matrix also written to `output.csv`) and `detail` (the
#' long-format table also written to `detail.csv`).
#'
#' @importFrom stats wilcox.test t.test p.adjust
#'
#' @export

test.unmix.comparison.channels <- function(
    results,
    reference.folder = "AutoSpectral",
    test.method       = c( "wilcox", "t.test" ),
    p.adjust.method   = "BH",
    min.pairs         = 4,
    log.transform     = FALSE,
    output.csv        = "unmix_comparison_stats_channels.csv",
    detail.csv        = "unmix_comparison_stats_channels_detail.csv",
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
    "folder", "on.target.fluorophore", "off.target.fluorophore", "metric", "value"
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
  
  # --- restrict to rows with both an on-target and off-target fluorophore --
  # (drops Unstained.rSD, which has no on-target fluorophore, and
  # Mahalanobis, which has no off-target fluorophore)
  
  channel.rows <- results.df[
    !is.na( results.df$on.target.fluorophore ) &
      !is.na( results.df$off.target.fluorophore ),
  ]
  
  if ( nrow( channel.rows ) == 0 ) {
    stop(
      "No rows in `results` have both an on-target and off-target fluorophore.",
      call. = FALSE
    )
  }
  
  pairing.df <- data.frame(
    folder = channel.rows$folder,
    pair.id = paste(
      channel.rows$on.target.fluorophore, channel.rows$off.target.fluorophore,
      sep = " -> "
    ),
    metric = channel.rows$metric,
    value  = channel.rows$value,
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