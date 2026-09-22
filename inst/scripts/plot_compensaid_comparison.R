# plot_compensaid_comparison.R
#
# Re-plotting helper for compare.compensaid.folders(): reproduces the
# barchart(s) that function writes directly from its `results`/`summary`
# output, without re-running CompensAID on anything.
# compare.compensaid.folders() calls this function internally for its own
# figures, so the two never drift out of sync.
#
# Depends on .resolve.df() (defined in test_unmix_comparison.R -- despite
# the filename, this is a production dependency, not test-only code). Both
# files must be loaded into the package namespace together.

#' @title Plot a Previously Run CompensAID Comparison
#'
#' @description
#' Produces (or reproduces) the barchart(s) [compare.compensaid.folders()]
#' writes, directly from its `results`/`summary` output, without re-running
#' CompensAID. Useful for re-plotting with different cosmetics, a different
#' folder order, or after manually editing `summary.csv`.
#'
#' Always produces one figure: the mean number of CompensAID-flagged marker
#' combinations per folder (across that folder's replicates), as a stacked
#' barchart with one bar segment per severity band (most severe at the
#' base). Error bars on the top of each stack show the dispersion of the
#' total flagged count across replicates (`se.n.flagged`; SEM or SD
#' according to how `summary` was built), and are omitted for folders with a
#' single replicate. Set `normalize = TRUE` to additionally produce a second
#' figure of the mean percentage of tested pairs flagged per folder -- useful
#' when folders don't all have the same number of usable fluorophores, so a
#' raw count isn't directly comparable.
#'
#' @param results Either the long-format `results` data frame returned by
#' `compare.compensaid.folders()` (or read back in from its `output.csv`),
#' or a character path to that CSV. Currently unused by the plot itself
#' (the barchart is built entirely from `summary`) but accepted for
#' interface symmetry with `compare.compensaid.folders()` and to allow
#' future per-sample figures without a signature change.
#' @param summary Either the per-folder `summary` data frame returned by
#' `compare.compensaid.folders()` (or read back in from its `summary.csv`),
#' or a character path to that CSV. Must contain `folder`,
#' `mean.n.flagged`, `mean.frac.flagged`, and one `mean.<label>` column per
#' entry of `severity.labels`. The matching `se.n.flagged` and
#' `se.frac.flagged` columns are used for error bars when present.
#' @param folder.levels Optional character vector of folder labels. When
#' supplied, the x-axis follows this order rather than the default
#' alphabetical ordering. Default `NULL`.
#' @param severity.labels Character vector of severity band labels, most
#' severe first, matching the `mean.<label>` columns present in `summary`.
#' Default `c("Severe", "Moderate", "Mild")`.
#' @param severity.colors Named character vector of hex colors, names
#' matching (a subset of) `severity.labels`. Default `NULL` uses a built-in
#' red-shade ramp (darkest = most severe); any label not covered by the
#' default is assigned an additional shade automatically.
#' @param error.label Character, default `"SEM"`. What the `se.*` columns of
#' `summary` represent (`"SEM"` or `"SD"`); used only for the error-bar
#' caption. It is not recorded in `summary`, so it must match the
#' `error.type` the comparison was run with.
#' @param plot.dir Character. Directory for output figures. Created if
#' absent. Default `"./figure_compensaid_comparison"`.
#' @param normalize Logical, default `FALSE`. Adds the percentage-flagged
#' figure described above.
#' @param plot.width,plot.height Numeric, defaults `7` and `5` (inches).
#' @param base.font.size Numeric, default `11`.
#' @param title.size Numeric, default `NULL`. `NULL` leaves the
#' `theme_classic()` default in place.
#' @param text.angle Numeric, default `45`. Rotation angle (degrees) of the
#' x-axis (folder) labels.
#' @param verbose Logical, default `TRUE`.
#'
#' @return Invisibly, a named list with `count.plot` (always) and
#' `fraction.plot` (only when `normalize = TRUE`), the ggplot objects
#' produced. Figures are written to `plot.dir` and printed to the active
#' graphics device.
#'
#' @importFrom ggplot2 ggplot aes geom_col geom_errorbar position_stack
#' @importFrom ggplot2 scale_fill_manual labs
#' @importFrom ggplot2 theme_classic theme element_text ggsave
#' @importFrom ragg agg_jpeg
#' @importFrom stats setNames
#' @importFrom utils read.csv stack
#' @importFrom grDevices hcl.colors
#'
#' @seealso [compare.compensaid.folders()]
#'
#' @export

plot.compensaid.comparison <- function(
    results          = NULL,
    summary,
    folder.levels    = NULL,
    severity.labels  = c( "Severe", "Moderate", "Mild" ),
    severity.colors  = NULL,
    error.label      = "SEM",
    plot.dir         = "./figure_compensaid_comparison",
    normalize        = FALSE,
    plot.width       = 7,
    plot.height      = 5,
    base.font.size   = 11,
    title.size       = NULL,
    text.angle       = 45,
    verbose          = TRUE
) {

  summary.df <- .resolve.df( summary, "summary" )

  required.cols <- c( "folder", "mean.n.flagged", "mean.frac.flagged" )
  missing.cols  <- setdiff( required.cols, colnames( summary.df ) )
  if ( length( missing.cols ) > 0 ) {
    stop(
      paste0(
        "`summary` is missing required column(s): ",
        paste( missing.cols, collapse = ", " )
      ),
      call. = FALSE
    )
  }

  if ( !is.null( folder.levels ) ) {
    summary.df$folder <- factor( summary.df$folder, levels = folder.levels )
  }

  present.severity <- severity.labels[
    paste0( "mean.", severity.labels ) %in% colnames( summary.df )
  ]
  if ( length( present.severity ) == 0 ) {
    stop(
      paste0(
        "`summary` does not contain any of the expected severity columns (",
        paste( paste0( "mean.", severity.labels ), collapse = ", " ), ")."
      ),
      call. = FALSE
    )
  }

  long.df <- do.call( rbind, lapply( present.severity, function( s ) {
    data.frame(
      folder   = summary.df$folder,
      severity = s,
      n        = summary.df[[ paste0( "mean.", s ) ]],
      stringsAsFactors = FALSE
    )
  } ) )
  long.df$severity <- factor( long.df$severity, levels = present.severity )

  se.count <- if ( "se.n.flagged" %in% colnames( summary.df ) ) {
    summary.df$se.n.flagged
  } else {
    rep( NA_real_, nrow( summary.df ) )
  }

  count.err <- data.frame(
    folder = summary.df$folder,
    top    = summary.df$mean.n.flagged,
    se     = se.count,
    stringsAsFactors = FALSE
  )
  count.err <- count.err[ !is.na( count.err$se ), , drop = FALSE ]

  if ( is.null( severity.colors ) ) {

    default.colors <- c(
      Severe   = "#B2182B",
      Moderate = "#F4A582",
      Mild     = "#FDDBC7"
    )

    severity.colors <- default.colors[ present.severity ]
    unmapped        <- present.severity[ !present.severity %in% names( default.colors ) ]

    if ( length( unmapped ) > 0 ) {
      extra <- grDevices::hcl.colors( length( unmapped ), palette = "Reds 2" )
      severity.colors[ unmapped ] <- extra
    }
  }

  if ( !dir.exists( plot.dir ) ) dir.create( plot.dir, recursive = TRUE )

  p.count <- ggplot2::ggplot(
    long.df, ggplot2::aes( x = folder, y = n, fill = severity )
  ) +
    ggplot2::geom_col( position = ggplot2::position_stack( reverse = TRUE ) ) +
    ggplot2::scale_fill_manual( values = severity.colors, name = "Severity" ) +
    ggplot2::labs(
      title = "CompensAID-flagged marker combinations per operator",
      x     = "Operator / method",
      y     = "Mean number of flagged marker combinations"
    ) +
    ggplot2::theme_classic( base_size = base.font.size ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text( angle = text.angle, hjust = 1 )
    )

  if ( nrow( count.err ) > 0 ) {
    p.count <- p.count +
      ggplot2::geom_errorbar(
        data = count.err,
        ggplot2::aes( x = folder, ymin = top - se, ymax = top + se ),
        width = 0.25,
        inherit.aes = FALSE
      ) +
      ggplot2::labs( caption = paste0( "Error bars: ", error.label, " across replicates" ) )
  }

  if ( !is.null( title.size ) ) {
    p.count <- p.count +
      ggplot2::theme( plot.title = ggplot2::element_text( size = title.size ) )
  }

  count.file <- file.path( plot.dir, "CompensAID_flagged_count.jpg" )
  ggplot2::ggsave(
    count.file, p.count, device = ragg::agg_jpeg, width = plot.width, height = plot.height
  )
  print( p.count )
  if ( verbose ) message( sprintf( "\033[32mWrote: %s\033[0m", count.file ) )

  out <- list( count.plot = p.count )

  if ( normalize ) {

    summary.df$pct.flagged <- 100 * summary.df$mean.frac.flagged

    se.frac <- if ( "se.frac.flagged" %in% colnames( summary.df ) ) {
      100 * summary.df$se.frac.flagged
    } else {
      rep( NA_real_, nrow( summary.df ) )
    }

    frac.err <- data.frame(
      folder = summary.df$folder,
      top    = summary.df$pct.flagged,
      se     = se.frac,
      stringsAsFactors = FALSE
    )
    frac.err <- frac.err[ !is.na( frac.err$se ), , drop = FALSE ]

    p.frac <- ggplot2::ggplot(
      summary.df, ggplot2::aes( x = folder, y = pct.flagged )
    ) +
      ggplot2::geom_col( fill = "grey40" ) +
      ggplot2::labs(
        title = "CompensAID-flagged fraction of tested marker combinations",
        x     = "Operator / method",
        y     = "Mean flagged combinations (%)"
      ) +
      ggplot2::theme_classic( base_size = base.font.size ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text( angle = text.angle, hjust = 1 )
      )

    if ( nrow( frac.err ) > 0 ) {
      p.frac <- p.frac +
        ggplot2::geom_errorbar(
          data = frac.err,
          ggplot2::aes( x = folder, ymin = top - se, ymax = top + se ),
          width = 0.25,
          inherit.aes = FALSE
        ) +
        ggplot2::labs( caption = paste0( "Error bars: ", error.label, " across replicates" ) )
    }

    if ( !is.null( title.size ) ) {
      p.frac <- p.frac +
        ggplot2::theme( plot.title = ggplot2::element_text( size = title.size ) )
    }

    frac.file <- file.path( plot.dir, "CompensAID_flagged_fraction.jpg" )
    ggplot2::ggsave(
      frac.file, p.frac, device = ragg::agg_jpeg, width = plot.width, height = plot.height
    )
    print( p.frac )
    if ( verbose ) message( sprintf( "\033[32mWrote: %s\033[0m", frac.file ) )

    out$fraction.plot <- p.frac
  }

  invisible( out )
}
