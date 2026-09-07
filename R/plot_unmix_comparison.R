# plot_unmix_comparison.R

## Internal helper. Loads the bundled fluorophore_database.csv when the
## caller did not supply one, matching the loading convention used
## elsewhere in the package.
##
## @keywords internal
.load.default.fluorophore.database <- function() {
  fluor.data.path <- system.file(
    "extdata", "fluorophore_database.csv", package = "AutoSpectral"
  )
  fluorophore.database <- utils::read.csv( fluor.data.path )
  fluorophore.database[ fluorophore.database == "" ] <- NA
  fluorophore.database
}


#' @title Plot a Previously Run Unmixing Comparison
#'
#' @description
#' Produces (or reproduces) every figure `compare.unmix.folders()` writes,
#' directly from its `results`/`summary` output, without re-scanning
#' folders or re-unmixing anything. Useful for re-plotting with different
#' cosmetics (colors, sizes, a log axis) or after editing the CSVs by hand,
#' without paying the cost of the unmixing pipeline again.
#'
#' `compare.unmix.folders()` calls this internally for its own figures, so
#' the two never drift out of sync.
#'
#' @param results Either the long-format `results` data frame returned by
#' `compare.unmix.folders()` (or read back in from its `output.csv`), or a
#' character path to that CSV.
#' @param summary Either the per-fluorophore `summary` data frame returned
#' by `compare.unmix.folders()` (or read back in from its `summary.csv`), or
#' a character path to that CSV.
#' @param setup.files Optional named character vector, in the same form
#' passed to `compare.unmix.folders()` (names matching the `folder` values
#' in `results`/`summary`, values ignored here). When supplied, the x-axis
#' of every figure follows `names(setup.files)` in order rather than the
#' default alphabetical ordering. Default `NULL`.
#' @param fluorophore.database Data frame of fluorophore names and metadata
#' (must include `fluorophore`, `excitation.laser`, and `nominal.wavelength`
#' columns), used to build the fluorophore color palette. Default `NULL`
#' loads the bundled `fluorophore_database.csv`.
#' @param plot.dir Character. Directory for output figures. Created if
#' absent. Default `"./figure_unmix_comparison"`.
#' @param plot.per.fluorophore Logical, default `FALSE`. The overview
#' figures (one point per fluorophore, one column per folder) are always
#' produced. Set `TRUE` to additionally produce the more granular
#' one-figure-per-fluorophore plots (one point per off-target channel).
#' @param log.scale Logical, default `FALSE`. Plot the y-axis on a log10
#' scale. Every metric here is non-negative (sum of absolute values or a distance),
#' so this is generally safe, but any exact-zero values are dropped from a
#' log-scale plot (log10(0) is undefined) -- a warning names how many
#' points that affects, per figure, when it happens.
#' @param plot.width,plot.height Numeric, defaults `7` and `5` (inches).
#' Dimensions passed to `ggplot2::ggsave()` for every figure.
#' @param base.font.size Numeric, default `11`. Base font size (points)
#' passed to `ggplot2::theme_classic()`.
#' @param title.size Numeric, default `NULL`. Explicit plot-title font size
#' (points). `NULL` leaves the `theme_classic()` default in place.
#' @param point.size Numeric, default `1.5`. Size of the jittered points.
#' @param point.alpha Numeric, default `0.6`. Opacity of the jittered
#' points, in `[0, 1]`.
#' @param text.angle Numeric, default `45`. Rotation angle (degrees) of the
#' x-axis (folder) labels.
#' @param legend.max.rows Numeric, default `25`. When a figure's legend
#' (fluorophore color key) would need more than this many entries in a
#' single column, it wraps into additional columns rather than
#' overflowing the top of the plot. Ignored when `legend.ncol` is set.
#' @param legend.ncol Numeric, default `NULL`. Fixes the number of legend
#' columns explicitly, overriding `legend.max.rows`. Useful for keeping a
#' consistent legend layout across a batch of figures with varying
#' fluorophore counts.
#' @param legend.font.size Numeric, default `NULL`. Font size (points) for
#' legend text. `NULL` scales it from `base.font.size` (70%).
#' @param legend.key.size Numeric, default `0.8`. Size (lines) of the
#' legend color swatches; smaller values fit more entries per column.
#' @param legend.width.per.col Numeric, default `1.1`. Extra figure width
#' (inches) added for each legend column beyond the first, so a wrapped
#' legend never crowds the plot panel.
#' @param verbose Logical, default `TRUE`.
#'
#' @return Invisibly, the fluorophore color map (named character vector,
#' fluorophore -> hex color) used for the figures.
#'
#' @importFrom utils read.csv
#'
#' @export

unmix.comparison.plot <- function(
    results,
    summary,
    setup.files           = NULL,
    fluorophore.database  = NULL,
    plot.dir              = "./figure_unmix_comparison",
    plot.per.fluorophore  = FALSE,
    log.scale             = FALSE,
    plot.width            = 7,
    plot.height            = 5,
    base.font.size        = 11,
    title.size            = NULL,
    point.size            = 1.5,
    point.alpha           = 0.6,
    text.angle            = 45,
    legend.max.rows       = 25,
    legend.ncol           = NULL,
    legend.font.size      = NULL,
    legend.key.size       = 0.8,
    legend.width.per.col  = 1.1,
    verbose               = TRUE
) {

  results.df <- .resolve.df( results, "results" )
  summary.df <- .resolve.df( summary, "summary" )

  if ( !is.null( setup.files ) ) {

    folder.levels <- names( setup.files )

    missing.levels <- setdiff( unique( results.df$folder ), folder.levels )
    if ( length( missing.levels ) > 0 ) {
      stop(
        paste0(
          "`setup.files` is missing folder(s) present in `results`: ",
          paste( missing.levels, collapse = ", " )
        ),
        call. = FALSE
      )
    }

    results.df$folder <- factor( results.df$folder, levels = folder.levels )
    summary.df$folder <- factor( summary.df$folder, levels = folder.levels )
  }

  if ( is.null( fluorophore.database ) ) {
    fluorophore.database <- .load.default.fluorophore.database()
  }

  fluorophore.color.map <- .build.fluorophore.palette(
    fluorophore.names = c(
      results.df$on.target.fluorophore, results.df$off.target.fluorophore
    ),
    fluorophore.database = fluorophore.database
  )

  if ( !dir.exists( plot.dir ) ) dir.create( plot.dir, recursive = TRUE )

  mad.df <- results.df[ results.df$metric == "Unstained.rSD", ]
  if ( nrow( mad.df ) > 0 ) {
    .plot.metric.boxplot(
      mad.df,
      title          = "Unstained robust SD (MAD)",
      y.label        = "rSD (MAD), unmixed fluorescence channels",
      file.path.out  = file.path( plot.dir, "Unstained_rSD.jpg" ),
      log.scale      = log.scale,
      plot.width     = plot.width,
      plot.height    = plot.height,
      base.font.size = base.font.size,
      title.size     = title.size,
      point.size     = point.size,
      point.alpha    = point.alpha,
      text.angle     = text.angle
    )
  }

  summary.specs <- list(
    SSI             = "Sum of SSI",
    Delta.MFI       = "Sum of spillover errors",
    Spillover.ratio = "Sum of absolute spillover ratio",
    FPR             = "Sum of false-positive rate errors",
    Mahalanobis     = "Median Mahalanobis distance"
  )

  for ( m in names( summary.specs ) ) {

    m.rows <- summary.df[ summary.df$metric == m, ]
    if ( nrow( m.rows ) == 0 ) next

    .plot.metric.boxplot(
      m.rows,
      title                = paste0( summary.specs[[ m ]], " (all fluorophores)" ),
      y.label              = summary.specs[[ m ]],
      file.path.out        = file.path( plot.dir, paste0( "Summary_", m, ".jpg" ) ),
      color.col            = "fluorophore",
      color.map            = fluorophore.color.map,
      log.scale            = log.scale,
      plot.width           = plot.width,
      plot.height          = plot.height,
      base.font.size       = base.font.size,
      title.size           = title.size,
      point.size           = point.size,
      point.alpha          = point.alpha,
      text.angle           = text.angle,
      legend.max.rows      = legend.max.rows,
      legend.ncol          = legend.ncol,
      legend.font.size     = legend.font.size,
      legend.key.size      = legend.key.size,
      legend.width.per.col = legend.width.per.col
    )
  }

  # --- optional per-fluorophore detail plots (off by default; the summary
  # figures above are the intended final-figure output) ---

  if ( plot.per.fluorophore ) {

    boxplot.specs <- list(
      SSI             = "Secondary Stain Index",
      Delta.MFI       = "Spillover (delta MFI)",
      Spillover.ratio = "Spillover ratio ((pos - neg) / neg)",
      FPR             = "False-positive rate"
    )

    for ( m in names( boxplot.specs ) ) {

      metric.rows <- results.df[ results.df$metric == m, ]
      if ( nrow( metric.rows ) == 0 ) next

      for ( fl in unique( metric.rows$on.target.fluorophore ) ) {

        fl.rows <- metric.rows[ metric.rows$on.target.fluorophore == fl, ]

        file.stem <- paste0( m, "_", gsub( "[^A-Za-z0-9_-]", "_", fl ) )

        .plot.metric.boxplot(
          fl.rows,
          title                = paste0( boxplot.specs[[ m ]], ": ", fl ),
          y.label              = boxplot.specs[[ m ]],
          file.path.out        = file.path( plot.dir, paste0( file.stem, ".jpg" ) ),
          color.col            = "off.target.fluorophore",
          color.map            = fluorophore.color.map,
          log.scale            = log.scale,
          plot.width           = plot.width,
          plot.height          = plot.height,
          base.font.size       = base.font.size,
          title.size           = title.size,
          point.size           = point.size,
          point.alpha          = point.alpha,
          text.angle           = text.angle,
          legend.max.rows      = legend.max.rows,
          legend.ncol          = legend.ncol,
          legend.font.size     = legend.font.size,
          legend.key.size      = legend.key.size,
          legend.width.per.col = legend.width.per.col
        )
      }
    }

    mahal.rows <- results.df[ results.df$metric == "Mahalanobis", ]

    for ( fl in unique( mahal.rows$on.target.fluorophore ) ) {

      fl.rows <- mahal.rows[ mahal.rows$on.target.fluorophore == fl, ]

      file.stem <- paste0( "Mahalanobis_", gsub( "[^A-Za-z0-9_-]", "_", fl ) )

      .plot.metric.boxplot(
        fl.rows,
        title                = paste0( "Median Mahalanobis distance: ", fl ),
        y.label              = "Median Mahalanobis distance (off-target channels)",
        file.path.out        = file.path( plot.dir, paste0( file.stem, ".jpg" ) ),
        color.col            = "on.target.fluorophore",
        color.map            = fluorophore.color.map,
        log.scale            = log.scale,
        plot.width           = plot.width,
        plot.height          = plot.height,
        base.font.size       = base.font.size,
        title.size           = title.size,
        point.size           = point.size,
        point.alpha          = point.alpha,
        text.angle           = text.angle,
        legend.max.rows      = legend.max.rows,
        legend.ncol          = legend.ncol,
        legend.font.size     = legend.font.size,
        legend.key.size      = legend.key.size,
        legend.width.per.col = legend.width.per.col
      )
    }
  }

  invisible( fluorophore.color.map )
}
