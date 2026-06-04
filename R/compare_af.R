#' Compare Autofluorescence Spectra Sets
#'
#' Evaluates the quality of multiple AF spectra matrices -- for example,
#' produced by different parameter settings of \code{\link{get.af.spectra}} --
#' against the same unstained FCS file. A single AF assignment function is
#' applied consistently across every entry in \code{af.spectra.list}, so that
#' differences in cosine similarity reflect the AF spectra themselves rather
#' than the assignment strategy.
#'
#' For each \code{af.spectra} matrix the function:
#' \enumerate{
#'   \item Assigns each cell to its best-matching AF spectrum row using
#'     \code{assign.fn}.
#'   \item Computes the cosine similarity between the cell's raw detector
#'     signal and its assigned AF spectrum.
#'   \item Summarises per-cell similarities into \code{Mean_Sim} and
#'     \code{rSD_Sim}.
#' }
#' A grand baseline using the first row of the first list entry (i.e. the mean
#' AF spectrum of the first candidate) is always prepended so every plot has a
#' common anchor.
#'
#' Results are returned as a list and, optionally, as a summary bar chart saved
#' to \code{plot.dir}.
#'
#' @param unstained.fcs Character scalar. Path to the unstained FCS file.
#' @param spectra Numeric matrix of fluorophore spectra (fluorophores x
#'   detectors). Row names must be fluorophore names; column names must match
#'   the detector channels in the FCS file. Any row named \code{"AF"} is
#'   removed automatically.
#' @param af.spectra.list Named list of AF spectra matrices. Each element must
#'   be a numeric matrix with columns matching \code{colnames(spectra)}. The
#'   first row of each matrix is treated as that candidate's mean AF spectrum.
#'   Names are used as labels throughout; if the list is unnamed, entries are
#'   labelled \code{"af1"}, \code{"af2"}, etc.
#' @param assign.fn Character scalar. Name of the assign-type function used to
#'   map each cell to an AF spectrum row. Must be available in the current
#'   search path and follow the assign-type calling convention:
#'   \code{fn(raw.data, spectra, af.spectra)} returning an integer vector of
#'   row indices into \code{af.spectra}. Default:
#'   \code{"assign.af.fluorophores"}.
#' @param n.downsample Integer scalar. Maximum number of events used from the
#'   FCS file. A random subsample is drawn when the file contains more events.
#'   Set to \code{Inf} to use all events. Default: \code{1000L}.
#' @param plot.dir Character scalar. Directory in which to save the summary
#'   plot PDF. Created recursively if it does not exist. Set to \code{NULL} to
#'   skip saving. Default: \code{"figure_af_accuracy"}.
#' @param title Character scalar. Stem used to name the output PDF
#'   (\code{<plot.dir>/<title>.pdf}). Default: \code{"compare_af"}.
#'
#' @return A named list with one entry per candidate (plus \code{"baseline"}).
#'   Each entry contains:
#'   \describe{
#'     \item{\code{Assignments}}{Integer vector of per-cell AF-spectrum row
#'       indices.}
#'     \item{\code{Similarity}}{Numeric vector of per-cell cosine similarities
#'       between the raw detector signal and the assigned AF spectrum.}
#'     \item{\code{Mean_Sim}}{Mean of \code{Similarity} (NAs excluded).}
#'     \item{\code{rSD_Sim}}{Standard deviation of \code{Similarity} (NAs
#'       excluded).}
#'     \item{\code{n.variants}}{Number of AF spectrum rows (variants) in this
#'       candidate's matrix.}
#'   }
#'
#' @seealso \code{\link{get.af.spectra}}, \code{\link{test.af.accuracy}},
#'   \code{\link{assign.af.fluorophores}}
#'
#' @importFrom ggplot2 ggplot aes geom_col geom_errorbar geom_point
#' @importFrom ggplot2 scale_fill_manual scale_y_continuous theme_bw theme element_blank
#' @importFrom ggplot2 element_text element_line labs coord_cartesian
#'
#' @export
compare.af <- function(
  unstained.fcs,
  spectra,
  af.spectra.list,
  assign.fn    = "assign.af.fluorophores",
  n.downsample = 1000L,
  plot.dir     = "figure_af_accuracy",
  title        = "compare_af"
) {

  # ---- input validation ------------------------------------------------------

  stopifnot(
    is.character( unstained.fcs ), length( unstained.fcs ) == 1L,
    file.exists( unstained.fcs ),
    is.matrix( spectra ), is.numeric( spectra ),
    is.list( af.spectra.list ), length( af.spectra.list ) >= 1L,
    all( vapply( af.spectra.list, function( m )
      is.matrix( m ) && is.numeric( m ) && ncol( m ) == ncol( spectra ),
      logical( 1L ) ) ),
    is.character( assign.fn ), length( assign.fn ) == 1L,
    exists( assign.fn, mode = "function" ),
    is.numeric( n.downsample ) || is.integer( n.downsample ),
    length( n.downsample ) == 1L, n.downsample >= 1L
  )

  # ---- name the list if unnamed ----------------------------------------------

  if ( is.null( names( af.spectra.list ) ) )
    names( af.spectra.list ) <- paste0( "af", seq_along( af.spectra.list ) )

  # replace any empty names
  empty <- names( af.spectra.list ) == ""
  if ( any( empty ) )
    names( af.spectra.list )[ empty ] <-
      paste0( "af", seq_along( af.spectra.list ) )[ empty ]

  candidate.names <- names( af.spectra.list )

  # ---- remove AF row from spectra if present ---------------------------------

  if ( "AF" %in% rownames( spectra ) )
    spectra <- spectra[ rownames( spectra ) != "AF", , drop = FALSE ]

  fluorophore.n <- nrow( spectra )

  # ---- read FCS and subset to spectral columns -------------------------------

  raw.data <- readFCS( unstained.fcs )[ , colnames( spectra ) ]

  # ---- downsample ------------------------------------------------------------

  if ( is.finite( n.downsample ) && nrow( raw.data ) > n.downsample )
    raw.data <- raw.data[ sample( nrow( raw.data ), n.downsample ), ]

  n.cells  <- nrow( raw.data )
  fn       <- get( assign.fn )

  # ---- grand baseline: mean AF from the first candidate ----------------------
  # Every cell is assigned to row 1 of af.spectra.list[[1]] regardless of
  # which variant it most resembles. This gives a fixed anchor for the plot.

  baseline.vec <- af.spectra.list[[ 1L ]][ 1L, ]

  baseline.sim <- vapply( seq_len( n.cells ), function( cell ) {
    cell.vec <- raw.data[ cell, ]
    denom    <- sqrt( sum( cell.vec^2 ) ) * sqrt( sum( baseline.vec^2 ) )
    if ( denom == 0 ) return( NA_real_ )
    sum( cell.vec * baseline.vec ) / denom
  }, numeric( 1L ) )

  baseline.result <- list(
    Assignments = rep( 1L, n.cells ),
    Similarity  = baseline.sim,
    Mean_Sim    = mean( baseline.sim, na.rm = TRUE ),
    rSD_Sim     = stats::mad(   baseline.sim, na.rm = TRUE ),
    n.variants  = 1L
  )

  # ---- per-candidate evaluation ----------------------------------------------

  candidate.results <- lapply( candidate.names, function( cand ) {

    af.spectra <- af.spectra.list[[ cand ]]
    af.n       <- nrow( af.spectra )

    message( sprintf( "Evaluating '%s'  (%d AF variant%s)",
                      cand, af.n, if ( af.n == 1L ) "" else "s" ) )

    # assign each cell to its best AF spectrum row
    af.idx <- fn(
      raw.data   = raw.data,
      spectra    = spectra,
      af.spectra = af.spectra
    )

    # cosine similarity: raw signal vs assigned AF spectrum
    sim <- vapply( seq_len( n.cells ), function( cell ) {
      cell.vec <- raw.data[ cell, ]
      af.vec   <- af.spectra[ af.idx[ cell ], ]
      denom    <- sqrt( sum( cell.vec^2 ) ) * sqrt( sum( af.vec^2 ) )
      if ( denom == 0 ) return( NA_real_ )
      sum( cell.vec * af.vec ) / denom
    }, numeric( 1L ) )

    list(
      Assignments = af.idx,
      Similarity  = sim,
      Mean_Sim    = mean( sim, na.rm = TRUE ),
      rSD_Sim     = stats::mad(   sim, na.rm = TRUE ),
      n.variants  = af.n
    )
  } )

  names( candidate.results ) <- candidate.names

  # prepend the grand baseline
  all.results <- c( list( baseline = baseline.result ), candidate.results )

  # ---- build summary data frame for plotting ---------------------------------

  all.labels <- names( all.results )

  summary.df <- data.frame(
    label      = factor( all.labels, levels = all.labels ),
    Mean_Sim   = vapply( all.results, `[[`, numeric( 1L ), "Mean_Sim" ),
    rSD_Sim    = vapply( all.results, `[[`, numeric( 1L ), "rSD_Sim"   ),
    n.variants = vapply( all.results, `[[`, integer( 1L ), "n.variants" ),
    is.baseline = c( TRUE, rep( FALSE, length( candidate.names ) ) ),
    stringsAsFactors = FALSE
  )

  # ---- plot ------------------------------------------------------------------

  if ( !is.null( plot.dir ) ) {

    if ( !dir.exists( plot.dir ) ) dir.create( plot.dir, recursive = TRUE )

    n.bars    <- nrow( summary.df )
    bar.cols  <- grDevices::hcl.colors( n.bars, palette = "Dark 3" )
    bar.cols[ summary.df$is.baseline ] <- "#AAAAAA"   # grey for baseline
    names( bar.cols ) <- summary.df$label

    # secondary label: append variant count to x-axis tick
    x.labels <- ifelse(
      summary.df$is.baseline,
      as.character( summary.df$label ),
      sprintf( "%s\n(%d var.)", summary.df$label, summary.df$n.variants )
    )
    names( x.labels ) <- summary.df$label

    # y range: clip to [floor - 0.02, 1] so small differences are visible
    y.lo <- max( 0, min( summary.df$Mean_Sim - summary.df$rSD_Sim,
                         na.rm = TRUE ) - 0.02 )

    compare.plot <- ggplot2::ggplot(
      summary.df,
      ggplot2::aes( x = label, y = Mean_Sim, fill = label )
    ) +
      ggplot2::geom_col( width = 0.65, colour = "white", linewidth = 0.3 ) +
      ggplot2::geom_errorbar(
        ggplot2::aes( ymin = Mean_Sim - rSD_Sim, ymax = Mean_Sim + rSD_Sim ),
        width     = 0.25,
        linewidth = 0.6,
        colour    = "grey30"
      ) +
      ggplot2::geom_point(
        size   = 2,
        shape  = 21,
        fill   = "white",
        colour = "grey30"
      ) +
      ggplot2::scale_fill_manual( values = bar.cols, guide = "none" ) +
      ggplot2::scale_x_discrete( labels = x.labels ) +
      ggplot2::scale_y_continuous(
        name   = "Mean cosine similarity (raw vs assigned AF spectrum)",
        expand = ggplot2::expansion( mult = c( 0, 0.03 ) )
      ) +
      ggplot2::coord_cartesian( ylim = c( y.lo, 1 ) ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.title.x     = ggplot2::element_blank(),
        axis.text.x      = ggplot2::element_text( size = 9 ),
        axis.text.y      = ggplot2::element_text( size = 9 ),
        axis.title.y     = ggplot2::element_text( size = 10 ),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor   = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_line( colour = "grey90" )
      ) +
      ggplot2::labs(
        caption = sprintf(
          "Assignment function: %s  |  n = %d cells  |  error bars = +/- 1 rSD",
          assign.fn, n.cells
        )
      )

    pdf.path <- file.path( plot.dir, paste0( title, ".pdf" ) )
    grDevices::pdf( pdf.path, width = 2 + 1.2 * n.bars, height = 5 )
    print( compare.plot )
    grDevices::dev.off()
    message( "Saved comparison plot to: ", pdf.path )

    all.results[[ ".plot" ]] <- compare.plot
  }

  return( all.results )
}
