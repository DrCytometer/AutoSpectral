#' Benchmark AF Assignment Accuracy Against Spectral Panel Size
#'
#' Repeatedly subsamples the fluorophore spectra matrix to a range of panel
#' sizes, runs \code{\link{test.af.accuracy}} for each subsample, and
#' summarises how cosine similarity between raw detector signals and the
#' assigned AF spectra changes with the number of fluorophores. Results are
#' saved as a line-plot PDF and returned as a list.
#'
#' @param unstained.fcs Character scalar. Path to the unstained FCS file
#'   passed through to \code{\link{test.af.accuracy}}.
#' @param spectra Numeric matrix of fluorophore spectra (fluorophores x
#'   detectors). Row names must be fluorophore names; column names must match
#'   the detector channels in the FCS file. Any row named \code{"AF"} is
#'   removed automatically before processing.
#'   Column names must match those of \code{spectra}.
#' @param af.spectra Numeric matrix of AF spectra (AF variants x detectors).
#'   The first row is treated as the mean AF spectrum and is used for the
#'   baseline comparison. Column names must match those of \code{spectra}.
#' @param asp Aspect-ratio value passed through to
#'   \code{\link{test.af.accuracy}}.
#' @param functions Character vector of AF-function names to benchmark. See
#'   \code{\link{test.af.accuracy}} for the required calling conventions.
#' @param n.fluors Integer vector of panel sizes (number of fluorophores) to
#'   evaluate. Values outside \code{[1, nrow(spectra)]} are silently dropped.
#'   Default: \code{c(5, 10, 15, 20, 25, 30, 35, 40)}.
#' @param n.draws Integer scalar. Number of random draws (unique fluorophore
#'   subsets) per panel size. A unique seed is derived for each
#'   \code{(n.fluors, draw)} combination from \code{seed}, so results are
#'   fully reproducible. Default: \code{5L}.
#' @param seed Integer scalar. Base seed for reproducible subsampling.
#'   Default: \code{42L}.
#' @param n.downsample Integer scalar. Maximum number of events passed to each
#'   \code{\link{test.af.accuracy}} call. A random subsample is drawn when the
#'   FCS file contains more events. Set to \code{Inf} to use all events.
#'   Default: \code{1000L}.
#' @param plot.dir Character scalar. Directory for output files. Created
#'   recursively if it does not exist. Default: \code{"figure_af_accuracy"}.
#' @param filename Character scalar. Stem used to name the summary PDF (the
#'   file will be \code{<plot.dir>/<filename>.pdf}).
#'   Default: \code{"af_accuracy_spectra_benchmark"}.
#'
#' @return A named list with three elements:
#'   \describe{
#'     \item{\code{raw}}{Data frame of per-draw results with columns
#'       \code{n.fluors}, \code{draw}, \code{method}, \code{Mean_Sim},
#'       and \code{SD_Sim}.}
#'     \item{\code{summary}}{Data frame summarising \code{raw} across draws
#'       for each \code{(n.fluors, method)} combination, with columns
#'       \code{n.fluors}, \code{method}, \code{mean_Mean_Sim},
#'       \code{sd_Mean_Sim} (variability across draws), and
#'       \code{mean_SD_Sim} (average within-draw spread).}
#'     \item{\code{plot}}{The \code{ggplot2} object for the summary line plot.}
#'   }
#'
#' @seealso \code{\link{test.af.accuracy}}
#'
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_point labs
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous scale_colour_manual
#' @importFrom ggplot2 scale_fill_manual theme_bw theme element_blank element_text
#'
#' @export
benchmark.af.spectra <- function(
  unstained.fcs,
  spectra,
  af.spectra,
  asp,
  functions    = c( "assign.af.fluorophores", "assign.af.residuals",
                    "assign.af.joint.cov" ),
  n.fluors     = c( 5, 10, 15, 20, 25, 30, 35, 40 ),
  n.draws      = 5L,
  seed         = 42L,
  n.downsample = 1000L,
  plot.dir     = "figure_af_accuracy",
  filename     = "af_accuracy_spectra_benchmark"
) {

  # ---- input validation ------------------------------------------------------

  stopifnot(
    is.character( unstained.fcs ), length( unstained.fcs ) == 1L,
    file.exists( unstained.fcs ),
    is.matrix( spectra ), is.numeric( spectra ),
    is.matrix( af.spectra ), is.numeric( af.spectra ),
    ncol( spectra ) == ncol( af.spectra ),
    is.character( functions ), length( functions ) >= 1L,
    is.numeric( n.draws ) || is.integer( n.draws ),
    length( n.draws ) == 1L, n.draws >= 1L,
    is.numeric( n.downsample ) || is.integer( n.downsample ),
    length( n.downsample ) == 1L, n.downsample >= 1L
  )

  # ---- remove AF row from spectra if present ---------------------------------

  if ( "AF" %in% rownames( spectra ) )
    spectra <- spectra[ rownames( spectra ) != "AF", , drop = FALSE ]

  total.fluors <- nrow( spectra )

  # ---- validate and clamp n.fluors -------------------------------------------

  n.fluors <- sort( unique( as.integer( n.fluors ) ) )
  n.fluors <- n.fluors[ n.fluors >= 1L & n.fluors <= total.fluors ]
  if ( length( n.fluors ) == 0L )
    stop(
      "No valid values in n.fluors after filtering against nrow(spectra) = ",
      total.fluors, call. = FALSE
    )

  if ( !dir.exists( plot.dir ) ) dir.create( plot.dir, recursive = TRUE )

  # ---- main loop: each (n, draw) combination ---------------------------------

  all.rows <- vector( "list", length( n.fluors ) * n.draws )
  row.idx  <- 0L

  for ( n in n.fluors ) {
    message( sprintf( "--- Benchmarking n.fluors = %d ---", n ) )

    for ( draw in seq_len( n.draws ) ) {

      # reproducible, unique seed per (n, draw)
      set.seed( seed + n * 100L + draw )
      sub.idx     <- sample( total.fluors, n )
      sub.spectra <- spectra[ sub.idx, , drop = FALSE ]

      # direct biplot output to a temporary subdirectory
      tmp.dir <- file.path( plot.dir, "tmp_benchmark" )
      if ( !dir.exists( tmp.dir ) ) dir.create( tmp.dir, recursive = TRUE )

      result <- tryCatch(
        test.af.accuracy(
          unstained.fcs = unstained.fcs,
          spectra       = sub.spectra,
          af.spectra    = af.spectra,
          asp           = asp,
          functions     = functions,
          n.downsample  = n.downsample,
          plot.dir      = tmp.dir,
          title         = sprintf( "n%02d_draw%02d", n, draw )
        ),
        error = function( e ) {
          message( sprintf( "  Draw %d failed: %s", draw, e$message ) )
          NULL
        }
      )

      if ( is.null( result ) ) next

      for ( fn.name in names( result ) ) {
        row.idx <- row.idx + 1L
        all.rows[[ row.idx ]] <- data.frame(
          n.fluors         = n,
          draw             = draw,
          method           = fn.name,
          Mean_Sim         = result[[ fn.name ]]$Mean_Sim,
          SD_Sim           = result[[ fn.name ]]$SD_Sim,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  # ---- collate results -------------------------------------------------------

  all.rows <- all.rows[ !vapply( all.rows, is.null, logical( 1L ) ) ]
  if ( length( all.rows ) == 0L )
    stop( "All benchmark runs failed; nothing to plot.", call. = FALSE )

  results.df           <- do.call( rbind, all.rows )
  rownames( results.df ) <- NULL

  # ---- summarise across draws ------------------------------------------------

  summary.df <- do.call( rbind, lapply(
    split( results.df,
           list( results.df$n.fluors, results.df$method ),
           drop = TRUE ),
    function( grp ) {
      data.frame(
        n.fluors         = grp$n.fluors[ 1L ],
        method           = grp$method[ 1L ],
        mean_Mean_Sim    = mean( grp$Mean_Sim, na.rm = TRUE ),
        sd_Mean_Sim      = stats::sd(   grp$Mean_Sim, na.rm = TRUE ),
        mean_SD_Sim      = mean( grp$SD_Sim,   na.rm = TRUE ),
        stringsAsFactors = FALSE
      )
    }
  ) )
  rownames( summary.df ) <- NULL

  # ---- plot ------------------------------------------------------------------

  summary.df$method <- factor(
    summary.df$method,
    levels = c( "mean.af", functions )
  )

  n.methods     <- length( levels( summary.df$method ) )
  method.colors <- grDevices::hcl.colors( n.methods, palette = "Dark 3" )
  names( method.colors ) <- levels( summary.df$method )

  ribbon.alpha <- 0.15

  bench.plot <- ggplot2::ggplot(
    summary.df,
    ggplot2::aes(
      x      = n.fluors,
      y      = mean_Mean_Sim,
      colour = method,
      fill   = method,
      group  = method
    )
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = mean_Mean_Sim - sd_Mean_Sim,
        ymax = mean_Mean_Sim + sd_Mean_Sim
      ),
      alpha  = ribbon.alpha,
      colour = NA
    ) +
    ggplot2::geom_line(  linewidth = 0.8 ) +
    ggplot2::geom_point( size      = 2.5 ) +
    ggplot2::scale_x_continuous(
      name   = "Number of fluorophores in spectra",
      breaks = n.fluors
    ) +
    ggplot2::scale_y_continuous(
      name   = "Mean cosine similarity (raw vs AF spectrum)",
      limits = c( NA, 1 )
    ) +
    ggplot2::scale_colour_manual( values = method.colors, name = "Method" ) +
    ggplot2::scale_fill_manual(   values = method.colors, name = "Method" ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position  = "right",
      panel.grid.minor = ggplot2::element_blank(),
      axis.text        = ggplot2::element_text( size = 10 ),
      axis.title       = ggplot2::element_text( size = 11 ),
      legend.text      = ggplot2::element_text( size =  9 )
    ) +
    ggplot2::labs(
      caption = sprintf(
        "Ribbon = +/- SD of Mean_Sim across %d random draws per n; seed = %d",
        n.draws, seed
      )
    )

  pdf.path <- file.path( plot.dir, paste0( filename, ".pdf" ) )
  grDevices::pdf( pdf.path, width = 7, height = 5 )
  print( bench.plot )
  grDevices::dev.off()
  message( "Saved benchmark plot to: ", pdf.path )

  return(
    list(
      raw     = results.df,
      summary = summary.df,
      plot    = bench.plot
    )
  )
}
