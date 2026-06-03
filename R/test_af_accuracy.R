#' Test Autofluorescence Assignment Accuracy
#'
#' Benchmarks autofluorescence (AF) assignment and unmixing functions against
#' an unstained FCS file. For each function supplied, cells are assigned to an
#' AF spectrum, unmixed, and evaluated by cosine similarity between the raw
#' detector signal and the assigned AF spectrum. A mean-AF baseline (every cell
#' assigned to the mean AF spectrum) is always prepended to the results. Biplot
#' panels for each method are saved as a single PDF.
#'
#' @param unstained.fcs Character scalar. Path to the unstained FCS file used
#'   as the reference data set.
#' @param spectra Numeric matrix of fluorophore spectra (fluorophores x
#'   detectors). Row names must be fluorophore names; column names must match
#'   the detector channels in the FCS file. Any row named \code{"AF"} is
#'   removed automatically before processing.
#' @param af.spectra Numeric matrix of AF spectra (AF variants x detectors).
#'   The first row is treated as the mean AF spectrum and is used for the
#'   baseline comparison. Column names must match those of \code{spectra}.
#' @param asp Aspect-ratio value passed to \code{\link{create.biplot}}.
#' @param functions Character vector of AF-function names to benchmark. Each
#'   name must resolve to a function in the current search path. Functions
#'   whose names start with \code{"fit."} are called with the signature
#'   \code{fn(raw.data, unmixed, unmixing.matrix, spectra, af.spectra)} and
#'   must return a list with elements \code{$unmixed} (cells x fluorophores,
#'   no AF column) and \code{$af.idx} (integer vector of per-cell AF-spectrum
#'   indices). All other functions are treated as assign-type and called with
#'   \code{fn(raw.data, spectra, af.spectra)}, returning an integer vector of
#'   AF-spectrum indices.
#' @param n.downsample Integer scalar. Maximum number of events read from the
#'   FCS file. A random subsample of this size is drawn when the file contains
#'   more events. Set to \code{Inf} to use all events. Default: \code{1000L}.
#' @param plot.dir Character scalar. Directory in which to save the biplot PDF.
#'   Created recursively if it does not exist. Default: \code{"figure_af_accuracy"}.
#' @param title Character scalar. Stem used to name the output PDF (the file
#'   will be \code{<plot.dir>/<title>_biplots.pdf}). Default:
#'   \code{"af_accuracy"}.
#'
#' @return A named list with one entry per tested method (including the
#'   \code{"mean.af"} baseline). Each entry is itself a list with elements:
#'   \describe{
#'     \item{\code{Assignments}}{Integer vector of per-cell AF-spectrum indices
#'       (all \code{1L} for \code{"mean.af"}).}
#'     \item{\code{Unmixed}}{Numeric matrix of unmixed fluorophore values
#'       (cells x fluorophores, no AF column).}
#'     \item{\code{Similarity}}{Numeric vector of per-cell cosine similarities
#'       between the raw detector signal and the assigned AF spectrum.}
#'     \item{\code{Mean_Sim}}{Mean of \code{Similarity} (NAs excluded).}
#'     \item{\code{SD_Sim}}{Standard deviation of \code{Similarity} (NAs
#'       excluded).}
#'   }
#'
#' @seealso \code{\link{benchmark.af.spectra}}, \code{\link{unmix.ols}},
#'   \code{\link{unmix.ols.fast}}, \code{\link{create.biplot}}
#'
#' @importFrom ggplot2 ggsave
#'
#' @export
test.af.accuracy <- function(
  unstained.fcs,
  spectra,
  af.spectra,
  asp,
  functions    = c( "assign.af.fluorophores", "assign.af.residuals",
                    "assign.af.joint.cov" ),
  n.downsample = 1000L,
  plot.dir     = "figure_af_accuracy",
  title        = "af_accuracy"
) {

  # ---- input validation ------------------------------------------------------

  stopifnot(
    is.character( unstained.fcs ), length( unstained.fcs ) == 1L,
    file.exists( unstained.fcs ),
    is.matrix( spectra ), is.numeric( spectra ),
    is.matrix( af.spectra ), is.numeric( af.spectra ),
    ncol( spectra ) == ncol( af.spectra ),
    is.character( functions ), length( functions ) >= 1L,
    is.numeric( n.downsample ) || is.integer( n.downsample ),
    length( n.downsample ) == 1L, n.downsample >= 1L
  )

  # ---- remove AF row from spectra if present ---------------------------------

  if ( "AF" %in% rownames( spectra ) )
    spectra <- spectra[ rownames( spectra ) != "AF", , drop = FALSE ]

  # ---- read FCS and subset to spectral columns -------------------------------

  raw.data <- readFCS( unstained.fcs )[ , colnames( spectra ) ]

  # ---- downsample ------------------------------------------------------------

  if ( is.finite( n.downsample ) && nrow( raw.data ) > n.downsample )
    raw.data <- raw.data[ sample( nrow( raw.data ), n.downsample ), ]

  # ---- baseline unmixings (no AF; mean AF) -----------------------------------

  unmixed.no.af <- unmix.ols( raw.data, spectra )

  combined.mean <- rbind( spectra, af.spectra[ 1, , drop = FALSE ] )
  rownames( combined.mean )[ nrow( combined.mean ) ] <- "AF"
  unmixed.mean.af <- unmix.ols( raw.data, combined.mean )

  # biplot axes: two channels with the highest SD in the no-AF unmix
  channel.sd     <- apply( unmixed.no.af, 2, stats::sd )
  worst.channels <- colnames( unmixed.no.af )[
    order( channel.sd, decreasing = TRUE )[ 1:2 ] ]

  # ---- pre-compute OLS quantities used by fit-type functions -----------------

  XtX             <- tcrossprod( spectra )
  unmixing.matrix <- solve.default( XtX, spectra )   # fluorophores x detectors
  unmixed.no.af.m <- raw.data %*% t( unmixing.matrix )

  # ---- combined spectra template for the assign -> unmix loop ----------------

  af.n          <- nrow( af.spectra )
  fluorophore.n <- nrow( spectra )
  fluors.af     <- c( rownames( spectra ), "AF" )

  combined.spectra <- matrix(
    0,
    nrow     = fluorophore.n + 1L,
    ncol     = ncol( spectra ),
    dimnames = list( fluors.af, colnames( spectra ) )
  )
  combined.spectra[ seq_len( fluorophore.n ), ] <- spectra

  # ---- mean-AF baseline ------------------------------------------------------

  mean.af.vec <- af.spectra[ 1, ]

  mean.af.sim <- vapply( seq_len( nrow( raw.data ) ), function( cell ) {
    cell.vec <- raw.data[ cell, ]
    denom    <- sqrt( sum( cell.vec^2 ) ) * sqrt( sum( mean.af.vec^2 ) )
    if ( denom == 0 ) return( NA_real_ )
    sum( cell.vec * mean.af.vec ) / denom
  }, numeric( 1 ) )

  unmixed.mean.af.fluors <- unmixed.mean.af[
    , colnames( unmixed.mean.af ) != "AF", drop = FALSE ]

  mean.af.result <- list(
    Assignments = rep( 1L, nrow( raw.data ) ),
    Unmixed     = unmixed.mean.af.fluors,
    Similarity  = mean.af.sim,
    Mean_Sim    = mean( mean.af.sim, na.rm = TRUE ),
    SD_Sim      = stats::sd(   mean.af.sim, na.rm = TRUE )
  )

  # ---- main per-function loop ------------------------------------------------

  af.results <- lapply( functions, function( f ) {

    message( paste( "Extracting AF using", f ) )

    is.fit.fn <- grepl( "^fit\\.", f )
    fn        <- get( f )

    if ( is.fit.fn ) {

      # fit-type: returns list( unmixed, AF, af.idx, fitted.af )
      fit.out <- fn(
        raw.data        = raw.data,
        unmixed         = unmixed.no.af.m,
        unmixing.matrix = unmixing.matrix,
        spectra         = spectra,
        af.spectra      = af.spectra
      )
      af.idx  <- fit.out$af.idx
      unmixed <- fit.out$unmixed     # cells x fluorophores (no AF column)

    } else {

      # assign-type: returns integer vector of af.idx
      af.idx <- fn(
        raw.data   = raw.data,
        spectra    = spectra,
        af.spectra = af.spectra
      )

      unmixed <- matrix(
        0,
        nrow     = nrow( raw.data ),
        ncol     = fluorophore.n + 1L,
        dimnames = list( NULL, fluors.af )
      )

      for ( af in seq_len( af.n ) ) {
        combined.spectra[ fluorophore.n + 1L, ] <- af.spectra[ af, ]
        cell.idx <- which( af.idx == af )
        if ( length( cell.idx ) > 0L ) {
          unmixed[ cell.idx, ] <- unmix.ols.fast(
            raw.data[ cell.idx, , drop = FALSE ],
            combined.spectra
          )
        }
      }

      # drop AF column to match fit-type output shape
      unmixed <- unmixed[ , fluors.af[ fluors.af != "AF" ], drop = FALSE ]
    }

    # cosine similarity: raw signal vs assigned AF spectrum
    sim.results <- vapply( seq_len( nrow( raw.data ) ), function( cell ) {
      cell.vec <- raw.data[ cell, ]
      af.vec   <- af.spectra[ af.idx[ cell ], ]
      denom    <- sqrt( sum( cell.vec^2 ) ) * sqrt( sum( af.vec^2 ) )
      if ( denom == 0 ) return( NA_real_ )
      sum( cell.vec * af.vec ) / denom
    }, numeric( 1 ) )

    list(
      Assignments = af.idx,
      Unmixed     = unmixed,
      Similarity  = sim.results,
      Mean_Sim    = mean( sim.results, na.rm = TRUE ),
      SD_Sim      = stats::sd(   sim.results, na.rm = TRUE )
    )
  } )

  names( af.results ) <- functions

  # prepend the mean-AF baseline
  af.results <- c( list( mean.af = mean.af.result ), af.results )

  # ---- plotting --------------------------------------------------------------

  if ( !dir.exists( plot.dir ) ) dir.create( plot.dir, recursive = TRUE )

  plot.list   <- list()
  plot.labels <- character( 0 )

  make.bp <- function( data, label ) {
    tryCatch(
      create.biplot(
        plot.data  = data,
        x.dim      = worst.channels[ 1 ],
        y.dim      = worst.channels[ 2 ],
        asp        = asp,
        title      = label,
        output.dir = plot.dir,
        save       = FALSE
      ),
      error = function( e ) {
        message( "Biplot failed for '", label, "': ", e$message )
        NULL
      }
    )
  }

  p.no.af <- make.bp( unmixed.no.af, "No AF" )
  if ( !is.null( p.no.af ) ) {
    plot.list[[ length( plot.list ) + 1L ]] <- p.no.af
    plot.labels <- c( plot.labels, "No AF" )
  }

  for ( f in names( af.results ) ) {
    p <- make.bp( af.results[[ f ]]$Unmixed, f )
    if ( !is.null( p ) ) {
      plot.list[[ length( plot.list ) + 1L ]] <- p
      plot.labels <- c( plot.labels, f )
    }
  }

  if ( length( plot.list ) > 0L &&
       requireNamespace( "cowplot", quietly = TRUE ) ) {

    labelled.plots <- mapply(
      function( p, lbl ) {
        cowplot::plot_grid(
          cowplot::ggdraw() +
            cowplot::draw_label( lbl, fontface = "bold", size = 10 ),
          p,
          ncol        = 1,
          rel_heights = c( 0.08, 1 )
        )
      },
      plot.list, plot.labels,
      SIMPLIFY = FALSE
    )

    n.cols   <- min( 3L, length( labelled.plots ) )
    grid     <- cowplot::plot_grid( plotlist = labelled.plots, ncol = n.cols )
    pdf.path <- file.path( plot.dir, paste0( title, "_biplots.pdf" ) )

    grDevices::pdf(
      pdf.path,
      width  = 5 * n.cols,
      height = 5 * ceiling( length( labelled.plots ) / n.cols )
    )
    print( grid )
    grDevices::dev.off()
    message( "Saved biplot grid to: ", pdf.path )

  } else if ( length( plot.list ) > 0L ) {

    warning( "cowplot not available; skipping PDF grid. ",
             "Individual plots saved to plot.dir." )

    for ( i in seq_along( plot.list ) ) {
      ggplot2::ggsave(
        filename = file.path(
          plot.dir,
          paste0( gsub( "[^A-Za-z0-9_]", "_", plot.labels[ i ] ), ".jpg" )
        ),
        plot   = plot.list[[ i ]],
        device = ragg::agg_jpeg,
        width  = 5, height = 5
      )
    }
  }

  return( af.results )
}
