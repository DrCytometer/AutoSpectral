# get_af_spectra.r

#' @title Get Autofluorescence Spectra
#'
#' @description
#' Extracts autofluorescence spectra from an unstained sample. Intended for use
#' with `unmix.autospectral`. Uses SOM clustering for rapid identification of
#' cells with similar AF profiles.
#'
#' Optionally deduplicates the resulting spectra by cosine similarity
#' (`deduplicate = TRUE`, default) to remove near-identical profiles that cause
#' spurious over-correction of near-zero events in fully stained samples. When
#' `refine = TRUE`, a second round of targeted modulation is performed on cells
#' that remain far from zero after the first-pass correction; modulated spectra
#' are screened for redundancy against each other and against the base library
#' before being appended.
#'
#' @importFrom parallelly availableCores
#' @importFrom FNN knnx.index
#' @importFrom cowplot plot_grid
#' @importFrom ggplot2 ggsave ggtitle
#' @importFrom ragg agg_jpeg
#'
#' @param unstained.sample Path and file name for an unstained sample FCS file.
#'   The sample type and processing (protocol) method should match the fully
#'   stained samples to which the AF will be applied, ideally.
#' @param asp The AutoSpectral parameter list. Prepare using
#'   `get.autospectral.param`.
#' @param spectra Spectral signatures of fluorophores, normalized between 0 and
#'   1, with fluorophores in rows and detectors in columns.
#' @param som.dim Number of x and y dimensions for the SOM. Default is `10`.
#' @param dist Integer 1:4, distance function (1 manhattan, 2 euclidean,
#'   3 chebyshev, 4 cosine). Default `4`.
#' @param figures Logical, whether to plot the spectral traces and heatmap for
#'   the AF signatures. Default is `TRUE`.
#' @param save Logical, whether to save the CSV file for the AF signatures.
#'   Default is `TRUE`.
#' @param plot.dir Directory (folder) where the plots will be saved. Default is
#'   `NULL`, which inherits from `asp$figure.af.dir`.
#' @param table.dir Directory (folder) where the spectra csv file will be saved.
#'   Default is `NULL`, which inherits from `asp$table.af.dir`.
#' @param title Title for the output spectral plots and csv file. Default is
#'   `"Autofluorescence spectra"`.
#' @param verbose Logical, controls messaging. Default is `TRUE`.
#' @param deduplicate Logical, default `FALSE`. Whether to deduplicate AF spectra
#'   by cosine similarity after the base clustering stage and again after the
#'   refinement stage. Deduplication removes near-identical spectral profiles
#'   that can cause overzealous matching of near-zero events in fully stained
#'   samples, reducing apparent "squishing". Deduplication is slightly less
#'   accurate. Set to `TRUE` to us it.
#' @param duplication.threshold Numeric, default `0.99`. The cosine similarity
#'   threshold used for deduplication. A spectrum is dropped if its cosine
#'   similarity to any already-retained spectrum meets or exceeds this value.
#'   Only used when `deduplicate = TRUE`.
#' @param use.unmixed Logical, default `TRUE`. Whether to include an OLS
#'   unmixed-space projection (`unmix.ols.fast(unstained.exprs, spectra)`)
#'   alongside the raw detector data as SOM clustering input. Set to `FALSE`
#'   to cluster on raw detector space only, which is appropriate when
#'   `spectra` contains several similar or collinear fluorophores (e.g. a
#'   bead-cell comparison panel), where an OLS unmix is itself unstable and
#'   would corrupt the clustering features rather than enrich them.
#'   `use.unmixed = FALSE` also forces `refine = FALSE` and
#'   `plot.unmixed = FALSE`, since both rely on the same per-cell OLS
#'   unmixing residuals and are subject to the same instability.
#' @param af.basis.components Integer, default `NULL`. When supplied, appends
#'   this many components of the panel-oblique autofluorescence structure of
#'   the unstained sample to the SOM training features, alongside the raw
#'   detector data and (when `use.unmixed = TRUE`) the unmixed coefficients.
#'   Computed directly from the panel residual of the unstained events
#'   (trimming the top 1% by residual norm, subsampling to 50000 events for
#'   the decomposition, then projecting every event), reusing the unmixing
#'   matrix already computed for `unmixed.no.af` rather than re-deriving it.
#'   Training-time use carries none of the collinearity risk a per-cell
#'   regression design would, since clustering is not a regression; this
#'   only gives the SOM extra shape-discriminating features to place nodes
#'   with. Forced to `NULL` when `use.unmixed = FALSE`, since it depends on
#'   the same unmixing matrix. Default `NULL` disables it, matching prior
#'   behaviour.
#' @param raw.pca.components Integer, default `NULL`. When supplied, appends
#'   this many raw-event principal components of the unstained sample
#'   (ordinary singular vectors of the raw detector data, not restricted to
#'   the part the panel cannot explain, subsampled to 50000 events for the
#'   decomposition) to the SOM training features. Unlike
#'   `af.basis.components`, these components are free to cross into the
#'   panel's span; that carries no collinearity risk here for the same
#'   reason as above. Not affected by `use.unmixed`, since it does not use
#'   the unmixing matrix at all. Default `NULL` disables it.
#' @param refine Logical, default `FALSE`. Controls whether to perform a second
#'   round of autofluorescence measurement on "problem cells": those with the
#'   highest residual fluorophore signal after the first-pass per-cell AF
#'   extraction, as defined by `problem.quantile`. When `FALSE`, behavior is
#'   identical to versions of AutoSpectral prior to 1.0.0. If you are working
#'   with samples containing complex autofluorescence, e.g. tissues or tumors,
#'   using `refine = TRUE` will improve autofluorescence extraction at the cost
#'   of an increase in unmixing time.
#' @param problem.quantile Numeric, default `0.99`. The quantile for determining
#'   which cells are "problematic" after first-pass per-cell AF extraction. Cells
#'   at or above this quantile with respect to the L2 norm of their unmixed
#'   fluorophore channels (i.e. still furthest from zero) are selected for the
#'   second-round modulation. A value of `0.99` means the top 1% of cells.
#' @param plot.unmixed Logical, default `FALSE`. Whether to unmix the
#'   unstained sample before and after AF extraction and plot the comparison
#'   as a single side-by-side biplot (`unmixed.no.af`, `unmixed`, and, when
#'   `refine = TRUE` and modulation succeeds, `unmixed.second`). This runs a
#'   full per-cell AF unmixing pass purely for diagnostic plotting, so it
#'   defaults off. When `refine = FALSE`, the comparison plot is drawn
#'   immediately from the first-pass unmixing. When `refine = TRUE`, plotting
#'   is deferred until the refinement loop finishes, so the plot always
#'   reflects the final (possibly modulated) AF spectra rather than an
#'   intermediate state.
#' @param plot.unmixed.n Integer, default `30000`. When `plot.unmixed = TRUE`
#'   and `refine = FALSE`, the unstained sample is subsampled to this many
#'   events before the diagnostic unmixing pass, since the full refinement
#'   population isn't otherwise needed. Ignored when `refine = TRUE`, since
#'   the full population is already required for problem-cell identification.
#' @param remove.contaminants Logical, default `TRUE`. A QC check is performed
#'   to exclude any autofluorescence spectrum that is nearly identical to a
#'   fluorophore signature in `spectra`. This guards against low-level
#'   contamination of the unstained sample by single-stained controls.
#' @param contaminant.threshold Numeric, default `0.99`. When
#'   `remove.contaminants = TRUE`, events in the unstained sample whose cosine
#'   similarity to any fluorophore spectrum in `spectra` meets or exceeds this
#'   value are removed **before** SOM construction. This per-event filter is more
#'   sensitive than the post-SOM centroid check because contaminating events are
#'   unlikely to dominate an entire SOM node. Lower values are more aggressive;
#'   the practical range is roughly 0.98–0.999.
#' @param parallel Logical, default `TRUE`, which enables parallel processing
#'   for per-cell AF identification. Used when `refine = TRUE`.
#' @param threads Numeric, defaults to a single thread for sequential processing
#'   (`parallel = FALSE`) or all available cores if `parallel = TRUE`. Used when
#'   `refine = TRUE`.
#' @param return.model Logical. When `TRUE`, attaches an `"af.model"` attribute
#'   to the returned spectra containing per-node covariance, occupancy priors,
#'   abundance priors and scatter statistics, for use by `unmix.af.gls()`. The
#'   return value is still a matrix, so existing callers are unaffected.
#'   Default `FALSE`. Requires refine = FALSE for well-populated per-node
#'   covariances.
#' @param model.rank Integer, maximum rank retained for each node's spectral
#'   covariance. Default `6`.
#' @param model.var.explained Numeric, fraction of within-node variance to
#'   retain. Default `0.95`.
#' @param model.min.events Integer, minimum events for a node to receive a
#'   covariance estimate. Nodes below this get the pooled covariance.
#'   Default `50`.
#' @param model.shrinkage Numeric in `[0, 1]`, shrinkage of each node covariance
#'   toward the pooled covariance. Guards nodes with few events. Default `0.10`.
#' @param heatmap.color.palette Optional character string defining the viridis
#'   color palette for the fluorophore heatmap. Default is `"viridis"`. Options:
#'   `"magma"`, `"inferno"`, `"plasma"`, `"viridis"`, `"cividis"`, `"rocket"`,
#'   `"mako"`, `"turbo"`.
#' @param spectral.trace.color.palette Optional character string defining the
#'   color palette for the AF traces. Default is `NULL` (default R Brewer
#'   colors). Options: same as `heatmap.color.palette`.
#' @param af.fill.color Color for the shaded region indicating the range of
#'   autofluorescence variation in the variant plot. Default is `"red"`.
#' @param af.line.color Color for the median autofluorescence line in the
#'   variant plot. Default is `"black"`.
#'
#' @return A matrix of autofluorescence spectra (spectra in rows, detectors in
#'   columns). Row 1 is the population mean of the base spectra; subsequent rows
#'   are the deduplicated base spectra and, if `refine = TRUE`, modulated spectra
#'   for problem cells.
#'
#' @export
#'
#' @references
#' Van Gassen S et al. (2015). "FlowSOM: Using self-organizing maps for
#' 87(7), 636-645. \doi{10.1002/cyto.a.22625}
#' Wehrens R, Kruisselbrink J (2018). "Flexible Self-Organizing Maps in kohonen
#' 3.0." \emph{Journal of Statistical Software}, \emph{87}(7), 1-18.
#' \doi{10.18637/jss.v087.i07}

get.af.spectra <- function(
    unstained.sample,
    asp,
    spectra,
    som.dim              = 10,
    dist                 = 4L,
    figures              = TRUE,
    save                 = TRUE,
    plot.dir             = NULL,
    table.dir            = NULL,
    title                = "Autofluorescence spectra",
    verbose              = TRUE,
    deduplicate          = FALSE,
    duplication.threshold = 0.99,
    use.unmixed          = TRUE,
    af.basis.components  = NULL,
    raw.pca.components   = NULL,
    refine               = FALSE,
    problem.quantile     = 0.99,
    plot.unmixed         = FALSE,
    plot.unmixed.n       = 30000,
    remove.contaminants  = TRUE,
    contaminant.threshold = 0.99,
    parallel             = TRUE,
    threads              = if ( parallel ) 0 else 1,
    return.model         = FALSE,
    model.rank           = 6L,
    model.var.explained  = 0.95,
    model.min.events     = 50L,
    model.shrinkage      = 0.10,
    heatmap.color.palette         = "viridis",
    spectral.trace.color.palette  = NULL,
    af.fill.color        = "red",
    af.line.color        = "black"
) {

  # ---------------------------------------------------------------------------
  # Setup
  # ---------------------------------------------------------------------------

  if ( is.null( plot.dir ) )  plot.dir  <- asp$figure.af.dir
  if ( !dir.exists( plot.dir ) ) dir.create( plot.dir )
  if ( is.null( table.dir ) ) table.dir <- asp$table.spectra.dir
  if ( !dir.exists( table.dir ) ) dir.create( table.dir )
  if ( is.null( title ) ) title <- asp$af.file.name

  if ( is.null( threads ) ) threads <- asp$worker.process.n
  if ( parallel & threads == 0 ) threads <- parallelly::availableCores()

  if ( !use.unmixed && refine ) {
    warning(
      "`use.unmixed = FALSE` forces `refine = FALSE`: the refinement stage ",
      "identifies problem cells from per-cell OLS unmixing residuals against ",
      "`spectra`, which is exactly the instability `use.unmixed = FALSE` is ",
      "meant to avoid.",
      call. = FALSE
    )
    refine <- FALSE
  }

  if ( !use.unmixed && plot.unmixed ) {
    warning(
      "`use.unmixed = FALSE` forces `plot.unmixed = FALSE`: the before/after ",
      "diagnostic plot requires the same OLS `unmixed.no.af` baseline against ",
      "`spectra`, which is exactly the instability `use.unmixed = FALSE` is ",
      "meant to avoid.",
      call. = FALSE
    )
    plot.unmixed <- FALSE
  }

  if ( !use.unmixed && !is.null( af.basis.components ) ) {
    warning(
      "`use.unmixed = FALSE` forces `af.basis.components = NULL`: it ",
      "requires unmixing against `spectra` to find the panel-oblique ",
      "directions, which is exactly the instability `use.unmixed = FALSE` ",
      "is meant to avoid. `raw.pca.components` does not have this problem, ",
      "since it never unmixes against `spectra`.",
      call. = FALSE
    )
    af.basis.components <- NULL
  }

  # spectral reference matrix must contain exactly one row per fluorophore
  # before any unmixing method runs. Checked before the "AF" row is stripped
  # below, since subsetting a matrix does not reliably preserve the
  # "fluorophore" attribute check.spectra.duplicates() relies on.
  check.spectra.duplicates( spectra )

  # remove any pre-existing AF row from the fluorophore spectra
  if ( "AF" %in% rownames( spectra ) )
    spectra <- spectra[ rownames( spectra ) != "AF", , drop = FALSE ]

  spectral.channels <- colnames( spectra )

  # ---------------------------------------------------------------------------
  # Import and prepare unstained sample
  # ---------------------------------------------------------------------------

  unstained.ff   <- readFCS( unstained.sample, return.keywords = TRUE )
  file.name      <- resolve.file.name( unstained.sample, unstained.ff$keywords[[ "$FIL" ]], verbose = verbose )

  # Shared timestamp for this run's final spectral-profile outputs (CSV,
  # trace, heatmap, variant plots), so repeated calls against the same
  # `plot.dir` / `table.dir` never overwrite a prior run's saved spectra.
  # Processing/QC outputs (contamination report, AF-extraction comparison)
  # intentionally keep static names and are expected to overwrite.
  run.timestamp <- format( Sys.time(), "%Y%m%d_%H%M%S" )

  # retain scatter (and, where configured, imaging) columns alongside the
  # spectral data so that per-node scatter statistics can be computed when
  # return.model = TRUE. Rows are kept aligned with unstained.exprs through
  # every subsetting step below.
  scatter.channels  <- intersect( asp$default.scatter.parameter,
                                  colnames( unstained.ff$data ) )
  unstained.scatter <- if ( length( scatter.channels ) > 0 )
    unstained.ff$data[ , scatter.channels, drop = FALSE ] else NULL

  unstained.exprs <- unstained.ff$data[ , spectral.channels ]

  # ------------------------------------------------------------------
  # Per-event contaminant filter (before SOM)
  # Remove events that are too similar to any fluorophore spectrum.
  # This is more sensitive than post-SOM centroid QC because a small
  # number of contaminating events will not dominate an entire SOM node.
  # ------------------------------------------------------------------
  if ( remove.contaminants ) {
    # subtract mean background
    sample.mean <- colMeans( unstained.exprs )
    unstained.exprs.orth <- sweep( unstained.exprs, 2, sample.mean, "-" )
    # check for fluor contamination on background-subtracted data
    keep.events <- .filter.contaminant.events(
      unstained.exprs.orth,
      spectra,
      contaminant.threshold
    )
    n.removed <- sum( !keep.events )
    if ( n.removed > 0 ) {
      if ( verbose ) {
        message(
          n.removed, " event(s) in the unstained sample were removed prior to",
          " SOM construction due to cosine similarity \u2265 ",
          contaminant.threshold,
          " with a fluorophore spectrum (possible stained-sample contamination)."
        )
      }
      unstained.exprs <- unstained.exprs[ keep.events, , drop = FALSE ]
      if ( !is.null( unstained.scatter ) )
        unstained.scatter <- unstained.scatter[ keep.events, , drop = FALSE ]
    }
  }

  # OLS unmix without AF - combined with raw data for richer clustering
  # features. Skipped entirely when `use.unmixed = FALSE`, since an OLS
  # unmix against a collinear `spectra` (e.g. several similar fluorophores
  # in a bead-cell comparison panel) is itself unstable and would corrupt
  # rather than enrich the clustering features.
  #
  # The unmixing matrix is computed directly here, rather than through
  # unmix.ols.fast(), so af.basis.components below can reuse it instead of
  # re-deriving the same solve a second time.
  if ( use.unmixed || !is.null( af.basis.components ) ) {
    unmixing.matrix <- solve.default( tcrossprod( spectra ), spectra )
  }

  if ( use.unmixed ) {
    unmixed.no.af <- unstained.exprs %*% t( unmixing.matrix )
    cluster.data  <- cbind( unstained.exprs, unmixed.no.af )
  } else {
    cluster.data  <- unstained.exprs
  }

  # Optional extra shape-discriminating features for SOM training. Neither
  # carries the collinearity risk a per-cell regression design would, since
  # clustering is not a regression; both only change which directions the
  # SOM has to differentiate nodes along.
  if ( !is.null( af.basis.components ) ) {

    # Panel residual of every event -- the part the fluorophores cannot
    # explain, which is where the leading out-of-span directions live.
    # Trimmed and subsampled for the decomposition only; the resulting
    # directions are then applied to every event.
    resid      <- unstained.exprs - ( unstained.exprs %*% t( unmixing.matrix ) ) %*% spectra
    resid.norm <- sqrt( rowSums( resid^2 ) )
    keep.basis <- resid.norm <= stats::quantile( resid.norm, 0.99, na.rm = TRUE )

    events.basis <- unstained.exprs[ keep.basis, , drop = FALSE ]
    if ( nrow( events.basis ) > 5e4 ) {
      set.seed( asp$bird.seed )
      events.basis <- events.basis[ sample( nrow( events.basis ), 5e4 ), , drop = FALSE ]
    }
    resid.basis <- events.basis - ( events.basis %*% t( unmixing.matrix ) ) %*% spectra

    n.pc <- min( as.integer( af.basis.components ), nrow( resid.basis ) - 1L,
                 length( spectral.channels ) )
    sv   <- svd( t( resid.basis ), nu = n.pc, nv = 0L )

    emp.scores <- unstained.exprs %*% sv$u[ , seq_len( n.pc ), drop = FALSE ]
    colnames( emp.scores ) <- paste0( "EmpAF", seq_len( n.pc ) )
    cluster.data <- cbind( cluster.data, emp.scores )
  }

  if ( !is.null( raw.pca.components ) ) {

    fit.data <- unstained.exprs
    if ( nrow( fit.data ) > 5e4 ) {
      set.seed( asp$bird.seed )
      fit.data <- fit.data[ sample( nrow( fit.data ), 5e4 ), , drop = FALSE ]
    }

    n.pc <- min( as.integer( raw.pca.components ), nrow( fit.data ) - 1L,
                 length( spectral.channels ) )
    sv   <- svd( fit.data, nu = 0L, nv = n.pc )

    pca.scores <- unstained.exprs %*% sv$v[ , seq_len( n.pc ), drop = FALSE ]
    colnames( pca.scores ) <- paste0( "RawPC", seq_len( n.pc ) )
    cluster.data <- cbind( cluster.data, pca.scores )
  }

  # ---------------------------------------------------------------------------
  # Stage 1: Base AF spectra via SOM
  # ---------------------------------------------------------------------------

  if ( verbose ) message( "Creating a self-organizing map of the autofluorescence" )

  cell.n <- nrow( cluster.data )

  if ( cell.n < 100 ) {
    stop( paste( "Inadequate cell numbers provided:", cell.n ) )
  } else if ( cell.n < 500 ) {
    som.dim <- max( 2, floor( sqrt( cell.n / 3 ) ) )
  }

  map <- get.som.codes(
    data    = cluster.data,
    som.dim = som.dim,
    dist    = dist,
    seed    = asp$bird.seed,
    threads = if ( parallel ) threads else 1L
  )

  # L-infinity normalise SOM node codes
  af.spectra <- t( apply( map$codes[ , spectral.channels ], 1,
                          function( x ) x / max( abs( x ) ) ) )
  af.spectra <- as.matrix( stats::na.omit( af.spectra ) )

  # Prepend population mean
  mean.af    <- colMeans( af.spectra )
  af.spectra <- rbind( mean.af, af.spectra )
  rownames( af.spectra ) <- paste0( "AF", seq_len( nrow( af.spectra ) ) )

  # Contamination QC: remove spectra resembling fluorophores
  af.spectra <- qc.af.spectra( af.spectra, spectra, plot.dir, remove.contaminants, pass = 2,
                               sample.label = file.name )

  # Deduplication of base spectra
  if ( deduplicate ) {
    if ( verbose ) message( "Deduplicating base AF spectra by cosine similarity" )
    n.before   <- nrow( af.spectra )
    af.spectra <- deduplicate.spectra( af.spectra, threshold = duplication.threshold )
    n.after    <- nrow( af.spectra )
    if ( verbose )
      message(
        sprintf(
          "  %d base spectra retained after deduplication (dropped %d)",
          n.after, n.before - n.after
        )
      )
  }

  # Refresh row names after any removals
  rownames( af.spectra ) <- paste0( "AF", seq_len( nrow( af.spectra ) ) )

  # ---------------------------------------------------------------------------
  # Stage 1 figures: base spectra
  # ---------------------------------------------------------------------------

  if ( figures ) {
    if ( verbose ) message( "Plotting autofluorescence spectra" )

    # clamp negatively vectored AF spectra for plotting only (occurs on S8, A8)
    af.spectra.plot <- af.spectra
    af.spectra.plot[ af.spectra.plot < 0 ] <- 0

    tryCatch(
      expr = {
        spectral.trace(
          spectral.matrix      = af.spectra.plot,
          asp                  = asp,
          title                = paste( file.name, title, run.timestamp ),
          plot.dir             = plot.dir,
          split.lasers         = FALSE,
          color.palette        = spectral.trace.color.palette
        )
        spectral.heatmap(
          spectra              = af.spectra.plot,
          title                = paste( file.name, title, run.timestamp ),
          plot.dir             = plot.dir,
          color.palette        = heatmap.color.palette
        )
        spectral.variant.plot.dens(
          spectra.variants   = af.spectra.plot,
          median.spectrum    = mean.af,
          title              = paste0( file.name, " ", title, " density ", run.timestamp ),
          save               = TRUE,
          plot.dir           = plot.dir,
          variant.color = af.fill.color,
          median.line.color  = af.line.color
        )
      },
      error = function( e ) {
        message( "Error in plotting AF spectra: ", e$message )
        return( NULL )
      }
    )
  }

  # ---------------------------------------------------------------------------
  # Stage 2: first-pass per-cell AF unmixing
  # ---------------------------------------------------------------------------
  # Runs whenever the refinement loop needs problem cells (`refine = TRUE`) or
  # a before/after diagnostic plot has been requested (`plot.unmixed = TRUE`).
  # When only the diagnostic plot is wanted, the unstained sample is
  # subsampled to `plot.unmixed.n` events first, since the full refinement
  # population is not otherwise needed and per-cell unmixing is the slow part
  # of this function.

  unmixed.second <- NULL

  if ( refine || plot.unmixed ) {

    if ( refine && return.model ) {
      warning( "`return.model = TRUE` with `refine = TRUE`: refined AF spectra ",
               "are synthesised rather than drawn from SOM node populations, ",
               "so several dictionary entries may attract few events and will ",
               "fall back to the pooled covariance. Use `refine = FALSE` when ",
               "building an AF model.", call. = FALSE )
    }

    if ( verbose ) message( "Identifying best-fitting AF - first pass" )

    # Full population when refining (problem-cell quantiles need the real
    # population); a speed subsample when only plotting.
    if ( refine ) {
      plot.idx <- seq_len( nrow( unstained.exprs ) )
    } else if ( nrow( unstained.exprs ) > plot.unmixed.n ) {
      set.seed( asp$bird.seed )
      plot.idx <- sample( nrow( unstained.exprs ), plot.unmixed.n )
    } else {
      plot.idx <- seq_len( nrow( unstained.exprs ) )
    }

    unstained.exprs.pass <- unstained.exprs[ plot.idx, , drop = FALSE ]
    unmixed.no.af.pass   <- unmixed.no.af[ plot.idx, , drop = FALSE ]

    # Per-cell AF assignment on the unstained sample using the base spectra
    if ( requireNamespace( "AutoSpectralRcpp", quietly = TRUE ) &&
         "assign.af.fluor.fast" %in% ls( getNamespace( "AutoSpectralRcpp" ) ) ) {
      af.assignments <- AutoSpectralRcpp::assign.af.fluor.fast(
        raw.data  = unstained.exprs.pass,
        spectra   = spectra,
        af.spectra = af.spectra,
        threads   = asp$worker.process.n
      )
    } else {
      af.assignments <- assign.af.fluorophores(
        raw.data   = unstained.exprs.pass,
        spectra    = spectra,
        af.spectra = af.spectra
      )
    }

    # Unmix each cell with its assigned AF spectrum, tracking residuals and
    # projected fluorophore signal so we can compute a detector-space error
    fluor.idx    <- 2:( nrow( spectra ) + 1 )
    af.abundance <- rep( 0, nrow( unstained.exprs.pass ) )
    unmixed      <- cbind( af.abundance, unmixed.no.af.pass )

    af.fit <- unmix.af.fwl(
      raw.data         = unstained.exprs.pass,
      spectra          = spectra,
      af.spectra       = af.spectra,
      af.index         = af.assignments,
      unmixed.no.af    = unmixed.no.af.pass,
      return.fitted.af = TRUE
    )

    unmixed[ , 1 ]          <- af.fit$af
    unmixed[ , fluor.idx ]  <- af.fit$fluorophores

    proj.fluor <- af.fit$fluorophores %*% spectra
    residuals  <- unstained.exprs.pass - proj.fluor - af.fit$fitted.af

    # detector-space error = fluorophore projection + raw residuals
    error <- residuals + proj.fluor

    # -------------------------------------------------------------------
    # Stage 2 Refine: targeted modulation for problem cells
    # -------------------------------------------------------------------

    if ( refine ) {

      # ---- Identify problem cells (those still furthest from zero) ------

      if ( verbose ) message( "Refine: calculating error magnitude for problem cell selection" )

      if ( length( fluor.idx ) > 1 ) {
        error.magnitude <- sqrt( rowSums( unmixed[ , fluor.idx ]^2 ) )
      } else {
        error.magnitude <- abs( unmixed[ , fluor.idx ] )
      }

      # Step the quantile down in 5 % increments until we have >= 500 cells
      while ( TRUE ) {
        threshold      <- stats::quantile( error.magnitude, problem.quantile )
        problem.idx    <- which( error.magnitude > threshold )
        problem.cell.n <- length( problem.idx )

        if ( problem.cell.n >= 500 ) break

        problem.quantile <- problem.quantile - 0.05

        if ( problem.quantile < 0.5 ) {
          threshold      <- stats::quantile( error.magnitude, problem.quantile )
          problem.idx    <- which( error.magnitude > threshold )
          problem.cell.n <- length( problem.idx )
          break
        }
      }

      if ( verbose )
        message(
          sprintf(
            "Refine: %d problem cells selected (quantile = %.2f, threshold = %.2f)",
            problem.cell.n, problem.quantile, threshold
          )
        )

      # ---- Modulate base spectra using error clusters --------------------

      if ( problem.cell.n > 10 ) {

        # AF abundance for the problem cells (normalisation denominator)
        af.abundance.problem <- unmixed[ problem.idx, 1 ]
        af.abundance.problem[ af.abundance.problem == 0 ] <- 1e-6

        # Spill ratios: per-channel error normalised by AF abundance,
        # giving a dimensionless signature of how the current AF estimate is wrong
        spill.ratios <- error[ problem.idx, ] / af.abundance.problem

        if ( verbose )
          message(
            paste( "Refine: clustering", problem.cell.n, "problem cells by spillover error pattern" )
          )

        som.dim.error <- max( 2, floor( sqrt( problem.cell.n / 3 ) ) )

        colnames( spill.ratios ) <- colnames( spectra )
        map.error <- get.som.codes(
          data    = spill.ratios,
          som.dim = som.dim.error,
          dist    = dist,
          seed    = asp$bird.seed,
          threads = if ( parallel ) threads else 1L
        )

        error.assign <- as.integer(
          FNN::knnx.index( data = map.error$codes, query = spill.ratios, k = 1 )
        )

        cluster.ids <- unique( error.assign )

        modulated.list <- lapply( cluster.ids, function( cl ) {
          cl.sub.idx <- which( error.assign == cl )
          global.idx <- problem.idx[ cl.sub.idx ]

          # median correction pattern for this error cluster
          median.ratio <- apply(
            spill.ratios[ cl.sub.idx, , drop = FALSE ],
            2,
            stats::median
          )

          # which base AF spectra were assigned to cells in this cluster?
          contributing.af.ids <- unique( af.assignments[ global.idx ] )

          # modulate each contributing base spectrum
          new.specs <- lapply( contributing.af.ids, function( id ) {
            base.spec <- af.spectra[ id, ]
            updated   <- base.spec * ( 1 + median.ratio )
            peak      <- max( abs( updated ) )
            if ( peak > 1e-12 ) updated <- updated / peak
            return( updated )
          } )

          return( do.call( rbind, new.specs ) )
        } )

        modulated.af.spectra <- do.call( rbind, modulated.list )
        modulated.af.spectra <- as.matrix( stats::na.omit( modulated.af.spectra ) )

        if ( nrow( modulated.af.spectra ) > 0 && deduplicate ) {

          # Step 1: deduplicate modulated spectra against each other
          n.mod.before         <- nrow( modulated.af.spectra )
          modulated.af.spectra <- deduplicate.spectra(
            modulated.af.spectra,
            threshold = duplication.threshold
          )

          # Step 2: drop any modulated spectrum too similar to an already-kept
          # base spectrum (cross-deduplication)
          cross.sim  <- cosine.similarity.cross( modulated.af.spectra, af.spectra )
          # cross.sim is (n_modulated x n_existing); keep rows where max sim < threshold
          novel.mask <- apply( cross.sim, 1, max ) < duplication.threshold
          modulated.af.spectra <- modulated.af.spectra[ novel.mask, , drop = FALSE ]

          n.novel <- nrow( modulated.af.spectra )
          if ( verbose )
            message(
              sprintf(
                "Refine: %d novel modulated spectra retained after deduplication (dropped %d)",
                n.novel, n.mod.before - n.novel
              )
            )
        }

        if ( nrow( modulated.af.spectra ) > 0 ) {

          af.spectra <- rbind( af.spectra, modulated.af.spectra )
          af.spectra <- as.matrix( stats::na.omit( af.spectra ) )

          # Contamination QC on the expanded set
          af.spectra <- qc.af.spectra( af.spectra, spectra, plot.dir, remove.contaminants,
                                       sample.label = file.name )

          rownames( af.spectra ) <- paste0( "AF", seq_len( nrow( af.spectra ) ) )

          if ( verbose )
            message(
              sprintf(
                "Refine: %d total AF spectra after modulation and QC",
                nrow( af.spectra )
              )
            )

        } else {
          if ( verbose )
            message( "Refine: all modulated spectra were redundant with base spectra - nothing appended." )
        }

        # ---- Second-pass unmixing for the diagnostic plot -----------------

        if ( plot.unmixed ) {
          if ( verbose ) message( "Refine: identifying best-fitting AF - second pass" )

          if ( requireNamespace( "AutoSpectralRcpp", quietly = TRUE ) &&
               "unmix.autospectral.rcpp" %in% ls( getNamespace( "AutoSpectralRcpp" ) ) ) {
            unmixed.second <- AutoSpectralRcpp::unmix.autospectral.rcpp(
              raw.data   = unstained.exprs.pass,
              spectra    = spectra,
              af.spectra = af.spectra,
              verbose    = FALSE,
              parallel   = TRUE,
              threads    = threads
            )
          } else {
            af.assignments.second <- assign.af.fluorophores(
              raw.data   = unstained.exprs.pass,
              spectra    = spectra,
              af.spectra = af.spectra
            )

            af.fit.second <- unmix.af.fwl(
              raw.data   = unstained.exprs.pass,
              spectra    = spectra,
              af.spectra = af.spectra,
              af.index   = af.assignments.second
            )

            unmixed.second <- unmixed
            unmixed.second[ , 1 ]         <- af.fit.second$af
            unmixed.second[ , fluor.idx ] <- af.fit.second$fluorophores
          }
        }

      } else {
        message( "Refine: insufficient problem cells found - skipping modulation." )
      }

    }   # end refine

    # -------------------------------------------------------------------
    # Before/after diagnostic plot
    # -------------------------------------------------------------------
    # Fires once, after the refinement loop (if any) has finished, so the
    # plot always reflects the final AF spectra rather than an intermediate
    # state. Falls back to two panels when refine found too few problem
    # cells to produce a second pass.

    if ( plot.unmixed && figures && ncol( unmixed.no.af.pass ) > 1 ) {

      if ( verbose ) message( "Plotting impact of AF extraction" )

      channel.sd     <- apply( unmixed.no.af.pass, 2, stats::sd )
      worst.channels <- colnames( unmixed.no.af.pass )[ order( channel.sd, decreasing = TRUE )[ 1:2 ] ]

      tryCatch(
        expr = {
          panel.list <- list(
            create.biplot(
              unmixed.no.af.pass,
              x.dim      = worst.channels[ 1 ],
              y.dim      = worst.channels[ 2 ],
              asp        = asp,
              title      = paste0( file.name, "_", title, "_No_AF_Extraction" ),
              output.dir = plot.dir,
              save       = FALSE
            ) + ggplot2::ggtitle( "No AF extraction" ),
            create.biplot(
              unmixed,
              x.dim      = worst.channels[ 1 ],
              y.dim      = worst.channels[ 2 ],
              asp        = asp,
              title      = paste0( file.name, "_", title, "_PerCell_AF_Extraction_First_Pass" ),
              output.dir = plot.dir,
              save       = FALSE
            ) + ggplot2::ggtitle( "First pass" )
          )

          if ( !is.null( unmixed.second ) ) {
            panel.list <- c( panel.list, list(
              create.biplot(
                unmixed.second,
                x.dim      = worst.channels[ 1 ],
                y.dim      = worst.channels[ 2 ],
                asp        = asp,
                title      = paste0( file.name, "_", title, "_PerCell_AF_Extraction_Second_Pass" ),
                output.dir = plot.dir,
                save       = FALSE
              ) + ggplot2::ggtitle( "Second pass" )
            ) )
          }

          combined.plot <- cowplot::plot_grid( plotlist = panel.list, nrow = 1 )

          ggplot2::ggsave(
            filename  = file.path(
              plot.dir,
              paste0( file.name, "_", title, "_AF_Extraction_Comparison.jpg" )
            ),
            plot      = combined.plot,
            device    = ragg::agg_jpeg,
            width     = 5 * length( panel.list ),
            height    = 5,
            limitsize = FALSE
          )
        },
        error = function( e ) {
          message( "Error in plotting AF extraction: ", e$message )
          return( NULL )
        }
      )
    }

  }   # end first pass / refine / plot.unmixed

  # ---------------------------------------------------------------------------
  # Save and final figures
  # ---------------------------------------------------------------------------

  if ( save ) {
    af.file.name <- paste0( file.name, " ", title, " ", run.timestamp, ".csv" )
    utils::write.csv( af.spectra, file = file.path( table.dir, af.file.name ) )
  }

  if ( figures ) {
    if ( verbose ) message( "Plotting autofluorescence variation" )

    tryCatch(
      expr = {
        spectral.variant.plot(
          af.spectra.plot,
          mean.af,
          title               = paste( file.name, "Autofluorescence variation", run.timestamp ),
          save                = TRUE,
          plot.dir            = plot.dir,
          variant.fill.color  = af.fill.color,
          median.line.color   = af.line.color
        )
      },
      error = function( e ) e
    )
  }

  if ( return.model ) {

    if ( verbose ) message( "Computing per-node AF covariance model" )

    attr( af.spectra, "af.model" ) <- .compute.af.node.model(
      event.spectral = cluster.data[ , spectral.channels, drop = FALSE ],
      event.scatter  = unstained.scatter,
      af.spectra     = af.spectra,
      rank           = model.rank,
      var.explained  = model.var.explained,
      min.events     = model.min.events,
      shrinkage      = model.shrinkage
    )

    if ( verbose ) {
      mdl <- attr( af.spectra, "af.model" )
      shape.frac <- vapply( mdl$nodes, function( z ) {
        v <- z$shape.frac
        if ( is.null( v ) || length( v ) != 1L ) NA_real_ else as.numeric( v )
      }, numeric( 1 ) )
      message( sprintf(
        "  %d nodes, median rank %d, median occupancy %d events",
        length( mdl$nodes ),
        stats::median( vapply( mdl$nodes, function( z ) z$rank, numeric( 1 ) ) ),
        stats::median( vapply( mdl$nodes, function( z ) z$n,    numeric( 1 ) ) ) ) )
      message( sprintf(
        "  shape.frac: median %.2f, %d / %d nodes scored (rest below min.events)",
        stats::median( shape.frac, na.rm = TRUE ),
        sum( !is.na( shape.frac ) ), length( shape.frac ) ) )
    }
  }

  return( af.spectra )

}


# ---------------------------------------------------------------------------
# Helper: greedy cosine-similarity deduplication
# ---------------------------------------------------------------------------

#' @title Deduplicate Spectra by Cosine Similarity
#'
#' @description
#' Iterates through rows of a spectral matrix in order, retaining a row only
#' if its cosine similarity to every already-retained row is strictly below
#' `threshold`. Intended for internal use by `get.af.spectra`.
#'
#' @param spectra Numeric matrix, spectra in rows and detectors in columns.
#'   Assumed to be L-infinity normalised.
#' @param threshold Numeric scalar in (0, 1]. Rows at or above this similarity
#'   to any retained row are dropped. Default `0.99`.
#'
#' @return Numeric matrix with redundant rows removed.
#'
#' @keywords internal

deduplicate.spectra <- function( spectra, threshold = 0.99 ) {

  if ( nrow( spectra ) <= 1 ) return( spectra )

  # Row-normalise to unit length for cosine similarity via dot product
  norms  <- sqrt( rowSums( spectra^2 ) )
  norms  <- ifelse( norms < 1e-12, 1, norms )
  s.norm <- spectra / norms

  kept <- integer( 0 )

  for ( i in seq_len( nrow( spectra ) ) ) {
    if ( length( kept ) == 0 ) {
      kept <- c( kept, i )
    } else {
      # cosine similarities of row i against all kept rows
      sims <- s.norm[ kept, , drop = FALSE ] %*% s.norm[ i, ]
      if ( all( sims < threshold ) )
        kept <- c( kept, i )
    }
  }

  return( spectra[ kept, , drop = FALSE ] )

}

# -----------------------------------------------------------------------------
# Per-node AF model: within-node spectral covariance, occupancy and abundance
# priors, and scatter statistics.
#
# Assignment is recomputed against the FINAL dictionary rather than tracked
# through the SOM, because rows are added (population mean), removed (QC,
# deduplication) and modified (refinement) after mapping. Recomputing also
# matches how unmix.af.gls() will assign at run time.
# -----------------------------------------------------------------------------

.compute.af.node.model <- function(
    event.spectral,
    event.scatter,
    af.spectra,
    rank          = 6L,
    var.explained = 0.95,
    min.events    = 50L,
    shrinkage     = 0.10
) {

  event.spectral <- as.matrix( event.spectral )
  node.n <- nrow( af.spectra )
  det.n  <- ncol( af.spectra )

  # scalar least-squares fit of each event to each dictionary entry:
  # alpha = <y, a> / <a, a>, residual sum of squares = |y|^2 - alpha^2 |a|^2
  aa    <- rowSums( af.spectra^2 )
  ya    <- event.spectral %*% t( af.spectra )              # events x nodes
  alpha <- sweep( ya, 2L, aa, "/" )
  rss   <- rowSums( event.spectral^2 ) - sweep( alpha^2, 2L, aa, "*" )

  assign.k <- max.col( -rss, ties.method = "first" )
  alpha.k  <- alpha[ cbind( seq_len( nrow( alpha ) ), assign.k ) ]

  # shape residual: event rescaled to unit abundance, minus its dictionary entry
  usable <- alpha.k > 0 & is.finite( alpha.k )

  delta.all <- matrix( NA_real_, nrow( event.spectral ), det.n )
  delta.all[ usable, ] <-
    event.spectral[ usable, , drop = FALSE ] / alpha.k[ usable ] -
    af.spectra[ assign.k[ usable ], , drop = FALSE ]

  pooled.cov <- stats::cov( delta.all[ usable, , drop = FALSE ] )

  # Pooled shape fraction: the same intercept regression as per node, over
  # all usable events. Deltas mix true shape variation with measurement
  # noise scaled by 1 / alpha^2; the GLS covariance adds detector noise on
  # its own diagonal, so any noise left inside the AF covariance would be
  # counted twice. Deflating by the estimated shape fraction removes the
  # first-order double count. Also serves as the fallback for nodes too
  # small to estimate their own fraction.
  pooled.shape.frac <- 1
  d2.all  <- rowSums( delta.all[ usable, , drop = FALSE ]^2 )
  ia2.all <- 1 / alpha.k[ usable ]^2
  if ( stats::var( ia2.all ) > 0 && mean( d2.all ) > 0 ) {
    cf.all <- stats::coef( stats::lm( d2.all ~ ia2.all ) )
    pooled.shape.frac <- min( max( cf.all[ 1L ], 0 ) / mean( d2.all ), 1 )
  }
  pooled.cov <- pooled.shape.frac * pooled.cov

  nodes <- lapply( seq_len( node.n ), function( k ) {

    idx <- which( assign.k == k & usable )
    n.k <- length( idx )

    # Decompose the within-node spread into shape variation and measurement
    # noise. Noise enters delta as (photon noise / alpha), so it scales as
    # 1/alpha^2 while genuine shape variation does not. The intercept of
    # |delta|^2 ~ 1/alpha^2 is the shape component.
    shape.frac <- NA_real_
    if ( n.k >= min.events ) {
      d2  <- rowSums( delta.all[ idx, , drop = FALSE ]^2 )
      ia2 <- 1 / alpha.k[ idx ]^2
      if ( stats::var( ia2 ) > 0 && mean( d2 ) > 0 ) {
        cf <- stats::coef( stats::lm( d2 ~ ia2 ) )
        shape.frac <- min( max( cf[ 1L ], 0 ) / mean( d2 ), 1 )
      }
    }

    node.shape.frac <- if ( is.na( shape.frac ) ) pooled.shape.frac else shape.frac

    cov.k <- if ( n.k >= min.events ) {
      ( 1 - shrinkage ) * node.shape.frac *
        stats::cov( delta.all[ idx, , drop = FALSE ] ) +
        shrinkage * pooled.cov
    } else pooled.cov

    ev  <- eigen( cov.k, symmetric = TRUE )
    lam <- pmax( ev$values, 0 )

    keep <- seq_len( min( rank, sum( lam > 0 ) ) )
    if ( length( keep ) > 1L ) {
      cum <- cumsum( lam[ keep ] ) / sum( lam )
      hit <- which( cum >= var.explained )
      if ( length( hit ) > 0L ) keep <- keep[ seq_len( hit[ 1L ] ) ]
    }

    basis <- t( ev$vectors[ , keep, drop = FALSE ] )
    colnames( basis ) <- colnames( af.spectra )

    log.alpha <- if ( n.k >= 3L ) log( alpha.k[ idx ] ) else NA_real_

    list(
      n            = n.k,
      basis        = basis,
      lambda       = lam[ keep ],
      rank         = length( keep ),
      total.var    = sum( lam ),
      shape.frac   = if ( length( shape.frac ) == 1L )
        as.numeric( shape.frac ) else NA_real_,
      log.alpha.mu = if ( n.k >= 3L ) mean( log.alpha ) else NA_real_,
      log.alpha.sd = if ( n.k >= 3L ) stats::sd( log.alpha ) else NA_real_,
      scatter.mean = if ( !is.null( event.scatter ) && n.k >= 3L )
        colMeans( event.scatter[ idx, , drop = FALSE ] ) else NULL,
      scatter.cov  = if ( !is.null( event.scatter ) && n.k >= min.events )
        stats::cov( event.scatter[ idx, , drop = FALSE ] ) else NULL
    )
  } )

  names( nodes ) <- rownames( af.spectra )

  counts <- vapply( nodes, function( z ) z$n, numeric( 1 ) )

  list(
    nodes      = nodes,
    prior      = ( counts + 1 ) / sum( counts + 1 ),   # Laplace-smoothed
    pooled.cov = pooled.cov,
    n.assigned = sum( usable ),
    detectors  = colnames( af.spectra )
  )
}

# ---------------------------------------------------------------------------
# Private helper: per-event contaminant filter
# ---------------------------------------------------------------------------
# Uses the Rcpp version from AutoSpectralRcpp when available; falls back to a
# vectorised pure-R implementation that is fast enough for typical FCS files.

.filter.contaminant.events <- function( event.mat, spectra.mat, threshold ) {

  # Use Rcpp helper if available
  if ( requireNamespace( "AutoSpectralRcpp", quietly = TRUE ) &&
       "filter_contaminant_events_cpp" %in% ls( getNamespace( "AutoSpectralRcpp" ) ) ) {
    return(
      AutoSpectralRcpp::filter_contaminant_events_cpp(
        event_mat   = event.mat,
        spectra_mat = spectra.mat,
        threshold   = threshold
      )
    )
  }

  # Pure-R fallback: vectorised cosine similarity, spectrum by spectrum
  # event norms (length n.events)
  event.norms <- sqrt( rowSums( event.mat^2 ) ) + 1e-9

  # start with all events kept
  keep <- rep( TRUE, nrow( event.mat ) )

  for ( f in seq_len( nrow( spectra.mat ) ) ) {
    spec.vec  <- spectra.mat[ f, ]
    spec.norm <- sqrt( sum( spec.vec^2 ) ) + 1e-9
    # dot products via matrix-vector multiply
    dots <- as.numeric( event.mat %*% spec.vec )
    cs   <- dots / ( event.norms * spec.norm )
    keep <- keep & ( cs < threshold )
    # early exit if everything already flagged
    if ( !any( keep ) ) break
  }

  keep
}

# ---------------------------------------------------------------------------
# Helper: cross cosine similarity matrix
# ---------------------------------------------------------------------------

#' @title Cross Cosine Similarity
#'
#' @description
#' Computes pairwise cosine similarity between rows of two matrices `a` and
#' `b`. Returns an (nrow(a) x nrow(b)) matrix. Intended for internal use by
#' `get.af.spectra`.
#'
#' @param a Numeric matrix.
#' @param b Numeric matrix with the same number of columns as `a`.
#'
#' @return Numeric matrix of cosine similarities, shape (nrow(a), nrow(b)).
#'
#' @keywords internal

cosine.similarity.cross <- function( a, b ) {
  t( apply( b, 1, function( brow ) .cosine.sim.rows( a, brow ) ) )
}
