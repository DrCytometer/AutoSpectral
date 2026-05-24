# get_fluor_variants.r

#' @title Get Fluorophore Variants
#'
#' @description
#' Assesses variation in the spectral signature of a single-stained flow
#' cytometry control sample using SOM clustering on scatter-matched,
#' background-corrected positive events.
#'
#' Autofluorescence is characterised \strong{in situ} from the paired
#' universal-negative file specified in the control table (or from the lower
#' 25\% of events within the control file itself when no universal negative is
#' available). The AF mean vector is unit-normalised and projected out of each
#' event to identify the empirical peak detector, mirroring the approach in
#' \code{get.spectra.automated}. All positive events above the raw threshold
#' in the (empirical) peak channel are selected, up to \code{n.cells}
#' (randomly downsampled when more are present). For each selected event,
#' \eqn{k} scatter-space nearest neighbours are found in the unstained pool
#' and their spectral values are averaged to form a per-event background
#' estimate, which is then subtracted. SOM clustering on the resulting
#' background-corrected matrix recovers the population-level distribution of
#' spectral shapes. A cosine-similarity QC step retains only SOM centroids
#' sufficiently similar to the reference spectrum, followed by off-peak
#' smoothing.
#'
#' @importFrom FlowSOM SOM
#' @importFrom FNN knnx.index
#'
#' @param fluor Character. Name of the fluorophore.
#' @param file.name Named character vector of control FCS filenames, named by
#'   fluorophore.
#' @param control.dir Character. Directory containing the control FCS files.
#' @param asp The AutoSpectral parameter list from \code{get.autospectral.param()}.
#' @param spectra Numeric matrix. Reference spectra; fluorophores in rows,
#'   detectors in columns.
#' @param figures Logical. Whether to save a spectral-variant plot. Default
#'   \code{TRUE}.
#' @param output.dir Character. Directory for figures.
#' @param verbose Logical. Whether to print progress messages. Default
#'   \code{TRUE}.
#' @param spectral.channel Character vector of spectral detector channel names.
#' @param scatter.channel Character vector of scatter parameter names
#'   (e.g. \code{"FSC-A"}, \code{"SSC-A"}) used for KNN scatter matching
#'   against the unstained pool.
#' @param universal.negative Named character vector mapping fluorophore names
#'   to their paired unstained FCS filename, or \code{"FALSE"} / \code{NA}
#'   when none is available.
#' @param raw.thresholds Named numeric vector of per-channel positivity
#'   thresholds (typically the 99.5th percentile of the unstained sample).
#' @param flow.channel Named character vector of expected peak raw channels,
#'   one per fluorophore.
#' @param n.cells Integer, default \code{10000}. Maximum number of positive
#'   events used for SOM clustering. Files with more events above threshold
#'   are randomly downsampled to this number.
#' @param som.dim Integer, default \code{10}. Side length of the square SOM
#'   grid. Produces up to \code{som.dim^2} candidate variant spectra before
#'   cosine QC.
#' @param k.neighbors Integer, default \code{3}. Number of scatter-space
#'   nearest neighbours from the unstained pool used to form the per-event
#'   background estimate.
#' @param sim.threshold Numeric, default \code{0.99}. Minimum cosine
#'   similarity between a SOM centroid and the reference spectrum for the
#'   centroid to be retained as a variant.
#' @param variant.fill.color Color for the shaded ribbon in the variant plot.
#'   Default \code{"red"}.
#' @param variant.fill.alpha Alpha for \code{variant.fill.color}. Default
#'   \code{0.7}.
#' @param median.line.color Color for the reference-spectrum line. Default
#'   \code{"black"}.
#' @param median.linewidth Width of the reference-spectrum line. Default
#'   \code{1}.
#'
#' @return A numeric matrix; variants in rows, detectors in columns, values
#'   normalised to \eqn{[0, 1]}. When no centroids survive cosine QC the
#'   single reference spectrum is returned (one row).
#'
#' @references
#' Van Gassen S et al. (2015). FlowSOM. \emph{Cytometry Part A}, 87(7),
#' 636-645. \doi{10.1002/cyto.a.22625}

get.fluor.variants.test <- function(
    fluor,
    file.name,
    control.dir,
    asp,
    spectra,
    figures,
    output.dir,
    verbose,
    spectral.channel,
    scatter.channel,
    universal.negative,
    raw.thresholds,
    flow.channel,
    n.cells            = 10000L,
    som.dim            = 10L,
    k.neighbors        = 3L,
    sim.threshold      = 0.99,
    variant.fill.color = "red",
    variant.fill.alpha = 0.7,
    median.line.color  = "black",
    median.linewidth   = 1
) {

  if ( verbose )
    message( paste0( "\033[34m", "Getting spectral variants for ", fluor, "\033[0m" ) )

  # reference spectrum for this fluorophore (1 x D)
  original.spectrum <- spectra[ fluor, , drop = FALSE ]
  orig.vec          <- as.numeric( original.spectrum )

  cols.read <- unique( c( scatter.channel, spectral.channel ) )

  # -------------------------------------------------------------------------
  # 1. Read fluorophore FCS file (scatter + spectral channels)
  # -------------------------------------------------------------------------

  pos.full    <- readFCS( file.path( control.dir, file.name[ fluor ] ) )
  pos.present <- intersect( cols.read, colnames( pos.full ) )
  pos.full    <- pos.full[ , pos.present, drop = FALSE ]

  scat.present <- intersect( scatter.channel, colnames( pos.full ) )
  spec.present <- intersect( spectral.channel, colnames( pos.full ) )

  if ( length( spec.present ) == 0 ) {
    warning( paste0( "No spectral channels in file for ", fluor,
                     ". Returning reference spectrum." ) )
    return( original.spectrum )
  }

  spec.data <- pos.full[ , spec.present, drop = FALSE ]   # all events x D

  # -------------------------------------------------------------------------
  # 2. Derive AF reference in situ from the unstained file
  # -------------------------------------------------------------------------

  is.valid.neg <- !is.na( universal.negative[ fluor ] ) &&
    universal.negative[ fluor ] != "FALSE"              &&
    grepl( "\\.fcs$", universal.negative[ fluor ], ignore.case = TRUE )

  if ( is.valid.neg ) {

    neg.path    <- file.path( control.dir, universal.negative[ fluor ] )
    neg.full    <- readFCS( neg.path )
    neg.present <- intersect( cols.read, colnames( neg.full ) )
    neg.full    <- neg.full[ , neg.present, drop = FALSE ]

    neg.spec    <- intersect( spec.present, colnames( neg.full ) )
    neg.scat    <- intersect( scat.present, colnames( neg.full ) )

    af.spectral  <- neg.full[ , neg.spec, drop = FALSE ]   # M x D
    af.scatter   <- neg.full[ , neg.scat, drop = FALSE ]   # M x S
    af.mean      <- colMeans( af.spectral )
    common.scat  <- neg.scat

    if ( verbose )
      message( paste0( "\033[32m  Universal negative: ",
                       nrow( af.spectral ), " events\033[0m" ) )

  } else {

    # Internal-negative: lower 25% of events by expected (or mean) signal
    exp.peak <- flow.channel[ fluor ]
    ord.vals <- if ( exp.peak %in% spec.present ) spec.data[ , exp.peak ] else
                  rowMeans( spec.data )

    n.neg       <- max( 2L, floor( nrow( spec.data ) * 0.25 ) )
    i.neg       <- order( ord.vals )[ seq_len( n.neg ) ]

    af.spectral <- spec.data[ i.neg, , drop = FALSE ]
    af.scatter  <- pos.full[  i.neg, scat.present, drop = FALSE ]
    af.mean     <- colMeans( af.spectral )
    common.scat <- scat.present

    if ( verbose )
      message( paste0( "\033[33m  Internal negative: lower 25% (",
                       n.neg, " events)\033[0m" ) )
  }

  # unit-normalised AF mean vector for projection
  af.unit <- af.mean / ( sqrt( sum( af.mean^2 ) ) + 1e-9 )

  # -------------------------------------------------------------------------
  # 3. Select positive events via AF projection + peak threshold
  # -------------------------------------------------------------------------

  # Project AF out of all events to reveal the empirical peak detector
  # (mirrors get.spectra.automated step 4b)
  proj.vals <- spec.data %*% af.unit
  mat.orth  <- spec.data - proj.vals %*% t( af.unit )
  emp.peak  <- names( which.max( colMeans( mat.orth ) ) )

  exp.peak   <- flow.channel[ fluor ]
  thresh.col <- if ( emp.peak %in% spec.present )  emp.peak  else
                if ( exp.peak %in% spec.present )  exp.peak  else
                spec.present[ 1L ]

  thresh.val <- if ( thresh.col %in% names( raw.thresholds ) )
    raw.thresholds[ thresh.col ] else 0

  pos.idx <- which( spec.data[ , thresh.col ] > thresh.val )

  if ( length( pos.idx ) < 20L ) {
    warning( paste0( "Insufficient positive events for ", fluor,
                     " (", length( pos.idx ), "). Returning reference spectrum." ) )
    return( original.spectrum )
  }

  # random downsample to n.cells cap
  if ( length( pos.idx ) > n.cells ) {
    set.seed( 42L )
    pos.idx <- sample( pos.idx, n.cells )
  }

  if ( verbose )
    message( paste0( "\033[32m  ", length( pos.idx ),
                     " positive events selected (peak: ", thresh.col, ")\033[0m" ) )

  pos.spec <- spec.data[ pos.idx, spec.present, drop = FALSE ]  # N x D
  pos.scat <- pos.full[  pos.idx, common.scat,  drop = FALSE ]  # N x S

  # -------------------------------------------------------------------------
  # 4. Per-event scatter-matched background subtraction
  # -------------------------------------------------------------------------

  has.scatter <- length( common.scat ) >= 1L && nrow( af.scatter ) > 0L

  if ( has.scatter ) {

    knn.idx <- FNN::knnx.index(
      data  = as.matrix( af.scatter[ , common.scat, drop = FALSE ] ),
      query = as.matrix( pos.scat ),
      k     = k.neighbors
    )

    # average spectral values of the k nearest unstained neighbours per event
    n.ev       <- nrow( knn.idx )
    bg.matched <- matrix( 0, n.ev, length( spec.present ) )
    colnames( bg.matched ) <- spec.present

    for ( ki in seq_len( k.neighbors ) ) {
      bg.matched <- bg.matched +
        af.spectral[ knn.idx[ , ki ], spec.present, drop = FALSE ]
    }
    bg.matched    <- bg.matched / k.neighbors
    pos.corrected <- pos.spec - bg.matched

  } else {

    # Fallback: subtract global AF median when no scatter channels are present
    af.med        <- apply( af.spectral[ , spec.present, drop = FALSE ], 2, stats::median )
    pos.corrected <- sweep( pos.spec, 2, af.med, "-" )

    if ( verbose )
      message( paste0( "\033[33m  No scatter channels; ",
                       "falling back to global AF median subtraction\033[0m" ) )
  }

  # clip physical impossibilities
  pos.corrected[ pos.corrected < 0 ] <- 0

  # -------------------------------------------------------------------------
  # 5. SOM clustering on background-corrected events
  # -------------------------------------------------------------------------

  event.n     <- nrow( pos.corrected )
  som.dim.use <- if ( event.n < 500L )
    max( 2L, floor( sqrt( event.n / 3 ) ) ) else som.dim

  set.seed( 42L )
  som.map <- FlowSOM::SOM(
    pos.corrected,
    xdim   = som.dim.use,
    ydim   = som.dim.use,
    silent = TRUE
  )

  # L-inf normalise SOM centroids -> candidate variant spectra
  variant.spectra <- t( apply(
    som.map$codes, 1,
    function( x ) { mx <- max( x ); if ( mx > 0 ) x / mx else x }
  ) )
  colnames( variant.spectra ) <- spec.present

  # pad to full spectral.channel width if any channels were absent in this file
  if ( !identical( spec.present, spectral.channel ) ) {
    full.mat <- matrix(
      0, nrow( variant.spectra ), length( spectral.channel ),
      dimnames = list( NULL, spectral.channel )
    )
    full.mat[ , spec.present ] <- variant.spectra
    variant.spectra <- full.mat
  }

  # drop all-zero centroids (can arise with very sparse data)
  variant.spectra <- variant.spectra[ rowSums( variant.spectra ) > 0, , drop = FALSE ]

  if ( nrow( variant.spectra ) == 0 ) {
    warning( paste0( "All SOM centroids zero for ", fluor,
                     ". Returning reference spectrum." ) )
    return( original.spectrum )
  }

  # -------------------------------------------------------------------------
  # 6. Cosine similarity QC
  # -------------------------------------------------------------------------

  similar <- apply( variant.spectra, 1, function( v ) {
    sim.mat <- cosine.similarity( rbind( orig.vec, v ) )
    sim.mat[ lower.tri( sim.mat ) ] > sim.threshold
  } )

  if ( !any( similar ) ) {
    warning( paste0(
      "\033[31mNo SOM centroids passed cosine QC (threshold = ",
      sim.threshold, ") for ", fluor,
      ". Returning reference spectrum.\033[0m"
    ) )
    return( original.spectrum )
  }

  variant.spectra <- variant.spectra[ similar, , drop = FALSE ]

  # -------------------------------------------------------------------------
  # 7. Off-peak smoothing toward the reference spectrum
  # -------------------------------------------------------------------------

  # Shrink variant values toward the reference in channels where the
  # fluorophore contributes negligible signal, to avoid inflating cross-talk.
  peak.idx <- orig.vec > 0.05

  variant.spectra <- t( apply( variant.spectra, 1, function( x ) {
    y              <- x
    y[ !peak.idx ] <- 0.5 * x[ !peak.idx ] + 0.5 * orig.vec[ !peak.idx ]
    y
  } ) )

  rownames( variant.spectra ) <- paste0( fluor, "_", seq_len( nrow( variant.spectra ) ) )

  # -------------------------------------------------------------------------
  # 8. Plot
  # -------------------------------------------------------------------------

  if ( figures ) {
    if ( verbose )
      message( paste0( "\033[32m  Plotting spectral variation for ",
                       fluor, "\033[0m" ) )
    spectral.variant.plot(
      spectra.variants   = variant.spectra,
      median.spectrum    = orig.vec,
      title              = paste0( fluor, "_variants" ),
      save               = TRUE,
      plot.dir           = output.dir,
      variant.fill.color = variant.fill.color,
      variant.fill.alpha = variant.fill.alpha,
      median.line.color  = median.line.color,
      median.linewidth   = median.linewidth
    )
    spectral.variant.plot.dens(
      spectra.variants   = variant.spectra,
      median.spectrum    = orig.vec,
      title              = paste0( fluor, "_variants_density" ),
      save               = TRUE,
      plot.dir           = output.dir,
      variant.color = variant.fill.color,
      variant.alpha = variant.fill.alpha,
      median.line.color  = median.line.color,
      median.linewidth   = median.linewidth
    )
  }

  return( variant.spectra )
}
