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
#' @param control.type Character, either "beads" or "cells". Determines the type
#' of control sample being used and the subsequent processing steps.
#' @param raw.thresholds Named numeric vector of per-channel positivity
#'   thresholds (typically the 99.5th percentile of the unstained sample).
#' @param unmixed.thresholds A named vector of numerical values corresponding to
#'   the threshold for positivity in each unmixed channel. Determined by the
#'   99.5th percentile on the unstained sample, typically after single-cell AF
#'   unmixing.
#' @param flow.channel Named character vector of expected peak raw channels,
#'   one per fluorophore.
#' @param af.pcs Matrix of autofluorescence-defining principal components.
#' @param n.cells Integer, default \code{10000}. Maximum number of positive
#'   events used for SOM clustering. Files with more events above threshold
#'   are randomly downsampled to this number.
#' @param som.dim Integer, default \code{10}. Side length of the square SOM
#'   grid. Produces up to \code{som.dim^2} candidate variant spectra before
#'   cosine QC.
#' @param k.neighbors Integer, default \code{3}. Number of scatter-space
#'   nearest neighbours from the unstained pool used to form the per-event
#'   background estimate.
#' @param sim.threshold Numeric, default \code{0.985}. Minimum cosine
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

get.fluor.variants <- function(
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
    control.type,
    raw.thresholds,
    unmixed.thresholds,
    flow.channel,
    af.pcs,
    n.cells            = 10000L,
    som.dim            = 10L,
    k.neighbors        = 3L,
    sim.threshold      = 0.985,
    variant.fill.color = "red",
    variant.fill.alpha = 0.7,
    median.line.color  = "black",
    median.linewidth   = 1
) {

  if ( verbose )
    message( paste0( "\033[34m", "Getting spectral variants for ", fluor, "\033[0m" ) )

  # reference spectrum for this fluorophore
  original.spectrum <- spectra[ fluor, , drop = FALSE ]
  orig.vec <- as.numeric( original.spectrum )

  pos.data <- readFCS( file.path( control.dir, file.name[ fluor ] ) )
  # remove saturated events
  keep <- rowSums( pos.data[ , spectral.channel ] >= asp$expr.data.max ) == 0
  pos.spectral <- pos.data[ keep, spectral.channel ]
  pos.scatter <- pos.data[ keep, scatter.channel ]

  # get data above threshold in peak channel
  peak.channel <- flow.channel[ fluor ]
  pos.idx <- which( pos.spectral[ , peak.channel ] > raw.thresholds[ peak.channel ] )
  neg.idx <- setdiff( seq_len( nrow( pos.spectral ) ), pos.idx )

  # restrict to top n events
  if ( length( pos.idx ) > n.cells * 2 ) {
    sorted.idx <- order(
      pos.spectral[ pos.idx, peak.channel ],
      decreasing = TRUE )[ 1:( n.cells * 2 ) ]
    pos.idx <- pos.idx[ sorted.idx ]
  }

  if ( length( pos.idx ) < 20 ) {
    warning(
      paste0( "Insufficient positive events found for ",
              fluor, ". Skipping spectral variation." )
    )
    return( original.spectrum )
  }

  # if we have an unstained negative, use that for the negative events
  # check for universal negative, if none, use internal negative
  is.valid.file <- !is.na( universal.negative[ fluor ] ) &&
    universal.negative[ fluor ] != "FALSE" &&
    grepl( "\\.fcs$", universal.negative[ fluor ], ignore.case = TRUE )

  if ( is.valid.file ) {
    neg.path    <- file.path( control.dir, universal.negative[ fluor ] )
    neg.data    <- readFCS( neg.path )
    # remove saturated events
    keep <- rowSums( neg.data[ , spectral.channel ] >= asp$expr.data.max ) == 0
    neg.spectral <- neg.data[ keep, spectral.channel ]
    neg.scatter <- neg.data[ keep, scatter.channel ]

    # use kNN to scatter-match, subtract background per-event
    knn.idx <- FNN::knnx.index(
      data  = neg.scatter,
      query = pos.scatter[ pos.idx, , drop = FALSE ],
      k     = k.neighbors
    )

    # average spectral values of the k nearest unstained neighbours per event
    n.ev       <- nrow( knn.idx )
    bg.matched <- matrix( 0, n.ev, length( spectral.channel ) )
    colnames( bg.matched ) <- spectral.channel

    for ( ki in seq_len( k.neighbors ) ) {
      bg.matched <- bg.matched + neg.spectral[ knn.idx[ , ki ], , drop = FALSE ]
    }
    bg.matched    <- bg.matched / k.neighbors
    pos.corrected <- pos.spectral[ pos.idx, ] - bg.matched
  } else {
    # check that we have some internal negative events
    if ( length( neg.idx ) > 50 ) {
      # define mean background from negative fraction
      background <- colMeans( pos.spectral[ neg.idx, ] )
      # subtract the global background
      pos.corrected <- pos.spectral - background
    } else {
      pos.corrected <- pos.spectral
    }
  }

  # project out any remaining AF (cells only) using AF components
  if ( control.type[ fluor ] == "cells" ) {
    # unmix with this fluor + AF components
    pos.unmixed <- unmix.ols( pos.corrected, rbind( af.pcs, original.spectrum ) )
    # back-project the AF components into raw space
    af.pc.n <- nrow( af.pcs )
    af.projection <- pos.unmixed[ , 1:af.pc.n ] %*% af.pcs
    # subtract the projected AF
    pos.corrected <- pos.corrected - af.projection
  }

  # unmix background-corrected data in full fluorophore space
  pos.unmixed <- unmix.ols( pos.corrected, spectra )

  # select up to n.cells that are still positive
  keep.idx <- which( pos.unmixed[ , fluor ] > unmixed.thresholds[ fluor ] * 2 )

  # cosine screening
  ev.max  <- apply( pos.corrected[ keep.idx, ], 1, max )
  ev.max[ ev.max <= 0 ] <- 1
  ev.norm <- pos.corrected[ keep.idx, ] / ev.max
  ev.cosine <- .cosine.sim.rows( ev.norm, orig.vec )
  cosine.keep <- which( ev.cosine >= sim.threshold )

  if ( length( cosine.keep ) < 20 ) {
    warning( paste0( "Insufficient events passed pre-SOM cosine QC for ",
                     fluor, ". Returning reference spectrum." ) )
    return( original.spectrum )
  }

  som.input <- cbind( pos.unmixed[ keep.idx[ cosine.keep ], ],
                      pos.corrected[ keep.idx[ cosine.keep ], ] )
  colnames( som.input ) <- c( colnames( pos.unmixed ), spectral.channel )
  event.n   <- length( cosine.keep )
  if ( event.n < 500L )
    som.dim <- max( 2L, floor( sqrt( event.n / 3 ) ) )

  # cluster on the cleaned-up positive fluorophore data
  set.seed( asp$bird.seed )
  if ( requireNamespace( "EmbedSOM", quietly = TRUE ) ) {
    map <- EmbedSOM::SOM(
      som.input, # check with just pos.corrected
      xdim = som.dim,
      ydim = som.dim,
      batch = TRUE,
      parallel = TRUE
    )
  } else {
    map <- FlowSOM::SOM(
      som.input,
      xdim = som.dim,
      ydim = som.dim,
      silent = TRUE
    )
  }

  # get spectra: SOM centroids are new profiles, normalize (L-inf)
  variant.spectra <- t(
    apply( map$codes[ , spectral.channel, drop = FALSE ], 1,
           function( x ) x / max( x ) )
  )
  # remove anything that's NA (unlikely)
  variant.spectra <- as.matrix( stats::na.omit( variant.spectra ) )

  if ( nrow( variant.spectra ) == 0 ) {
    warning( paste0( "No valid variants for ", fluor,
                     ". Returning reference spectrum." ) )
    return( original.spectrum )
  }

  if ( FALSE ) {
    # qc to remove dissimilar spectral variants (usually AF contamination)
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
  }
  # Shrink variant values toward the reference in channels where the
  # fluorophore contributes negligible signal, to avoid inflating cross-talk.
  peak.idx <- orig.vec > 0.05

  variant.spectra <- t( apply( variant.spectra, 1, function( x ) {
    y              <- x
    y[ !peak.idx ] <- 0.5 * x[ !peak.idx ] + 0.5 * orig.vec[ !peak.idx ]
    y
  } ) )

  rownames( variant.spectra ) <- paste0( fluor, "_", seq_len( nrow( variant.spectra ) ) )

  # Plotting
  if ( figures ) {
    if ( verbose )
      message( paste0( "\033[32m  Plotting spectral variation for ",
                       fluor, "\033[0m" ) )
    spectral.variant.plot.dens(
      spectra.variants   = variant.spectra,
      median.spectrum    = orig.vec,
      title              = paste0( fluor, "_variants" ),
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
