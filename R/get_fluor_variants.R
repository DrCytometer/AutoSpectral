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
#' The Spillover Spreading Matrix inputs are now computed differently: instead of
#' re-unmixing this control's positive population against the full panel and
#' taking an empirical MAD, this fits a well-conditioned, low-rank model
#' (autofluorescence components plus this fluorophore's own reference row --
#' or this fluorophore alone for beads, or a cell control too collinear with
#' AF to separate safely) and derives spread from the residual of that fit.
#'
#' The low-rank fit is exactly the AF-projection step \code{get.fluor.variants()}
#' already performs; the difference is that its residual and this control's
#' own recovered abundance are kept rather than discarded. Because the
#' residual is, by construction, orthogonal to every column of this small
#' design, its magnitude is not amplified by collinearity elsewhere in the
#' panel the way a full \code{n}-parameter OLS solve of the raw positive
#' population is. Projecting that residual through the full-panel
#' pseudoinverse (\code{unmix.ols()}) and robustly regressing its squared
#' value against this control's own recovered abundance -- via
#' \code{.fix.huber.slope()} (falling back to \code{fit.robust.linear.model()}
#' if \code{AutoSpectralRcpp} is unavailable) -- gives, for every target
#' fluorophore at once, the same intercept-plus-slope relationship the
#' Residual Model of Cai et al. (bioRxiv 2026, USERM) derives analytically:
#' \code{Var(unmixed_k) = intercept_k + slope_k * b_h}, where \code{b_h} is
#' this control's own on-channel abundance. \code{slope_k} is exactly
#' \code{SS(fluor, k)} in the units used throughout this package (added
#' variance per unit of on-channel abundance), with no separate
#' unstained-population baseline required -- the regression's own intercept
#' already is that baseline, estimated across this control's own dim-to-
#' bright range rather than from a different sample and a different
#' estimator entirely.
#'
#' This uses every event that clears the raw peak-channel threshold
#' (\code{pos.idx}), not just the narrower \code{keep.idx}/cosine-QC'd subset
#' used for SOM clustering, so the regression sees the widest available
#' range of on-channel brightness and does not depend on cosine QC
#' succeeding -- a control that fails cosine QC (and would otherwise fall
#' back to a cruder spread estimate) gets the identical Residual Model
#' treatment here.
#'
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
#' @param af.pcs Named list of autofluorescence-defining principal component
#'   matrices, one per unique unstained FCS file. Names are FCS filenames
#'   matching entries in \code{universal.negative}.
#' @param neg.cache Optional named list, one entry per unique unstained FCS
#'   file (names are filenames), each holding `$spectral` and `$scatter`
#'   matrices already read from disk. Built once in `get.spectral.variants()`
#'   and reused across all fluorophores sharing the same negative, avoiding a
#'   repeat disk read/decode per fluorophore. Default `NULL` falls back to
#'   reading `universal.negative[fluor]` directly.
#' @param use.unmixed Logical, default \code{TRUE}. Whether to unmix
#'   background-corrected positive events against the full \code{spectra}
#'   matrix for positivity selection and SOM clustering features, and whether
#'   the Residual Model's residual-projection step runs at all. Set to
#'   \code{FALSE} when \code{spectra} contains several similar or collinear
#'   fluorophores (e.g. a bead-cell comparison panel), where the full-spectra
#'   unmix is itself unstable or unsolvable. When \code{FALSE}, clustering
#'   falls back to raw detector space only, and the \code{"spillover.spread"}
#'   family of attributes is not computed (see Value).
#' @param n.cells Integer, default \code{10000}. Maximum number of positive
#'   events used for SOM clustering. Files with more events above threshold
#'   are randomly downsampled to this number. Does not limit the event set
#'   used for the Residual Model regression, which uses every event in
#'   \code{pos.idx}.
#' @param som.dim Integer, default \code{10}. Side length of the square SOM
#'   grid. Produces up to \code{som.dim^2} candidate variant spectra before
#'   cosine QC.
#' @param k.neighbors Integer, default \code{3}. Number of scatter-space
#'   nearest neighbours from the unstained pool used to form the per-event
#'   background estimate.
#' @param sim.threshold Numeric, default \code{0.985}. Minimum cosine
#'   similarity between a SOM centroid and the reference spectrum for the
#'   centroid to be retained as a variant.
#' @param sim.threshold.floor Numeric, default \code{0.90}. Lower bound for
#'   adaptive relaxation of \code{sim.threshold} when the initial cutoff
#'   retains fewer than 20 events. Relaxation is logged via \code{warning()}
#'   and the threshold actually used is returned as the
#'   \code{"cosine.threshold.used"} attribute.
#' @param af.collinear.threshold Numeric, default \code{0.95}. Minimum
#'   cosine similarity between \code{fluor}'s reference spectrum and any of
#'   its paired unstained file's AF principal directions (\code{af.pcs}) at
#'   or above which the AF-component projection step -- and the low-rank
#'   fit feeding the Residual Model -- drops the AF-PC terms and uses this
#'   fluorophore's reference row alone, since a joint OLS fit against
#'   near-collinear AF and fluorophore directions can push real fluorophore
#'   signal into the AF term. Recorded as the \code{"af.collinear"} attribute.
#' @param noise.floor.tail.fraction Numeric in (0, 1), default \code{0.20}.
#'   Per-detector noise floor is the MAD (scaled to a Gaussian-equivalent
#'   variance) of the lowest fraction of raw values in that detector's
#'   column, using every event in the control file. Lower values isolate a
#'   purer background tail but with fewer events to estimate from; higher
#'   values are more stable but risk pulling in dim positive events.
#' @param noise.mask.threshold Numeric in (0, 1), default \code{0.05}. A
#'   detector is excluded from this control's contribution to the pooled
#'   noise-model regression (see \code{"noise.mask"} below) when this
#'   fluorophore's own reference spectrum exceeds this fraction of its own
#'   peak there -- the channels where unmodelled spectral-variant wobble,
#'   not photon noise, would otherwise dominate the residual.
#' @param noise.n.cells Integer, default \code{2000L}. Maximum events used
#'   for the noise-model residual (see \code{"noise.resid"} below), sampled
#'   from \code{keep.idx} -- the same unambiguously-positive events already
#'   selected for SOM input. The negative/dim majority of a single-stained
#'   control carries no information the unstained sample doesn't already
#'   supply, and pools to a very large, AF-dominated mass across many
#'   controls; restricting to the positive gate avoids re-fitting that
#'   problem at high cost.
#' @param huber.k Numeric, default \code{1.345}. Huber tuning constant passed
#'   to \code{.fix.huber.slope()} for the spillover-spread regression.
#' @param huber.max.iter Integer, default \code{100L}. Maximum IRLS
#'   iterations passed to \code{.fix.huber.slope()}.
#' @param variant.fill.color Color for the shaded ribbon in the variant plot.
#'   Default \code{"red"}.
#' @param variant.fill.alpha Alpha for \code{variant.fill.color}. Default
#'   \code{0.7}.
#' @param median.line.color Color for the reference-spectrum line. Default
#'   \code{"black"}.
#' @param median.linewidth Width of the reference-spectrum line. Default
#'   \code{1}.
#' @param parallel Logical, default `TRUE`. Enable OpenMP multi-threading
#'   for this call's batch SOM (`get.som.codes()`).
#' @param threads Numeric or `NULL`. OpenMP threads for the batch SOM.
#'   `NULL` defaults to `0` (all available cores) when `parallel = TRUE`.
#'
#' @return A numeric matrix; variants in rows, detectors in columns, values
#' normalised to \eqn{[0, 1]}. Row 1 is always the library reference
#' spectrum for \code{fluor}; subsequent rows are SOM-derived variants that
#' passed cosine QC. When too few positive events are available, or no
#' centroids survive cosine QC, the single reference spectrum is returned
#' (one row), still carrying the full set of attributes below. Carries:
#' \code{"noise.floor"} (per-detector background SD, described above),
#' \code{"noise.events"} (up to \code{noise.n.cells} events x detectors
#' matrix, sampled from \code{keep.idx}, background-corrected -- and
#' AF-corrected where the projection above ran -- but not otherwise fit;
#' pooled and fit jointly against the full panel in
#' \code{get.spectral.variants()}), \code{"noise.mask"} (logical vector,
#' length \code{ncol(spectra)}, \code{TRUE} at detectors this fluorophore's
#' own spectrum dominates), \code{"spillover.spread"} (named numeric vector,
#' one Huber-robust slope per target fluorophore -- \code{SS(fluor, .)} --
#' or \code{NULL} if fewer than 20 events were available or
#' \code{use.unmixed = FALSE}), \code{"spillover.spread.intercept"} (named
#' numeric vector, the matching intercept -- this control's own estimate of
#' each channel's baseline unmixed variance at zero abundance), \code{"on.channel.mfi"}
#' (this control's own median recovered abundance from the low-rank fit,
#' \code{NA} under the same conditions as `"spillover.spread"`),
#' \code{"spillover.spread.n"} (integer, the number of events behind the
#' regression), \code{"spillover.spread.range"} (numeric, the range of this
#' control's own recovered abundance used in the regression -- a control
#' whose events barely vary in brightness gives the slope little to fit
#' against, regardless of event count), and \code{"spillover.spread.source"}
#' (character, \code{"residual"} when computed, \code{NA} otherwise).
#'
#' @references
#' Van Gassen S et al. (2015). FlowSOM. \emph{Cytometry Part A}, 87(7),
#' 636-645. \doi{10.1002/cyto.a.22625}
#'
#' Cai X et al. (2026). Residual Model for unmixed spread prediction in
#' spectral flow cytometry. \emph{bioRxiv} 2026.01.27.701929.

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
    neg.cache          = NULL,
    use.unmixed        = TRUE,
    n.cells            = 10000L,
    som.dim            = 10L,
    k.neighbors        = 3L,
    sim.threshold      = 0.985,
    sim.threshold.floor    = 0.90,
    af.collinear.threshold = 0.95,
    noise.floor.tail.fraction = 0.20,
    noise.mask.threshold = 0.05,
    noise.n.cells = 2000L,
    huber.k         = 1.345,
    huber.max.iter  = 100L,
    variant.fill.color = "red",
    variant.fill.alpha = 0.7,
    median.line.color  = "black",
    median.linewidth   = 1,
    parallel           = TRUE,
    threads            = NULL
) {

  if ( verbose )
    message( paste0( "\033[34m", "Getting spectral variants for ", fluor, "\033[0m" ) )

  # reference spectrum for this fluorophore
  original.spectrum <- spectra[ fluor, , drop = FALSE ]
  orig.vec <- as.numeric( original.spectrum )
  af.collinear <- FALSE

  pos.data <- readFCS(
    file.path( control.dir, file.name[ fluor ] ),
    columns = union( spectral.channel, scatter.channel )
  )
  # remove saturated events
  keep <- rowSums( pos.data[ , spectral.channel ] >= asp$expr.data.max ) == 0
  pos.spectral <- pos.data[ keep, spectral.channel ]
  pos.scatter <- pos.data[ keep, scatter.channel ]

  # get data above threshold in peak channel
  peak.channel <- flow.channel[ fluor ]
  pos.idx <- which( pos.spectral[ , peak.channel ] > raw.thresholds[ peak.channel ] )
  neg.idx <- setdiff( seq_len( nrow( pos.spectral ) ), pos.idx )

  # ---------------------------------------------------------------------------
  # Noise floor upper bound (per detector)
  # ---------------------------------------------------------------------------
  # Estimated from the lowest noise.floor.tail.fraction of raw values in each
  # detector column, using every event in the control file rather than a
  # row-wise positive/negative split. Unchanged from get.fluor.variants().

  noise.floor.est <- stats::setNames(
    rep( NA_real_, length( spectral.channel ) ), spectral.channel )

  if ( nrow( pos.spectral ) >= 200 ) {
    noise.floor.est <- apply(
      pos.spectral, 2,
      function( v ) {
        v <- v[ is.finite( v ) ]
        if ( length( v ) < 200 ) return( NA_real_ )
        cutoff   <- stats::quantile( v, noise.floor.tail.fraction, names = FALSE )
        low.tail <- v[ v <= cutoff ]
        if ( length( low.tail ) < 20 ) return( NA_real_ )
        stats::mad( low.tail, constant = 1.4826 )
      } )
    names( noise.floor.est ) <- spectral.channel
  }

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
    attr( original.spectrum, "noise.floor" ) <- noise.floor.est
    return( original.spectrum )
  }

  # if we have an unstained negative, use that for the negative events
  # check for universal negative, if none, use internal negative
  is.valid.file <- !is.na( universal.negative[ fluor ] ) &&
    universal.negative[ fluor ] != "FALSE" &&
    grepl( "\\.fcs$", universal.negative[ fluor ], ignore.case = TRUE )

  if ( is.valid.file ) {
    neg.file <- universal.negative[ fluor ]

    if ( !is.null( neg.cache ) && neg.file %in% names( neg.cache ) ) {
      neg.data <- cbind( neg.cache[[ neg.file ]]$spectral, neg.cache[[ neg.file ]]$scatter )
    } else {
      neg.path <- file.path( control.dir, neg.file )
      neg.data <- readFCS( neg.path, columns = union( spectral.channel, scatter.channel ) )
    }

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
    pos.corrected <- pos.spectral[ pos.idx, , drop = FALSE ] - bg.matched
  } else {
    # check that we have some internal negative events
    if ( length( neg.idx ) > 50 ) {
      # define mean background from negative fraction
      background <- colMeans( pos.spectral[ neg.idx, , drop = FALSE ] )
      # subtract the global background
      pos.corrected <- pos.spectral - background
    } else {
      pos.corrected <- pos.spectral
    }
  }

  # ---------------------------------------------------------------------------
  # Low-rank AF + fluorophore fit (Residual Model input)
  # ---------------------------------------------------------------------------
  # This control's own recovered abundance (b.h) and the residual of a
  # well-conditioned, low-rank fit: AF PCs plus this fluorophore's own
  # reference row for cells where AF projection is safe, this fluorophore
  # alone otherwise (beads, or a cell control too collinear with AF to
  # separate). The residual is, by construction, orthogonal to every column
  # of this small design -- unlike a full-panel OLS solve of the raw
  # positive population, its magnitude is not amplified by collinearity
  # elsewhere in the panel. Both are reused below to predict spillover
  # spread, instead of re-unmixing this control's events against the full
  # panel and taking an empirical MAD.

  low.rank.fit.done <- FALSE

  # project out any remaining AF (cells only) using AF components
  if ( control.type[ fluor ] == "cells" ) {
    # select the AF PCs derived from this control's paired unstained file;
    # fall back to the first entry if the filename is not found
    neg.fn <- universal.negative[ fluor ]
    if ( is.list( af.pcs ) ) {
      af.pcs.mat <- if ( !is.na( neg.fn ) && neg.fn %in% names( af.pcs ) )
        af.pcs[[ neg.fn ]] else af.pcs[[ 1L ]]
    } else {
      af.pcs.mat <- af.pcs
    }

    # a fluorophore whose reference spectrum is nearly collinear with one of
    # the AF principal directions cannot have AF safely partitioned out by a
    # joint OLS fit against both: the fit is free to push real fluorophore
    # signal into the AF term and vice versa, and the more collinear the two
    # are the more of the fluorophore's own on-peak signal gets absorbed by
    # the AF term instead. Detect that case up front and skip the
    # projection rather than risk silently distorting the reference shape.
    af.pc.cosine <- max( abs( .cosine.sim.rows( af.pcs.mat, orig.vec ) ) )
    af.collinear <- af.pc.cosine >= af.collinear.threshold

    if ( af.collinear ) {

      if ( verbose )
        message( paste0(
          "  ", fluor, " is highly collinear with autofluorescence (cosine = ",
          round( af.pc.cosine, 3 ), " >= af.collinear.threshold = ",
          af.collinear.threshold, "). Skipping AF-component projection; ",
          "relying on scatter-matched background subtraction only."
        ) )

    } else {

      # unmix with this fluor + AF components -- the low-rank design the
      # Residual Model needs. Kept (b.h, resid.lr) before AF is subtracted
      # out of pos.corrected for the SOM/cosine-QC steps that follow.
      low.rank.design <- rbind( af.pcs.mat, original.spectrum )
      pos.unmixed.lr  <- unmix.ols( pos.corrected, low.rank.design )
      b.h             <- as.numeric( pos.unmixed.lr[ , ncol( pos.unmixed.lr ) ] )
      resid.lr        <- pos.corrected - pos.unmixed.lr %*% low.rank.design
      low.rank.fit.done <- TRUE

      # back-project the AF components into raw space
      af.pc.n <- nrow( af.pcs.mat )
      af.projection <- pos.unmixed.lr[ , 1:af.pc.n, drop = FALSE ] %*% af.pcs.mat
      # subtract the projected AF
      pos.corrected <- pos.corrected - af.projection

    }
  }

  # rank-1 fallback (fluorophore alone): beads, or a cell control too
  # collinear with AF for the joint fit above to run.
  if ( !low.rank.fit.done ) {
    pos.unmixed.lr <- unmix.ols( pos.corrected, original.spectrum )
    b.h            <- as.numeric( pos.unmixed.lr[ , 1 ] )
    resid.lr       <- pos.corrected - pos.unmixed.lr %*% original.spectrum
  }

  # unmix background-corrected data in full fluorophore space. Skipped when
  # `use.unmixed = FALSE`, since a full-spectra OLS unmix against several
  # similar or collinear fluorophores (e.g. a bead-cell comparison panel) is
  # itself unstable and would corrupt rather than inform the selection and
  # clustering steps that follow.
  if ( use.unmixed ) {

    pos.unmixed <- unmix.ols( pos.corrected, spectra )

    # select up to n.cells that are still positive
    keep.idx <- which( pos.unmixed[ , fluor ] > unmixed.thresholds[ fluor ] * 2 )

  } else {

    pos.unmixed <- NULL
    # already thresholded on the raw peak channel above (`pos.idx`); retain
    # every background-corrected event
    keep.idx <- seq_len( nrow( pos.corrected ) )

  }

  # ---------------------------------------------------------------------------
  # Noise-model events (downsampled positive population)
  # ---------------------------------------------------------------------------
  # Unchanged from get.fluor.variants(): pooled and fit jointly against the
  # full panel in get.spectral.variants(), independent of the Residual
  # Model spread estimate below.

  noise.idx <- keep.idx
  if ( length( noise.idx ) > noise.n.cells ) {
    set.seed( asp$bird.seed )
    noise.idx <- sample( noise.idx, noise.n.cells )
  }

  if ( length( noise.idx ) < 20 ) {

    noise.events <- NULL
    noise.mask   <- rep( NA, length( spectral.channel ) )

  } else {

    noise.events <- pos.corrected[ noise.idx, , drop = FALSE ]
    noise.mask   <- orig.vec > noise.mask.threshold * max( orig.vec )
  }

  # ---------------------------------------------------------------------------
  # Spillover spread (Residual Model)
  # ---------------------------------------------------------------------------
  # For each target fluorophore k, project this control's low-rank residual
  # through the full-panel pseudoinverse via unmix.ols() -- proj[, k] is
  # exactly a_k' resid.lr, the same linear functional unmix.ols() applies to
  # every event, where a_k is the k-th row of the pseudoinverse of `spectra`.
  # Var(proj[, k] | b.h) = a_k' Sigma_epsilon(b.h) a_k, which the Residual
  # Model takes to be linear in b.h: Sigma_0 + b.h * Sigma_1. Robustly
  # regressing proj[, k]^2 on b.h therefore recovers a_k' Sigma_1 a_k
  # directly as its slope -- exactly SS(fluor, k), in the same
  # variance-per-unit-abundance units used throughout this package -- with
  # a_k' Sigma_0 a_k as its intercept, this control's own estimate of
  # channel k's baseline unmixed variance at zero abundance. No separate
  # unstained-population baseline is subtracted; the regression's own
  # intercept already is that baseline, estimated from this control's own
  # events rather than a different sample and a different estimator.
  #
  # Computed once, here, on every event in pos.idx (not the narrower
  # keep.idx/cosine-QC'd subset used for SOM input below), so a control that
  # ultimately fails cosine QC still gets the identical treatment -- unlike
  # get.fluor.variants(), there is no separate pre-QC fallback estimate.

  spread.mad       <- NULL
  spread.intercept <- NULL
  on.channel.mfi   <- NA_real_
  spread.n         <- NA_integer_
  spread.range     <- NA_real_
  spread.source    <- NA_character_

  if ( use.unmixed && length( b.h ) >= 20 ) {

    proj <- unmix.ols( resid.lr, spectra )

    spread.coef <- apply( proj, 2, function( p ) {
      tryCatch(
        .fix.huber.slope( b.h, p^2, k = huber.k, max.iter = huber.max.iter ),
        error = function( e ) fit.robust.linear.model(
          b.h, p^2, x.name = "b.h", y.name = "proj2", fix.unmix = TRUE
        )
      )
    } )

    spread.intercept <- spread.coef[ 1, ]
    spread.mad       <- spread.coef[ 2, ]
    on.channel.mfi    <- stats::median( b.h )
    spread.n          <- length( b.h )
    spread.range      <- diff( range( b.h ) )
    spread.source     <- "residual"
  }

  # cosine screening. Cosine similarity is scale-invariant, so per-event
  # rescaling before the comparison has no effect on the result -- compare
  # directly against the reference spectrum.
  pos.corrected.keep <- pos.corrected[ keep.idx, , drop = FALSE ]
  ev.cosine <- .cosine.sim.rows( pos.corrected.keep, orig.vec )
  cosine.keep    <- which( ev.cosine >= sim.threshold )
  threshold.used <- sim.threshold

  # Fall back to a relaxed threshold when too few events pass. A fluorophore
  # collinear with AF can have its whole event population sit just under
  # sim.threshold rather than a handful of true outliers -- in that regime a
  # fixed cutoff throws away a real, if noisier, population instead of
  # screening genuine contaminants. Relax in fixed steps down to
  # sim.threshold.floor, logging the threshold actually used so this is
  # visible rather than silent.
  if ( length( cosine.keep ) < 20 && sim.threshold > sim.threshold.floor ) {

    relaxed.grid <- seq( sim.threshold - 0.01, sim.threshold.floor, by = -0.01 )

    for ( relaxed in relaxed.grid ) {
      candidate.keep <- which( ev.cosine >= relaxed )
      if ( length( candidate.keep ) >= 20 ) {
        cosine.keep    <- candidate.keep
        threshold.used <- relaxed
        warning( paste0(
          "Cosine QC for ", fluor, " retained fewer than 20 events at ",
          "sim.threshold = ", sim.threshold, "; relaxed to ",
          round( relaxed, 3 ), " to retain ", length( cosine.keep ),
          " event(s). Inspect this fluorophore's variant plot and ",
          "spillover spread -- it is likely collinear with autofluorescence,",
          "very dim, or otherwise unreliable."
        ), call. = FALSE )
        break
      }
    }
  }

  if ( length( cosine.keep ) < 20 ) {
    warning( paste0( "Insufficient events passed cosine QC for ", fluor,
                     " even after relaxing to sim.threshold.floor = ",
                     sim.threshold.floor, ". Returning reference spectrum." ) )
    attr( original.spectrum, "noise.floor" )             <- noise.floor.est
    attr( original.spectrum, "noise.events" )            <- noise.events
    attr( original.spectrum, "noise.mask" )              <- noise.mask
    attr( original.spectrum, "spillover.spread" )        <- spread.mad
    attr( original.spectrum, "spillover.spread.intercept" ) <- spread.intercept
    attr( original.spectrum, "on.channel.mfi" )          <- on.channel.mfi
    attr( original.spectrum, "spillover.spread.n" )      <- spread.n
    attr( original.spectrum, "spillover.spread.range" )  <- spread.range
    attr( original.spectrum, "spillover.spread.source" ) <- spread.source
    attr( original.spectrum, "cosine.threshold.used" )   <- NA_real_
    attr( original.spectrum, "af.collinear" )            <- af.collinear
    return( original.spectrum )
  }

  if ( use.unmixed ) {
    som.input <- cbind( pos.unmixed[ keep.idx[ cosine.keep ], , drop = FALSE ],
                        pos.corrected[ keep.idx[ cosine.keep ], , drop = FALSE ] )
    colnames( som.input ) <- c( colnames( pos.unmixed ), spectral.channel )
  } else {
    som.input <- pos.corrected[ keep.idx[ cosine.keep ], , drop = FALSE ]
    colnames( som.input ) <- spectral.channel
  }
  event.n   <- length( cosine.keep )

  if ( event.n < 500L )
    som.dim <- max( 2L, floor( sqrt( event.n / 3 ) ) )

  # cluster on the cleaned-up positive fluorophore data
  som.threads <- if ( isTRUE( parallel ) ) {
    if ( is.null( threads ) ) 0L else as.integer( threads )
  } else {
    1L
  }

  map <- get.som.codes(
    data    = som.input,
    som.dim = som.dim,
    seed    = asp$bird.seed,
    threads = som.threads
  )

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

  # Shrink variant values toward the reference in channels where the
  # fluorophore contributes negligible signal, to avoid inflating cross-talk.
  peak.idx <- orig.vec > 0.05

  variant.spectra <- t( apply( variant.spectra, 1, function( x ) {
    y              <- x
    y[ !peak.idx ] <- 0.5 * x[ !peak.idx ] + 0.5 * orig.vec[ !peak.idx ]
    y
  } ) )

  variant.spectra <- rbind(
    original.spectrum,
    variant.spectra
  )

  rownames( variant.spectra ) <- paste0( fluor, "_", seq_len( nrow( variant.spectra ) ) )

  # Plotting
  if ( figures ) {
    if ( verbose )
      message( paste0( "\033[32m  Plotting spectral variation for ",
                       fluor, "\033[0m" ) )
    tryCatch(
      expr = {
        spectral.variant.plot.dens(
          spectra.variants   = variant.spectra,
          median.spectrum    = orig.vec,
          title              = paste0( fluor, "_variants" ),
          save               = TRUE,
          plot.dir           = output.dir,
          variant.color      = variant.fill.color,
          variant.alpha      = variant.fill.alpha,
          median.line.color  = median.line.color,
          median.linewidth   = median.linewidth
        )
      },
      error = function( e ) {
        warning( paste0( "Spectral variant plot failed for ", fluor,
                         ": ", conditionMessage( e ) ) )
      }
    )
  }

  attr( variant.spectra, "noise.floor" )             <- noise.floor.est
  attr( variant.spectra, "noise.events" )            <- noise.events
  attr( variant.spectra, "noise.mask" )              <- noise.mask
  attr( variant.spectra, "spillover.spread" )        <- spread.mad
  attr( variant.spectra, "spillover.spread.intercept" ) <- spread.intercept
  attr( variant.spectra, "on.channel.mfi" )          <- on.channel.mfi
  attr( variant.spectra, "spillover.spread.n" )      <- spread.n
  attr( variant.spectra, "spillover.spread.range" )  <- spread.range
  attr( variant.spectra, "spillover.spread.source" ) <- spread.source
  attr( variant.spectra, "cosine.threshold.used" )   <- threshold.used
  attr( variant.spectra, "af.collinear" )            <- af.collinear
  return( variant.spectra )
}
