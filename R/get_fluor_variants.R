# get_fluor_variants.r

#' @title Get Fluorophore Variants
#'
#' @description
#' Assesses variation in the spectral signature of a single-stained flow
#' cytometry control sample. Uses SOM-based clustering on the brightest positive
#' events in the file.
#'
#' @importFrom FlowSOM SOM
#'
#' @param fluor The name of the fluorophore.
#' @param file.name A named vector of file names for the samples.
#' @param control.dir The directory containing the control files.
#' @param asp The AutoSpectral parameter list.
#' @param spectra A matrix containing the spectral data. Fluorophores in rows,
#' detectors in columns.
#' @param af.spectra Spectral signatures of autofluorescences, normalized
#' between 0 and 1, with fluorophores in rows and detectors in columns. Prepare
#' using `get.af.spectra`.
#' @param n.cells Numeric. Number of cells to use for defining the variation in
#' spectra. Up to `n.cells` cells will be selected as positive events in the peak
#' channel for each fluorophore, above the 99.5th percentile level in the
#' unstained sample.
#' @param som.dim Numeric. Number of x and y dimensions to use in the SOM for
#' clustering the spectral variation.
#' @param figures Logical, controls whether the variation in spectra for each
#' fluorophore is plotted in `output.dir`. Default is `TRUE`.
#' @param output.dir File path to whether the figures and .rds data file will be
#' saved. Default is `NULL`, in which case `asp$variant.dir` will be used.
#' @param verbose Logical, default is `TRUE`. Set to `FALSE` to suppress messages.
#' @param spectral.channel A vector of spectral channels.
#' @param universal.negative A named vector of unstained negative samples, with
#' names corresponding to the fluorophores.
#' @param control.type Character, either "beads" or "cells". Determines the type
#' of control sample being used and the subsequent processing steps.
#' @param raw.thresholds A named vector of numerical values corresponding to
#' the threshold for positivity in each raw detector channel. Determined by the
#' 99.5th percentile on the unstained sample, typically.
#' @param unmixed.thresholds A named vector of numerical values corresponding to
#' the threshold for positivity in each unmixed channel. Determined by the
#' 99.5th percentile on the unstained sample, typically after single-cell AF
#' unmixing.
#' @param flow.channel A named vector of peak raw channels, one per fluorophore.
#' @param refine Logical, default is `TRUE`. Controls whether to perform a second
#' round of variation measurement on "problem cells", which are those with the
#' highest spillover, as defined by `problem.quantile`.
#' @param problem.quantile Numeric, default `0.95`. The quantile for determining
#' which cells will be considered "problematic" after unmixing with per-cell AF
#' extraction. Cells in the `problem.quantile` or above with respect to total
#' signal in the fluorophore (non-AF) channels after per-cell AF extraction will
#' be used to determine additional autofluorescence spectra, using a second round
#' of clustering and modulation of the previously selected autofluorescence
#' spectra. A value of `0.95` means the top 5% of cells, those farthest from zero,
#' will be selected for further investigation.
#'
#' @return A matrix with the flow expression data.
#'
#' @references
#' Van Gassen S et al. (2015). "FlowSOM: Using self-organizing maps for
#' visualization and interpretation of cytometry data." \emph{Cytometry Part A},
#' 87(7), 636-645. \doi{10.1002/cyto.a.22625}
#' Wehrens R, Kruisselbrink J (2018). “Flexible Self-Organizing Maps in kohonen
#' 3.0.” \emph{Journal of Statistical Software}, \emph{87}(7), 1-18.
#' \doi{10.18637/jss.v087.i07}

get.fluor.variants <- function(
    fluor,
    file.name,
    control.dir,
    asp,
    spectra,
    af.spectra,
    n.cells,
    som.dim,
    figures,
    output.dir,
    verbose,
    spectral.channel,
    universal.negative,
    control.type,
    raw.thresholds,
    unmixed.thresholds,
    flow.channel,
    refine = TRUE,
    problem.quantile = 0.95
) {

  if ( verbose )
    message( paste0( "\033[34m", "Getting spectral variants for ", fluor, "\033[0m" ) )

  # get the original "best" spectrum for this fluorophore
  original.spectrum <- spectra[ fluor, , drop = FALSE ]

  # read in the FCS data, spectral channels only
  pos.data <- readFCS( file.path( control.dir, file.name[ fluor ] ) )[ , spectral.channel ]

  # how many detectors do we have?
  detector.n <- ncol( spectra )

  # get data above threshold in peak channel
  peak.channel <- flow.channel[ fluor ]
  raw.idx <- which( pos.data[ , peak.channel ] > raw.thresholds[ peak.channel ] )
  neg.idx <- setdiff( seq_len( nrow( pos.data ) ), raw.idx )

  # restrict to top n events
  if ( length( raw.idx ) > n.cells * 2 ) {
    sorted.idx <- order(
      pos.data[ raw.idx, peak.channel ],
      decreasing = TRUE )[ 1:( n.cells * 2 ) ]
    raw.idx <- raw.idx[ sorted.idx ]
  }

  if ( length( raw.idx ) < 20 ) {
    warning(
      paste0(
        "Insufficient positive events found for ",
        fluor,
        ". Skipping spectral variation."
      )
    )
    return( original.spectrum )
  }

  # remove the AF channel in spectra if present
  if ( "AF" %in% rownames( spectra ) )
    no.af.spectra <- spectra[ rownames( spectra ) != "AF", , drop = FALSE ]
  else
    no.af.spectra <- spectra

  fluorophores <- rownames( no.af.spectra )
  fluorophore.n <- nrow( no.af.spectra )

  ### Processing method depends on whether we use cells or beads ###

  if ( control.type[ fluor ] == "cells" ) {
    #################################################
    ### Identify spectral variation: cell control ###
    #################################################

    ### Cells: per-cell AF extraction is required ###

    # score the AF spectra per cell
    # use residual alignment in the context of the only fluorophore as the metric
    af.assignments <- assign.af.residuals(
      pos.data[ raw.idx, ],
      original.spectrum,
      af.spectra
    )

    # set up for per-cell AF extraction
    af.n <- nrow( af.spectra )
    af.idx <- fluorophore.n + 1
    combined.spectra <- matrix(
      0,
      nrow = af.idx,
      ncol = detector.n
    )
    colnames( combined.spectra ) <- colnames( no.af.spectra )
    fluors.af <- c( fluorophores, "AF" )
    rownames( combined.spectra ) <- fluors.af
    combined.spectra[ 1:fluorophore.n, ] <- no.af.spectra
    fitted.af <- matrix(
      0,
      nrow = length( raw.idx ),
      ncol = detector.n
    )

    # create empty matrix for collection of unmixed data (fluor channel + AF)
    pos.unmixed <- matrix(
      0,
      nrow = length( raw.idx ),
      ncol = af.idx
    )
    colnames( pos.unmixed ) <- fluors.af

    # unmix cells using assigned AF spectra
    for ( af in seq_len( af.n ) ) {
      # set this AF as the spectrum to use
      combined.spectra[ af.idx, ] <- af.spectra[ af, ]

      # get the cells using this AF
      cell.idx <- which( af.assignments == af )
      if( length( cell.idx ) == 0 ) next

      if ( length( cell.idx ) > 0 ) {
        # unmix with this AF
        pos.unmixed[ cell.idx, ] <- unmix.ols(
          pos.data[ raw.idx[ cell.idx ], , drop = FALSE ],
          combined.spectra
        )
        # track AF component
        fitted.af[ cell.idx, ] <- pos.unmixed[ cell.idx, af.idx, drop = FALSE ] %*%
          af.spectra[ af, , drop = FALSE ]
      }
    }

    # check for data above threshold in unmixed fluor channel
    pos.idx <- which( pos.unmixed[ , fluor ] > unmixed.thresholds[ fluor ]*2 )
    event.n <- length( pos.idx )

    # check that we still have data; if not, return original spectrum
    if ( event.n < 20 ) {
      warning(
        paste0(
          "Insufficient positive events found for ",
          fluor,
          ". Skipping spectral variation."
        )
      )

      return( spectra[ fluor, , drop = FALSE ] )
    }

    # subtract fitted af component from raw data
    raw.input <- pos.data[ raw.idx[ pos.idx ], ] - fitted.af[ pos.idx, ]

    # corresponding unmixed data
    unmixed.input <- pos.unmixed[ pos.idx, ]

    # prepare for clustering
    som.input <- cbind( unmixed.input, raw.input )

    # at this point, we can merge the cell and bead pipelines

  } else {
    #################################################
    ### Identify spectral variation: bead control ###
    #################################################

    ### Beads: we just need to subtract the background, no AF ###

    # check for universal negative, if none, use internal negative
    is.valid.file <- !is.na( universal.negative[ fluor ] ) &&
      universal.negative[ fluor ] != "FALSE" &&
      grepl( "\\.fcs$", universal.negative[ fluor ], ignore.case = TRUE )

    if ( is.valid.file ) {
      neg.data <- readFCS(
        file.path( control.dir, universal.negative[ fluor ] )
      )[ , spectral.channel ]

      # get background on up to 10k events
      if ( nrow( neg.data ) > asp$gate.downsample.n.beads ) {
        set.seed( 42 )
        neg.idx <- sample( nrow( neg.data ), asp$gate.downsample.n.beads )
        background <- apply( neg.data[ neg.idx, ], 2, stats::median )
      } else {
        background <- apply( neg.data, 2, stats::median )
      }

    } else {
      # select events below unstained raw threshold

      if ( length( neg.idx ) > asp$gate.downsample.n.beads ) {
        # downsample if lots of events
        set.seed( 42 )
        neg.idx <- sample( neg.idx, asp$gate.downsample.n.beads )
        background <- apply( pos.data[ neg.idx, spectral.channel ], 2, stats::median )

      } else if ( length( neg.idx ) > asp$min.cell.warning.n ) {
        # use selected data below threshold if moderate numbers of events
        background <- apply( pos.data[ neg.idx, spectral.channel ], 2, stats::median )
      } else if ( nrow( pos.data ) < asp$min.cell.stop.n ) {
        warning(
          paste0(
            "Minimal data present in sample: ",
            fluor,
            ". Variation assessment not possible for this fluorophore."
          ),
          call. = FALSE
        )
        # return the original spectrum only
        return( spectra[ fluor, , drop = FALSE ] )

      } else {
        # low event sample--take lower half of distribution
        idx.low <- order(
          pos.data[ , peak.channel ], decreasing = FALSE
          )[ 1:( nrow( pos.data ) / 2 ) ]
        neg.selected <- pos.data[ idx.low, spectral.channel, drop = FALSE ]
        background <- apply( neg.selected, 2, stats::median )
      }
    }

    # unmix without AF extraction because these are beads
    pos.unmixed <- unmix.ols( pos.data[ raw.idx, ], no.af.spectra )

    # prepare for clustering
    unmixed.input <- pos.unmixed
    raw.input <- pos.data[ raw.idx, ]
    som.input <- cbind( unmixed.input, raw.input )
  }


  #############################################
  ### Initial clustering on positive events ###
  #############################################

  # safeguard SOM dimensions based on number of events
  event.n <- nrow( raw.input )
  if ( event.n < 500 )
    som.dim <- max( 2, floor( sqrt( event.n / 3 ) ) )

  # cluster using both unmixed and raw data as input for better discrimination
  set.seed( 42 )
  map <- FlowSOM::SOM(
    som.input,
    xdim = som.dim,
    ydim = som.dim,
    silent = TRUE
  )

  # get spectra: SOM centroids are new profiles, normalize (L-inf)
  variant.spectra <- t(
    apply(
      map$codes[ , spectral.channel ],
      1,
      function( x ) x / max( x )
    )
  )
  # remove anything that's NA (unlikely)
  variant.spectra <- as.matrix( stats::na.omit( variant.spectra ) )

  # qc to remove dissimilar spectral variants (usually AF contamination)
  similar <- sapply( seq_len( nrow( variant.spectra ) ), function( sp ) {
    sim <- cosine.similarity( rbind( original.spectrum, variant.spectra[ sp, ] ) )
    sim <- sim[ lower.tri( sim ) ]
    sim > asp$sim.threshold
  } )

  # revert to the original spectrum only if all the variants look odd
  if ( ! any( similar ) ) {
    warning(
      paste0(
        "\033[31m",
        "No similar spectral variants found for ",
        fluor,
        "\033[0m"
      )
    )
    return( original.spectrum )
  } else {
    # otherwise restrict to the similar variants only
    variant.spectra <- variant.spectra[ similar, , drop = FALSE ]
  }

  # identify true signature areas in original spectrum
  peak.idx <- original.spectrum > 0.05

  # smooth variation towards original outside original spectrum peaks
  variant.spectra <- t(
    apply( variant.spectra, 1, function( x ) {
      y <- x
      # shrink only off-peak channels
      y[ !peak.idx ] <- 0.5 * x[ !peak.idx ] + 0.5 * original.spectrum[ !peak.idx ]
      y
    } )
  )

  # add names for tracking
  rownames( variant.spectra ) <- paste0( fluor, "_", 1:nrow( variant.spectra ) )

  if ( refine ) {
    ##########################################################
    ### Refine spectral variation: identify problem events ###
    ##########################################################

    # To perform the refinement, we first need to run one pass of AutoSpectral's
    # per-cell fluorophore optimization for this fluorophore. This will tell us
    # which cells are still problematic (creating spillover into other channels)
    # when only this first batch of variants is used. We can then focus on these
    # cells to perform a second round of variant selection.

    # Note: At some point, this should probably include plotting of the changes as
    # in the updated get.af.spectra.

    # since this is a single-color control, we know only the single color is present
    unmixed.fl <- unmix.ols( raw.input, original.spectrum )
    resid <- raw.input - ( unmixed.fl %*% original.spectrum )
    initial.error <- rowSums( abs( resid ) )

    # now we pick which of the first batch of variants best aligns with the residual per cell
    # calculate deltas (change and normalized distance)
    delta.fl <- variant.spectra - matrix(
      original.spectrum,
      nrow = nrow( variant.spectra ),
      ncol = detector.n,
      byrow = TRUE
    )
    delta.norm <- sqrt( rowSums( delta.fl^2 ) )

    # delta.fl %*% t(resid) calculates the dot product for every row of resid
    dot.products <- delta.fl %*% t( resid )

    numerator <- sweep( dot.products, 2, unmixed.fl, "*" )
    # row-wise norms
    resid.norms <- sqrt( rowSums( resid^2 ) )
    denom.matrix <- outer( delta.norm, resid.norms, "*" )

    # score for alignment
    joint.score.matrix <- numerator / denom.matrix
    # find which factor (column) is the "winner" for each cell
    best.variant.per.cell <- max.col( t( joint.score.matrix ), ties.method = "first" )
    # track actual assignments (starting with 1 for all cells--original spectrum)
    var.assignments <- rep( 0, nrow( raw.input ) )

    # set up a spectral matrix for testing variants
    spectra.curr <- no.af.spectra

    # store new unmixed data, starting with original spectrum unmixing
    unmixed.first.pass <- unmixed.input[ , fluorophores ]

    # only do this if the best variant has changed for at least one cell
    if ( ! all( best.variant.per.cell == 1 ) ) {
      # unmix in full space, cycling through variants and assigning cells
      for ( var in 1:nrow( variant.spectra ) ) {
        # which cells have been assigned this variant?
        cell.idx <- which( best.variant.per.cell == var )

        if ( length( cell.idx ) > 0 ) {
          # supplant the base spectrum with this variant
          spectra.curr[ fluor, ] <- variant.spectra[ var, ]

          # re-unmix with this variant
          trial.unmix <- unmix.ols( raw.input[ cell.idx, , drop = FALSE ], spectra.curr )

          # assess the residual error with this variant
          trial.resid <- raw.input[ cell.idx, , drop = FALSE ] -
            ( trial.unmix %*% spectra.curr )
          trial.error <- rowSums( abs( trial.resid ) )

          # assess which cells are actually better (lower residual)
          improved <- which( trial.error < initial.error[ cell.idx ] )

          # update any improved cells
          if ( length( improved ) > 0 ) {
            unmixed.first.pass[ cell.idx[ improved ], ] <- trial.unmix[ improved, ]
            var.assignments[ cell.idx[ improved ] ] <- var
          }
        }
      }
    }

    # We can define the error (spillover spread) as the signal in other fluorophore
    # channels. Technically, this might be better as the signal above the unstained
    # (say 99.5th percentile), but we can just pick the most problematic events
    # (those farthest from zero in other channels) as the error, which avoids any
    # need for thresholds. Thresholds can be problematic for various reasons.

    # we want to assess this error on the cells after they have undergone one round
    # of unmixing improvement using the variant spectral optimization above, so we
    # use unmixed.first.pass

    # which are the other channels?
    other.channels <- which( rownames( no.af.spectra ) != fluor )

    # project signal in other channels back into raw space
    spillover.spread.proj <- unmixed.first.pass[ , other.channels, drop = FALSE ] %*%
      no.af.spectra[ other.channels, , drop = FALSE ]

    ### run a second round of clustering on cells that are far from zero ###
    if ( verbose )
      message(
        paste0(
          "\033[32m",
          "Calculating target spectra for problematic cells",
          "\033[0m"
        )
      )

    # in a single-stained sample, any other fluorophore signal is error.
    if ( fluorophore.n > 2 ) {
      error.magnitude <- sqrt( rowSums( unmixed.first.pass[, other.channels ]^2 ) )
    } else {
      error.magnitude <- abs( unmixed.first.pass[ , other.channels ] )
    }

    # filter for problematic cells
    threshold <- stats::quantile( error.magnitude, problem.quantile )
    problem.idx <- which( error.magnitude > threshold )
    # with an initial input of 2000 cells, we can expect no more than 100 in the top 5%
    problem.cell.n <- length( problem.idx )

    # if we have fewer than 100 cells, adjust problem threshold until we have more
    while ( problem.cell.n < 100 ) {
      problem.quantile <- problem.quantile - 0.05

      # if we hit 50%, adjust SOM dimensions and exit
      if ( problem.quantile < 0.5 ) break

      # update for the next iteration check
      threshold <- stats::quantile( error.magnitude, problem.quantile )
      problem.idx <- which( error.magnitude > threshold )
      problem.cell.n <- length( problem.idx )
    }

    if ( problem.cell.n < 100 ) {
      if ( verbose ) {
        message(
          paste0(
            "\033[33m",
            "Insufficient error-prone events found. Skipping modulation.",
            "\033[0m"
          )
        )
      }

      if ( figures ) {
        if ( verbose )
          message(
            paste0(
              "\033[32m",
              "Plotting spectral variation for ",
              fluor,
              "\033[0m"
            )
          )
        spectral.variant.plot(
          spectra.variants = variant.spectra,
          median.spectrum = as.numeric( original.spectrum ),
          title = paste0( fluor, "_variants" ),
          save = TRUE,
          plot.dir = output.dir
        )
      }

      return( variant.spectra )
    }

    # current fluorophore estimate
    fl.abundance <- unmixed.first.pass[ problem.idx, fluor ]
    fl.abundance[ fl.abundance == 0 ] <- 1 # safeguard 0 division

    # what are the error signatures?
    spill.ratios <- spillover.spread.proj[ problem.idx, ] / fl.abundance

    if ( verbose ) {
      message(
        paste0(
          "\033[32m",
          "Clustering ",
          length( problem.idx ),
          " events with highest spillover error",
          "\033[0m"
        )
      )
    }

    # use a minimum of 2 dimensions, or around 3 cells per node
    som.dim <- max( 2, floor( sqrt( problem.cell.n / 3 ) ) )

    # cluster only the problematic data
    set.seed( 42 )
    map.error <- FlowSOM::SOM(
      spill.ratios,
      xdim = som.dim,
      ydim = som.dim,
      silent = TRUE
    )

    # for each cluster, we find which base variants were assigned and update them
    cluster.ids <- unique( map.error$mapping[ , 1 ] )

    modulated.list <- lapply( cluster.ids, function( cl ) {
      cl.sub.idx <- which( map.error$mapping[ , 1 ] == cl )
      global.idx <- problem.idx[ cl.sub.idx ]

      # get the median correction pattern for this cluster
      median.ratio <- apply(
        spill.ratios[ cl.sub.idx, , drop = FALSE ],
        2,
        stats::median
      )

      # identify all unique base fluor spectra that fell into this error cluster
      contributing.var.ids <- unique( var.assignments[ global.idx ] )

      # create a new version of each contributing AF using this cluster's ratio
      new.specs <- lapply( contributing.var.ids, function( id ) {
        if ( id == 0 )
          base.spec <- original.spectrum
        else
          base.spec <- variant.spectra[ id, ]
        updated <- base.spec * ( 1 + median.ratio )
        return( updated / max( updated ) ) # Re-normalize
      } )

      return( do.call( rbind, new.specs ) )
    } )

    modulated.var.spectra <- do.call( rbind, modulated.list )
    variant.spectra <- rbind( variant.spectra, modulated.var.spectra )

    # unlikely, but check for NAs
    variant.spectra <- as.matrix( stats::na.omit( variant.spectra ) )

    # add names
    rownames( variant.spectra ) <- paste0( fluor, "_", 1:nrow( variant.spectra ) )

    #####################
    ### Similarity QC ###
    #####################

    # qc to remove dissimilar spectral variants (usually AF contamination)
    similar <- sapply( seq_len( nrow( variant.spectra ) ), function( sp ) {
      sim <- cosine.similarity( rbind( original.spectrum, variant.spectra[ sp, ] ) )
      sim <- sim[ lower.tri( sim ) ]
      sim > 0.99
    } )

    # revert to the original spectrum only if all the variants look odd
    # otherwise restrict to the similar variants only
    if ( ! any( similar ) )
      variant.spectra <- spectra[ fluor, , drop = FALSE ]
    else
      variant.spectra <- variant.spectra[ similar, , drop = FALSE ]

    # smooth variation towards original outside original spectrum peaks
    variant.spectra <- t(
      apply( variant.spectra, 1, function( x ) {
        y <- x
        # shrink only off-peak channels
        y[ !peak.idx ] <- 0.5 * x[ !peak.idx ] + 0.5 * original.spectrum[ !peak.idx ]
        y
      } )
    )
  }


  ################
  ### Plotting ###
  ################

  # plot the variation in the spectrum for this fluorophore
  if ( figures ) {
    if ( verbose )
      message(
        paste0(
          "\033[32m",
          "Plotting spectral variation for ",
          fluor,
          "\033[0m"
        )
      )
    spectral.variant.plot(
      spectra.variants = variant.spectra,
      median.spectrum = as.numeric( original.spectrum ),
      title = paste0( fluor, "_variants" ),
      save = TRUE,
      plot.dir = output.dir
    )
  }

  return( variant.spectra )
}

