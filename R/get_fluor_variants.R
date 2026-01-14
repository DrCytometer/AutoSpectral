# get_fluor_variants.r

#' @title Get Fluorophore Variants
#'
#' @description
#' Assesses variation in the spectral signature of a single-stained flow
#' cytometry control sample. Uses SOM-based clustering on the brightest positive
#' events in the file.
#'
#' @importFrom flowCore read.FCS exprs
#' @importFrom EmbedSOM SOM
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
#'
#' @return A matrix with the flow expression data.

get.fluor.variants <- function( fluor,
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
                                flow.channel ) {

  if ( verbose )
    message( paste0( "\033[34m", "Getting spectral variants for ", fluor, "\033[0m" ) )

  pos.data <- suppressWarnings(
    flowCore::read.FCS(
      file.path( control.dir, file.name[ fluor ] ),
      transformation = NULL,
      truncate_max_range = FALSE,
      emptyValue = FALSE
      )
    )

  # read exprs for spectral channels only
  pos.data <- flowCore::exprs( pos.data )[ , spectral.channel ]

  # get data above threshold in peak channel
  # restrict to top n events
  peak.channel <- flow.channel[ fluor ]

  raw.idx <- which( pos.data[ , peak.channel ] > raw.thresholds[ peak.channel ] )
  neg.idx <- setdiff( seq_len( nrow( pos.data ) ), raw.idx )

  if ( length( raw.idx ) > n.cells * 2 ) {
    sorted.idx <- order(
      pos.data[ raw.idx, peak.channel ],
      decreasing = TRUE )[ 1:( n.cells * 2 ) ]
    raw.idx <- raw.idx[ sorted.idx ]
  }

  if ( "AF" %in% rownames( spectra ) )
    no.af.spectra <- spectra[ rownames( spectra ) != "AF", , drop = FALSE ]
  else
    no.af.spectra <- spectra

  if ( control.type[ fluor ] == "cells" ) {

    # subset to this fluor only
    fluor.spectrum <- no.af.spectra[ fluor, , drop = FALSE ]

    # autospectral per-cell AF unmixing
    pos.unmixed <- unmix.ols( pos.data[ raw.idx, ], fluor.spectrum )

    af.n <- nrow( af.spectra )
    detector.n <- ncol( fluor.spectrum )
    combined.spectra <- matrix( NA_real_, nrow = 2, ncol = detector.n )
    colnames( combined.spectra ) <- colnames( fluor.spectrum )
    fluor.af <- c( fluor, "AF" )
    rownames( combined.spectra ) <- fluor.af
    combined.spectra[ 1, ] <- fluor.spectrum

    # initial residual error
    error <- rowSums( abs( pos.data[ raw.idx, ] - pos.unmixed %*% fluor.spectrum ) )
    # initial AF (0 value)
    initial.af <- matrix( 0, nrow = length( raw.idx ), ncol = 1 )
    colnames( initial.af ) <- c( "AF" )
    pos.unmixed <- cbind( pos.unmixed, initial.af )
    fitted.af <- matrix( 0, nrow = length( raw.idx ), ncol = detector.n )

    for ( af in seq_len( af.n ) ) {
      combined.spectra[ 2, ] <- af.spectra[ af, , drop = FALSE ]
      unmixed.af <- unmix.ols( pos.data[ raw.idx, ], combined.spectra )

      # residual error
      error.af <- rowSums(
        abs(
          pos.data[ raw.idx, ] - ( unmixed.af %*% combined.spectra )
          )
        )
      improved <- which( error.af < error )

      error[ improved ] <- error.af[ improved ]
      pos.unmixed[ improved, fluor.af ] <- unmixed.af[ improved, ]

      fitted.af[ improved, ] <- unmixed.af[ improved, "AF", drop = FALSE ] %*%
        af.spectra[ af, , drop = FALSE ]
    }

    # check for data above threshold in unmixed fluor channel
    pos.idx <- which( pos.unmixed[ , fluor ] > unmixed.thresholds[ fluor ]*2 )

    # check that we still have data; if not, return original spectrum
    if ( length( pos.idx ) < asp$min.cell.stop.n )
      return( spectra[ fluor, , drop = FALSE ] )

    # subtract fitted af component from raw data
    remaining.raw <- pos.data[ raw.idx[ pos.idx ], ] - fitted.af[ pos.idx, ]

    # re-unmix af-subtracted data in full fluorophore space
    pos.unmixed <- unmix.ols( remaining.raw, no.af.spectra )

    # cluster
    som.input <- cbind( pos.unmixed, remaining.raw )
    set.seed( 42 )
    map <- EmbedSOM::SOM( som.input, xdim = som.dim, ydim = som.dim )

    # get spectra
    variant.spectra <- t(
      apply(
        map$codes[ , spectral.channel ], 1, function( x ) x / max( x )
        )
      )
    variant.spectra <- as.matrix( na.omit( variant.spectra ) )
    rownames( variant.spectra ) <- paste0( fluor, "_", 1:nrow( variant.spectra ) )

  } else { # get bead control data and spectral variation
    # get negative background to subtract
    # check for universal negative, if none, use internal negative
    if ( universal.negative[ fluor ] != FALSE ) {
      neg.data <- suppressWarnings(
        flowCore::read.FCS(
          file.path( control.dir, universal.negative[ fluor ] ),
          transformation = NULL,
          truncate_max_range = FALSE,
          emptyValue = FALSE
          )
        )
      neg.data <- flowCore::exprs( neg.data )[ , spectral.channel ]

      # get background on up to 10k events
      if ( nrow( neg.data ) > asp$gate.downsample.n.beads ) {
        set.seed( 42 )
        neg.idx <- sample( nrow( neg.data ), asp$gate.downsample.n.beads )
        background <- apply( neg.data[ neg.idx, spectral.channel ], 2, median )
      } else {
        background <- apply( neg.data[ , spectral.channel ], 2, median )
      }

    } else {
      # select events below unstained raw threshold

      if ( length( neg.idx ) > asp$gate.downsample.n.beads ) {
        # downsample if lots of events
        set.seed( 42 )
        neg.idx <- sample( neg.idx, asp$gate.downsample.n.beads )
        background <- apply( pos.data[ neg.idx, spectral.channel ], 2, median )

      } else if ( length( neg.idx ) > asp$min.cell.warning.n ) {
        # use selected data below threshold if moderate numbers of events
        background <- apply( pos.data[ neg.idx, spectral.channel ], 2, median )
      } else if ( nrow( pos.data ) < asp$min.cell.stop.n ) {
        warning( paste0( "Minimal data present in sample: ", fluor,
                         "Variation assessment not possible for this fluorophore." ) )
        return( spectra[ fluor, , drop = FALSE ] )
      } else {
        # low event sample--take lower half of distribution
        idx.low <- order(
          pos.data[ , peak.channel ], decreasing = FALSE
          )[ 1:( nrow( pos.data ) / 2 ) ]
        neg.selected <- pos.data[ idx.low, spectral.channel, drop = FALSE ]
        background <- apply( neg.selected, 2, median )
      }
    }

    # unmix without AF extraction because beads
    pos.unmixed <- unmix.ols( pos.data[ raw.idx, ], no.af.spectra )

    # cluster
    som.input <- cbind( pos.unmixed, pos.data[ raw.idx, ] )
    set.seed( 42 )
    map <- EmbedSOM::SOM(
      som.input,
      xdim = som.dim,
      ydim = som.dim
      )

    # get spectra, subtracting background from unstained/negative
    variant.spectra <- sweep(
      map$codes[ , spectral.channel ], 2, background, FUN = "-"
      )
    variant.spectra <- t(
      apply(
        variant.spectra, 1, function( x ) x / max( x )
        )
      )
    variant.spectra <- as.matrix( na.omit( variant.spectra ) )
    rownames( variant.spectra ) <- paste0( fluor, "_", 1:nrow( variant.spectra ) )
  }

  # qc to remove dissimilar spectral variants (usually AF contamination)
  original.spectrum <- spectra[ fluor, ]

  similar <- sapply( seq_len( nrow( variant.spectra ) ), function( sp ) {
    sim <- cosine.similarity( rbind( original.spectrum, variant.spectra[ sp, ] ) )
    sim <- sim[ lower.tri( sim ) ]
    sim > 0.99
  } )

  if ( ! any( similar ) )
    variant.spectra <- spectra[ fluor, , drop = FALSE ]
  else
    variant.spectra <- variant.spectra[ similar, , drop = FALSE ]

  # identify true signature areas in original spectrum
  peak.idx <- original.spectrum > 0.05

  # smooth variation towards original outside original spectrum peaks
  variant.spectra.denoised <- t(
    apply( variant.spectra, 1, function( x ) {
    y <- x
    # shrink only off-peak channels
    y[ !peak.idx ] <- 0.5 * x[ !peak.idx ] + 0.5 * original.spectrum[ !peak.idx ]
    y
  } )
  )

  if ( figures )
    spectral.variant.plot(
      variant.spectra.denoised,
      original.spectrum,
      title = paste0( fluor, "_variants" ),
      save = TRUE,
      plot.dir = output.dir
      )

  return( variant.spectra.denoised )
}

