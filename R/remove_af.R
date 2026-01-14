# remove_af.r

#' @title Remove Autofluorescence Contamination
#'
#' @description
#' This function removes autofluorescence contamination from a sample, using
#' the specified parameters and settings.
#'
#' @importFrom sp point.in.polygon
#' @importFrom flowWorkspace flowjo_biexp
#' @importFrom stats prcomp sd
#'
#' @param samp Sample identifier.
#' @param clean.expr List containing cleaned expression data.
#' @param spectral.channel Vector of spectral channel names.
#' @param peak.channel Vector of peak detection channels for fluorophores.
#' @param universal.negative Named vector mapping samples to their matching
#' negatives.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#' @param scatter.param Vector of scatter parameters.
#' @param negative.n Integer. Number of events to include in the downsampled
#' negative population. Default is `500`.
#' @param positive.n Integer. Number of events to include in the downsampled
#' positive population. Default is `1000`.
#' @param scatter.match Logical, default is `TRUE`. Whether to select negative
#' events based on scatter profiles matching the positive events. Defines a
#' region of FSC and SSC based on the distribution of selected positive events.
#' @param main.figures Logical, if `TRUE` creates the main figures to show the
#' impact of intrusive autofluorescent event removal and scatter-matching for
#' the negatives.
#' @param intermediate.figures Logical, if `TRUE` returns additional figures to
#' show the inner workings of the cleaning, including definition of low-AF cell
#' gates on the PCA-unmixed unstained and spectral ribbon plots of the AF
#' exclusion from the unstained.
#' @param verbose Logical, default is `TRUE`. Set to `FALSE` to suppress messages.
#'
#' @return A matrix containing the expression data with autofluorescent events
#' removed for the sample.

remove.af <- function( samp,
                       clean.expr,
                       spectral.channel,
                       peak.channel,
                       universal.negative,
                       asp,
                       scatter.param,
                       negative.n = 500,
                       positive.n = 1000,
                       scatter.match = TRUE,
                       main.figures = TRUE,
                       intermediate.figures = FALSE,
                       verbose = TRUE ) {

  if ( verbose ) {
    message(
      paste0(
        "\033[34m",
        "Identifying autofluorescence contamination in ",
        samp,
        "\033[0m"
      )
    )
  }

  # match universal negative
  matching.negative <- universal.negative[[ samp ]]

  # get expr data for negative
  expr.data.neg <- clean.expr[[ matching.negative ]]

  if ( nrow( expr.data.neg ) > asp$af.gate.downsample.n.cells ) {
    set.seed( asp$gate.downsample.seed )
    downsample.idx <- sample( nrow( expr.data.neg ), asp$af.gate.downsample.n.cells )
    expr.data.neg <- expr.data.neg[ downsample.idx, ]
  }

  # use PCA to project autofluorescence
  expr.data.center <- apply( expr.data.neg[ , spectral.channel ], 2, median )
  expr.data.scale <- apply( expr.data.neg[ , spectral.channel ], 2, mad )

  scaled.data <- scale( expr.data.neg[ , spectral.channel ], center = expr.data.center,
                        scale = expr.data.scale )

  af.pca <- prcomp( scaled.data, center = FALSE, scale. = FALSE )

  # use the first 2 components to identify main AF signatures to be removed
  af.components <- apply( af.pca$rotation[ , 1:2 ], 2, function( x ) x / max( x ) )
  af.components <- t( apply( af.components, 2, function( x ) x / max( x ) ) )

  # unmix using first two components plus crude fluorophore spectrum
  unmixed.neg <- unmix.ols( expr.data.neg[ , spectral.channel ], af.components )

  biexp.transform <- flowWorkspace::flowjo_biexp(
    channelRange = asp$default.transformation.param$length,
    maxValue = asp$default.transformation.param$max.range,
    pos = asp$default.transformation.param$pos,
    neg = asp$default.transformation.param$neg,
    widthBasis = asp$default.transformation.param$width,
    inverse = FALSE )
  zero.point <- biexp.transform( 0 )

  # transform for easier visualization and faster gating
  unmixed.neg <- apply( unmixed.neg, 2, biexp.transform )
  unmixed.neg <- unmixed.neg - zero.point

  # determine which component has the most variation (intrusive AF)
  af.sd <- apply( unmixed.neg[ , 1:2 ], 2, sd )
  af.order <- order( af.sd, decreasing = TRUE )

  # get gate defining low AF region
  af.gate.idx <- do.gate.af( unmixed.neg[ , af.order ], matching.negative, asp,
                             intermediate.figures )

  # get corresponding raw data (non-PCA)
  af.cells <- expr.data.neg[ -af.gate.idx, , drop = FALSE ]

  # include up to 500 non-AF cells
  if ( length( af.gate.idx ) > 500 ) {
    set.seed( asp$gate.downsample.seed )
    sample.idx <- sample( af.gate.idx, 500 )
    non.af.cells <- expr.data.neg[ sample.idx, ]
  } else {
    non.af.cells <- expr.data.neg[ af.gate.idx, ]
  }

  # get median channel difference as spectrum of intrusive AF
  af.median <- apply( af.cells[ , spectral.channel ], 2, median )
  non.af.median <- apply( non.af.cells[ , spectral.channel ], 2, median )
  af.spectrum <- af.median - non.af.median
  af.spectrum <- af.spectrum / max( abs( af.spectrum ) )

  # flip if the non-AF and AF cells have gotten inverted
  if ( max( abs( af.spectrum ) ) > max( af.spectrum ) )
    af.spectrum <- af.spectrum * -1

  # find peak, subpeak
  af.peak <- which.max( af.spectrum )
  af.subpeak <- which.max( af.spectrum[ af.spectrum < 0.5 ] )

  # get peak channel for fluorophore
  fluor.peak <- peak.channel[[ samp ]]

  # get fluorophore-stained sample data
  expr.data.pos <- clean.expr[[ samp ]]

  if ( fluor.peak == "other" ) {
    fluor.means <- colMeans( expr.data.pos[ , spectral.channel ] )
    fluor.peak <- names( which.max( fluor.means ) )
  }

  # use AF peak if not the same as the fluorophore peak
  if ( names( af.peak ) == fluor.peak )
    af.peak <- af.subpeak

  af.data <- data.frame(
    x = af.cells[ , names( af.peak ) ],
    y = af.cells[ , fluor.peak ]
  )
  non.af.data <- data.frame(
    x = non.af.cells[ , names( af.peak ) ],
    y = non.af.cells[ , fluor.peak ]
  )

  # optionally later pass af.peak, fluor.peak as variable names
  af.boundaries <- fit.af.spline( af.data, non.af.data, asp )

  if ( verbose ) {
    message(
      paste0(
        "\033[34m",
        "Removing autofluorescence contamination in ",
        samp,
        "\033[0m"
      )
    )
  }

  # find events in this bound in the stained sample
  gate.data.pos <- expr.data.pos[ , c( names( af.peak ), fluor.peak ) ]

  gate.population.pip <- sp::point.in.polygon(
    gate.data.pos[ , 1 ], gate.data.pos[ , 2 ],
    af.boundaries$upper$x, af.boundaries$upper$y )

  gate.population.idx <- which( gate.population.pip == 0 )

  # define negative clean-up for plotting and threshold
  gate.data.neg <- expr.data.neg[ , c( names( af.peak ), fluor.peak ) ]

  gate.population.pip <- sp::point.in.polygon(
    gate.data.neg[ , 1 ], gate.data.neg[ , 2 ],
    af.boundaries$upper$x, af.boundaries$upper$y )

  gate.neg.idx <- which( gate.population.pip == 0 )

  # plot data pre/post-removal
  if ( main.figures ) {

    if ( verbose ) {
      message(
        paste0(
          "\033[34m",
          "Plotting AF removal for ",
          samp,
          "\033[0m"
        )
      )
    }

    # set limit for plotting of gate
    af.boundary.ggp <- data.frame(
      x = c( af.boundaries$upper$x,
             af.boundaries$upper$x[ 1 ] ),
      y = c( af.boundaries$upper$y,
             af.boundaries$upper$y[ 1 ] )
    )

    af.boundary.ggp[ af.boundary.ggp > asp$expr.data.max ] <- asp$expr.data.max
    af.boundary.ggp[ af.boundary.ggp < asp$expr.data.min ] <- asp$expr.data.min

    # plot stained control clean-up
    # handle cases of no removal
    if ( length( gate.population.idx ) == nrow( expr.data.pos ) ) {
      removed.data <- matrix(
        0,
        nrow = nrow( expr.data.pos ),
        ncol = ncol( expr.data.pos ),
        dimnames = dimnames( expr.data.pos )
        )
    } else {
      removed.data <- expr.data.pos[ -gate.population.idx, , drop = FALSE ]
    }

    spectral.ribbon.plot(
      pos.expr.data = expr.data.pos,
      neg.expr.data = expr.data.pos[ gate.population.idx, , drop = FALSE ],
      spectral.channel = spectral.channel,
      asp = asp, fluor.name = samp,
      title = asp$af.plot.filename,
      af = TRUE,
      removed.data = removed.data
      )

    # plot AF removal gating on stained control
    gate.af.sample.plot( gate.data.pos, samp, af.boundary.ggp, asp )

    if ( intermediate.figures ) {

      negative.label <- paste( samp, "negative", matching.negative )

      # plot negative clean-up
      if ( length( gate.neg.idx ) == nrow( expr.data.neg ) ) {
        removed.data <- matrix(
          0,
          nrow = nrow( expr.data.neg ),
          ncol = ncol( expr.data.neg ),
          dimnames = dimnames( expr.data.neg )
          )
      } else {
        removed.data <- expr.data.neg[ -gate.neg.idx, , drop = FALSE ]
      }

      spectral.ribbon.plot(
        pos.expr.data = expr.data.neg,
        neg.expr.data = expr.data.neg[ gate.neg.idx, , drop = FALSE ],
        spectral.channel = spectral.channel,
        asp = asp,
        fluor.name = negative.label,
        title = asp$af.plot.filename,
        af = TRUE,
        removed.data = removed.data
        )

      # plot AF removal gating on negative control
      gate.af.sample.plot( gate.data.neg, negative.label, af.boundary.ggp, asp )
    }
  }

  if ( scatter.match ) {

    if ( verbose ) {
      message(
        paste0(
          "\033[34m",
          "Getting scatter-matched negatives for ",
          samp,
          "\033[0m"
        )
      )
    }

    # define positive events as those above a threshold (default 99.5%) in the negative
    if ( samp == "AF" )
      threshold <- asp$positivity.threshold.af
    else
      threshold <- asp$positivity.threshold

    # peak channel data
    pos.peak.channel <- expr.data.pos[ gate.population.idx, fluor.peak ]
    neg.peak.channel <- expr.data.neg[ gate.neg.idx, fluor.peak ]

    # determine threshold for positivity based on cleaned negative
    positivity.threshold <- quantile( neg.peak.channel, threshold )
    pos.above.threshold <- pos.peak.channel[ pos.peak.channel > positivity.threshold ]

    # warn if few events in positive
    if ( length( pos.above.threshold ) < asp$min.cell.warning.n )
      warning(
        paste0(
          "\033[31m",
          "Warning! Fewer than ",
          asp$min.cell.warning.n,
          " positive events in ",
          samp,
          "\033[0m",
          "\n"
          )
        )

    # stop if fewer than minimum acceptable events, returning original data
    if ( length( pos.above.threshold ) < asp$min.cell.stop.n ) {
      warning(
        paste0(
          "\033[31m",
          "Warning! Fewer than ",
          asp$min.cell.stop.n,
          " positive events in ",
          samp,
          "\n",
          "Returning original data",
          "\033[0m",
          "\n"
          )
        )
      return( clean.expr[[ samp ]][ gate.population.idx, ] )
    }

    # select only brightest positive.n events
    if ( length( pos.above.threshold ) >= positive.n )
      pos.selected <- sort(
        pos.above.threshold, decreasing = TRUE )[ 1:positive.n ]
    else
      pos.selected <- pos.above.threshold

    # scatter-match based on positive scatter profile
    pos.selected.expr <- expr.data.pos[ names( pos.selected ), ]
    pos.scatter.coord <- unique( pos.selected.expr[ , scatter.param ] )

    pos.scatter.gate <- suppressWarnings(
      tripack::convex.hull(
        tripack::tri.mesh(
        pos.scatter.coord[ , 1 ],
        pos.scatter.coord[ , 2 ]
      ) ) )

    neg.scatter.matched.pip <- sp::point.in.polygon(
      expr.data.neg[ , scatter.param[ 1 ] ],
      expr.data.neg[ , scatter.param[ 2 ] ],
      pos.scatter.gate$x, pos.scatter.gate$y )

    neg.population.idx <- which( neg.scatter.matched.pip != 0 )

    # warn if few events in negative
    if ( length( neg.population.idx ) < asp$min.cell.warning.n )
      warning(
        paste0(
          "\033[31m",
          "Warning! Fewer than ",
          asp$min.cell.warning.n,
          " scatter-matched negative events for ",
          samp,
          "\033[0m",
          "\n"
        )
      )

    # stop if fewer than minimum acceptable events, returning original negative
    if ( length( neg.population.idx ) < asp$min.cell.stop.n ) {
      warning(
        paste0(
          "\033[31m", "Warning! Fewer than ",
          asp$min.cell.stop.n,
          " scatter-matched negative events for ",
          samp,
          "\n",
          "Reverting to original negative.",
          "\n",
          "\033[0m"
          )
        )

      return( clean.expr[[ samp ]][ gate.population.idx, ] )
    }

    if ( length( neg.population.idx ) > negative.n )
      neg.population.idx <- sample( neg.population.idx, negative.n )

    neg.scatter.matched <- expr.data.neg[ neg.population.idx, ]

    if ( main.figures ) {
      scatter.match.plot(
        pos.expr.data = pos.selected.expr,
        neg.expr.data = neg.scatter.matched,
        fluor.name = samp,
        scatter.param = scatter.param,
        asp = asp
        )

      if ( intermediate.figures )
        spectral.ribbon.plot(
          pos.expr.data = pos.selected.expr,
          neg.expr.data = neg.scatter.matched,
          spectral.channel = spectral.channel,
          asp = asp,
          fluor.name = samp
          )

    }

    return( rbind( pos.selected.expr, neg.scatter.matched ) )
  }

  return( clean.expr[[ samp ]][ gate.population.idx, ] )
}














