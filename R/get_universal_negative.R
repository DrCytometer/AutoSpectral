# get_universal_negative.r

#' @title Get Universal Negative Control
#'
#' @description
#' This function identifies and processes the universal negative
#' control for a given sample, including scatter matching and plotting.
#'
#' @importFrom FNN knnx.index
#'
#' @param clean.expr.data List containing cleaned expression data.
#' @param samp Sample identifier.
#' @param universal.negatives Named vector mapping samples to their universal
#' negatives.
#' @param scatter.param Vector of scatter parameters.
#' @param peak.channels Named vector mapping samples to their peak channels.
#' @param downsample Logical indicating whether to downsample the data.
#' @param negative.n Number of negative events to select.
#' @param positive.n Number of positive events to select.
#' @param spectral.channel Vector of spectral channel names.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#' @param control.type Named vector mapping samples to their control types.
#' @param scatter.match Logical indicating whether to perform scatter matching.
#' Default is `TRUE`.
#' @param k.neighbors Numeric, number of scatter-matched unstained events to
#' pair with every positive event for background determination. Default is `3`.
#' @param intermediate.figures Logical, if `TRUE` returns additional figures to
#' show the inner workings of the cleaning, including definition of low-AF cell
#' gates on the PCA-unmixed unstained and spectral ribbon plots of the AF
#' exclusion from the unstained. Default is `FALSE` to speed up processing.
#' @param main.figures Logical, if `TRUE` creates the main figures to show the
#' impact of intrusive autofluorescent event removal and scatter-matching for
#' the negatives.
#' @param verbose Logical, default is `TRUE`. Set to `FALSE` to suppress messages.
#'
#' @return A data frame containing the selected positive and scatter-matched
#' negative events.

get.universal.negative <- function(
    clean.expr.data,
    samp,
    universal.negatives,
    scatter.param,
    peak.channels,
    downsample,
    negative.n,
    positive.n,
    spectral.channel,
    asp,
    control.type,
    scatter.match = TRUE,
    k.neighbors = 3L,
    intermediate.figures = FALSE,
    main.figures = TRUE,
    verbose = TRUE
  ) {

  if ( verbose )
    message( paste0( "\033[34m", "Getting universal negative for: ", samp, "\033[0m" ) )

  pos.control.expr <- clean.expr.data[[ samp ]]
  neg.control.expr <- clean.expr.data[[ universal.negatives[ samp ] ]]

  peak.channel <- peak.channels[ samp ]
  pos.peak.channel <- pos.control.expr[ , peak.channel ]
  neg.peak.channel <- neg.control.expr[ , peak.channel ]

  # define positive events as those above a threshold (default 99.5%) in the negative
  if ( samp == "AF" )
    threshold <- asp$positivity.threshold.af
  else
    threshold <- asp$positivity.threshold

  positivity.threshold <- stats::quantile( neg.peak.channel, threshold )
  pos.above.threshold <- pos.peak.channel[ pos.peak.channel > positivity.threshold ]

  # warn if few events in positive
  if ( length( pos.above.threshold ) < asp$min.cell.warning.n ) {
    warning(
      paste0( "\033[31mWarning! Fewer than ", asp$min.cell.warning.n,
              " positive events in ", samp, "\033[0m", "\n" ),
      call. = FALSE
    )
  }

  # stop if fewer than minimum acceptable events, returning original data
  if ( length( pos.above.threshold ) < asp$min.cell.stop.n ) {
    warning(
      paste0( "\033[31mWarning! Fewer than ", asp$min.cell.stop.n,
              " positive events in ", samp, "\n",
              "Returning original data.\033[0m", "\n" ),
      call. = FALSE
    )
    return( rbind( pos.control.expr, neg.control.expr ) )
  }

  # select only brightest positive.n events
  if ( length( pos.above.threshold ) >= positive.n )
    pos.selected <- sort( pos.above.threshold, decreasing = TRUE )[ 1:positive.n ]
  else
    pos.selected <- pos.above.threshold

  ## scatter-match negative -----
  # recover full data
  pos.selected.expr <- pos.control.expr[ names( pos.selected ), ]

  # find scatter-matched events in the universal negative
  sample.control.type <- control.type[[ samp ]]

  if ( scatter.match ) {
    # get the scatter coordinates of the brightest positive events in the stained control
    pos.scatter.coord <- pos.selected.expr[ , scatter.param ]
    neg.scatter.coord <- neg.control.expr[ , scatter.param ]

    # kNN: for each positive event find its k nearest neighbors in the negative
    knn.idx <- tryCatch(
      FNN::knnx.index(
        data  = neg.scatter.coord,
        query = pos.scatter.coord,
        k     = min( k.neighbors, nrow( neg.scatter.coord ) )
      ),
      error = function( e ) {
        warning( paste( "kNN scatter match failed for", samp, "-",
                        e$message, "- reverting to random sample." ),
                 call. = FALSE )
        NULL
      }
    )

    if ( is.null( knn.idx ) ) {
      # fallback: random subsample of negative
      neg.population.idx <- sample( nrow( neg.control.expr ),
                                    min( negative.n, nrow( neg.control.expr ) ) )
      neg.scatter.matched <- neg.control.expr[ neg.population.idx, ]

    } else {
      neg.population.idx <- unique( as.vector( knn.idx ) )

      # warn / stop checks
      min.neg.n <- ceiling( asp$min.cell.warning.n / 2 )
      if ( length( neg.population.idx ) < min.neg.n ) {
        warning(
          paste0( "\033[31mWarning! Fewer than ", min.neg.n,
                  " kNN-matched negative events for ", samp, "\033[0m\n" ),
          call. = FALSE
        )
      }

      if ( length( neg.population.idx ) < asp$min.cell.stop.n ) {
        warning(
          paste0( "\033[31mWarning! Fewer than ", asp$min.cell.stop.n,
                  " kNN-matched negative events for ", samp,
                  "\nReverting to original negative.\n\033[0m" ),
          call. = FALSE
        )
        if ( nrow( neg.control.expr ) > negative.n )
          neg.control.expr <- neg.control.expr[
            sample( nrow( neg.control.expr ), negative.n ), ]
        return( rbind( pos.selected.expr, neg.control.expr ) )
      }

      neg.scatter.matched <- neg.control.expr[ neg.population.idx, ]
    }

  } else {
    if ( nrow( neg.control.expr ) > negative.n ) {
      neg.population.idx  <- sample( nrow( neg.control.expr ), negative.n )
      neg.scatter.matched <- neg.control.expr[ neg.population.idx, ]
    } else {
      neg.scatter.matched <- neg.control.expr
    }
  }


  if ( main.figures ) {
    asp$ribbon.bins <- 150
    tryCatch(
      expr = {
        scatter.match.plot(
          pos.expr.data = pos.selected.expr,
          neg.expr.data = neg.scatter.matched,
          fluor.name = samp,
          scatter.param = scatter.param,
          asp = asp
        )

        if ( intermediate.figures ) {
          spectral.ribbon.plot(
            pos.expr.data = pos.selected.expr,
            neg.expr.data = neg.scatter.matched,
            spectral.channel = spectral.channel,
            asp = asp,
            fluor.name = samp,
            max.points = 2e5
          )
        }

      },
      error = function( e ) {
        message( "Error in plotting universal negative: ", e$message )
        return( NULL )
      }
    )

  }

  return( rbind( pos.selected.expr, neg.scatter.matched ) )

}
