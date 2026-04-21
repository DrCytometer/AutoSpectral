# calculate_secondary_spillover.R

#' @title Calculate Secondary Stain Index and Spillover
#'
#' @description
#' Calculates the Secondary Stain Index (SSI) and secondary spillover for every
#' fluorophore channel in an unmixed dataset, given a single-stained control and
#' a matched unstained control. SSI is defined as:
#'
#' \deqn{SSI = \frac{median(pos) - median(neg)}{2 \times rSD(neg)}}
#'
#' where *rSD* is the robust standard deviation (median absolute deviation, MAD).
#' Secondary spillover is defined as:
#'
#' \deqn{Spillover = \frac{\Delta MFI_{channel}}{\Delta MFI_{on\text{-}channel}}}
#'
#' Positive events are identified as single-stained cells whose intensity in the
#' on-target channel exceeds the 99th percentile of the unstained control. Up to
#' 500 of the brightest positive events and up to 500 randomly sampled unstained
#' events are used for robustness and speed.
#'
#' @param unstained.unmixed A numeric matrix of unmixed fluorescence values for
#' the unstained control. Rows are cells; columns are fluorophore channels.
#' @param single.stained.unmixed A numeric matrix of unmixed fluorescence values
#' for the single-stained control. Must have the same column names as
#' `unstained.unmixed`.
#' @param fluor.name Character string. The column name of the target fluorophore
#' (the "on-channel") used to define positive events and calculate
#' \eqn{\Delta MFI_{on\text{-}channel}}.
#'
#' @return A data frame with one row per fluorophore channel (matching the column
#' names of `single.stained.unmixed`) and two columns:
#' \describe{
#'   \item{`Spillover`}{Fractional secondary spillover into each channel.}
#'   \item{`SSI`}{Secondary Stain Index for each channel. The SSI for
#'   `fluor.name` itself is set to `0` (on-channel reference).}
#' }
#'
#' @export

calculate.ssi <- function(
    unstained.unmixed,
    single.stained.unmixed,
    fluor.name
  ) {

  # --- input validation ---
  if ( !is.matrix( unstained.unmixed ) && !is.data.frame( unstained.unmixed ) ) {
    stop( "`unstained.unmixed` must be a matrix or data frame.", call. = FALSE )
  }
  if ( !is.matrix( single.stained.unmixed ) && !is.data.frame( single.stained.unmixed ) ) {
    stop( "`single.stained.unmixed` must be a matrix or data frame.", call. = FALSE )
  }
  if ( !fluor.name %in% colnames( single.stained.unmixed ) ) {
    stop(
      paste0( "`fluor.name` ('", fluor.name, "') not found in column names of ",
              "`single.stained.unmixed`." ),
      call. = FALSE
    )
  }
  if ( !fluor.name %in% colnames( unstained.unmixed ) ) {
    stop(
      paste0( "`fluor.name` ('", fluor.name, "') not found in column names of ",
              "`unstained.unmixed`." ),
      call. = FALSE
    )
  }

  # --- positive event selection ---
  # identify positives as events above the 99th percentile of the unstained
  unstained.threshold <- stats::quantile(
    unstained.unmixed[ , fluor.name ], 0.99
  )

  pos.mask <- single.stained.unmixed[ , fluor.name ] > unstained.threshold
  stained.sample <- single.stained.unmixed[ pos.mask, , drop = FALSE ]

  # warn if very few positive events were found; SSI will be unreliable
  if ( nrow( stained.sample ) < 10 ) {
    warning(
      paste0(
        "Fewer than 10 positive events found for '", fluor.name, "'. ",
        "SSI estimates may be unreliable."
      ),
      call. = FALSE
    )
  }

  # take top 500 brightest positive events
  if ( nrow( stained.sample ) > 500 ) {
    top.idx <- order(
      stained.sample[ , fluor.name ], decreasing = TRUE
    )[ 1:500 ]
    stained.sample <- stained.sample[ top.idx, , drop = FALSE ]
  }

  # downsample unstained to 500 events
  if ( nrow( unstained.unmixed ) > 500 ) {
    set.seed( 42 )
    neg.idx <- sample( nrow( unstained.unmixed ), 500 )
    unstained.sample <- unstained.unmixed[ neg.idx, , drop = FALSE ]
  } else {
    unstained.sample <- unstained.unmixed
  }

  # --- SSI and spillover calculation ---
  fluor.median.neg <- stats::median( unstained.sample[ , fluor.name ] )
  fluor.median.pos <- stats::median( stained.sample[ , fluor.name ] )
  delta.mfi.fluor  <- fluor.median.pos - fluor.median.neg

  result <- sapply( colnames( single.stained.unmixed ), function( fl ) {
    unstained.rsd <- stats::mad( unstained.sample[ , fl ] )
    median.neg    <- stats::median( unstained.sample[ , fl ] )
    median.pos    <- stats::median( stained.sample[ , fl ] )
    delta.mfi     <- median.pos - median.neg
    spill         <- delta.mfi / delta.mfi.fluor
    ssi           <- delta.mfi / ( 2 * unstained.rsd )
    c( spill, ssi )
  } )

  result <- data.frame(
    Spillover = result[ 1, ],
    SSI       = result[ 2, ]
  )

  # on-channel SSI is defined as zero (self-reference)
  result[ fluor.name, "SSI" ] <- 0

  return( result )
}
