#' @title Estimate Residual Spillover From a Known-Negative Mask
#' @description
#' Thin, in-memory wrapper over the batched pair estimator, for callers that
#' already know which events are target-negative and do not need the
#' estimator to infer it.
#' @export
estimate.residual.spillover <- function(
    unmixed, source, targets, negative.mask,
    threshold.source,
    threshold.target = NULL,
    spread.var       = NULL,
    neg.var          = NULL,
    ...
) {
  x.source <- unmixed[ , source ]
  X.target <- unmixed[ , targets, drop = FALSE ]
  BIG <- 1e12

  Threshold.target <- if ( is.null( threshold.target ) ) {
    matrix( -BIG, nrow = nrow( X.target ), ncol = ncol( X.target ) )
  } else {
    threshold.target[ , targets, drop = FALSE ]
  }
  Threshold.target[ which( negative.mask ) ]  <-  BIG
  Threshold.target[ which( !negative.mask ) ] <- -BIG

  if ( is.null( spread.var ) ) spread.var <- rep( 0, length( targets ) )
  if ( is.null( neg.var ) )
    neg.var <- apply( X.target, 2, function( v ) stats::mad( v )^2 )

  .fix.envelope.slope.batch(
    x.source = x.source, X.target = X.target,
    threshold.source = threshold.source, Threshold.target = Threshold.target,
    spread.var = spread.var, neg.var = neg.var, ...
  )
}
