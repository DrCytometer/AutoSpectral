#' @title Estimate Residual Spillover From a Known-Negative Mask
#' @description
#' Thin, in-memory wrapper over the batched pair estimator, for callers that
#' already know which events are target-negative and do not need the
#' estimator to infer it.
#'
#' @param unmixed Numeric matrix, cells x fluorophores, already-unmixed
#'   abundances. Must contain a column named `source` and one column per
#'   entry of `targets`.
#' @param source Character scalar, the name of the source fluorophore's
#'   column in `unmixed`.
#' @param targets Character vector, the names of the target fluorophores'
#'   columns in `unmixed` to estimate residual spillover into.
#' @param negative.mask Logical vector, length `nrow(unmixed)`, `TRUE` for
#'   events already known to be target-negative. Used directly in place of
#'   the batched estimator's own negative-event inference.
#' @param threshold.source Numeric scalar or vector of length
#'   `nrow(unmixed)`, the source fluorophore's own per-event positivity
#'   boundary.
#' @param threshold.target Optional numeric matrix, cells x fluorophores,
#'   containing at least the columns in `targets`, giving each target's own
#'   per-event positivity boundary. When `NULL`, every event is treated as
#'   below threshold before `negative.mask` is applied.
#' @param spread.var Optional numeric vector, length `length(targets)`, the
#'   source's contribution to each target's spillover-spread variance.
#'   Defaults to zero for every target when `NULL`.
#' @param neg.var Optional numeric vector, length `length(targets)`, each
#'   target's negative-population variance. When `NULL`, computed as
#'   `stats::mad()^2` on each target's column of `unmixed`.
#' @param ... Additional arguments passed through to
#'   `.fix.envelope.slope.batch()`.
#'
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
