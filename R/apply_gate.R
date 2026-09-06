# apply_gate.r

#' @title Apply Gate
#'
#' @description
#' Applies a previously-defined scatter gate boundary to a set of flow
#' cytometry expression data and returns only the events falling inside it.
#' This is a lightweight, standalone counterpart to the gating step used
#' internally by `get.gated.flow.expression.data()`, intended for interactive
#' use and for testing other `AutoSpectral` functions against a gated subset
#' of data without going through the full control-file/FCS-reading pipeline.
#'
#' @importFrom sp point.in.polygon
#'
#' @param flow.data A matrix or data frame of flow cytometry data (for
#' example, unmixed or raw expression data) with named columns, including the
#' two scatter parameters named in `scatter.param`.
#' @param gate.boundary A gate boundary, as returned by
#' `define.gate.landmarks()`, `define.gate.density()`, or `do.gate()` — a list
#' containing at least numeric `x` and `y` components describing the polygon
#' vertices.
#' @param scatter.param Character vector of length 2 giving the names of the
#' two scatter columns in `flow.data` to gate on. Default is
#' `asp$default.scatter.parameter`.
#' @param asp The AutoSpectral parameter list, prepared using
#' `get.autospectral.param()`. Only used to supply the default for
#' `scatter.param`; not required if `scatter.param` is supplied directly.
#' @param min.fraction Numeric between `0` and `1`, default `0.01`. If the
#' fraction of events retained by the gate falls below this value, a warning
#' is issued (the function still returns the gated data). Set to `0` to
#' disable this check.
#'
#' @seealso
#' * [define.gate.landmarks()]
#' * [define.gate.density()]
#' * [do.gate()]
#' * [check.gates()]
#' * [get.gated.flow.expression.data()]
#'
#' @return `flow.data`, subset to only those events falling inside
#' `gate.boundary`.
#'
#' @export

apply.gate <- function(
    flow.data,
    gate.boundary,
    scatter.param = asp$default.scatter.parameter,
    asp = NULL,
    min.fraction = 0.01
) {
  
  # flow.data must be present and have named columns
  if ( missing( flow.data ) || is.null( flow.data ) ) {
    stop( "flow.data must be supplied.", call. = FALSE )
  }
  if ( is.null( dim( flow.data ) ) ) {
    stop( "flow.data must be a matrix or data frame.", call. = FALSE )
  }
  if ( is.null( colnames( flow.data ) ) ) {
    stop( "flow.data must have named columns.", call. = FALSE )
  }
  
  # scatter.param must resolve to exactly two names
  if ( is.null( scatter.param ) ) {
    stop(
      "scatter.param was not supplied and could not be determined from asp$default.scatter.parameter. ",
      "Either supply scatter.param directly or supply asp.",
      call. = FALSE
    )
  }
  if ( !is.character( scatter.param ) || length( scatter.param ) != 2 ) {
    stop( "scatter.param must be a character vector of length 2.", call. = FALSE )
  }
  
  # scatter.param columns must exist in flow.data
  missing.param <- setdiff( scatter.param, colnames( flow.data ) )
  if ( length( missing.param ) > 0 ) {
    stop(
      paste0(
        "The following scatter.param column(s) were not found in flow.data: ",
        paste( missing.param, collapse = ", " )
      ),
      call. = FALSE
    )
  }
  
  # gate.boundary must be a well-formed polygon
  if ( !is.list( gate.boundary ) || !all( c( "x", "y" ) %in% names( gate.boundary ) ) ) {
    stop( "gate.boundary must be a list containing x and y components.", call. = FALSE )
  }
  if ( !is.numeric( gate.boundary$x ) || !is.numeric( gate.boundary$y ) ) {
    stop( "gate.boundary$x and gate.boundary$y must be numeric.", call. = FALSE )
  }
  if ( length( gate.boundary$x ) != length( gate.boundary$y ) ) {
    stop( "gate.boundary$x and gate.boundary$y must be the same length.", call. = FALSE )
  }
  if ( length( gate.boundary$x ) < 3 ) {
    stop( "gate.boundary must contain at least 3 vertices to describe a polygon.", call. = FALSE )
  }
  
  n.events <- nrow( flow.data )
  if ( n.events == 0 ) {
    stop( "flow.data contains no events.", call. = FALSE )
  }
  
  gate.data <- flow.data[ , scatter.param ]
  
  flow.population.pip <- sp::point.in.polygon(
    gate.data[ , 1 ], gate.data[ , 2 ],
    gate.boundary$x, gate.boundary$y
  )
  
  gate.population.idx <- which( flow.population.pip != 0 )
  
  # don't allow a gate that excludes every event to pass silently
  if ( length( gate.population.idx ) == 0 ) {
    stop(
      paste0(
        "Gate excludes all ", n.events, " events using scatter.param = ",
        paste( scatter.param, collapse = ", " ), ". ",
        "Check that flow.data and gate.boundary are on the same scale ",
        "and refer to the same cytometer/channels."
      ),
      call. = FALSE
    )
  }
  
  # warn (but don't stop) if the retained fraction looks suspiciously small
  retained.fraction <- length( gate.population.idx ) / n.events
  if ( min.fraction > 0 && retained.fraction < min.fraction ) {
    warning(
      paste0(
        "Gate retained only ", length( gate.population.idx ), " of ", n.events,
        " events (", sprintf( "%.2f%%", 100 * retained.fraction ), "). ",
        "Check flow.data and gate.boundary correspond to the same sample/scale."
      ),
      call. = FALSE
    )
  }
  
  return( flow.data[ gate.population.idx, ] )
}