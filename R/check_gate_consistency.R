# check_gate_consistency.r

#' @title Check Gate Consistency
#'
#' @description
#' A small helper function to check for non-unique input factors when setting
#' gate definitions.
#'
#' @param vec Vector of a factor, such as `large.gate` settings.
#' @param var.name Name of the factor (variable) being checked for multiplicity
#'
#' @seealso
#' * [tune.gate()]
#' * [define.gate.landmarks()]
#' * [define.gate.density()]
#' * [do.gate()]
#' * [define.flow.control()]
#'
#' @return The unique element in the factor `vec`, if all elements of `vec` are
#' identical.
#'
#' @export

check.consistency <- function( vec, var.name ) {
  u <- unique( vec )
  if ( length( u ) > 1 ) {
    stop( paste( "Inconsistent values for", var.name, "found in selected samples." ), call. = FALSE )
  }
  return( u )
}
