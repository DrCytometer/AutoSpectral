# check_gates.r

#' @title Check Gates For Errors
#'
#' @description
#' A small helper function to check for non-standard or crash-inducing input when
#' users supply gates to `define.flow.control()`.
#'
#' @param gate.list Named list of gates. To use this, pre-define the
#' gates using `define.gate.landmarks()` and/or `define.gate.density()`, ensure
#' that the names of the gates correspond to the names in the `control.def.file`,
#' and ensure that the `gate.name` column has been filled in for the
#' `control.def.file`.
#' @param control.table Dataframe or table of the control file, read in via
#' `define.flow.control()` and cleaned up by that function.
#' @param asp The AutoSpectral parameter list defined using `get.autospectral.param`.
#'
#' @seealso
#' * [tune.gate()]
#' * [define.gate.landmarks()]
#' * [define.gate.density()]
#' * [do.gate()]
#' * [define.flow.control()]
#'
#' @return Silently returns `TRUE` if all checks pass. If any check fails, the
#' pipeline halts.
#'
#' @export

check.gates <- function( gate.list, control.table, asp ) {

  # must be a list
  if ( !is.list(gate.list) ) {
    stop("Pre-defined gates must be supplied as a list object.",
         call. = FALSE)
  }

  # must have names
  if ( is.null(names(gate.list)) || any(names(gate.list) == "") ) {
    stop("All elements in gate.list must be named.",
         call. = FALSE)
  }

  # names must match names in the control.table
  valid.names <- control.table$gate.name
  invalid.supplied <- setdiff( names(gate.list), valid.names )
  if ( length(invalid.supplied) > 0 ) {
    stop(
      paste0(
        "The following gate names do not match those in the control file:",
        "\n",
        paste(invalid.supplied, collapse = "\n")
      ),
      call. = FALSE
    )
  }

  # if excessive gates, stop
  if ( length(gate.list) > length(valid.names) ) {
    stop("gate.list contains more gates than defined in the control table.",
         call. = FALSE)
  }

  # if not enough gates, warn
  missing.gates <- setdiff(valid.names, names(gate.list))
  if ( length(missing.gates) > 0 ) {
    warning(
      paste0(
        "The following expected gates are missing:",
        "\n",
        paste(missing.gates, collapse = "\n ")
      ),
      call. = FALSE
    )
  }

  # structure, empty check, and boundary validation
  lapply(names(gate.list), function(gn) {
    g <- gate.list[[gn]]

    # Check for required components
    if ( !all(c("x", "y", "i") %in% names(g)) ) {
      stop(sprintf("Gate '%s' is missing required components (x, y, or i).", gn),
           call. = FALSE)
    }

    # Check for empty gates
    if ( length(g$x) == 0 ) {
      stop(sprintf("Gate '%s' is empty (zero coordinates provided).", gn),
           call. = FALSE)
    }

    # Check for matching lengths
    if ( length(unique(c(length(g$x), length(g$y), length(g$i)))) > 1 ) {
      stop(sprintf("Lengths of x, y, and i in gate '%s' do not match.", gn),
           call. = FALSE)
    }

    # Check boundaries for X
    if ( any(g$x < asp$scatter.data.min.x) || any(g$x > asp$scatter.data.max.x) ) {
      stop(sprintf("Gate '%s' has X-coordinates outside the allowed range [%f, %f].",
                   gn, asp$scatter.data.min.x, asp$scatter.data.max.x),
           call. = FALSE)
    }

    # Check boundaries for Y
    if ( any(g$y < asp$scatter.data.min.y) || any(g$y > asp$scatter.data.max.y) ) {
      stop(sprintf("Gate '%s' has Y-coordinates outside the allowed range [%f, %f].",
                   gn, asp$scatter.data.min.y, asp$scatter.data.max.y),
           call. = FALSE)
    }
  })

  return( invisible( TRUE ) )
}
