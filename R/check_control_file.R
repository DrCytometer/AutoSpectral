# check_control_file.r

#' @title Check Control File
#' @description
#' Attempts to find potential failure points and input errors in the control file
#' `control.def.file` used to define the single-stained control setup for
#' AutoSpectral.
#'
#' @param control.dir File path to the single-stained control FCS files.
#' @param control.def.file CSV file defining the single-color control file names,
#' fluorophores they represent, marker names, peak channels, and gating requirements.
#' @param asp The AutoSpectral parameter list defined
#' using `get.autospectral.param`.
#' @param strict Logical. Controls whether the function triggers a break or
#' continues and outputs a list of errors. Default is `FALSE`.
#'
#' @return A dataframe of errors and warnings intended to help the user fix
#' problems with the `control.def.file`.
#'
#' @export

check.control.file <- function(
    control.dir,
    control.def.file,
    asp,
    strict = FALSE
) {

  issues <- validate.control.file( control.dir, control.def.file, asp )

  if ( nrow( issues ) == 0 & strict ) {
    return( invisible( TRUE ) )
  } else if ( nrow( issues ) == 0 & ! strict ) {
    message( "\033[34m No critical errors found in control file.\033[0m" )
    return( issues )
  }

  ## ---- group and summarize ----
  error.rules   <- unique( issues$rule[ issues$severity == "error" ] )
  warning.rules <- unique( issues$rule[ issues$severity == "warning" ] )

  ## ---- long-form messages ----
  long.messages <- list(

    missing_af = function( x ) {
      paste(
        "\033[31m",
        "No autofluorescence (AF) control was detected.",
        "AF controls are strongly recommended for spectral unmixing,",
        "particularly when using primary cell samples,",
        "and are currently required for the AutoSpectral workflow.",
        "\033[0m"
      )
    },

    parameter_mismatch = function( x ) {
      paste(
        "\033[31m",
        "The FCS files do not share a common parameter set.",
        "This usually indicates:",
        "- files exported from different cytometer configurations, or",
        "- post-acquisition parameter filtering.",
        "All control files must have identical parameters in identical order.",
        "\033[0m"
      )
    },

    multiple_cytometers = function( x ) {
      paste(
        "\033[31m",
        "Control files were acquired on multiple cytometers.",
        "Spectral unmixing assumes a single optical configuration.",
        "Mixing cytometers is not supported.",
        "\033[0m"
      )
    }
  )

  ## ---- send out warnings ----
  for ( rule in warning.rules ) {
    if ( !is.null( long.messages[[ rule ]] ) ) {
      warning( long.messages[[ rule ]]( issues[ issues$rule == rule, ] ),
               call. = FALSE )
    }
  }

  ## ---- enforce strict mode ----
  if ( length( error.rules ) > 0 && strict ) {

    msg <- paste(
      "Control file validation failed for the following reasons:\n\n",
      paste(
        vapply(
          error.rules,
          function( r ) {
            if ( !is.null( long.messages[[ r ]] ) ) {
              paste0("- ", long.messages[[ r ]]( issues[ issues$rule == r, ] ) )
            } else {
              paste0( "- ", r )
            }
          },
          character( 1 )
        ),
        collapse = "\n\n"
      )
    )

    stop( msg, call. = FALSE )
  }

  ## ---- non-strict mode: return issues ----
  return( issues )
}
