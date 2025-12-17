# new_issue.r

# helper function to log issues

#' @title New Issue
#' @description
#' Helper function to log issues in control file
#'
#' @param severity Character, severity of the issue, e.g., "error", "warning".
#' @param rule Character, rule being broken.
#' @param filename Character, name of the FCS file linked to the issue.
#' @param column Character, column of the control file linked to the issue.
#' @param message Character, output message to guide the user in fixing the error.
#'
#' @return A dataframe of issues

.new_issue <- function( severity, rule, filename = NA, column = NA, message ) {
  data.frame(
    severity = severity,
    rule     = rule,
    filename = filename,
    column   = column,
    message  = message,
    stringsAsFactors = FALSE
  )
}
