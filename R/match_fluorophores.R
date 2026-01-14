# match_fluorophores.r

#' @title Match Fluorophores
#'
#' @description
#' This function matches control filenames to fluorophores in the fluorophore
#' database, including synonyms, and returns the matched fluorophores.
#'
#' @param control.filenames Vector of control filenames.
#' @param fluorophore.database Data frame containing fluorophore information,
#' including synonyms.
#'
#' @return A named vector of matched fluorophores for each control filename.
#'
#' @export

match.fluorophores <- function( control.filenames, fluorophore.database ) {

  delim.start <- "(?<![A-Za-z0-9-])"
  delim.end   <- "(?![A-Za-z0-9-])"

  fluorophore.matches <- list()

  for ( filename in control.filenames ) {
    fluorophore <- character( 0 )

    # Columns to check in order of priority
    fluor.cols <- c( "fluorophore", "synonym1", "synonym2", "synonym3" )

    for ( col in fluor.cols ) {
      vals <- fluorophore.database[[ col ]]

      for ( i in seq_along( vals ) ) {
        fluor <- vals[ i ]
        if ( is.na( fluor ) || fluor == "" ) next

        # Escape regex metacharacters
        fluor.escaped <- gsub( "([][{}()^$.|*+?\\\\])", "\\\\\\1", fluor )
        fluor.escaped <- gsub( " ", "\\\\s*", fluor.escaped)

        pattern <- paste0( delim.start, fluor.escaped, delim.end )

        if ( grepl( pattern, filename, ignore.case = TRUE, perl = TRUE ) ) {
          fluorophore <- fluorophore.database$fluorophore[ i ]
          message(
            paste0(
              "\033[32mMatch: ",
              fluor, " to ", fluorophore, " in ", filename,
              "\033[0m"
              )
            )
          break
        }
      }

      if ( length( fluorophore ) > 0 ) break
    }

    if ( length( fluorophore ) == 0) {
      fluorophore <- "No match"
      message( paste0( "\033[31m", "No Match for: ", filename, "\033[0m" ) )
    }

    fluorophore.matches[ filename ] <- fluorophore
  }

  fluorophore.matches <- unlist( fluorophore.matches )

  return( fluorophore.matches )
}
