# match_markers.r

#' @title Match Markers
#'
#' @description
#' This function matches control filenames to markers in the marker database,
#' including synonyms, and returns the matched markers
#'
#' @param control.filenames Vector of control filenames.
#' @param marker.database Data frame containing marker information.
#'
#' @return A named vector of matched markers for each control filename.
#'
#' @export

match.markers <- function( control.filenames, marker.database ) {

  delim.start <- "(?<![A-Za-z0-9-])"
  delim.end   <- "(?![A-Za-z0-9-])"

  marker.matches <- list()

  # columns in order of importance
  antigen.cols <- c(
    "marker", paste0( "synonym", 1:9 )
  )

  for ( filename in control.filenames ) {

    all.matches <- list()

    for ( col in antigen.cols ) {
      vals <- marker.database[[ col ]]

      for ( i in seq_along( vals ) ) {
        antigen <- vals[ i ]
        if ( is.na( antigen ) || antigen == "" ) next

        # escape regex
        antigen.escaped <- gsub( "([][{}()^$.|*+?\\\\])", "\\\\\\1", antigen )
        antigen.escaped <- gsub( " ", "\\\\s*", antigen.escaped )

        pattern <- paste0( delim.start, antigen.escaped, delim.end )

        if ( grepl( pattern, filename, ignore.case = TRUE, perl = TRUE ) ) {
          all.matches[[ length( all.matches ) + 1 ]] <- list(
            antigen = antigen,
            marker  = marker.database$marker[ i ],
            nchar   = nchar( antigen )
          )
        }
      }
    }

    # decide best match
    if ( length( all.matches ) == 0 ) {
      marker.matches[[filename]] <- ""
      message( sprintf( "\033[31mNo matching marker for: %s\033[0m", filename ) )

    } else {
      # choose longest antigen match
      best <- all.matches[[ which.max( sapply( all.matches, `[[`, "nchar" ) ) ]]

      marker.matches[[ filename ]] <- best$marker
      message( sprintf(
        "\033[32mMarker match: %s -> %s in %s\033[0m",
        best$antigen, best$marker, filename
      ) )
    }
  }

  return( unlist( marker.matches ) )
}


