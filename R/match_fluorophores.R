# match_fluorophores.r

#' @title Match Fluorophores
#'
#' @description
#' This function matches control filenames to fluorophores in the fluorophore
#' database, including synonyms, and returns the matched fluorophores based on 
#' the longest character match to ensure specificity.
#'
#' @param control.filenames Vector of control filenames.
#' @param fluorophore.database Data frame containing fluorophore information.
#'
#' @return A named vector of matched fluorophores for each control filename.
#'
#' @export

match.fluorophores <- function( control.filenames, fluorophore.database ) {
  
  delim.start <- "(?<![A-Za-z0-9-])"
  delim.end   <- "(?![A-Za-z0-9-])"
  
  fluorophore.matches <- list()
  
  # Columns to check
  fluor.cols <- c( "fluorophore", paste0("synonym", 1:4) )
  
  for ( filename in control.filenames ) {
    
    all.matches <- list()
    
    for ( col in fluor.cols ) {
      vals <- fluorophore.database[[ col ]]
      
      for ( i in seq_along( vals ) ) {
        fluor <- vals[ i ]
        if ( is.na( fluor ) || fluor == "" ) next
        
        # Escape regex metacharacters and handle spaces
        fluor.escaped <- gsub( "([][{}()^$.|*+?\\\\])", "\\\\\\1", fluor )
        fluor.escaped <- gsub( " ", "\\\\s*", fluor.escaped )
        
        pattern <- paste0( delim.start, fluor.escaped, delim.end )
        
        if ( grepl( pattern, filename, ignore.case = TRUE, perl = TRUE ) ) {
          # Store match details to compare length later
          all.matches[[ length( all.matches ) + 1 ]] <- list(
            matched.text = fluor,
            fluorophore  = fluorophore.database$fluorophore[ i ],
            nchar        = nchar( fluor )
          )
        }
      }
    }
    
    # Decide best match based on length
    if ( length( all.matches ) == 0 ) {
      fluorophore.matches[[ filename ]] <- "No match"
      message( sprintf( "\033[31mNo matching fluorophore for: %s\033[0m", filename ) )
      
    } else {
      # Choose the match with the longest string length
      # (e.g., "PE Fire 640" wins over "PE")
      best <- all.matches[[ which.max( sapply( all.matches, `[[`, "nchar" ) ) ]]
      
      fluorophore.matches[[ filename ]] <- best$fluorophore
      message( sprintf(
        "\033[32mFluorophore match: %s -> %s in %s\033[0m",
        best$matched.text, best$fluorophore, filename
      ) )
    }
  }
  
  return( unlist( fluorophore.matches ) )
}