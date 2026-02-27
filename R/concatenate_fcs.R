# concatenate_fcs.r

#' @title Concatenate Multiple FCS Files
#'
#' @description
#' Takes a string/list of input .fcs file paths, reads in the data, combines it
#' and writes a new .fcs file with all the data together.
#'
#' @param fcs.paths A character vector of full file paths to the input .fcs files.
#' @param output.name A character string for the name of the new .fcs file.
#' Default is `Concatenated.fcs`.
#' @param output.dir A character string specifying the directory to save the new
#' file. Default is `./concatenated_fcs`
#'
#' @return Path to the newly created FCS file.
#'
#' @export

concatenateFCS <- function(
    fcs.paths,
    output.name = "Concatenated.fcs",
    output.dir = "./concatenated_fcs"
  ) {

  if ( length( fcs.paths ) < 1 ) stop( "No file paths provided." )

  message( paste( "Concatenating", length( fcs.paths ), "FCS files..." ) )

  # read the first file to establish the column names and keywords
  first.file <- readFCS( fcs.paths[ 1 ], return.keywords = TRUE )
  all.data <- list( first.file$data )
  final.keywords <- first.file$keywords

  col.names.template <- colnames( first.file$data )

  # read subsequent files and collect all data
  if ( length( fcs.paths ) > 1 ) {
    for ( i in 2:length( fcs.paths ) ) {
      current.data <- readFCS( fcs.paths[ i ], return.keywords = FALSE )

      # ensure column counts and names match
      if ( !identical( colnames( current.data ), col.names.template ) ) {
        stop( paste( "Column mismatch in file:", fcs.paths[ i ] ) )
      }

      all.data[[ i ]] <- current.data
    }
  }

  # combine the data
  combined.mat <- do.call( rbind, all.data )

  # update $TOT (Total events) and $FIL (Filename)
  final.keywords[[ "$TOT" ]] <- as.character( nrow( combined.mat ) )
  final.keywords[[ "$FIL" ]] <- output.name

  if ( !dir.exists( output.dir ) ) dir.create( output.dir )

  # Write the file
  writeFCS(
    mat = combined.mat,
    keys = final.keywords,
    file.name = output.name,
    output.dir = output.dir
  )

  return( file.path( output.dir, output.name ) )
}
