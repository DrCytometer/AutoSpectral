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

  # -- Pass 1: headers only -- event counts and column-consistency check,
  # without reading any binary DATA segment
  headers <- lapply( fcs.paths, function( p ) readFCSheader( p )[[ 1 ]] )

  get.col.names <- function( kw ) {
    n.par <- as.numeric( kw[[ "$PAR" ]] )
    unname( vapply( seq_len( n.par ), function( i ) {
      val <- kw[[ paste0( "$P", i, "N" ) ]]
      if ( is.null( val ) ) paste0( "Channel_", i ) else val
    }, character( 1 ) ) )
  }

  col.names.template <- get.col.names( headers[[ 1 ]] )

  if ( length( fcs.paths ) > 1 ) {
    for ( i in 2:length( fcs.paths ) ) {
      if ( !identical( get.col.names( headers[[ i ]] ), col.names.template ) ) {
        stop( paste( "Column mismatch in file:", fcs.paths[ i ] ) )
      }
    }
  }

  event.counts <- vapply( headers, function( kw ) as.numeric( kw[[ "$TOT" ]] ), numeric( 1 ) )
  total.events <- sum( event.counts )

  # -- Pass 2: pre-allocate the combined matrix once, then fill it file by
  # file -- avoids holding every input file plus a second, rbind-ed copy in
  # memory at the same time
  combined.mat <- matrix( 0, nrow = total.events, ncol = length( col.names.template ) )
  colnames( combined.mat ) <- col.names.template

  row.offset <- 0L
  for ( i in seq_along( fcs.paths ) ) {
    current.data <- readFCS( fcs.paths[ i ], return.keywords = FALSE )
    n.i <- nrow( current.data )
    combined.mat[ ( row.offset + 1L ):( row.offset + n.i ), ] <- current.data
    row.offset <- row.offset + n.i
    rm( current.data )
    if ( i %% 5 == 0 ) gc()
  }

  # readFCSheader() returns a named character vector rather than a list;
  # convert so downstream key assignment/writeFCS() behave identically to
  # the previous readFCS(return.keywords = TRUE)$keywords
  final.keywords <- as.list( headers[[ 1 ]] )

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
