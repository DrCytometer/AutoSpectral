# resolve_file_name.R

#' @title Resolve File Name
#'
#' @description
#' Determines the file name to use for output naming and plot titles.
#' Prefers the name embedded in the FCS `$FIL` keyword, but falls back to
#' the actual file name on disk when the file has been renamed since
#' acquisition -- a common occurrence that otherwise leaves output files
#' and plot titles referring to a name the user no longer recognizes.
#'
#' @param path Character. Path to the FCS file as provided by the caller.
#' @param fil Character or `NULL`. The `$FIL` keyword value read from the
#'   file, if present.
#' @param verbose Logical, default `TRUE`. Whether to message when a
#'   mismatch is detected and resolved.
#'
#' @return Character, the resolved file name (without directory path).
#'
#' @keywords internal

resolve.file.name <- function( path, fil, verbose = TRUE ) {
  
  disk.name <- basename( path )
  
  # no $FIL keyword at all -- disk name is the only option
  if ( is.null( fil ) || !nzchar( fil ) ) return( disk.name )
  
  # nothing usable on disk (shouldn't normally happen, readFCS would have
  # already failed) -- fall back to $FIL
  if ( !nzchar( disk.name ) ) return( fil )
  
  # compare names ignoring extension, since some instruments write $FIL
  # without ".fcs" while the file on disk always has it
  fil.stem  <- tools::file_path_sans_ext( fil )
  disk.stem <- tools::file_path_sans_ext( disk.name )
  
  if ( identical( fil.stem, disk.stem ) ) return( fil )
  
  # mismatch: the file has been renamed since acquisition. Prefer the name
  # on disk, since that's what the user sees and expects in output.
  if ( verbose )
    message(
      sprintf(
        "File name mismatch: FCS keyword $FIL is \"%s\" but the file on disk is \"%s\". Using the name on disk.",
        fil, disk.name
      )
    )
  
  disk.name
}