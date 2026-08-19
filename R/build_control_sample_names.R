# build_control_sample_names.r

#' @title Build Unique Control Sample Names
#'
#' @description
#' Derives a unique `sample` identifier for every row of a control table,
#' allowing multiple single-stained controls to be supplied for the same
#' fluorophore. Rows whose `fluorophore` is not duplicated are returned
#' unchanged. For a duplicated fluorophore, uniqueness is resolved in stages,
#' each applied only to the rows still in conflict after the previous stage:
#' \enumerate{
#'   \item `"fluorophore (control.type)"`
#'   \item `"fluorophore (control.type) (marker)"`
#'   \item `"fluorophore (control.type) (marker)_1"`, `"..._2"`, etc.
#' }
#'
#' @param fluorophore Character vector of fluorophore names, one per control.
#' @param control.type Character vector of control types (`"cells"`/`"beads"`),
#' one per control.
#' @param marker Character vector of marker/antigen names, one per control.
#' May contain `NA`.
#'
#' @return A character vector the same length as `fluorophore`, guaranteed to
#' contain no duplicates.

.build.control.sample.names <- function( fluorophore, control.type, marker ) {
  
  sample.name <- fluorophore
  
  dup.fluor <- duplicated( fluorophore ) | duplicated( fluorophore, fromLast = TRUE )
  if ( !any( dup.fluor ) ) return( sample.name )
  
  marker.label <- marker
  marker.label[ is.na( marker.label ) ] <- "unspecified"
  
  # stage 1: fluorophore (control.type)
  candidate <- sample.name
  candidate[ dup.fluor ] <- paste0(
    fluorophore[ dup.fluor ], " (", control.type[ dup.fluor ], ")"
  )
  
  still.dup <- ( duplicated( candidate ) | duplicated( candidate, fromLast = TRUE ) ) & dup.fluor
  
  # stage 2: fluorophore (control.type) (marker)
  if ( any( still.dup ) ) {
    candidate[ still.dup ] <- paste0(
      fluorophore[ still.dup ], " (", control.type[ still.dup ], ") (",
      marker.label[ still.dup ], ")"
    )
    
    still.dup <- ( duplicated( candidate ) | duplicated( candidate, fromLast = TRUE ) ) & dup.fluor
  }
  
  # stage 3: trailing _1, _2, ... within each remaining conflict group
  if ( any( still.dup ) ) {
    for ( f in unique( fluorophore[ still.dup ] ) ) {
      grp.idx <- which( still.dup & fluorophore == f )
      candidate[ grp.idx ] <- paste0( candidate[ grp.idx ], "_", seq_along( grp.idx ) )
    }
  }
  
  candidate
}