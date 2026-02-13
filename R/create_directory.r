# create_directory.r

#' @title Create Directory
#'
#' @description
#' Creates figure and table directories.
#'
#' @param asp The AutoSpectral parameter list defined
#' using `get.autospectral.param`
#'
#' @return No returns. Creates directories in the current working directory.
#' @export

create.directory <- function( asp ) {

    figure.dir <- c(
      asp$figure.gate.dir,
      asp$figure.spectra.dir,
      asp$figure.similarity.heatmap.dir
    )

    for ( ftd in c( figure.dir, asp$table.spectra.dir, asp$unmixed.fcs.dir ) )
        if ( ! is.null( ftd ) && ! file.exists( ftd ) )
            dir.create( ftd, recursive = TRUE )
}

