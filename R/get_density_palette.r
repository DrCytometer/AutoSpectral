# get_density_palette.r

#' @title Get Density Color Palette
#'
#' @description
#' Returns a good color palette for plotting cell density distribution.
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom stats ecdf
#'
#' @param dens A numeric vector representing the density distribution.
#'
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#'
#' @return A vector of colors representing the density palette.
#'
#' @export

get.density.palette <- function( dens, asp )
{
    rainbow.palette <- colorRampPalette(
      asp$density.palette.base.color )( asp$density.palette.base.n )

    dens.range <- range( dens, na.rm = TRUE )

    dens.grid <- seq(
      dens.range[ 1 ], dens.range[ 2 ], length.out = asp$density.palette.n
      )

    density.palette.idx <- round(
      ecdf( dens )( dens.grid ) * asp$density.palette.base.n
      )

    rainbow.palette[ density.palette.idx ]
}
