# contour_polygons_from_grid.r

#' @title Contour Polygons From a Density Grid
#'
#' @description
#' Converts a pre-computed 2D density grid (as returned by
#' `AutoSpectralRcpp::fast_kde2d_cpp()` or `MASS::kde2d()`) directly into
#' filled contour polygons using `isoband::isobands()`, ready to hand to
#' `ggplot2::geom_polygon()`. This bypasses `ggplot2::stat_contour_filled()`,
#' which would otherwise re-derive the same grid from a long-format melt of
#' the density matrix before running the identical isobanding step.
#'
#' @importFrom isoband isobands
#'
#' @param x Numeric vector of grid coordinates along the matrix rows of `z`.
#' @param y Numeric vector of grid coordinates along the matrix columns of
#' `z`.
#' @param z Numeric matrix of density values, indexed `[x, y]` (rows
#' correspond to `x`, columns to `y`) -- the orientation returned by both
#' `fast_kde2d_cpp()` and `MASS::kde2d()`. Transposed internally to match
#' `isoband::isobands()`'s `[y, x]` convention.
#' @param breaks Numeric vector of band boundaries; consecutive pairs define
#' each filled band, matching the `breaks` argument of
#' `ggplot2::geom_contour_filled()`.
#'
#' @return A data frame with columns `x`, `y`, `subgroup` (isoband's own
#' per-piece/hole identifier, for the `subgroup` aesthetic) and `level` (an
#' ordered factor spanning every band implied by `breaks`, low to high, so
#' colour-to-band mapping stays fixed regardless of which bands are empty in
#' a given plot). Pass as:
#' `geom_polygon(data = ..., aes(x, y, group = level, subgroup = subgroup,
#' fill = level))`. Returns `NULL` if `breaks` has fewer than two elements or
#' every band is empty.

.contour.polygons.from.grid <- function( x, y, z, breaks ) {

  if ( length( breaks ) < 2 ) return( NULL )

  bands <- isoband::isobands(
    x = x, y = y, z = t( z ),
    levels_low  = breaks[ -length( breaks ) ],
    levels_high = breaks[ -1 ]
  )

  # isobands() always returns one list element per requested band, in
  # order, even when a band has no area (length-zero x/y). Keep the full
  # label set as factor levels so colour-to-band mapping stays fixed
  # regardless of which bands happen to be empty in a given plot.
  level.labels <- names( bands )
  has.area     <- vapply( bands, function( b ) length( b$x ) > 0, logical( 1 ) )
  if ( !any( has.area ) ) return( NULL )

  poly.df <- do.call( rbind, lapply( which( has.area ), function( i ) {
    data.frame(
      x        = bands[[ i ]]$x,
      y        = bands[[ i ]]$y,
      subgroup = bands[[ i ]]$id,
      level    = level.labels[ i ]
    )
  } ) )

  poly.df$level <- factor( poly.df$level, levels = level.labels, ordered = TRUE )
  poly.df
}
