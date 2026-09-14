# point_in_polygon.r

#' Point-in-Polygon Test (Ray Casting)
#'
#' Standard PNPOLY ray-casting test, vectorised over query points against a
#' single fixed polygon. Dependency-free replacement for
#' `sp::point.in.polygon()` -- returns a logical vector rather than `sp`'s
#' 0/1/2/3 status code, since every caller in this package only ever tests
#' for "inside or on the boundary" (`!= 0`) or "outside" (`== 0`).
#'
#' @param x Numeric vector, x-coordinates of the query points.
#' @param y Numeric vector, y-coordinates of the query points, same length as `x`.
#' @param poly.x Numeric vector, x-coordinates of the polygon vertices, in order.
#' @param poly.y Numeric vector, y-coordinates of the polygon vertices, same length as `poly.x`.
#'
#' @return Logical vector, same length as `x`/`y`: `TRUE` where the point
#'   falls inside or on the boundary of the polygon.
#'
#' @keywords internal
.point.in.polygon <- function( x, y, poly.x, poly.y ) {
  n <- length( poly.x )
  inside <- rep( FALSE, length( x ) )
  j <- n
  for ( i in seq_len( n ) ) {
    cond <- ( ( poly.y[ i ] > y ) != ( poly.y[ j ] > y ) ) &
      ( x < ( poly.x[ j ] - poly.x[ i ] ) * ( y - poly.y[ i ] ) /
          ( poly.y[ j ] - poly.y[ i ] ) + poly.x[ i ] )
    inside[ cond ] <- xor( inside[ cond ], TRUE )
    j <- i
  }
  inside
}