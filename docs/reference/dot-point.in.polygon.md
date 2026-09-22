# Point-in-Polygon Test (Ray Casting)

Standard PNPOLY ray-casting test, vectorised over query points against a
single fixed polygon. Dependency-free replacement for
[`sp::point.in.polygon()`](https://edzer.github.io/sp/reference/point.in.polygon.html)
– returns a logical vector rather than `sp`'s 0/1/2/3 status code, since
every caller in this package only ever tests for "inside or on the
boundary" (`!= 0`) or "outside" (`== 0`).

## Usage

``` r
.point.in.polygon(x, y, poly.x, poly.y)
```

## Arguments

- x:

  Numeric vector, x-coordinates of the query points.

- y:

  Numeric vector, y-coordinates of the query points, same length as `x`.

- poly.x:

  Numeric vector, x-coordinates of the polygon vertices, in order.

- poly.y:

  Numeric vector, y-coordinates of the polygon vertices, same length as
  `poly.x`.

## Value

Logical vector, same length as `x`/`y`: `TRUE` where the point falls
inside or on the boundary of the polygon.
