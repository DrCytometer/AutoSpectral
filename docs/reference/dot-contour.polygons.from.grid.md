# Contour Polygons From a Density Grid

Converts a pre-computed 2D density grid (as returned by
[`AutoSpectralRcpp::fast_kde2d_cpp()`](https://rdrr.io/pkg/AutoSpectralRcpp/man/fast_kde2d_cpp.html)
or [`MASS::kde2d()`](https://rdrr.io/pkg/MASS/man/kde2d.html)) directly
into filled contour polygons using
[`isoband::isobands()`](http://isoband.r-lib.org/reference/isobands.md),
ready to hand to
[`ggplot2::geom_polygon()`](https://ggplot2.tidyverse.org/reference/geom_polygon.html).
This bypasses
[`ggplot2::stat_contour_filled()`](https://ggplot2.tidyverse.org/reference/geom_contour.html),
which would otherwise re-derive the same grid from a long-format melt of
the density matrix before running the identical isobanding step.

## Usage

``` r
.contour.polygons.from.grid(x, y, z, breaks)
```

## Arguments

- x:

  Numeric vector of grid coordinates along the matrix rows of `z`.

- y:

  Numeric vector of grid coordinates along the matrix columns of `z`.

- z:

  Numeric matrix of density values, indexed `[x, y]` (rows correspond to
  `x`, columns to `y`) – the orientation returned by both
  `fast_kde2d_cpp()` and
  [`MASS::kde2d()`](https://rdrr.io/pkg/MASS/man/kde2d.html). Transposed
  internally to match
  [`isoband::isobands()`](http://isoband.r-lib.org/reference/isobands.md)'s
  `[y, x]` convention.

- breaks:

  Numeric vector of band boundaries; consecutive pairs define each
  filled band, matching the `breaks` argument of
  [`ggplot2::geom_contour_filled()`](https://ggplot2.tidyverse.org/reference/geom_contour.html).

## Value

A data frame with columns `x`, `y`, `subgroup` (isoband's own
per-piece/hole identifier, for the `subgroup` aesthetic) and `level` (an
ordered factor spanning every band implied by `breaks`, low to high, so
colour-to-band mapping stays fixed regardless of which bands are empty
in a given plot). Pass as:
`geom_polygon(data = ..., aes(x, y, group = level, subgroup = subgroup, fill = level))`.
Returns `NULL` if `breaks` has fewer than two elements or every band is
empty.
