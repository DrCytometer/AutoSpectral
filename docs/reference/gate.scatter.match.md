# Gate Matching Scatter

Calculates the density of the input scatter coordinates, cuts at the
defined density percentile and calculates a convex hull boundary around
the points within the dense region. Permits selection of events with
similar forward and side scatter profiles in a paired sample, e.g., an
unstained sample.

## Usage

``` r
gate.scatter.match(scatter.data, percentile = 0.5, bw.factor = 1, grid.n = 10)
```

## Arguments

- scatter.data:

  Matrix of forward and side scatter data

- percentile:

  Numeric 0-1. Percentile of density to use for defining the contour for
  the gate boundary. Default is `0.5`.

- bw.factor:

  Numeric, scalar bandwidth multiplier. Larger numbers smooth the
  density, expanding the gate outwards. Default is `1`.

- grid.n:

  Numeric, size of the grid for defining the kernel density estimate.
  Default is `10`.

## Value

A trimesh object gate boundary
