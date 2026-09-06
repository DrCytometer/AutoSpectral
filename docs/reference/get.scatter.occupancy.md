# Get Scatter Occupancy

Determines the region of 2D scatter space (typically FSC vs SSC) that a
sample's events actually occupy, using a density-relative threshold
rather than a fixed instrument range. A 2D kernel density estimate is
built over a robust, quantile-trimmed span of the data, so that a
handful of extreme events cannot inflate the estimation grid and coarsen
its resolution over the region that actually matters. Events falling in
density bins below `density.threshold` of the peak density are treated
as sparse noise or outlier events and excluded; the remaining "occupied"
events define a data-driven working range for gating and plotting.

This identifies dense clusters generally, and does not by itself
distinguish a real debris population from the population of interest – a
distinct, well-populated debris cloud is typically just as dense as real
cells or beads and will not be excluded here. Separating debris from the
target population remains the job of the downstream target-maximum
search (e.g. the bottom-left exclusion logic in
[`do.gate()`](https://drcytometer.github.io/AutoSpectral/reference/do.gate.md));
this function only keeps sparse, spread-out outlier events from
distorting the working scale used to set that search up.

## Usage

``` r
get.scatter.occupancy(
  scatter.data,
  density.threshold = 0.02,
  grid.n = 128,
  max.events = 1e+05,
  trim.quantile = 0.001,
  bird.seed = NULL
)
```

## Arguments

- scatter.data:

  A matrix or data frame with (at least) two columns: Forward Scatter
  and Side Scatter, in that order.

- density.threshold:

  Numeric 0-1, default `0.02`. Bins with a kernel density below this
  fraction of the peak density are excluded as sparse.

- grid.n:

  Numeric, default `128`. Binning grid for the kernel density
  estimation.

- max.events:

  Numeric, default `1e5`. Maximum number of events used to fit the
  kernel density surface. Above this, a random subsample is used to fit
  the surface; every event is still evaluated against it.

- trim.quantile:

  Numeric 0-0.5, default `0.001`. Quantile trimmed from each end of the
  data range before setting the density grid limits, so a small number
  of extreme events cannot stretch the grid.

- bird.seed:

  Integer, seed for reproducible subsampling. Default `NULL` skips
  seeding.

## Value

A list with:

- `keep`: Logical vector, length `nrow(scatter.data)`, `TRUE` for events
  in an occupied (dense) region.

- `x.range`, `y.range`: Two-element numeric vectors giving the occupied
  range on each axis, taken from the retained events.

- `density`: The `x`/`y`/`z` kernel density surface used to determine
  occupancy (same structure as
  [`MASS::kde2d()`](https://rdrr.io/pkg/MASS/man/kde2d.html) output),
  useful for diagnostic plotting.

## See also

- [`do.gate()`](https://drcytometer.github.io/AutoSpectral/reference/do.gate.md)

- [`define.gate.density()`](https://drcytometer.github.io/AutoSpectral/reference/define.gate.density.md)

- [`define.gate.landmarks()`](https://drcytometer.github.io/AutoSpectral/reference/define.gate.landmarks.md)
