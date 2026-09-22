# Cluster Data with a Self-Organizing Map

Replaces [`FlowSOM::SOM()`](https://rdrr.io/pkg/FlowSOM/man/SOM.html) /
[`EmbedSOM::SOM()`](https://rdrr.io/pkg/EmbedSOM/man/SOM.html) as the
SOM clustering engine behind
[`get.af.spectra()`](https://drcytometer.github.io/AutoSpectral/reference/get.af.spectra.md)
and
[`get.fluor.variants()`](https://drcytometer.github.io/AutoSpectral/reference/get.fluor.variants.md).
Uses the OpenMP- accelerated batch SOM (`som_train_batch_cpp()`) from
AutoSpectralRcpp when available, falling back to
[`FlowSOM::SOM()`](https://rdrr.io/pkg/FlowSOM/man/SOM.html) (pure R, no
compiled dependency) otherwise – so a plain `AutoSpectral` install with
no `AutoSpectralRcpp` still works, just without the OpenMP speedup.

## Usage

``` r
get.som.codes(
  data,
  som.dim,
  rlen = 10L,
  radius = NULL,
  dist = 4L,
  seed = 1337L,
  threads = 0L,
  unit.norm = FALSE
)
```

## Arguments

- data:

  Numeric matrix, training events x features. Must have colnames.

- som.dim:

  Integer, side length of the square SOM grid.

- rlen:

  Integer, number of epochs (full passes over the data). Default `10`.

- radius:

  Length-2 numeric vector, start/end neighbourhood radius. Default
  `NULL`, which derives both ends from the grid's own neighbour
  distances (67th percentile down to a small non-zero floor). Only used
  on the AutoSpectralRcpp path –
  [`FlowSOM::SOM()`](https://rdrr.io/pkg/FlowSOM/man/SOM.html) manages
  its own radius schedule internally.

- dist:

  Integer 1:4, distance function (1 manhattan, 2 euclidean, 3 chebyshev,
  4 cosine). Default `4`.

- seed:

  Integer, RNG seed for the initial codebook sample. Callers should pass
  `asp$bird.seed`.

- threads:

  Integer, OpenMP threads for the accelerated path. Default `0` (all
  available cores). Ignored on the
  [`FlowSOM::SOM()`](https://rdrr.io/pkg/FlowSOM/man/SOM.html) fallback
  path, which is single-threaded.

- unit.norm:

  Logical, default `FALSE`. Scale every training row to unit L2 length
  before training. Appropriate with `dist = 4` (cosine), where the
  distortion-minimising code is the mean of unit event vectors rather
  than the mean of raw ones; without it, bright events dominate each
  node's code. Leaves cosine distances, and therefore node assignments,
  unchanged. Has no principled justification under `dist = 2` and is
  left off by default so existing callers are unaffected. Note that the
  returned codes are then on the unit-normalised scale.

## Value

A list with `codes` (matrix, SOM nodes x features), `grid`, and
`nNodes`, matching the subset of
[`FlowSOM::SOM()`](https://rdrr.io/pkg/FlowSOM/man/SOM.html)'s return
value actually consumed elsewhere in AutoSpectral (`map$codes`).
