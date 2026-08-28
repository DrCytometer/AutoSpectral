# Spread-Scaled Positivity Thresholds

Converts a flat, per-fluorophore positivity threshold into a per-event
(or per-cluster) threshold that widens with the spillover spread
contributed by every fluorophore present in that event.

A flat cut taken from an unstained control describes the width of the
negative population only. Once a dye is bright, its photon and unmixing
noise inflates the estimate in every channel it spills into, so a bright
event can clear a flat cut in a channel it carries no signal in at all.
The Spillover Spreading Matrix from
[`get.spectral.variants()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectral.variants.md)
measures exactly this: entry `[a, b]` is the *variance* added to channel
`b` per unit of fluorophore `a`'s own on-channel abundance. The boundary
therefore grows as the square root of abundance, not in proportion to
it:

\$\$t\_{c,b} = m \cdot t_b + \kappa \sqrt{ \sum_a S\\S\_{a,b}
\max(x\_{c,a}, 0) }\$\$

Contributions are summed over every source, since in a stained sample
several dyes spread into the same channel at once.

## Usage

``` r
get.spread.thresholds(
  unmixed,
  thresholds,
  spillover.spread = NULL,
  spread.kappa = 2,
  margin = 1,
  verbose = TRUE
)
```

## Arguments

- unmixed:

  Numeric matrix (events or clusters x fluorophores), unmixed
  abundances. Column names must be the fluorophore names.

- thresholds:

  Named numeric vector of flat positivity thresholds in unmixed space,
  covering every column of `unmixed` (typically the 99.5th percentile of
  an unstained control, or `get.spectral.variants()$thresholds`).

- spillover.spread:

  Matrix (source fluorophore x target channel), the Spillover Spreading
  Matrix from `get.spectral.variants()$spillover.spread`. `NULL`
  (default) reproduces a flat threshold, broadcast to matrix shape.

- spread.kappa:

  Numeric, how many spread standard deviations to allow above the flat
  threshold. Default `2`.

- margin:

  Numeric multiplier applied to the flat component only, matching the
  `unstained.margin` convention in
  [`fix.my.unmix()`](https://drcytometer.github.io/AutoSpectral/reference/fix.my.unmix.md).
  Default `1`.

- verbose:

  Logical, controls messaging. Default `TRUE`.

## Value

Numeric matrix with the same dimensions and dimnames as `unmixed`,
giving the positivity threshold for each row and fluorophore.
