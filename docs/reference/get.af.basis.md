# Get Autofluorescence Basis

Builds a low-dimensional autofluorescence basis from an unstained
control as the leading right singular vectors of the raw, uncentred
unstained matrix. The leading component is therefore the mean background
direction and the remainder describe how background varies from cell to
cell.

A single autofluorescence spectrum cannot represent that variation. In a
real sample autofluorescence is a family of profiles that tracks cell
type, and cell type also determines which markers are expressed, so
whatever a single row leaves behind is correlated with marker expression
and is readily mistaken for spillover by any method that fits one
channel against another.

## Usage

``` r
get.af.basis(
  unstained,
  n.pc = "auto",
  max.pc = 6L,
  read.var = 125^2,
  n.permutations = 3L,
  max.events = 50000L,
  verbose = TRUE
)
```

## Arguments

- unstained:

  Numeric matrix (events x detectors), raw unstained control, already
  gated to the population of interest.

- n.pc:

  Integer or `"auto"`. Number of components to retain. `"auto"` keeps
  every component whose singular value exceeds the largest value a
  matrix of the same shape containing only read noise would produce.
  Default `"auto"`.

- max.pc:

  Integer, hard cap on the number of components. Default `6`.

- read.var:

  Numeric, per-detector read noise variance in raw units, the floor of
  the retention threshold when `n.pc = "auto"`. Default `125^2`.

- n.permutations:

  Integer, permutations used to build the retention threshold when
  `n.pc = "auto"`. Each permutation costs one decomposition. Default
  `3`.

- max.events:

  Integer, maximum events used for the decomposition. Default `50000`.

- verbose:

  Logical, controls messaging. Default `TRUE`.

## Value

Numeric matrix, components x detectors, with unit L2 rows and the
singular values attached as the `singular.values` attribute.
