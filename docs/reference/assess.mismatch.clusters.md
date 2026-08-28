# Cluster-Based Permutation Test for Spectral Mismatch Regions

Identifies contiguous regions of the detector array where the mismatch
between two particle types' spectra (e.g. beads vs. cells) is larger
than expected by chance, while accounting for two structural features of
the data: (1) adjacent detectors are correlated because a fluorophore's
emission is a smooth curve, so per-detector tests cannot be treated as
independent, and (2) a fluorophore can only show mismatch at a detector
where it has appreciable signal in the first place, so raw per-detector
averages are confounded by how many/which fluorophores peak there.

Both issues are handled together: mismatch is only evaluated at
fluorophore/detector pairs where the fluorophore has on-peak signal
(`peak.threshold`), a per-detector statistic is computed only from those
on-peak observations, and significance is assessed with a sign-flip
permutation test using cluster mass (summed \|statistic\| over
contiguous detectors) as the test statistic. Because permutations flip
the sign of whole fluorophore traces rather than shuffling individual
detector values, the detector-to-detector correlation structure already
present in the data is preserved in the null distribution, so no
independence assumption across detectors is required.

## Usage

``` r
assess.mismatch.clusters(
  dist.mat,
  ref.spectra,
  peak.threshold = 0.1,
  min.n = 3L,
  cluster.threshold = 1.96,
  n.perm = 2000L,
  bird.seed = NULL
)
```

## Arguments

- dist.mat:

  Numeric matrix of signed per-detector differences (fluorophore x
  detector), e.g. `comparison[[b]]$Distance` from
  [`run.bead.cell.comparison()`](https://drcytometer.github.io/AutoSpectral/reference/run.bead.cell.comparison.md)
  or the output of
  [`bead.cell.dist()`](https://drcytometer.github.io/AutoSpectral/reference/bead.cell.dist.md).
  Column order is assumed to be true spectral/detector order (as in the
  package's per-cytometer reference libraries); contiguity for
  clustering is defined by that order, not by column name.

- ref.spectra:

  Numeric matrix (fluorophore x detector) giving the reference particle
  type's spectral intensity, used only to build the on-peak mask. Must
  share row and column names with `dist.mat` (extra rows/columns are
  dropped; a warning is issued if the shared set is smaller than either
  input).

- peak.threshold:

  Numeric scalar in `[0, 1]`. A fluorophore/detector pair is included in
  the test if
  `ref.spectra[f, d] >= peak.threshold * max(ref.spectra[f, ])`.
  Defaults to `0.1` (10\\ of that fluorophore's peak height).

- min.n:

  Integer scalar. Minimum number of on-peak fluorophores required at a
  detector for it to receive a test statistic; detectors with fewer
  contributing fluorophores are set to `NA` and cannot seed or extend a
  cluster. Defaults to `3L`.

- cluster.threshold:

  Numeric scalar. Primary (per-detector) threshold on `abs(t-statistic)`
  used to form candidate clusters, in the usual
  cluster-based-permutation sense: this sets what counts as "locally
  unusual" before cluster mass is compared across permutations. Defaults
  to `1.96` (approx. two-sided alpha = 0.05 for a single detector,
  before cluster-level correction).

- n.perm:

  Integer scalar. Number of sign-flip permutations used to build the
  null distribution of the maximum cluster mass. Defaults to `2000L`.

- bird.seed:

  Integer scalar or `NULL`. Seed for the permutation draws, for
  reproducibility. Defaults to `NULL` (no seeding).

## Value

A list with components:

- stat:

  Named numeric vector (length `ncol(dist.mat)`) of the observed
  per-detector t-statistic, `NA` where `min.n` was not met.

- mask:

  Logical matrix (fluorophore x detector), the on-peak mask used.

- n.detector:

  Named integer vector of the number of on-peak fluorophores
  contributing to each detector's statistic.

- clusters:

  Data frame, one row per observed cluster, with columns `start`, `end`
  (detector names), `n.detectors`, `mass`, and `p.value`. Empty (0-row)
  if no detector met `cluster.threshold`.

- null.max.mass:

  Numeric vector of length `n.perm`, the permutation null distribution
  of the maximum cluster mass.

## Details

At each detector `d`, let `x` be the on-peak entries of `dist.mat[, d]`
(those `f` with `mask[f, d] == TRUE`). The per-detector statistic is the
one-sample t-statistic `mean(x) / (sd(x) / sqrt(length(x)))`, testing
the null hypothesis that the reference and comparison particle types
produce the same spectrum at that detector, restricted to fluorophores
that actually fluoresce there. Detectors are then thresholded at
`cluster.threshold` and contiguous runs of exceedance are collapsed into
clusters, each summarized by its cluster mass (sum of `abs(t)` across
its member detectors).

The null distribution is built by, for each of `n.perm` iterations,
independently flipping the sign of every fluorophore's entire on-peak
row, recomputing the per-detector statistic and cluster masses under
that permutation, and retaining only the single largest cluster mass
observed. This max-statistic approach controls the family-wise error
rate across the whole detector array without assuming detectors are
independent.

Each observed cluster's p-value is the proportion of permutation draws
whose maximum cluster mass met or exceeded that cluster's mass.
