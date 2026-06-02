# Calculate Optimization Necessity Scores for Spectral Variants

Determines, for each fluorophore with computed spectral variants,
whether per-cell spectral optimisation is likely to produce a meaningful
improvement in unmixing quality. Fluorophores whose spectral variation
does not overlap with the detector channels occupied by other panel
members cannot improve unmixing of those channels regardless of how many
variants are tested, so optimising them wastes computation without
benefit.

The score is computed in the same geometric space used by the C++ joint
pipeline (`unmix_autospectral_joint_cpp`):

1.  For each fluorophore `fl`, the pseudoinverse `U_nof` of the
    remaining fluorophores' spectra is computed. This is the same matrix
    used inside `FluorPrecomp` to build `pc.w_leakage`.

2.  The empirical covariance of `fl`'s delta spectra (variant spread) is
    propagated through `U_nof`:
    `leakage_cov = U_nof %*% delta_cov %*% t(U_nof)`. The diagonal of
    `leakage_cov` gives the variance of the unmixing error induced in
    each other fluorophore channel by `fl`'s spectral uncertainty.

3.  The score for `fl` is the sum of the standard deviations of that
    leakage: `score(fl) = sum( sqrt( abs( diag( leakage_cov ) ) ) )`.
    Scores are normalised to `[0, 1]` relative to the highest-scoring
    fluorophore in the panel.

Optionally, if a representative stained sample is provided (as a matrix
of per-fluorophore MFIs, typically the median positive signal), scores
are additionally weighted by `mu[fl]`, the fluorophore's own brightness.
This down-weights fluorophores that are geometrically capable of
cross-channel leakage but are always dim in the actual experiment being
unmixed. The MFI weighting is applied *after* normalisation so that the
geometric score remains interpretable on its own.

## Usage

``` r
calculate.optimize.necessity(
  spectra,
  delta.list,
  mu = NULL,
  threshold = 0.01,
  ridge = 1e-04,
  verbose = TRUE
)
```

## Arguments

- spectra:

  Numeric matrix of spectral signatures (fluorophores x detectors,
  values normalised 0-1). Must not include an `"AF"` row.

- delta.list:

  Named list of delta matrices, one per fluorophore. Each matrix is
  (n_variants x n_detectors) and contains the element-wise difference
  between each variant spectrum and the base spectrum for that
  fluorophore. Produced by
  [`get.spectral.variants()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectral.variants.md)
  as `spectra.variants$delta.list`.

- mu:

  Named numeric vector of per-fluorophore MFI values from a
  representative stained sample (one value per fluorophore in
  `spectra`). Values should be positive; negative or zero values are
  clamped to zero. Supply `NULL` (the default) to use the geometric
  score alone without brightness weighting.

- threshold:

  Numeric scalar in `[0, 1]`, default `0.01`. Normalised score below
  which optimisation is considered unnecessary and the fluorophore is
  excluded from the `optimize.recommended` output. A value of `0.01`
  means fluorophores scoring below 1% of the top-scoring fluorophore are
  skipped.

- ridge:

  Numeric scalar, default `1e-4`. Ridge regularisation added to the
  diagonal of `delta_cov` before propagation, matching the treatment
  inside the C++ precomputation. Prevents degenerate covariance
  estimates when a fluorophore has very few variants.

- verbose:

  Logical, default `TRUE`. When `TRUE`, prints a table of normalised
  scores and the resulting recommendation for each fluorophore.

## Value

A named list with three elements:

- `scores.raw`:

  Named numeric vector of raw (unnormalised) leakage propagation scores,
  one per fluorophore in `delta.list`.

- `scores.norm`:

  Named numeric vector of scores normalised to `[0, 1]`.

- `optimize.recommended`:

  Named logical vector: `TRUE` if the fluorophore's normalised score
  meets or exceeds `threshold`, `FALSE` otherwise. Fluorophores absent
  from `delta.list` are not included.

## See also

[`get.spectral.variants()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectral.variants.md),
[`sanitize.optimization.inputs()`](https://drcytometer.github.io/AutoSpectral/reference/sanitize.optimization.inputs.md)
