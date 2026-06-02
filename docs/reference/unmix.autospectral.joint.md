# Joint AutoSpectral Unmixing (Pure R)

Pure-R implementation of the joint covariance-weighted AutoSpectral
pipeline. Mirrors the algorithm in `unmix_autospectral_joint_cpp`
without requiring `AutoSpectralRcpp`. For each cell, AF variant
selection and fluorophore spectral variant optimisation are performed in
a single jointly-scored pass rather than sequentially.

The AF selection score is a multiplicative product of the residual ratio
and a covariance-propagated fluorophore leakage ratio, matching Section
A of the C++ pipeline. Fluorophore variant selection uses the same
composite score across all optimisable fluorophores before any swaps are
committed, with conflict resolution by residual-delta cosine similarity
(Section C).

Parallelisation follows the same pattern as
[`unmix.autospectral()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.autospectral.md),
using
[`create.parallel.lapply()`](https://drcytometer.github.io/AutoSpectral/reference/create.parallel.lapply.md)
to handle platform differences between Windows (PSOCK cluster) and Linux
(forked `mclapply`).

## Usage

``` r
unmix.autospectral.joint(
  raw.data,
  spectra,
  af.spectra,
  asp,
  spectra.variants = NULL,
  n.passes = 2L,
  n.variants = 1L,
  parallel = TRUE,
  threads = NULL,
  verbose = TRUE
)
```

## Arguments

- raw.data:

  Numeric matrix (cells x detectors). Columns must match those of
  `spectra`.

- spectra:

  Numeric matrix (fluorophores x detectors, normalised 0-1). Must not
  contain an `"AF"` row.

- af.spectra:

  Numeric matrix (n_af x detectors, normalised 0-1). Prepare using
  [`get.af.spectra()`](https://drcytometer.github.io/AutoSpectral/reference/get.af.spectra.md).
  At least two rows are required.

- asp:

  The AutoSpectral parameter list, prepared using
  [`get.autospectral.param()`](https://drcytometer.github.io/AutoSpectral/reference/get.autospectral.param.md).
  Required for parallel backend setup.

- spectra.variants:

  Named list as returned by
  [`get.spectral.variants()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectral.variants.md),
  containing at minimum `$variants`, `$delta.list`, and `$thresholds`.
  The optional `$optimize.recommended` slot is respected: fluorophores
  flagged `FALSE` are excluded from variant optimisation. Pass `NULL`
  for AF-only mode (no variant optimisation).

- n.passes:

  Integer, default `2L`. Number of joint optimisation passes per cell.
  Matches `n_passes` in `unmix_autospectral_joint_cpp`. Higher values
  allow more variant swaps to be committed per cell at the cost of
  additional computation.

- n.variants:

  Integer, default `1L`. Maximum number of top-scoring variants to
  evaluate per fluorophore per pass. Corresponds to `k_opt` in the
  legacy pipeline. Set via `speed` in the calling function, or supply
  directly.

- parallel:

  Logical, default `TRUE`. Whether to use parallel processing across
  cells. Uses
  [`create.parallel.lapply()`](https://drcytometer.github.io/AutoSpectral/reference/create.parallel.lapply.md)
  to handle platform differences.

- threads:

  Numeric, default `NULL`. Number of worker threads. When `NULL`,
  `asp$worker.process.n` is used. When `threads = 0`,
  [`parallelly::availableCores()`](https://parallelly.futureverse.org/reference/availableCores.html)
  is used. Ignored when `parallel = FALSE`.

- verbose:

  Logical, default `TRUE`.

## Value

Numeric matrix (cells x (n_fluorophores + 2)) with column names matching
the row names of `spectra` plus `"AF"` (abundance of the selected AF
component) and `"AF Index"` (1-based index into `af.spectra` of the
selected AF spectrum).

## See also

[`unmix.autospectral()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.autospectral.md),
[`unmix.autospectral.rcpp()`](https://rdrr.io/pkg/AutoSpectralRcpp/man/unmix.autospectral.rcpp.html),
[`get.spectral.variants()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectral.variants.md),
[`create.parallel.lapply()`](https://drcytometer.github.io/AutoSpectral/reference/create.parallel.lapply.md)
