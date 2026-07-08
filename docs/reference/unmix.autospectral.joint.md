# Joint AutoSpectral Unmixing (Pure R)

Pure-R implementation of the joint covariance-weighted AutoSpectral
pipeline. Mirrors `unmix_autospectral_joint_cpp` (in
`unmix_autospectral_joint_pipeline.cpp`) instruction-for-instruction,
without requiring `AutoSpectralRcpp`.

The pipeline is staged identically to the C++ version:

1.  Global weighted pseudoinverse + AF helper pre-computation (Section
    1).

2.  Per-fluorophore variant pre-computation, including leakage weights
    and rank-1 residual-update helpers (Section 2), plus a structural
    collinearity table used later for joint-pair conflict retries
    (Section 2B).

3.  Per-cell AF selection, with optional multi-pass AF refinement
    against the residual (`n.af.passes`).

4.  Per-cell fluorophore solve against the AF-subtracted residual (with
    optional per-cell Poisson-style weighting, `cell.weight`), followed
    by joint variant selection: candidates are scored with an
    alpha-weighted composite of residual ratio and covariance-propagated
    leakage ratio, conflicting swaps are resolved by residual-delta
    cosine similarity, and swaps flagged as a structurally collinear
    pair are queued for a retry once the winning partner's own swap has
    been applied. Each commit is individually verified against the
    currently accepted RSS and reverted if it does not improve it.

## Usage

``` r
unmix.autospectral.joint(
  raw.data,
  spectra,
  af.spectra,
  asp,
  spectra.variants = NULL,
  n.passes = 1L,
  parallel = TRUE,
  threads = NULL,
  cell.weight = FALSE,
  noise.floor = NULL,
  alpha = 0.5,
  collinear.thresh = 0.5,
  joint.pair.resolution = TRUE,
  n.af.passes = 1L,
  refine.af.quantile = 0.5,
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

  Integer, default `1L`. Number of joint optimisation passes per cell.
  Matches `n_passes` in `unmix_autospectral_joint_cpp`.

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

- cell.weight:

  Logical, default `FALSE`. Enables per-cell Poisson-style detector
  weighting (`1 / max(|y_hat|, noise.floor)`), matching `cell_weight` in
  the C++ pipeline. When `FALSE`, weighting collapses to the (all-ones,
  unless a global weight was requested) `w.global` used in Section 1.

- noise.floor:

  Numeric scalar or length-D vector, default `NULL`. Per-detector noise
  floor used by the weighting scheme; falls back to a constant `125.0`
  per detector when `NULL`, matching the C++ default. Only used when
  `cell.weight = TRUE`.

- alpha:

  Numeric in `[0, 1]`, default `0.5`. Weighting exponent for the joint
  candidate score: `resid.ratio^alpha * leakage.ratio^(1 - alpha)`.

- collinear.thresh:

  Numeric, default `0.5`. Cosine-similarity threshold (in the
  pseudoinverse row space) above which two optimisable fluorophores are
  flagged as structurally collinear for joint-pair conflict retries.

- joint.pair.resolution:

  Logical, default `TRUE`. When `TRUE`, a candidate swap that conflicts
  with an already-committed swap belonging to a structurally collinear
  partner is queued and retried once that partner's own swap has been
  applied (rather than simply dropped).

- n.af.passes:

  Integer, default `1L`. Number of AF selection passes. When `> 1`,
  cells with an initial AF abundance at or above the
  `refine.af.quantile` quantile are re-scored against their residual in
  subsequent passes, accumulating additional AF abundance where the
  refinement score is \< 1.0. The AF *index* used in the output is
  always the one selected on the first pass.

- refine.af.quantile:

  Numeric in `[0, 1]`, default `0.5`. Quantile (type 7) of the
  first-pass AF abundance used to select which cells are eligible for AF
  refinement passes.

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
