# Unmix by Generalised Least Squares

Per-cell GLS unmixing under a composite noise model: \$\$\Sigma_i =
diag(read.var) + diag(\mu_i)/\kappa + \frac{1}{\bar\kappa}\sum_f x\_{if}
m_f \[diag(p_f) - p_f p_f^\top\] + \sum_f x\_{if}^2 \Sigma_f^{var} +
\alpha_i^2 \Sigma\_{k_i}^{AF}\$\$

where \\m_f = \sum_d S\_{fd}\\ is the spectral mass of fluorophore f.
The spillover (multinomial), spectral-variant and AF terms are low rank,
so \\\Sigma_i^{-1}\\ is applied via the Woodbury identity without ever
forming a D x D inverse. Abundances are clamped at zero **only** when
constructing \\\Sigma_i\\ (a variance cannot be negative); the returned
estimate is unconstrained, preserving continuity through zero.

The spillover term and shot noise (`diag(mu_i)/kappa`) both derive from
photon partitioning; if `kappa` was estimated on cell or bead data by
[`estimate.noise.model()`](https://drcytometer.github.io/AutoSpectral/reference/estimate.noise.model.md)
(the normal case), the fitted slope has already absorbed the spillover
diagonal in full, and `include.spillover = TRUE` double-counts it. Leave
`include.spillover = FALSE` unless `kappa` came from a source measured
free of spillover (e.g. `BDSPECTRAL QSPE`).

## Usage

``` r
unmix.gls(
  raw.data,
  spectra,
  noise.model,
  variant.basis = NULL,
  af.spectra = NULL,
  af.index = NULL,
  af.basis = NULL,
  n.iter = 2L,
  method = c("woodbury", "dense"),
  include.spillover = FALSE,
  spillover.kappa = NULL,
  gain.cv = 0,
  active.threshold = 0.001,
  ridge = 1e-08,
  inner.ridge = 1e-10,
  return.variance = FALSE,
  return.fit.stats = FALSE,
  verbose = TRUE
)
```

## Arguments

- raw.data:

  Numeric matrix (events x detectors).

- spectra:

  Numeric matrix (fluorophores x detectors), no `"AF"` row.

- noise.model:

  List from
  [`estimate.noise.model()`](https://drcytometer.github.io/AutoSpectral/reference/estimate.noise.model.md).

- variant.basis:

  Optional list from
  [`build.variant.basis()`](https://drcytometer.github.io/AutoSpectral/reference/build.variant.basis.md).

- af.spectra:

  Optional AF dictionary. When supplied with `af.index`, each cell's
  assigned AF spectrum is appended to `spectra` for that cell.

- af.index:

  Optional integer vector, length `nrow(raw.data)`, giving the row of
  `af.spectra` assigned to each cell.

- af.basis:

  Optional list, one entry per AF dictionary node, each containing
  `basis` (components x detectors) and `lambda` (numeric). Supplies the
  within-node AF covariance term \\\alpha_i^2 \Sigma\_{k_i}^{AF}\\
  (section 2.3 of `CONTEXT_GLS_unmixing.md`). This does not need to be
  built separately: pass `attr(af.spectra, "af.model")$nodes` directly
  when `af.spectra` was produced by
  `get.af.spectra(..., return.model = TRUE)`, which already computes
  exactly this structure, keyed by `rownames(af.spectra)`. Ignored
  unless `af.spectra` and `af.index` are also supplied. Default `NULL`,
  which omits the term entirely – AF mismatch then appears as unmodelled
  residual rather than covariance.

- n.iter:

  Integer, number of GLS refinement iterations. Default `2`. The
  covariance depends on the abundances, so one iteration is a plug-in
  estimate from OLS and two is usually enough.

- method:

  `"woodbury"` (default) or `"dense"`. `"dense"` forms \\\Sigma_i\\
  explicitly and is intended as a correctness reference.

- include.spillover:

  Logical, whether to include the multinomial spillover term. Default
  `FALSE` (see description).

- spillover.kappa:

  Optional scalar \\\bar\kappa\\ for the multinomial term. Defaults to
  `noise.model$kappa.pooled`. A scalar is used deliberately: photon
  partitioning across detectors is a property of the emission spectrum,
  not of per-detector gain.

- gain.cv:

  Numeric, coefficient of variation of multiplicative gain fluctuation
  (illumination/transit-time variation), applied as a rank-1 term
  \\\code{gain.cv}^2\\ \mu_i \mu_i^\top\\. Distinct from shot noise and
  spillover, which are photon-counting effects; `gain.cv` captures gain
  variation at fixed photon count. If `curvature.coef` from
  [`estimate.noise.model()`](https://drcytometer.github.io/AutoSpectral/reference/estimate.noise.model.md)
  comes back positive and similar in magnitude across detectors on a
  bead run, that is this term's expected signature, and
  `sqrt(curvature.coef)` at a representative detector is a reasonable
  starting estimate. Default `0` (off): unvalidated on real data as of
  this writing.

- active.threshold:

  Numeric. Fluorophores below this fraction of the cell's largest
  abundance are excluded from the low-rank terms, keeping the Woodbury
  rank well below the detector count. Default `1e-3`.

- ridge:

  Numeric ridge added to the `XtX` solve. Default `1e-8`.

- inner.ridge:

  Numeric ridge added inside `.sigma.solve()`'s Woodbury inner matrix
  `M`. Kept separate from `ridge` because the two solves see different
  conditioning: `M` is `k x k` in the (typically small) active rank,
  `XtX` is `F x F` in the number of fluorophores. Default `1e-10`.

- return.variance:

  Logical. When `TRUE`, also returns per-cell estimated standard errors,
  i.e. \\\sqrt{diag((S\Sigma^{-1}S^\top)^{-1})}\\. Default `FALSE`.

- return.fit.stats:

  Logical. When `TRUE`, also returns per-cell goodness-of-fit: `chisq`
  (\\r^\top\Sigma^{-1}r\\, using the plug-in \\\Sigma\\ from the final
  iteration, evaluated at the converged \\\hat x\\ – the same plug-in
  convention `se` already uses), `df` (detectors minus active
  fluorophores, including AF), `logdet` (\\\log\|\Sigma_i\|\\, via the
  matrix determinant lemma at no extra solve), and `p` (upper-tail
  chi-square p-value). A small `p` flags a cell whose residual is
  inconsistent with the fitted noise model – the direct test proposed in
  `CONTEXT_GLS_unmixing.md` section 5 as a replacement for
  [`calculate.optimize.necessity()`](https://drcytometer.github.io/AutoSpectral/reference/calculate.optimize.necessity.md).
  Default `FALSE`.

- verbose:

  Logical. Default `TRUE`.

## Value

A matrix (cells x fluorophores), or when `return.variance` or
`return.fit.stats` is `TRUE`, a list with `unmixed` and whichever of
`se`, `chisq`, `df`, `logdet`, `p` were requested.
