# Unmix Autofluorescence by GLS with Per-Node Covariance

Selects an autofluorescence dictionary entry and abundance per cell by
maximising a Gaussian log-likelihood in which the AF spectrum itself is
uncertain: \$\$\Sigma_i = diag(read.var) + diag(\mu_i)/\kappa +
\alpha_i^2 \Sigma_k^{AF} + \sum_f x\_{if}^2 \Sigma_f^{var}\$\$

The AF covariance term softens the discreteness of the dictionary
without adding per-cell degrees of freedom: a cell whose true AF falls
between two nodes is no longer forced to attribute the mismatch to
photon noise. This is deliberately distinct from a continuous AF model,
which is known to inflate variance past the Cramer-Rao bound.

Applies the same low-rank Woodbury solve as
[`unmix.gls()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.gls.md)
(`.sigma.solve()`) rather than a dense Cholesky factorisation of the
full detector-space covariance. `af.model$nodes[[k]]` is already in the
same `basis`/`lambda` shape
[`build.variant.basis()`](https://drcytometer.github.io/AutoSpectral/reference/build.variant.basis.md)
produces, so the AF term for the candidate under test is folded into the
same low-rank accumulator as fluorophore variant uncertainty – no D x D
matrix is ever formed. Cost per candidate per iteration is
`O(D k^2 + k^3)` with `k` the active low-rank total (typically well
under 20), rather than `O(D^3)`.

## Usage

``` r
unmix.af.gls(
  raw.data,
  af.spectra,
  af.model = attr(af.spectra, "af.model"),
  spectra = NULL,
  variant.basis = NULL,
  noise.model = NULL,
  use.af.covariance = TRUE,
  use.node.prior = TRUE,
  use.abundance.prior = FALSE,
  include.spillover = FALSE,
  spillover.kappa = NULL,
  gain.cv = 0,
  active.threshold = 0.001,
  n.candidates = 5L,
  n.iter = 2L,
  ridge = 1e-08,
  inner.ridge = 1e-10,
  return.af.spectra = FALSE,
  bend.max.p = 1e-06,
  verbose = TRUE
)
```

## Arguments

- raw.data:

  Numeric matrix (events x detectors).

- af.spectra:

  AF dictionary (nodes x detectors), from
  `get.af.spectra(return.model = TRUE)`.

- af.model:

  Per-node model. Defaults to the `"af.model"` attribute of
  `af.spectra`. When supplied explicitly, must have one entry per row of
  `af.spectra`, in the same order.

- spectra:

  Optional fluorophore spectra (fluorophores x detectors). When
  supplied, fluorophores are fitted jointly with AF.

- variant.basis:

  Optional list from
  [`build.variant.basis()`](https://drcytometer.github.io/AutoSpectral/reference/build.variant.basis.md),
  applied to the fluorophore terms exactly as in
  [`unmix.gls()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.gls.md).
  Ignored unless `spectra` is supplied. Default `NULL`.

- noise.model:

  Optional list with `read.var` and `counts.per.unit`. When `NULL`, a
  flat diagonal is used and only the AF covariance shapes the fit.

- use.af.covariance:

  Logical, include the `alpha^2 Sigma_k^AF` term. Set `FALSE` for the
  A/B comparison. Default `TRUE`.

- use.node.prior:

  Logical, add `log(prior_k)` to the score. Default `TRUE`.

- use.abundance.prior:

  Logical, add the per-node lognormal prior on abundance. Default
  `FALSE`.

- include.spillover:

  Logical, include the multinomial spillover term for active
  fluorophores and the AF row itself. Default `FALSE` – see
  [`unmix.gls()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.gls.md)
  for why this double-counts shot noise under the usual `kappa`
  estimation route.

- spillover.kappa:

  Optional scalar spillover `kappa`. Defaults to
  `noise.model$kappa.pooled`.

- gain.cv:

  Numeric, multiplicative-gain coefficient of variation; see
  [`unmix.gls()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.gls.md).
  Default `0`.

- active.threshold:

  Numeric, fluorophores below this fraction of the cell's largest
  abundance are excluded from the low-rank terms. Default `1e-3`.

- n.candidates:

  Integer, number of dictionary entries scored per cell, shortlisted by
  the joint covariance-weighted score. `Inf` scores all. Default `5`.

- n.iter:

  Integer, GLS refinement iterations per candidate. Default `2`.

- ridge:

  Numeric ridge added to the `XtX` solve. Default `1e-8`.

- inner.ridge:

  Numeric ridge added inside the Woodbury inner matrix. Default `1e-10`.

- return.af.spectra:

  Logical, return the per-cell AF spectrum refined within the winning
  node's covariance basis (the best linear unbiased predictor of the
  shape deviation), rather than only the dictionary entry. Adds no free
  parameters: the deviation is determined by `Sigma` and the node's
  `lambda`, which the likelihood already uses. Costs an events x
  detectors matrix. Default `FALSE`.

- bend.max.p:

  Numeric in `(0, 1]`. The per-cell AF refinement is applied only to
  cells whose chi-square fit p-value falls below this, i.e. cells whose
  residual is inconsistent with the fitted noise model and therefore
  carry autofluorescence the dictionary cannot represent. Cells that fit
  are left on their discrete dictionary entry. The refinement can absorb
  genuine dye signal, so restricting it to misfit cells confines that
  exposure to the small population that gains from it. Requires a fitted
  `noise.model` to be meaningful, since the p-value is only calibrated
  against a real noise model. `1` bends every cell, which is not
  recommended: an ungated refinement absorbs dim fluorophore signal
  across the whole sample, while gating it to misfit cells retains the
  autofluorescence gain and removes the cost. Default `1e-6`.

- verbose:

  Logical. Default `TRUE`.

## Value

A list with `af.index`, `af.spectrum.name`, `af.abundance`, `loglik`,
`chisq`, `df`, `logdet`, `p`, `unmixed` (when `spectra` supplied), `se`,
and `af.spectrum` (when `return.af.spectra = TRUE`). `af.spectrum` rows
carry the same scale convention as the dictionary rows they refine, so
the reconstructed AF contribution for a cell is
`af.abundance * af.spectrum`; the refined rows are not re-normalised,
since re-normalising would require rescaling `af.abundance` to match.
