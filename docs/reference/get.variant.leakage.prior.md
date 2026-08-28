# Variant Leakage Prior

Propagates each fluorophore's measured spectral variability into unmixed
space, giving the covariance of the spillover a plausible reference
error could produce.

A reference spectrum is the centre of a distribution of possible
spectra, and the variants from
[`get.spectral.variants()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectral.variants.md)
measure that distribution. If the reference sits off-centre, the
resulting error must lie somewhere within it - a dye cannot be wrong in
a direction it never varies. Pushing the variant covariance through the
unmixing matrix converts "which ways can this spectrum move" into "which
channels can it leak into, and by how much".

For a cell carrying only fluorophore \\f\\ at abundance \\a\\, whose
true row is \\s_f + e_f\\, the unmixed values are \\a(\mathbf{1}\_f +
e_f U^\top)\\, so the induced spillover row is \\e_f U^\top\\ and its
covariance is \\U \Sigma_e U^\top\\. This is the same propagation used
for variant selection in the joint pipeline, but over the full design
rather than the design with the fluorophore removed, because here the
quantity of interest is the spillover actually induced under the
production unmix.

The result is a prior, not a constraint. It says where an error can
plausibly land, not which way it points; direction has to come from the
data. The ridge matters for the same reason: variants are measured on
controls, while the error being corrected is a control-versus-sample
difference, and the variant span covers only part of it. Without a floor
the prior would forbid corrections that genuinely need making.

## Usage

``` r
get.variant.leakage.prior(
  spectra,
  variants,
  extra.rows = NULL,
  af.name = "AF",
  span.fraction = 0.6,
  min.variants = 3L,
  verbose = TRUE
)
```

## Arguments

- spectra:

  Numeric matrix (fluorophores x detectors), the reference spectra used
  for unmixing.

- variants:

  The list returned by
  [`get.spectral.variants()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectral.variants.md),
  or a named list of delta matrices (variants x detectors) directly.

- extra.rows:

  Optional numeric matrix (rows x detectors) of background rows that
  form part of the production design, such as an autofluorescence basis
  from
  [`get.af.basis()`](https://drcytometer.github.io/AutoSpectral/reference/get.af.basis.md).
  Included in the unmixing matrix so the propagation matches how data
  are actually unmixed. Default `NULL`.

- af.name:

  Character or `NULL`, an autofluorescence row in `spectra` to exclude
  from the fluorophore set. Default `"AF"`.

- span.fraction:

  Numeric in (0, 1\], the fraction of real spectral error the variant
  family is expected to span. The ridge is sized from this, so the
  smaller it is the more freedom a correction has to lie outside the
  observed variant directions. Default `0.6`.

- min.variants:

  Integer, the fewest variants a fluorophore needs before its covariance
  is used. Fluorophores below this fall back to the panel median prior.
  Default `3`.

- verbose:

  Logical, controls messaging. Default `TRUE`.

## Value

A named list:

- `covariance`:

  Named list, one entry per fluorophore, each a fluorophores x
  fluorophores prior covariance of that row's induced spillover.

- `variance`:

  Fluorophores x fluorophores matrix of the diagonals, so
  `variance[f, c]` is the prior variance of the spillover from `f` into
  `c`. This is the quantity a per-coefficient trust weight needs.

- `sd`:

  Its square root, on the same scale as a spillover coefficient and
  suitable for plotting as an error propagation heatmap.

- `fallback`:

  Character vector of fluorophores given the median prior because they
  had too few variants.
