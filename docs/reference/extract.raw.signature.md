# Extract Raw Signature

Estimates one fluorophore's spectral signature directly from background
subtracted raw detector data, in the same way a signature is extracted
from a single stained control: bin the population by that fluorophore's
abundance, remove the contribution of every other fluorophore active in
the population at its current spectrum, and take the no-intercept slope
of the remaining detector signal against abundance.

The estimate is a positively weighted combination of measured detector
means, so it is a spectrum in the sense that matters: it has a peak
channel, it is non-negative up to noise, and it can be compared with a
control-derived row by cosine similarity. This is the difference from
back-solving a signature out of a compensation matrix, which returns
whichever linear combination of the existing rows reproduces the
required correction and need not resemble a spectrum at all.

Nothing is decided here. The function returns the candidate row together
with every quantity an acceptance gate needs, and the caller chooses
whether to adopt it.

## Usage

``` r
extract.raw.signature(
  raw.data,
  spectra,
  abundance,
  target,
  active = rownames(spectra),
  extra.abundance = NULL,
  intercept = TRUE,
  multivariate = TRUE,
  ridge = 1e-06,
  n.levels = 60L,
  min.bin.events = 50L,
  min.events = 200L,
  background.raw = NULL
)
```

## Arguments

- raw.data:

  Numeric matrix (events x detectors), already restricted to the
  population for this fluorophore and already background subtracted.

- spectra:

  Numeric matrix (fluorophores x detectors), the current reference
  spectra. Rows other than `target` are held fixed.

- abundance:

  Numeric matrix (events x fluorophores), row-matched to `raw.data`, the
  current best abundance estimates for the same population.

- target:

  Character, the fluorophore whose row is being estimated.

- active:

  Character vector, `target` plus every fluorophore whose contribution
  is to be removed at its current spectrum. Defaults to the whole panel,
  which is the correct block-coordinate step: a fluorophore whose
  abundance in this population is noise adds variance but no bias,
  whereas one that is left in is projected onto the target's abundance
  and absorbed into its row.

- extra.abundance:

  Optional numeric matrix (events x components), row-matched to
  `abundance`, additional regressors fitted jointly with the panel,
  typically the unclamped autofluorescence component coefficients.
  Background variation that tracks the abundances is an omitted variable
  in the joint fit, and with collinear abundances its bias is amplified
  onto every partial slope; giving it its own columns removes the
  transfer. Requires `multivariate`. Default `NULL`.

- intercept:

  Logical, whether to fit an intercept. The intercept absorbs any
  constant offset the removal left behind, such as the mean level of
  other markers multiplied by their own row errors. Default `TRUE`.

- multivariate:

  Logical, whether every active abundance is fitted jointly rather than
  the nuisance fluorophores being subtracted at their current spectra
  beforehand. Subtracting first books any error in a nuisance row
  against the target in proportion to how closely the two abundances
  track each other within this population, which is what `vif.target`
  measures. The joint fit removes that transfer at the cost of variance,
  and requires more abundance bins than active fluorophores. Default
  `TRUE`.

- ridge:

  Numeric, ridge penalty applied to the scaled joint design, keeping the
  solve stable when two fluorophores are nearly collinear within the
  population. The intercept is never penalised. Default `1e-6`.

- n.levels:

  Integer, the maximum number of abundance bins. The out-of-span signal
  is carried by the covariance between the binned residual and the
  binned abundances, so the bin count is the fit's sample size and the
  joint design spends one point per regressor; the count actually used
  is the largest that keeps `min.bin.events` events per bin. Default
  `60`.

- min.bin.events:

  Integer, the fewest events an abundance bin may contain. A population
  that cannot supply more bins than regressors returns `NULL`. Default
  `50`.

- min.events:

  Integer, minimum events in the population. Default `200`.

- background.raw:

  Optional numeric matrix (events x detectors) of background events
  after the same subtraction, used to report how much background the
  subtraction failed to remove. Default `NULL`.

## Value

`NULL` when the population cannot support a fit, otherwise a named list:

- `signature`:

  The candidate row, non-negative and L-infinity normalised.

- `signature.raw`:

  The unclamped, unnormalised slope.

- `stats`:

  One-row data frame of fit and plausibility quantities: `x.span`,
  `explained`, `explained.total`, `resid.rel`, `intercept.rel`,
  `bg.align`, `clamp.frac`, `anchor.rel`, `vif.target`, `vif.max`,
  `vif.partner`, `deg.change`, `peak.curr`, `peak.new`, `peak.new.rel`.
  `bg.align` is the cosine between the fitted intercept and the fitted
  signature; a value near minus one means the fit is trading a
  background floor against the slope, the signature of an unmodelled
  background confound. `peak.new.rel` is the current signature's own
  value at the candidate's new peak channel, relative to the current
  signature's own peak - 1 when the candidate does not propose a shift,
  falling toward 0 the less that channel has to do with the dye today.
