# Simulate Spectral Flow Cytometry Data

Generates synthetic spectral flow cytometry data for development and
validation of AutoSpectral pipelines. Cells are assigned random
fluorophore expression levels drawn from three discrete expression
layers (high, medium, low), mixed into detector space via provided
spectra (and optionally spectral variants), and corrupted by
configurable noise sources: binomial spillover sampling, Poisson shot
noise, and detector readout noise. Autofluorescence variation is
supported via a provided `af.spectra` matrix.

Ground-truth abundances, per-cell variant assignments, per-cell AF
spectra, and synthetic scatter are returned alongside the raw
detector-space matrix and a simple OLS unmixed result (no AF) for
immediate inspection.

## Usage

``` r
sim.flow.data(
  spectra,
  asp,
  n.cells = 10000L,
  complexity = 1/3,
  layer.centroids = c(0.5, 0.05, 0.005),
  layer.cv = 0.2,
  af.spectra = NULL,
  af.scale.range = c(0.001, 0.05),
  af.scatter.slope = 0.1,
  variants = NULL,
  scatter.data = NULL,
  fsc.mean.log = log(2e+05),
  fsc.sd.log = 0.4,
  ssc.mean.log = log(80000),
  ssc.sd.log = 0.5,
  fsc.ssc.cor = 0.6,
  shot.noise = TRUE,
  counts.per.unit = 0.5,
  spillover.noise = TRUE,
  detector.noise = TRUE,
  af.variation = TRUE,
  spectral.variation = TRUE,
  seed = 42L
)
```

## Arguments

- spectra:

  Numeric matrix of fluorophore spectral signatures, L-infinity
  normalised, fluorophores in rows and detectors in columns. Required.
  Defines both the panel and the detector layout.

- asp:

  The AutoSpectral parameter list from `get.autospectral.param`. Used to
  resolve `expr.data.max` and cytometer-specific detector noise presets.

- n.cells:

  Integer. Number of synthetic cells to generate. Default `10000`.

- complexity:

  Numeric in (0, 1\]. Fraction of fluorophores that are "on" (non-zero)
  per cell, drawn without replacement. Default `1/3`.

- layer.centroids:

  Numeric vector of length 3, giving the high, medium, and low
  expression centroids as fractions of `expr.data.max`. Default
  `c(0.5, 0.05, 0.005)`.

- layer.cv:

  Numeric scalar. Coefficient of variation (on the log scale) for
  expression within each layer. Default `0.2`.

- af.spectra:

  Numeric matrix of autofluorescence spectra, L-infinity normalised, AF
  components in rows and detectors in columns. As returned by
  `get.af.spectra`. If `NULL` (default), no AF is added.

- af.scale.range:

  Numeric vector of length 2. Range of the uniform distribution from
  which per-cell AF abundance is drawn, as a fraction of
  `expr.data.max`. Default `c(0.001, 0.05)`.

- af.scatter.slope:

  Numeric scalar. Controls how much SSC drives AF abundance. AF
  abundance is multiplied by
  `1 + af.scatter.slope * (ssc - median(ssc)) / mad(ssc)`. Default
  `0.1`.

- variants:

  List of spectral variant matrices as returned by
  `get.spectral.variants` or `get.fluor.variants`. Each element is named
  by fluorophore and contains a matrix of variant spectra (variants in
  rows, detectors in columns). If `NULL` (default), the base `spectra`
  rows are used for all cells.

- scatter.data:

  Optional numeric matrix of real scatter data (cells in rows, at least
  two columns for FSC and SSC). When provided, rows are sampled with
  replacement to generate synthetic scatter, bypassing the log-normal
  model. Default `NULL`.

- fsc.mean.log:

  Numeric. Mean of the log-normal FSC distribution (on the natural-log
  scale). Default `log(200000)`.

- fsc.sd.log:

  Numeric. SD of the log-normal FSC distribution. Default `0.4`.

- ssc.mean.log:

  Numeric. Mean of the log-normal SSC distribution. Default
  `log(80000)`.

- ssc.sd.log:

  Numeric. SD of the log-normal SSC distribution. Default `0.5`.

- fsc.ssc.cor:

  Numeric in (-1, 1). Pearson correlation between log(FSC) and log(SSC)
  in the synthetic scatter population. Default `0.6`.

- shot.noise:

  Logical. Whether to apply Poisson shot noise. Default `TRUE`.

- counts.per.unit:

  Numeric scalar for `shot.noise`. Default `0.5`. Works inversely;
  higher numbers will decrease shot noise-driven spread.

- spillover.noise:

  Logical. Whether to apply binomial spillover sampling. Default `TRUE`.

- detector.noise:

  Logical. Whether to apply detector readout noise (Gaussian,
  cytometer-specific). Default `TRUE`.

- af.variation:

  Logical. Whether to include autofluorescence variation (requires
  `af.spectra`). Default `TRUE`.

- spectral.variation:

  Logical. Whether to sample per-cell spectral variants (requires
  `variants`). Default `TRUE`.

- seed:

  Integer random seed for reproducibility. Default `42`.

## Value

A named list with elements:

- `raw`:

  Numeric matrix (cells x detectors) of simulated detector-space data,
  including all enabled noise sources.

- `unmixed.no.af`:

  Numeric matrix (cells x fluorophores) of OLS-unmixed data without AF,
  for immediate inspection.

- `truth`:

  Named list:

  `abundances`

  :   Matrix (cells x fluorophores+AF), true abundance per cell. AF
      column appended as the last column when `af.spectra` is provided.

  `expression.layer`

  :   Character matrix (cells x fluorophores), `"high"`, `"mid"`,
      `"low"`, or `NA` (off).

  `variant.index`

  :   Integer matrix (cells x fluorophores), index of the variant
      spectrum used per cell per fluorophore. `NA` when off or no
      variants provided.

  `af.row`

  :   Integer vector (length cells), which row of `af.spectra` was
      assigned to each cell. `NA` when no AF.

  `af.abundance`

  :   Numeric vector (length cells), true AF abundance per cell in
      instrument units. `0` when no AF.

  `scatter`

  :   Numeric matrix (cells x 2, columns FSC and SSC), synthetic or
      resampled scatter values.

- `params`:

  Named list of all simulation settings used, for reproducibility.
