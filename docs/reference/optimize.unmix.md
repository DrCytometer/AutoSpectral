# Optimize Spectral Unmixing

Parallel backend for per-cell spectral optimization in R.

## Usage

``` r
optimize.unmix(
  raw.data,
  unmixed,
  spectra,
  pos.thresholds,
  optimize.fluors,
  variants,
  delta.list,
  delta.norms,
  fluorophores,
  asp,
  k = 10L,
  nthreads = 1L,
  parallel = TRUE
)
```

## Arguments

- raw.data:

  Numeric matrix (cells x detectors)

- unmixed:

  Numeric matrix (cells x fluors)

- spectra:

  Numeric matrix (fluors x detectors)

- pos.thresholds:

  Numeric vector (n fluors)

- optimize.fluors:

  Character vector of fluorophores present in variants

- variants:

  List of variant matrices per fluorophore

- delta.list:

  List of delta matrices per fluorophore

- delta.norms:

  List of delta norms per fluorophore

- fluorophores:

  Character vector of fluorophore names

- asp:

  The AutoSpectral parameter list.

- k:

  Integer, number of variants to test

- nthreads:

  Integer, number of threads

- parallel:

  Logical, whether to use parallel processing

## Value

Unmixed data with cells in rows and fluorophores in columns.
