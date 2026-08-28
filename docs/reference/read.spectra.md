# Read In Saved Spectra

Reads in CSV files created by `AutoSpectral` or in the same format
(fluorophores in rows, detectors in columns, first row is detector
names, first column contains fluorophore names).

## Usage

``` r
read.spectra(
  spectra.file,
  spectra.dir = "./table_spectra",
  remove.af = FALSE,
  af.param = "AF",
  check.collinearity = TRUE,
  collinearity.threshold = 0.95
)
```

## Arguments

- spectra.file:

  File name for the spectra CSV file to be read.

- spectra.dir:

  File path to the folder containing `spectra.file`. Default is
  `table_spectra`, where `AutoSpectral` saves the spectra files.

- remove.af:

  Logical, default is `FALSE`. If `TRUE`, returns the spectral matrix
  without the default autofluorescence spectrum.

- af.param:

  Name of the autofluorescence parameter. Default is `AF`. Note that any
  fluorophores can be removed from the matrix by supplying a character
  vector, e.g., `c("BUV395", "PE")`, if desired.

- check.collinearity:

  Logical, default `TRUE`. CSV files cannot carry the `"fluorophore"`
  identity attribute that
  [`check.spectra.duplicates()`](https://drcytometer.github.io/AutoSpectral/reference/check.spectra.duplicates.md)
  checks exactly elsewhere, and a `sample`-disambiguated rowname (e.g.
  `"PE (cells)"` / `"PE (cells) (CD4)"`) is unique by construction, so a
  rownames-only check would never catch two controls for the same dye.
  When `TRUE`, screens for suspiciously similar rows as an early,
  load-time warning (the same screen
  [`check.spectra.duplicates()`](https://drcytometer.github.io/AutoSpectral/reference/check.spectra.duplicates.md)
  falls back to later).

- collinearity.threshold:

  Numeric, default `0.95`. Cosine similarity above which two rows
  trigger the warning. Matches the package default for
  `asp$similarity.warning.n`.

## Value

A matrix containing the fluorophore spectra (fluorophore x detectors).
