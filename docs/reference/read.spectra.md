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
  af.param = "AF"
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

## Value

A matrix containing the fluorophore spectra (fluorophore x detectors).
