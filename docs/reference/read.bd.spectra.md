# Extract Spectra From BD FCS File

Extracts the spectra (spillover) values from BD FCS files (e.g., A8, S8)

## Usage

``` r
read.bd.spectra(fcs.file)
```

## Arguments

- fcs.file:

  Path and filename to the FCS file to be read.

## Value

The spillover matrix as fluorophores x detectors, with fluorophore and
detector names applied.
