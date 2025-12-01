# Calculate Condition Number

Calculates the mixing matrix condition number ("Complexity IndexTM") for
a set of fluorophores in a spectral panel.

## Usage

``` r
calculate.condition.number(spectra)
```

## Arguments

- spectra:

  A matrix of spectral signatures of fluorophores, normalized between 0
  and 1. Fluorophores should be in rows and detectors (channels) in
  columns. More than one spectrum is required to perform the
  calculation.

## Value

A numeric value representing the condition number for the spectral flow
cytometry panel.
