# Rank AF Spectra by Spectral Angle

Reorders autofluorescence spectra by their spectral angle to a reference
(the mean AF spectrum), smallest angle first, and renames the rows
`AF1`, `AF2`, ... in the new order. The spectral angle is
`acos( cosine similarity )`, so the ordering is identical to ranking by
decreasing cosine similarity. Ties keep their original order. Spectra
with zero norm have an undefined angle and are placed last. Intended for
internal use by `get.af.spectra`.

## Usage

``` r
.rank.af.by.angle(af.spectra, reference.spectrum)
```

## Arguments

- af.spectra:

  Numeric matrix of AF spectra, spectra in rows, detectors in columns.

- reference.spectrum:

  Numeric vector, the reference (mean AF) spectrum, with one value per
  column of `af.spectra`.

## Value

`af.spectra` with rows reordered and renamed.
