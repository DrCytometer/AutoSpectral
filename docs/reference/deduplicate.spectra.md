# Deduplicate Spectra by Cosine Similarity

Iterates through rows of a spectral matrix in order, retaining a row
only if its cosine similarity to every already-retained row is strictly
below `threshold`. Intended for internal use by `get.af.spectra`.

## Usage

``` r
deduplicate.spectra(spectra, threshold = 0.99)
```

## Arguments

- spectra:

  Numeric matrix, spectra in rows and detectors in columns. Assumed to
  be L-infinity normalised.

- threshold:

  Numeric scalar in (0, 1\]. Rows at or above this similarity to any
  retained row are dropped. Default `0.99`.

## Value

Numeric matrix with redundant rows removed.
