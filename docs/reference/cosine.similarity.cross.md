# Cross Cosine Similarity

Computes pairwise cosine similarity between rows of two matrices `a` and
`b`. Returns an (nrow(a) x nrow(b)) matrix. Intended for internal use by
`get.af.spectra`.

## Usage

``` r
cosine.similarity.cross(a, b)
```

## Arguments

- a:

  Numeric matrix.

- b:

  Numeric matrix with the same number of columns as `a`.

## Value

Numeric matrix of cosine similarities, shape (nrow(a), nrow(b)).
