# dot Cosine Similarity Rows

Cosine similarity of each row of mat against a single reference vector.

## Usage

``` r
.cosine.sim.rows(mat, ref.vec)
```

## Arguments

- mat:

  The matrix (or dataframe), represented as events in rows and detectors
  in columns.

- ref.vec:

  The vector (numeric) of the reference spectrum to which the rows of
  the matrix will be compared.

## Value

Returns a numeric vector of length nrow(mat).
