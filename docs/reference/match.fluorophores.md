# Match Fluorophores

This function matches control filenames to fluorophores in the
fluorophore database, including synonyms, and returns the matched
fluorophores based on the longest character match to ensure specificity.

## Usage

``` r
match.fluorophores(control.filenames, fluorophore.database)
```

## Arguments

- control.filenames:

  Vector of control filenames.

- fluorophore.database:

  Data frame containing fluorophore information.

## Value

A named vector of matched fluorophores for each control filename.
