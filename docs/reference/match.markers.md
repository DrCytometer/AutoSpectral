# Match Markers

This function matches control filenames to markers in the marker
database, including synonyms, and returns the matched markers

## Usage

``` r
match.markers(control.filenames, marker.database)
```

## Arguments

- control.filenames:

  Vector of control filenames.

- marker.database:

  Data frame containing marker information.

## Value

A named vector of matched markers for each control filename.
