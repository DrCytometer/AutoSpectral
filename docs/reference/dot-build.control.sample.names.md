# Build Unique Control Sample Names

Derives a unique `sample` identifier for every row of a control table,
allowing multiple single-stained controls to be supplied for the same
fluorophore. Rows whose `fluorophore` is not duplicated are returned
unchanged. For a duplicated fluorophore, uniqueness is resolved in
stages, each applied only to the rows still in conflict after the
previous stage:

1.  `"fluorophore (control.type)"`

2.  `"fluorophore (control.type) (marker)"`

3.  `"fluorophore (control.type) (marker)_1"`, `"..._2"`, etc.

## Usage

``` r
.build.control.sample.names(fluorophore, control.type, marker)
```

## Arguments

- fluorophore:

  Character vector of fluorophore names, one per control.

- control.type:

  Character vector of control types (`"cells"`/`"beads"`), one per
  control.

- marker:

  Character vector of marker/antigen names, one per control. May contain
  `NA`.

## Value

A character vector the same length as `fluorophore`, guaranteed to
contain no duplicates.
