# Resolve File Name

Determines the file name to use for output naming and plot titles.
Prefers the name embedded in the FCS `$FIL` keyword, but falls back to
the actual file name on disk when the file has been renamed since
acquisition – a common occurrence that otherwise leaves output files and
plot titles referring to a name the user no longer recognizes.

## Usage

``` r
resolve.file.name(path, fil, verbose = TRUE)
```

## Arguments

- path:

  Character. Path to the FCS file as provided by the caller.

- fil:

  Character or `NULL`. The `$FIL` keyword value read from the file, if
  present.

- verbose:

  Logical, default `TRUE`. Whether to message when a mismatch is
  detected and resolved.

## Value

Character, the resolved file name (without directory path).
