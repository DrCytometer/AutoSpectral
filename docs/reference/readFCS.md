# Read FCS File

A light weight FCS file binary reader based on `flowstate`. Lower memory
usage and faster than `flowCore`.

## Usage

``` r
readFCS(fcs.path, return.keywords = FALSE, start.row = NULL, end.row = NULL)
```

## Arguments

- fcs.path:

  A character string specifying the file path (directory and file name)
  for the .fcs file to be read.

- return.keywords:

  Logical, default `FALSE`. Controls whether keywords are returned as
  well as the expression data matrix (useful for writing a new file).

- start.row:

  Optional numeric specifying the row to begin reading on. Can be useful
  for reading in just the metadata or for chunking files. Default is
  `NULL`, which will read the whole file.

- end.row:

  Optional numeric specifying the row to end reading on. Can be useful
  for reading in just the metadata or for chunking files. Default is
  `NULL`, which will read the whole file.

## Value

If `return.keywords = TRUE`, a list containing two elements:

1.  the expression data in a matrix, and 2) the keywords. If
    `return.keywords = FALSE`, only the expression data is returned as a
    matrix.

## References

Laniewski, Nathan. *flowstate*.
<https://github.com/nlaniewski/flowstate>

Ellis B, Haaland P, Hahne F, Le Meur N, Gopalakrishnan N, Spidlen J,
Jiang M, Finak G (2025). *flowCore: Basic structures for flow cytometry
data*.
[doi:10.18129/B9.bioc.flowCore](https://doi.org/10.18129/B9.bioc.flowCore)

## See also

[`read.FCS`](https://rdrr.io/pkg/flowCore/man/read.FCS.html)
