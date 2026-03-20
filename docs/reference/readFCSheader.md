# Read FCS Header

A light weight FCS file header reader based on `flowstate`, trying to
mimic `flowCore`.

## Usage

``` r
readFCSheader(fcs.path, keyword = NULL)
```

## Arguments

- fcs.path:

  A character string specifying the file path (directory and file name)
  for the .fcs file to be read.

- keyword:

  Optional argument specifying which keyword(s) to return if the whole
  set is not desired.

## Value

A named list of keywords (metadata)

## References

Laniewski, Nathan. *flowstate*.
<https://github.com/nlaniewski/flowstate>

Ellis B, Haaland P, Hahne F, Le Meur N, Gopalakrishnan N, Spidlen J,
Jiang M, Finak G (2025). *flowCore: Basic structures for flow cytometry
data*.
[doi:10.18129/B9.bioc.flowCore](https://doi.org/10.18129/B9.bioc.flowCore)
<https://bioconductor.org/packages/flowCore>

## See also

[`read.FCSheader`](https://rdrr.io/pkg/flowCore/man/read.FCSheader.html)
