# Write FCS File

A light weight FCS file binary writer based on `flowstate`. Lower memory
usage and faster than `flowCore`.

## Usage

``` r
writeFCS(mat, keys, file.name, output.dir)
```

## Arguments

- mat:

  The expression data, containing both unmixed data and any retained
  parameters such as scatter and time.

- keys:

  The keywords to be written to the file.

- file.name:

  The name of the FCS file to be written.

- output.dir:

  A character string specifying the directory to save the unmixed FCS
  file.

## References

Laniewski, Nathan. *flowstate*.
<https://github.com/nlaniewski/flowstate>

Ellis B, Haaland P, Hahne F, Le Meur N, Gopalakrishnan N, Spidlen J,
Jiang M, Finak G (2025). *flowCore: Basic structures for flow cytometry
data*.
[doi:10.18129/B9.bioc.flowCore](https://doi.org/10.18129/B9.bioc.flowCore)
<https://bioconductor.org/packages/flowCore>

## See also

[`write.FCS`](https://rdrr.io/pkg/flowCore/man/write.FCS.html)
