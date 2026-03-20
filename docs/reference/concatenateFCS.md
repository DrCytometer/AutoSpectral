# Concatenate Multiple FCS Files

Takes a string/list of input .fcs file paths, reads in the data,
combines it and writes a new .fcs file with all the data together.

## Usage

``` r
concatenateFCS(
  fcs.paths,
  output.name = "Concatenated.fcs",
  output.dir = "./concatenated_fcs"
)
```

## Arguments

- fcs.paths:

  A character vector of full file paths to the input .fcs files.

- output.name:

  A character string for the name of the new .fcs file. Default is
  `Concatenated.fcs`.

- output.dir:

  A character string specifying the directory to save the new file.
  Default is `./concatenated_fcs`

## Value

Path to the newly created FCS file.
