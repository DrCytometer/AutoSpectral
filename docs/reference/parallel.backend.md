# Parallel Backend

Sets up parallel processing backend using parallel and parallelly to
establish PSOCK clusters on Windows. Sets up clusters, defines
parLapply-based function.

## Usage

``` r
parallel.backend(
  asp,
  exports,
  threads,
  export.env = parent.frame(),
  dev.mode = FALSE,
  package.path = NULL
)
```

## Arguments

- asp:

  The AutoSpectral parameter list.

- exports:

  The vector of variables and functions to pass to the clusters.

- threads:

  Numeric, number of threads to use for parallel processing.

- export.env:

  The environment containing other functions and global variables
  potentially needed by the clusters.

- dev.mode:

  Logical, allows testing of function while in development. Default is
  `FALSE`.

- package.path:

  File.path to the R package files for AutoSpectral to permit loading of
  the functions via
  [`devtools::load_all()`](https://devtools.r-lib.org/reference/load_all.html)
  while in `dev.mode`.

## Value

An lapply function, either based on parLapply if parallel backend
initialization was successful or using sequential lapply if not.
