# Create Parallel Lapply

Sets up parallel processing lapply function for Windows, Mac OS or
Linux.

## Usage

``` r
create.parallel.lapply(
  asp,
  exports,
  parallel = TRUE,
  threads = NULL,
  export.env = parent.frame(),
  dev.mode = FALSE,
  package.path = NULL,
  allow.mclapply.mac = FALSE
)
```

## Arguments

- asp:

  The AutoSpectral parameter list.

- exports:

  The vector of variables and functions to pass to the clusters.

- parallel:

  Logical, controls whether parallel processing is used. Default is
  `TRUE`.

- threads:

  Numeric, number of threads to use for parallel processing. Default is
  `NULL` which will revert to `asp$worker.process.n` if `parallel=TRUE`.

- export.env:

  The environment containing other functions and global variables
  potentially needed by the clusters. Default is
  [`parent.frame()`](https://rdrr.io/r/base/sys.parent.html).

- dev.mode:

  Logical, allows testing of function while in development. Default is
  `FALSE`.

- package.path:

  File.path to the R package files for AutoSpectral to permit loading of
  the functions via
  [`devtools::load_all()`](https://devtools.r-lib.org/reference/load_all.html)
  while in `dev.mode`. Default is `NULL`.

- allow.mclapply.mac:

  Logical, if `TRUE` permits `mclapply()` forking on Mac OS. Default
  `FALSE` forces PSOCK cluster use on Mac to prevent multithreaded
  Accelerate BLAS from crashing parallels when matrix ops are used.

## Value

An lapply function, either based on `parLapply` for Windows or
`mcLapply` on Mac OS and Linux. If parallel backend initialization
fails, sequential `lapply` is returned.
