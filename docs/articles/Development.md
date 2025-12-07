# Development

Here’s some stuff that’s in progress. These changes will show up first
in the development `dev` branch, which you are welcome to install and
test.

I don’t have a particular timeframe or priority list for any of this, so
ping me if any of this is important to you.

- Parallelization improvements.
  - The parallelization of functions like `define.flow.control` can be
    problematic on Windows due to memory usage with `futures` and the
    tendency of Windows to usurp threads that R was using.
  - A new parallel backend using `parLapply` on Windows or `mclapply` on
    Mac/Linux is being implemented in `dev` branch. This should
    gracefully fall back to sequential `lapply` processing if it runs
    into any problems.
- More flexibility in negative/unstained samples.
- Allow multiple controls per fluorophore.
- An alternative to the automated gating.
  - This is now in progress.
- Cell-specific weighting for per-cell unmixing.
  - This is problematic for the ID7000, probably due to the PMT noise.
- Speed up per-cell fluorophore optimization.
  - Some progress on this has been implemented in `dev` branch. Should
    be ~4x faster now.
- Fix the issue causing discontinuities.
  - Some progress on this in `dev` branch.

To install the `dev` branch:

``` r
devtools::install_github("DrCytometer/AutoSpectral@dev")
```

To replace this with the (hopefully) stable version, run this:

``` r
remove.packages("AutoSpectral")
devtools::install_github("DrCytometer/AutoSpectral")
```
