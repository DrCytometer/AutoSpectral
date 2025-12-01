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
- More flexibility in negative/unstained samples.
- Allow multiple controls per fluorophore.
- An alternative to the automated gating.
- Cell-specific weighting for per-cell unmixing.
- Speed up per-cell fluorophore optimization.
- Fix the issue causing discontinuities.

To install the `dev` branch:

``` r
devtools::install_github("DrCytometer/AutoSpectral@dev")
```

To replace this with the (hopefully) stable version, run this:

``` r
remove.packages("AutoSpectral")
devtools::install_github("DrCytometer/AutoSpectral")
```
