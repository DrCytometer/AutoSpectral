# Development

Here’s some stuff that’s in progress. These changes will show up first
in the development `dev` branch, which you are welcome to install and
test.

If you would like to contribute to AutoSpectral, these are key areas I
have identified where it has shortcomings or needs improvement.
Suggestions for strategies, or better, code suggestions, would be
appreciated.

- Parallelization improvements.

  - The parallelization of functions like `define.flow.control` was
    problematic on Windows due to memory usage with `futures` and the
    tendency of Windows to usurp threads that R was using.
  - A new parallel backend using `parLapply` on Windows or `mclapply` on
    Mac/Linux is being implemented in `dev` branch. This should
    gracefully fall back to sequential `lapply` processing if it runs
    into any problems when being set-up.
  - There is currently a problem with the new backend in
    `get.spectral.variants` where in some cases it does not assess all
    fluorophores. Presumably there are some items missing from the
    exports to the clusters, but I don’t know.

- Code clean-up for `define.flow.control`.

  - This was the first part of the code I wrote, since it is required
    for everything downstream. It is, consequently, probably the worst.
  - I’ve made a little progress on this, but any changes here have
    ramifications throughout, so it is slow going.
  - Currently, we read in all of the FCS data into active memory in R,
    gating as we go to reduce data size. It might be better to
    restructure this to simply generate gates and retain indices of
    gated events in each file. FCS files would need to be read in on the
    fly for `clean.controls` and `get.fluorophore.spectra`. This would
    require a faster FCS reader, such as the `data.table`-based approach
    in Nathan Laniewski’s
    [flowState](https://github.com/nlaniewski/flowstate).

- More flexibility in negative/unstained samples.

- Allow multiple controls per fluorophore.

- An alternative to the automated gating.

  - This is now in progress.
  - This will be based on the “Cellular Landmarks” approach from Nathan
    Laniewski’s [flowState](https://github.com/nlaniewski/flowstate)
  - Initial definition of the gate will be based on the location of the
    brightest positive events in controls where the marker corresponds
    to known lymphocyte or monocyte populations. The code should fall
    back to the current gating approach when no known markers are
    provided.

- Cell-specific weighting for per-cell unmixing.

  - This is problematic for the ID7000, probably due to the PMT noise.
  - Benefits on other systems are surprisingly minimal, so this may not
    be implemented.

- Speed up per-cell AF extraction

  - This currently operates via a `for` loop, unmixing each AF signature
    sequentially. This is the best I have found. It is okay with fast
    BLAS.
  - Parallelization of this using `future` was not faster and fails with
    large data sets due to memory overrun.
  - I think it is unlikely that parallelization of this will help since
    computation of the unmixing will be faster than communication of the
    large expression matrix to different cores.
  - Approaches using `lapply` and `vapply` were slower and required more
    memory.

- Speed up per-cell fluorophore optimization.

  - Some progress on this has been implemented in `dev` branch. Should
    be ~4x faster now.
  - The primary improvement comes from generating fewer spectral
    variants to search through. My testing shows that a lower number is,
    if anything, better, presumably due to slightly higher quality of
    the spectra being generated while still covering the range of
    variation. This will help in R as well as Rcpp.
  - Secondary improvement comes from changing the solving strategy in
    C++ so that we aren’t recalculating the unmixing matrix (righthand
    solve) every time.
  - Tertiary improvment from changing the C++ compiler flags.
  - The pure R version will now benefit from parallelization,
    particularly on Mac/Linux.
  - A lot more work is needed on this. I need a smarter strategy.

- Fix the issue causing discontinuities.

  - Some progress on this in `dev` branch.
  - One improvement is better decontamination of AF from the spectral
    variants, which was causing cells to be pushed further away from the
    threshold for optimization.
  - A second is scaling the changes to cellular spectra based on the
    abundance of the fluorophore in question. This means cells close to
    the threshold will not change much, so the optimization will apply
    to higher expression cells, which is what is needed for reducing
    spillover spread and unmixing errors.
  - A third is the reduction in low-level noise introduced into the
    spectral variants due to fluctuations in electronic noise or
    autofluorescence. Any signal in the variant spectrum in a channel
    where the optimized single spectrum has less than 5% of the peak
    detector signal is now regularized towards the optimized spectrum.
    This focuses the variation on the “peaks”, which is what causes the
    unmixing errors and spread.
  - The primary remaining problem is, I think, with the crude
    approximation employed in the “fast” approach. The best solution
    would be to speed up the “slow” method to the point where the “fast”
    approach can be deprecated. Alternatively, some degree of smoothing
    based on kNN or regularization may be required.

- Better correction of unmixing errors

  - There are several strategies I’m investigating to handling the cases
    where the multi-colour samples contain obvious errors outside the
    range of variation seen in the single-colour controls. This is not
    stuff I’m going to put online, so get in touch if you’d like to work
    on this.

- Integration of Poisson IRLS

  - The Poisson IRLS in `AutoSpectralRcpp` is now working considerably
    better and faster. This may allow for integration of the Poisson
    optimization after identification of cellular autofluorescence and
    fluorophore signatures.
  - For this to be practical, per-cell fluorophore optimization needs to
    be faster.

To install the `dev` branch:

``` r
devtools::install_github("DrCytometer/AutoSpectral@dev")
```

To replace this with the (hopefully) stable version, run this:

``` r
remove.packages("AutoSpectral")
devtools::install_github("DrCytometer/AutoSpectral")
```
