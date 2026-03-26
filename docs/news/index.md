# Changelog

## AutoSpectral 1.5.2 (2026-03-26)

### Bug fixes

- A patch for
  [`get.af.spectra()`](https://drcytometer.github.io/AutoSpectral/reference/get.af.spectra.md)
  to better handle situations where a given AF spectrum is essentially
  inverted, with the greatest magnitude signal being negative. This
  occurs occasionally on the FACSDiscover S8 and A8 cytometers. The
  spectra will now be kept as negative, since this provides better
  unmixing, empirically. For the plotting, any negative spectra will be
  inverted to avoid issues with scaling.
- A patch for
  [`get.af.spectra()`](https://drcytometer.github.io/AutoSpectral/reference/get.af.spectra.md)
  to exclude plotting of the unmixed data on biplots when only a single
  fluorophore is being unmixed (so only one axis can be drawn).
- Exclusion of a problematic keyword from BD FACSDiscover files. This
  affects
  [`define.keywords()`](https://drcytometer.github.io/AutoSpectral/reference/define.keywords.md),
  which updates the keywords after unmixing. Files were being produced,
  but the keyword headers were corrupted.
- Rearrange the arguments for
  [`tune.gate()`](https://drcytometer.github.io/AutoSpectral/reference/tune.gate.md),
  [`define.gate.landmarks()`](https://drcytometer.github.io/AutoSpectral/reference/define.gate.landmarks.md)
  and also
  [`define.gate.density()`](https://drcytometer.github.io/AutoSpectral/reference/define.gate.density.md)
  to put the `gate.name` argument earlier. The default value for
  `gate.name` has been removed because the user should be providing
  this.

## AutoSpectral 1.5.1 (2026-03-26)

### Bug fixes

- A patch for
  [`get.spectral.variants()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectral.variants.md)
  to resolve an issue that occurred when only a single fluorophore was
  being assessed.

## AutoSpectral 1.5.0 (2026-03-25)

### Improvements

- New gating approach adapted from `flowstate`. This uses cellular
  landmarks to identify the position of key populations on forward and
  side-scatter. Essentially, we backgate. This means that you can use
  well-expressed markers on known cell populations (think CD3, CD19,
  CD14) to define the location of your cells. This is fast, appears to
  be pretty robust, and should be easy for you, the user, to change how
  it works. Control over this is provided via the CSV control file
  spreadsheet using two new columns: “gate.name” to define which
  controls should share the same gate, and “gate.define” (TRUE/FALSE) to
  specify which samples should be used to define the gate boundaries
  (e.g., we might use CD4 but not TIM-3 or IL-4).
- To assist with the new gating, there is a
  [`tune.gate()`](https://drcytometer.github.io/AutoSpectral/reference/tune.gate.md)
  function that allows you to put in a range of parameters to quickly
  see the impact on the gate boundary prior to running
  [`define.flow.control()`](https://drcytometer.github.io/AutoSpectral/reference/define.flow.control.md).
- Native, faster FCS read/write functionality adapted from `flowstate`.
- FCS file concatenation via
  [`concatenateFCS()`](https://drcytometer.github.io/AutoSpectral/reference/concatenateFCS.md).
- Faster gating by reducing
  [`MASS::kde2d`](https://rdrr.io/pkg/MASS/man/kde2d.html) calls and
  allowing C++ kernel density estimation if `AutoSpectralRcpp` is
  installed.
- Faster plotting along the same lines.
- Hopefully graceful error handling during
  [`clean.controls()`](https://drcytometer.github.io/AutoSpectral/reference/clean.controls.md).
- Hopefully graceful error handling with diagnostic plotting during gate
  definition, both with the gate.define functions and directly in
  [`define.flow.control()`](https://drcytometer.github.io/AutoSpectral/reference/define.flow.control.md).
- Success/failure reporting from
  [`clean.controls()`](https://drcytometer.github.io/AutoSpectral/reference/clean.controls.md).
- Reduced memory usage when unmixing.
- Chunking of files when unmixing to support unmixing of any size of
  file.
- Spectral signature QC when running
  [`get.fluorophore.spectra()`](https://drcytometer.github.io/AutoSpectral/reference/get.fluorophore.spectra.md).
- Autofluorescence profile QC when running
  [`get.af.spectra()`](https://drcytometer.github.io/AutoSpectral/reference/get.af.spectra.md).
- Faster processing in
  [`get.af.spectra()`](https://drcytometer.github.io/AutoSpectral/reference/get.af.spectra.md)
  through integration of `AutoSpectralRcpp`, when available.
- Additional checks on FCS files when setting up and unmixing:
  consistency in the spectral channels used for unmixing, both in the
  names and the voltages/gains.

## AutoSpectral 1.0.0 (2026-02-10)

### Improvements

- Version 1.0.0 brings a revamp to how AutoSpectral identifies the best
  spectra on a per-cell basis. The theory behind it remains the same–we
  are still trying to identify the variation in the autofluorescence and
  fluorophores that best reduces the residual on a per-cell basis. Now,
  however, we do not need to do that using brute force. Instead, we can
  search only through variants (or autofluorescences) that align with a
  given cell’s residual. Thus we can pre-screen the variants to a select
  few and then test just those. This means we can figure out the
  solution in way less time. It also means that a native R
  implementation of the algorithm is possible in R in a somewhat
  reasonable time frame. So, that may help for anyone struggling to use
  the fast C++ version in `AutoSpectralRcpp`. Specifics on this will be
  detailed in an article on GitHub and Colibri Cytometry.
- Since we can now quickly identify which variants are useful for a
  given cell, we can test more variants, allowing a finer-grained view
  of the variation, which may improve unmixing quality.
- Autofluorescence extraction and fluorophore variation extraction are
  now modified to allow searching for more variation, focusing on
  “problematic” cells that remain far from where they should be when the
  first batch of variation is applied. This is most helpful for
  extracting autofluorescence in complex tissue samples, where
  AutoSpectral previously struggled to deal with the last few messy
  cells. To enable this, set `refine=TRUE`.
- Speed in unmixing should be the biggest change, particularly if you
  run using `AutoSpectralRcpp`.
- When extracting autofluorescence using
  [`get.af.spectra()`](https://drcytometer.github.io/AutoSpectral/reference/get.af.spectra.md),
  you can now get a set of plots showing you the unmixed data for the
  channels most affected by the autofluorescence (“worst channels”). The
  same channels will be plotted after a single round of autofluorescence
  extraction per cell (as in AutoSpectral v0.9.2 and earlier) as well as
  after the second round, using data from more difficult cells. To see
  this, run with `refine = TRUE`.
- Autofluorescence is now assigned to each cell using a shortcut to
  “project” where the AF will impact on fluorophore or residual space.
  This is especially fast for residual-based assignment.
- Perhaps most importantly, discontinuities that sometimes appeared in
  the data after unmixing using per-cell-fluorophore optimization,
  particularly with the “fast” approximation, should now be gone or at
  least greatly diminished.

### Bug fixes

- Deprecation warnings in 0.9.1 were not done properly, causing errors
  when the deprecated arguments were specified. That should now be
  fixed.

## AutoSpectral 0.9.1 (2026-01-15)

### Improvements

- Faster OLS and WLS unmixing for per-cell optimization in R in
  [`unmix.autospectral()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.autospectral.md).
  Perhaps this should be classified as a bug fix. The use of singular
  value decomposition rolled out in 0.9.0 will remain for matrix
  unmixing, but for per-cell optimization loops where the unmixing
  matrix is recalculated multiple times, a faster version is needed.
  [`unmix.ols.fast()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.ols.fast.md)
  and
  [`unmix.wls.fast()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.wls.fast.md)
  use [`solve()`](https://rdrr.io/r/base/solve.html) for this and have
  been benchmarked as the best among variously tested options for base R
  unmixing.

### Lifecycle warnings

- The `calculate.error` option for calculation of root mean squared
  error (RMSE) has been deprecated as it slows down the unmixing and
  does not meaningfully measure the unmixing improvement.
- The `time.clean` option for
  [`clean.controls()`](https://drcytometer.github.io/AutoSpectral/reference/clean.controls.md)
  will be deprecated. This uses PeacoQC for time-based cleaning of
  single-stained control files. I’ve yet to see this have an impact.
- The `trim` option for
  [`clean.controls()`](https://drcytometer.github.io/AutoSpectral/reference/clean.controls.md)
  will be deprecated.

### Bug fixes

- Switch to FlowSOM for `SOM()` support.
  [`EmbedSOM::SOM()`](https://rdrr.io/pkg/EmbedSOM/man/SOM.html) appears
  to have a compilation error for Mac and has been removed from CRAN.
  Note that FlowSOM must be installed separately using BiocManager.
- Patch to writing of “-A” in the channel names of FCS files. This was
  implemented in 0.9.0 but was incorrectly applied to all channels
  rather than just the fluorescence parameters.

### Notes

- Dependencies have been slimmed down. `tidyr`, `dplyr` and `rlang` have
  all been removed in favor of base R. Base R packages `stats`, `utils`
  and `grDevices` are called via `::` rather than imported into the
  NAMESPACE.

## AutoSpectral 0.9.0 (2025-12-23)

### New features

- Unmixing matrix can be saved via save.unmixing.matrix()
- Weights can be calculated via calculate.weights()
- Plotting of unmixing matrix in get.fluorophore.spectra

### Improvements

- More stable, faster parallel backend allowing mclapply in Mac when
  appropriate.
- Changes to get.spectral.variants, including permanent fixing of
  previously user-modifiable parameters and low-level denoising of
  spectra.
- More checks in check.control.file.
- Faster AutoSpectral unmixing in base R.
- Adjustments to reduce any discontinuities produced during unmixing.
- See also updates in AutoSpectralRcpp, including a large speed up and
  general improvement to the Poisson IRLS unmixing.
- Calculation of the unmixing matrix (Moore-Penrose pseudoinverse) will
  now be done using singular value decomposition
  [`svd()`](https://rdrr.io/r/base/svd.html) for numerical stability for
  all approaches. Up to now, it has been done with normal equations via
  [`solve()`](https://rdrr.io/r/base/solve.html). This should be better
  in edge cases. In most cases, the only difference will be floating
  point error. Calculation time is equivalent because almost all of the
  computational effort is on projecting the raw data into the unmixed
  space via the unmixing matrix, not calculating the unmixing matrix.
- FCS files will now be written with “-A” in the channel names, e.g.,
  “PE-A” rather than just “PE”.

### Bug fixes

- Bug patch for situations with beads using internal negatives in
  get.fluor.variants
- Patch to
  [`reload.flow.control()`](https://drcytometer.github.io/AutoSpectral/reference/reload.flow.control.md)
  bug affecting ID7000 samples.
- Patch to
  [`define.flow.control()`](https://drcytometer.github.io/AutoSpectral/reference/define.flow.control.md)
  affecting universal negative definitions and impacting on
  [`clean.controls()`](https://drcytometer.github.io/AutoSpectral/reference/clean.controls.md).
- Patch to
  [`check.control.file()`](https://drcytometer.github.io/AutoSpectral/reference/check.control.file.md)
  affecting Opteon samples.

------------------------------------------------------------------------

## AutoSpectral 0.8.7 (2025-12-01)

### New features

- Support for Symphony A5 SE
- Support for Cytek Northern Lights
- Shiny app for control file setup via AutoSpectralHelper
- Marker names will now be added to the control file based on matches in
  the FCS file names, where possible.
- The Hotspot(TM) matrix will be calculated and plotted as per the
  pre-print by Peter Mage et al.

### Improvements

- More improvements to plotting.
