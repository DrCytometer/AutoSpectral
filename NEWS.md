# AutoSpectral 1.0.0 (2026-02-10)

## Improvements
- Version 1.0.0 brings a revamp to how AutoSpectral identifies the best spectra
on a per-cell basis. The theory behind it remains the same--we are still trying
to identify the variation in the autofluorescence and fluorophores that best
reduces the residual on a per-cell basis. Now, however, we do not need to do that
using brute force. Instead, we can search only through variants (or 
autofluorescences) that align with a given cell's residual. Thus we can pre-screen
the variants to a select few and then test just those. This means we can figure
out the solution in way less time. It also means that a native R implementation
of the algorithm is possible in R in a somewhat reasonable time frame. So, that
may help for anyone struggling to use the fast C++ version in `AutoSpectralRcpp`.
Specifics on this will be detailed in an article on GitHub and Colibri Cytometry.
- Since we can now quickly identify which variants are useful for a given cell,
we can test more variants, allowing a finer-grained view of the variation, which
may improve unmixing quality.

## Bug fixes

# AutoSpectral 0.9.1 (2026-01-15)

## Improvements
- Faster OLS and WLS unmixing for per-cell optimization in R in 
`unmix.autospectral()`. Perhaps this should be classified as a bug fix. The use
of singular value decomposition rolled out in 0.9.0 will remain for matrix
unmixing, but for per-cell optimization loops where the unmixing matrix is
recalculated multiple times, a faster version is needed. `unmix.ols.fast()` and
`unmix.wls.fast()` use `solve()` for this and have been benchmarked as the best
among variously tested options for base R unmixing.

## Lifecycle warnings
- The `calculate.error` option for calculation of root mean squared error (RMSE)
has been deprecated as it slows down the unmixing and does not meaningfully
measure the unmixing improvement.
- The `time.clean` option for `clean.controls()` will be deprecated. This uses
PeacoQC for time-based cleaning of single-stained control files. I've yet to see
this have an impact.
- The `trim` option for `clean.controls()` will be deprecated.

## Bug fixes
- Switch to FlowSOM for `SOM()` support. `EmbedSOM::SOM()` appears to have a
compilation error for Mac and has been removed from CRAN. Note that FlowSOM must
be installed separately using BiocManager.
- Patch to writing of "-A" in the channel names of FCS files. This was 
implemented in 0.9.0 but was incorrectly applied to all channels rather than
just the fluorescence parameters.

## Notes
- Dependencies have been slimmed down. `tidyr`, `dplyr` and `rlang` have all been
removed in favor of base R. Base R packages `stats`, `utils` and `grDevices` are
called via `::` rather than imported into the NAMESPACE.


# AutoSpectral 0.9.0 (2025-12-23)

## New features
- Unmixing matrix can be saved via save.unmixing.matrix()
- Weights can be calculated via calculate.weights()
- Plotting of unmixing matrix in get.fluorophore.spectra

## Improvements
- More stable, faster parallel backend allowing mclapply in Mac when appropriate.
- Changes to get.spectral.variants, including permanent fixing of previously 
user-modifiable parameters and low-level denoising of spectra.
- More checks in check.control.file.
- Faster AutoSpectral unmixing in base R.
- Adjustments to reduce any discontinuities produced during unmixing.
- See also updates in AutoSpectralRcpp, including a large speed up and general 
improvement to the Poisson IRLS unmixing.
- Calculation of the unmixing matrix (Moore-Penrose pseudoinverse) will now be
done using singular value decomposition `svd()` for numerical stability for all
approaches. Up to now, it has been done with normal equations via `solve()`.
This should be better in edge cases. In most cases, the only difference will be
floating point error. Calculation time is equivalent because almost all of the
computational effort is on projecting the raw data into the unmixed space via
the unmixing matrix, not calculating the unmixing matrix.
- FCS files will now be written with "-A" in the channel names, e.g., "PE-A"
rather than just "PE".

## Bug fixes
- Bug patch for situations with beads using internal negatives in 
get.fluor.variants
- Patch to `reload.flow.control()` bug affecting ID7000 samples.
- Patch to `define.flow.control()` affecting universal negative definitions and 
impacting on `clean.controls()`.
- Patch to `check.control.file()` affecting Opteon samples.


---
# AutoSpectral 0.8.7 (2025-12-01)

## New features
- Support for Symphony A5 SE
- Support for Cytek Northern Lights
- Shiny app for control file setup via AutoSpectralHelper
- Marker names will now be added to the control file based on matches in the 
FCS file names, where possible. 
- The Hotspot(TM) matrix will be calculated and plotted as per the pre-print by 
Peter Mage et al.

## Improvements
- More improvements to plotting.
