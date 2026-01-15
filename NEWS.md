# AutoSpectral 0.9.2 (2026-01-25)

## Improvements
- Faster base R per-cell optimization.

## Bug fixes
- Patch to writing of "-A" in the channel names of FCS files. This was 
implemented in 0.9.0 but was incorrectly applied to all channels rather than
just the fluorescence parameters.

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
will be deprecated as it slows down the unmixing and does not meaningfully
measure the unmixing improvement.
- The `time.clean` option for `clean.controls()` will be deprecated. This uses
PeacoQC for time-based cleaning of single-stained control files. I've yet to see
this have an impact.
- The `trim` option for `clean.controls()` will be deprecated.

## Bug fixes
- Switch to FlowSOM for `SOM()` support. `EmbedSOM::SOM()` appears to have a
compilation error for Mac and has been removed from CRAN. Note that FlowSOM must
be installed separately using BiocManager.


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
