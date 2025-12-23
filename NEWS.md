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
