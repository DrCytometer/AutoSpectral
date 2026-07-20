# AutoSpectral 1.6.4 (2026-07-20)

## Improvements

- Reduced default number of variants produced by `get.spectral.variants()` by
dropping `som.dim` from `10` to `5`, which speeds up the unmixing 1.5 - 3x. 
Profiling suggests the impact on unmixing quality is minimal.
- More markers.

## Bug fixes

- `get.bd.spectra()` was looking for a `SPILL` keyword, but the newer FCS files
from the A8 and S8 use `$SPILL`. Fixed.


# AutoSpectral 1.6.3 (2026-07-19)

## Improvements

- Installation should now be simplified by setting remote install sites.
- More fluorophores
- More markers
- More spectral references for the Aurora

## Bug fixes

- Add missing documentation for internal functions.


# AutoSpectral 1.6.2 (2026-07-13)

## Improvements

- Voltage consistency checks are now applied when reading in the single-stained
  control files to ensure all samples have been recorded on the same settings.
  Similarly, the voltages/gains read from the single-stained controls are kept 
  for checking sample FCS files when performing unmixing, to ensure that the
  samples were acquired on the same settings. Inconsistencies in this regard will
  prompt a warning, since AutoSpectral does not have the ability to access 
  individual instrument linearity or QC information, which would be required to 
  accurately adjust the unmixing for different acquisition settings.
- When loading in the data via `define.flow.control()`, there are now checks for
  low event counts after applying the gates. If any contrl sample has below 
  `asp$min.cell.warning.n` (default is 500) events left, a fallback approach is
  triggered, defining a new gate using `define.gate.landmarks()` for that control
  sample. This is designed to catch situations where users haven't understood
  how to use the gating and have not defined gating groups in the control file or
  provided pre-determined gates.


# AutoSpectral 1.6.1 (2026-07-12)

## Improvements

- Better peak channel assignment handling for the various versions of the Aurora.
  When an instrument lacks particular channels that would be assigned as the peak,
  the spectral reference library is used to infer the next best peak, either by
  looking at the reference for the fluorophore, if available, or by looking at
  the mean signature for other fluorophores with the same peak emission.

## Bug fixes

- 7-AAD peak channel had a typo for the Aurora.


# AutoSpectral 1.6.0 (2026-07-01)

## New Features

- Automated spectral profile extraction using approaches developed by Nathan
  Laniewski in `flowState`. This simplifies the workflow considerably, and
  eliminates the need for gating, which was a key cause of stress. The new 
  automated approach employs projection orthogonalization to determine the peak
  channel (which way is the control going that isn't like the autofluorescence?),
  and uses cosine filtering (which events are least contaminated by AF?) to find
  the clean spectral signatures.
- New joint unmixing pipeline (`unmix.autospectral.joint()`), implemented in C++,
  using covariance-weighted spillover and residual alignment to solve autofluorescence
  and fluorophore assignment together rather than sequentially. Structurally collinear
  fluorophore pairs (e.g., APC/BUV661) are handled via joint-pair conflict resolution
  so genuine double-positive signal isn't suppressed. Seven new tuning parameters
  (`n.af.passes`, `cell.weight`, `noise.floor`, `alpha`, `collinear.threshold`,
  `joint.pair.resolution`, `refine.af.quantile`) are now wired through both
  `unmix.fcs()` and `unmix.folder()`.
- `calculate.optimize.necessity()` to determine which cells actually require
  per-cell fluorophore optimization, avoiding unnecessary computation.
- Multipass autofluorescence refinement using a matching-pursuit approach, with
  fixed pass-1 baseline denominators to prevent noise from dominating the
  highest-abundance cells in later passes.
- `simulate.flow.data()` generates synthetic ground-truth data by summing
  single-stained controls, providing a benchmark for false-negative rates on
  double-positive populations.
- New AF autofluorescence comparison and benchmarking functions for evaluating
  accuracy against ground truth.
- A merging of the approaches for per-cell autofluorescence assignment using the
  covariance structure: `assign.af.joint.cov()` to assess how autofluorescence
  variation impacts fluorophore values, with joint scoring of residuals. This 
  improves on the existing residual- and projection-based AF assignment methods.
- Extension of the variance-covariance error propagation plus residuals scoring
  metric to per-cell fluorophore optimization. This avoids the previous need to
  optimize in the reduced space of only positive fluorophores in each cell, which
  was responsible for discontinuities appearing in several data sets around zero.
- Legacy mode support. All pre-existing workflows remain intact, but some must
  be accessed via a `legacy` argument in the function calls.
- AutoSpectral landing page. Type `?AutoSpectral` in your IDE to open a workflow
  overview with links to key functions and articles.
- FCS 3.1 compliance. Modified FCS files now carry the `$ORIGINALITY`,
  `$LAST_MODIFIER`, and `$LAST_MODIFIED` keywords per the FCS 3.1 standard.

## Improvements

- C++ accelerated FCS reading and writing via `AutoSpectralRcpp`. Adapted from
  work by Samuel Granjeaud.
- Much faster plotting throughout, using a combination of point downsampling,
  raster rendering (`ragg`/`scattermore`), and C++ acceleration when
  `AutoSpectralRcpp` is available.
- Faster variant and AF spectra determination via `EmbedSOM` when
  `AutoSpectralRcpp` is installed.
- `viridis` added as a dependency for improved colour scaling in density plots.
- Support for `patchwork` to allow manual faceting and combination of plots.
- Spectral reference library updated with additional entries.
- Scatter-matching via k-nearest neighbours for higher precision.
- Cosine filtering of unstained samples against the fluorophore spectral profiles
  to reduce potential errors when the "unstained" isn't actually completely 
  unstained.
- Cosine filtering and kNN background matching for `get.spectra.variants()` to
  reduce undesireable influence of autofluorescence.
- Non-negative clamping of spectral profiles to prevent artefactual negative
  spectral values propagating into unmixing.
- Type-checking and validation added at function entry (e.g.,
  `get.spectral.variants()`) to guard against deprecated positional arguments
  silently corrupting subsequent parameters, with explicit error messages.
- Cell weighting is now enabled by default for the ID7000, and detection
  thresholds have been relaxed for this platform.
- `NAMESPACE` injection for optional dependencies (`AutoSpectralRcpp`, `EmbedSOM`)
  replaced with `requireNamespace()` calls for a safer, more portable pattern.
- Removed the `cairo` and `methods` package dependencies.
- Defensive `tryCatch` wrapping around plotting calls, and guards against
  unintended matrix-to-vector dropping.
- Guards added for low cell count edge cases to prevent dimension errors.
- General function and code cleanup throughout.

## Bug fixes

- Keywords issue causing files from BD instruments to occasionally emit scrambled
  data should now be fixed. Only whitelisted channels are preserved from the
  input raw FCS files.
- Fixed handling of negative autofluorescence spectra in normalization and plotting.
- Fixed hexbin rendering in N x N plots.
- Fixed a typo (`apply` should have been `lapply`).
- Various fixes to get `unmix.fcs()` passing its test suite again.


# AutoSpectral 1.5.6 (2026-04-21)

## Improvements

- M x N and N x N plotting of unmixed data using `unmixed.mxn.plot()` and 
`unmixed.nxn.plot()`
- Calculate the secondary stain index per "Evaluating the performance of Slingshot 
SpectraComp particles as universal single stain controls in flow cytometry" by 
Oliveira et al. Call `calculate.ssi()`.
- Plot the mismatch between two spectral profiles (e.g., beads vs. cells) as in 
the pre-print by Konecny et al. on unmixing-dependent spread. Call `ppectral.mismatch.plot()`.
- Compare unmixing of the same data side-by-side with two different spectral mixing
matrices. Call `compare.unmix()`.
- Multiple file formats when creating plots: pdf, tiff, png, jpg, etc.

## Bug fixes

- Now passing random seed from `asp$downsample.seed` throughout.
- Fixed an issue where `spectral.variant.plot()` would not push the plot to the viewer.


# AutoSpectral 1.5.5 (2026-04-18)

## Improvements

- More fluorophores.
- More A5SE spectral signatures.

## Bug fixes

- Spectral normalization to the peak channel (L-inf) when reading spectral profiles
from BD S8/A8 FCS files.


# AutoSpectral 1.5.4 (2026-04-13)

## Improvements

- More fluorophores.
- More markers.

## Bug fixes

- A change to `validate.control.file()` to permit duplication of the unstained
samples in the control file. Necessary for integration with the new `AutoSpectralHelper`
app.


# AutoSpectral 1.5.3 (2026-03-31)

## Improvements

- Change to only create output folders on an as-needed basis.
- Improve fluorophore matching to take the longest fluorophore name found.
- Fix Xenith transformation settings.
- Add spectral references for the Xenith.

## Bug fixes

- A patch to allow per-cell fluorophore optimization when only a single fluorophore
is present in the data. This was being blocked by checks on the structure of the
spectral variants being passed.


# AutoSpectral 1.5.2 (2026-03-26)

## Bug fixes
- A patch for `get.af.spectra()` to better handle situations where a given AF
spectrum is essentially inverted, with the greatest magnitude signal being negative.
This occurs occasionally on the FACSDiscover S8 and A8 cytometers. The spectra
will now be kept as negative, since this provides better unmixing, empirically.
For the plotting, any negative spectra will be inverted to avoid issues with 
scaling.
- A patch for `get.af.spectra()` to exclude plotting of the unmixed data on biplots
when only a single fluorophore is being unmixed (so only one axis can be drawn).
- Exclusion of a problematic keyword from BD FACSDiscover files. This affects
`define.keywords()`, which updates the keywords after unmixing. Files were being
produced, but the keyword headers were corrupted.
- Patch for some annoying keywords issues affecting writing of FACSymphony A5SE files.
- Patch to remove extra raw channels from unmixed Xenith files by default.
- Rearrange the arguments for `tune.gate()`, `define.gate.landmarks()` and also
`define.gate.density()` to put the `gate.name` argument earlier. The default value
for `gate.name` has been removed because the user should be providing this.
- Patch to restore the functionality of `include.raw` (include the raw spectral
data in the unmixed FCS file) in `unmix.fcs()`.

# AutoSpectral 1.5.1 (2026-03-26)

## Bug fixes
- A patch for `get.spectral.variants()` to resolve an issue that occurred when 
only a single fluorophore was being assessed.


# AutoSpectral 1.5.0 (2026-03-25)

## Improvements

- New gating approach adapted from `flowstate`. This uses cellular landmarks to
identify the position of key populations on forward and side-scatter. Essentially,
we backgate. This means that you can use well-expressed markers on known cell
populations (think CD3, CD19, CD14) to define the location of your cells. This is
fast, appears to be pretty robust, and should be easy for you, the user, to change
how it works. Control over this is provided via the CSV control file spreadsheet
using two new columns: "gate.name" to define which controls should share the same
gate, and "gate.define" (TRUE/FALSE) to specify which samples should be used to
define the gate boundaries (e.g., we might use CD4 but not TIM-3 or IL-4).
- To assist with the new gating, there is a `tune.gate()` function that allows
you to put in a range of parameters to quickly see the impact on the gate boundary
prior to running `define.flow.control()`.
- Native, faster FCS read/write functionality adapted from `flowstate`.
- FCS file concatenation via `concatenateFCS()`.
- Faster gating by reducing `MASS::kde2d` calls and allowing C++ kernel density
estimation if `AutoSpectralRcpp` is installed.
- Faster plotting along the same lines.
- Hopefully graceful error handling during `clean.controls()`.
- Hopefully graceful error handling with diagnostic plotting during gate definition,
both with the gate.define functions and directly in `define.flow.control()`.
- Success/failure reporting from `clean.controls()`.
- Reduced memory usage when unmixing.
- Chunking of files when unmixing to support unmixing of any size of file.
- Spectral signature QC when running `get.fluorophore.spectra()`.
- Autofluorescence profile QC when running `get.af.spectra()`.
- Faster processing in `get.af.spectra()` through integration of `AutoSpectralRcpp`,
when available.
- Additional checks on FCS files when setting up and unmixing: consistency in the
spectral channels used for unmixing, both in the names and the voltages/gains.


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
- Autofluorescence extraction and fluorophore variation extraction are now
modified to allow searching for more variation, focusing on "problematic" cells that remain
far from where they should be when the first batch of variation is applied. This
is most helpful for extracting autofluorescence in complex tissue samples, where
AutoSpectral previously struggled to deal with the last few messy cells. To enable
this, set `refine=TRUE`.
- Speed in unmixing should be the biggest change, particularly if you run using
`AutoSpectralRcpp`.
- When extracting autofluorescence using `get.af.spectra()`, you can now get a
set of plots showing you the unmixed data for the channels most affected by the
autofluorescence ("worst channels"). The same channels will be plotted after a 
single round of autofluorescence extraction per cell (as in AutoSpectral v0.9.2
and earlier) as well as after the second round, using data from more difficult
cells. To see this, run with `refine = TRUE`.
- Autofluorescence is now assigned to each cell using a shortcut to "project"
where the AF will impact on fluorophore or residual space. This is especially fast
for residual-based assignment.
- Perhaps most importantly, discontinuities that sometimes appeared in the data
after unmixing using per-cell-fluorophore optimization, particularly with the
"fast" approximation, should now be gone or at least greatly diminished.

## Bug fixes
- Deprecation warnings in 0.9.1 were not done properly, causing errors when the
deprecated arguments were specified. That should now be fixed.


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
