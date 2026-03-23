# 13 Updates and Issues

## Updates and news

- Version 0.8.1: More fluorophores, rearranging detectors if needed

- Version 0.8.2: Support for Mosaic and Xenith cytometers

- Version 0.8.3: Patch for error introduced in 0.8.2

- Version 0.8.4: Changes to error messaging in check.control.file

- Version 0.8.5: Improvements to keyword handling in writing FCS files

- Version 0.8.6: Improvements to plotting, fluorophore matching

- Version 0.8.7: Support for Symphony A5 SE and Cytek Northern Lights.
  More improvements to plotting. Marker names will now be added to the
  control file based on matches in the FCS file names, where possible.
  The Hotspot(TM) matrix will be calculated and plotted as per the
  [preprint](https://www.biorxiv.org/content/10.1101/2025.04.17.649396v2.full.pdf)
  by Peter Mage et al.

- Version 0.9.0:

  - Changes to `get.spectral.variants`, including fixing of previously
    user-modifiable parameters, low-level denoising of spectra and a bug
    patch for situations with beads using internal negatives.
  - More checks in `check.control.file`.
  - New parallel backend.
  - Faster AutoSpectral unmixing in base R.
  - Adjustments to reduce any discontinuities produced during unmixing.
  - See also updates in `AutoSpectralRcpp`, including a large speed up
    and general improvement to the Poisson IRLS unmixing.
  - Patch to `reload.flow.control` bug affecting ID7000 samples.
  - Patch to
    [`define.flow.control()`](https://drcytometer.github.io/AutoSpectral/reference/define.flow.control.md)
    affecting universal negative definitions and impacting on
    [`clean.controls()`](https://drcytometer.github.io/AutoSpectral/reference/clean.controls.md).
  - Calculation of the unmixing matrix (Moore-Penrose pseudoinverse)
    will now be done using singular value decomposition
    [`svd()`](https://rdrr.io/r/base/svd.html) for numerical stability
    for all approaches. Up to now, it has been done with normal
    equations via [`solve()`](https://rdrr.io/r/base/solve.html). This
    should be better in edge cases. In most cases, the only difference
    will be floating point error. Calculation time is equivalent because
    almost all of the computational effort is on projecting the raw data
    into the unmixed space via the unmixing matrix, not calculating the
    unmixing matrix.
  - New functions to `save.unmixing.matrix` and `calculate.weights`
  - Patches to `define.flow.control` that were causing redundant gates
    to be created.
  - Code legibility formatting.

- Version 0.9.1:

- Switch to [`FlowSOM::SOM()`](https://rdrr.io/pkg/FlowSOM/man/SOM.html)
  from [`EmbedSOM::SOM()`](https://rdrr.io/pkg/EmbedSOM/man/SOM.html).

- Patch to appending “-A” suffix to parameter names.

- Version 0.9.2:

- Faster base R per-cell optimization.

- Patch to writing of “-A” in the channel names of FCS files.

- Version 1.0.0:

- Changes to
  [`unmix.autospectral()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.autospectral.md)
  to speed up processing and reduce discontinuities in the resulting
  unmixed data.

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
  now modified to search for more variation, focusing on “problematic”
  cells that remain far from where they should be when the first batch
  of variation is applied. This is most helpful for extracting
  autofluorescence in complex tissue samples, where AutoSpectral
  previously struggled to deal with the last few messy cells.

- Speed in unmixing should be the biggest change, particularly if you
  run using `AutoSpectralRcpp`.

- When extracting autofluorescence using
  [`get.af.spectra()`](https://drcytometer.github.io/AutoSpectral/reference/get.af.spectra.md),
  you will now get a set of plots showing you the unmixed data for the
  channels most affected by the autofluorescence (“worst channels”). The
  same channels will be plotted after a single round of autofluorescence
  extraction per cell (as in AutoSpectral v0.9.2 and earlier) as well as
  after the second round, using data from more difficult cells. To see
  this, run with `refine = TRUE`, which is the default setting now.

- Autofluorescence is now assigned to each cell using a shortcut to
  “project” where the AF will impact on fluorophore or residual space.
  This is especially fast for residual-based assignment.

- Perhaps most importantly, discontinuities that sometimes appeared in
  the data after unmixing using per-cell-fluorophore optimization,
  particularly with the “fast” approximation, should now be gone or at
  least greatly diminished.

- Another change to the base unmixing functions to provide a speed-up
  when unmixing in a per-cell loop. When performing standard OLS or WLS
  unmixing outside of per-cell optimization, singular value
  decomposition will still be used to calculate the unmixing matrix.

- Version 1.5.0:

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

- Native FCS read/write functionality adapted from `flowstate`.

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
