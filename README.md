
<!-- README.md is generated from README.Rmd. Please edit that file -->

# AutoSpectral

<!-- badges: start -->

<!-- badges: end -->

AutoSpectral is AutoSpill updated for the spectral flow era.

The goal of AutoSpectral is to provide you with the best possible
spectral signatures of the fluorophores in your single-stained controls.
Whether or not these accurately model your fully stained samples will
depend on what you’ve chosen to use for the controls, how they were
prepared and other factors such as machine condition and any divergence
in handling between samples and controls.

More to the point, AutoSpectral is intended to make working with messy
cell-based controls as easy as compensation beads. This should give you
better accuracy and precision in your spectral definition and thus in
your unmixing.

Plus, you can extract each cell’s individual autofluorescent background
in a manner specific to that cell, producing better unmixing with less
spread. Per-cell fluorescence spectral optimization can reduce unmixing
errors.

At the moment, the following cytometers are supported:

- Cytek Aurora (“aurora”)
- *New! Cytek Northern Lights (“auroraNL”)*
- Sony ID7000 (“id7000”)
- BD FACSDiscoverS8 (“s8”)
- BD FACSDiscoverA8 (“a8”)
- *New! BD FACSymphony A5 SE (“a5se”)*
- Agilent NovoCyte Opteon (“opteon”)
- Beckman Coulter CytoFLEX mosaic (“mosaic”)
- ThermoFisher Attune Xenith (“xenith”)

## Installation

[![Stable](https://img.shields.io/badge/stable-master-blue)](https://github.com/DrCytometer/AutoSpectral)
[![Dev](https://img.shields.io/badge/dev-branch-orange)](https://github.com/DrCytometer/AutoSpectral/tree/dev)

### Latest Release

**Version 0.9.0**

To install the latest, hopefully stable version, install using
`devtools` or `remotes`. You will need to install the Bioconductor
packages separately, I believe.

``` r
# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("flowWorkspace", "flowCore", "PeacoQC"))

# You'll need devtools or remotes to install from GitHub.
# install.packages("devtools")
devtools::install_github("DrCytometer/AutoSpectral")
```

As of version 0.8.7, there is a Shiny helper tool to assist you in
setting up your AutoSpectral control files. This is an interactive html
app that opens in RStudio. Hopefully this makes things easier. It is
new, so again, probably not perfect. To try it, visit
[AutoSpectralHelper](https://github.com/DrCytometer/AutoSpectralHelper).
If you update AutoSpectral, I recommend downloading a new version of the
app to ensure compatibility.

To install a specific release, e.g., the previous one, use the version
number:

``` r
remotes::install_github("DrCytometer/AutoSpectral@v0.8.7")
```

### Dev branch

If you’re feeling adventurous or simply want early access to the latest
features, you can try the `dev` branch. At any given point, this may not
be working well.

AutoSpectral is open source. If you are interested in contributing,
please visit
[Development](https://drcytometer.github.io/AutoSpectral/articles/Development.html)
for suggestions of where help is needed most.

You can install the development version of AutoSpectral from
[GitHub](https://github.com/) with:

``` r
devtools::install_github("DrCytometer/AutoSpectral@dev")
```

## Bug fixes and known issues

AutoSpectral is pretty complex and newly released, so there will be
bugs. Sorry. Thanks to all of you providing feedback. Please submit any
and all issues either using the Issues page or via email at
colibri-cytometry at gmail.

I will be taking a break from addressing issues over the holiday period.

To submit a bug report, go to
[Issues](https://github.com/DrCytometer/AutoSpectral/issues).

For more general problems, like not being clear on how to do things,
something doesn’t work well, or maybe you have an idea for something new
or better, visit the [Discussions
page](https://github.com/DrCytometer/AutoSpectral/discussions).

Since one of my recent updates broke things, I’ll be moving to using
tagged releases that should be easier to install if the latest version
has flaws. I’ve also set up a separate development branch, which will
get the updates first. Things probably should have been that way from
the start, but this is all new to me.

Please check the [help pages and
articles](https://drcytometer.github.io/AutoSpectral/) if you’re
struggling to understand how to do something. There’s a lot of info
there.

In particular, see the [Full
Workflow](https://drcytometer.github.io/AutoSpectral/articles/Full_AutoSpectral_Workflow.html).

### Known shortcomings

- Gating. The automated gating is not great. See the [help
  page](https://drcytometer.github.io/AutoSpectral/articles/Gating.html)
  for tips. I’m working on an alternative.
- Please note that FCS 3.2 files from the S8 and A8 cytometers are not
  fully supported in flowCore. You may receive warnings, but things
  should still work.
- More stuff in progress will appear in the [Development
  article](https://drcytometer.github.io/AutoSpectral/articles/Development.html).
- This is my first R package.

If you want to use data from another cytometer and are wiling to provide
files for establishing the workflow, contact the author/maintainer. See
existing information, which may also assist you in setting up your
control file, in the [cytometer
database](https://docs.google.com/spreadsheets/d/1wj7QPkgpsuPNeVKyt-WWdBu5R48aZTgEbH8-_bpKeBY/edit?usp=sharing).

AutoSpectral relies on a database of information containing fluorophore
emission details. If your fluorophore is not detected automatically by
`create.control.file()` and you want to add it, visit the Google sheet
for the [fluorophore
database](https://docs.google.com/spreadsheets/d/14j4lAQ6dkjDBKMborDv_MkSptyNBqZiBsq5jNNSCoiQ/edit?usp=sharing)
and add it there. New fluorophores will be incorporated into updates.

Similarly, there is a marker database to detect (and standardize) marker
names if they appear in the single-stained FCS control file names. Feel
free to add more markers or synonyms to the [marker
database](https://docs.google.com/spreadsheets/d/16FAinR_Nfnl00mpHvmQFJT_uJJY3VUWk29yAaQ7HHn8/edit?usp=sharing).

This work has received funding from the KU Leuven C1 program, the
European Union’s Horizon 2020 research and innovation programme under
grant agreement No 874707 (EXIMIOUS), Wellcome Investigator Award
222442/A/21/Z, and UKRI Proactive Vaccinology Award MR/Y004450/1
(IMMPROVE).

AutoSpectral is provided under an AGPL3 licence.

## Example

Below is a basic example of the workflow, using samples from the ID7000.
This only illustrates weighted least squares unmixing. For per-cell
autofluorescence extraction or per-cell fluorophore optimization, see
the articles on those topics.

Please see the [Full
Workflow](https://drcytometer.github.io/AutoSpectral/articles/Full_AutoSpectral_Workflow.html)

``` r
library( AutoSpectral )

# Define the location of the single-stained control files.
control.dir <- "./Cell controls"

# Get the parameters for your cytometer.
# Supported cytometers include "aurora", "id7000", "a8", "s8" and "opteon".
asp <- get.autospectral.param( cytometer = "id7000", figures = TRUE )

# Optionally, create a control file (see article on this).
# After creating it, manually edit to ensure it is correct.
create.control.file( control.dir, asp )

# Locate the edited control file.
control.def.file <- "fcs_control_file.csv"

# check the control file for errors
control.file.errors <- check.control.file( control.dir, control.def.file, asp )

## Adjustments to automated gating
# This will be covered more extensively elsewhere.
# in this case, the cells on scatter are very small, so we modify the target.
asp$gate.bound.density.max.target.cells <- 0

# Load the control files, gate, prepare for spectral extraction.
# This is usually the slowest step.
flow.control <- define.flow.control( control.dir, control.def.file, asp )

# For speed and accuracy, the recommended approach is to clean the controls
# using separate (universal) negative samples for both beads and cells.
# In doing so, we select cells with matching autofluorescent background and
# restrict the calculations to the brightest events.
flow.control <- clean.controls( flow.control, asp )

# Extract the clean spectra
univ.neg.spectra <- get.fluorophore.spectra( flow.control, asp,  
                                             use.clean.expr = TRUE,
                              title = "Universal Negative Cells" )

# By default AutoSpectral extracts autofluorescence as a spectrum.
# This is comparable to "Autofluorescence as a Fluorescent Tag" in SpectroFlo
# More on this later.
# To remove this AF parameter for unmixing:
no.af.spectra <- univ.neg.spectra[ ! rownames( univ.neg.spectra ) == "AF", ]

## Unmixing
# Set the location of the raw files to be unmixed.
sample.dir <- "./Raw samples"

# Perform unmixing, here we will unmix all fcs files in the folder using 
# Weighted least squares (WLS or in Sony lingo, WLSM).
unmix.folder( sample.dir, no.af.spectra, asp, flow.control, method = "WLS" )
```

## Installation and Runtime

Installation via GitHub should take only a minute or so. It takes less
than that on a Dell i7 core laptop.

Occasionally, the help gets corrupted. Just re-install if that happens.
If you know why this happens, let me know.

Installation of `AutoSpectralRcpp` will take a couple of minutes because
the code needs to compile. You will also first have to install Rtools to
have a compiler, and that will take longer, probably 10 minutes or so.
Upgrading the BLAS and LAPACK (see below) also takes several minutes if
you’re taking the time to read the instructions carefully.

Runtime will vary considerably depending on the size and quantity of
files you’re processing. For the benchmark 42-colour Aurora data set
using single-stained cell controls with plenty of data, the slow steps
in the pre-processing pipeline are `define.flow.control()` and
`clean.controls()`. For the same Aurora dataset on the aforementioned i7
Windows laptop, I get: \* `define.flow.control()` v0.8.7: 12min
sequential, 9min parallel \* `define.flow.control()` v0.9.0: 4min
sequential, 2min parallel \* `clean.controls()` v0.8.7: 11min
sequential, not possible parallel \* `clean.controls()` v0.9.0: 11min
sequential, 6.5min parallel, not fully optimized \*
`get.spectral.variants()` v.0.8.7: 2min sequential, 1min parallel \*
`get.spectral.variants()` v.0.9.0: 65sec sequential, 32sec parallel

Unmixing time depends on the following variable:

- file size (number of cells/events, including debris)
- number of detectors
- number of fluorophores
- unmixing algorithm. OLS and WLS are fast, per-cell AF extraction will
  be ~50x slower, per-cell fluorophore optimization will be 2-4x slower
  than per-cell AF when running the “fast” approximation using
  AutoSpectralRcpp, or about 20-50x slower when using the “slow” exact
  calculation, again using AutoSpectralRcpp. Using the “fast”
  approximation in pure R will likely be comparable to or slower than
  the “slow” method in C++.

Benchmarks for “C3 Lung_GFP_003_Samples.fcs”, again on the 8-core
laptop, using OpenBLAS and AutoSpectralRcpp, where applicable:

- `unmix.fcs()` WLS or OLS v0.8.7: 9sec
- `unmix.fcs()` WLS or OLS v0.9.0: 9sec
- `unmix.fcs()` perCell AF extraction v0.8.7: 2min
- `unmix.fcs()` perCell AF extraction v0.9.0: 1min
- `unmix.fcs()` perCell fluorophore optimization “fast” v0.8.7: 9min
- `unmix.fcs()` perCell fluorophore optimization “fast” v0.9.0: \<2min
- `unmix.fcs()` perCell fluorophore optimization “slow” v0.8.7: 62min
- `unmix.fcs()` perCell fluorophore optimization “slow” v0.9.0: 19min
- `unmix.folder()` WLS or OLS, 6 files, v0.8.7: 67sec sequential,
  interrupted parallel
- `unmix.folder()` WLS or OLS, 6 files, v0.9.0: 62sec sequential, 31sec
  parallel

Still some work to be done.

## Go Faster

R is not known for its speed.

For faster processing there are three things you can do, hopefully all
fairly easy.

First, upgrade the BLAS and LAPACK libraries used by R. These provide
algorithms for linear algebra, which is the heart of spectral unmixing.

On Windows, simply swapping out your .dll files as in this tutorial can
give speed ups of 5x. [Install
OpenBLAS](https://github.com/david-cortes/R-openblas-in-windows) All
this involves is downloading the files from the internet, placing them
in the right folder and doing a quick restart.

On Mac, you likely want to use the vecLib library BLAS that ships with
Mac OS. The following articles may be helpful in setting this as the
default BLAS for use in R: [BLAS for Mac in
R](https://cran.r-project.org/bin/macosx/RMacOSX-FAQ.html#Which-BLAS-is-used-and-how-can-it-be-changed_003f)
[Performance BLAS](https://csantill.github.io/RPerformanceWBLAS/)
Feedback from users on Mac who successfully upgrade their BLAS would be
appreciated.

Do not set multiple threads for the BLAS as this will conflict with
higher level parallelization, either in AutoSpectral or other packages.

After upgrading the BLAS, you probably need to reinstall `Rcpp` and
`RcppArmadillo` to have this compiled with the BLAS.

``` r
install.packages( c( "Rcpp", "RcppArmadillo" ) )
```

Now, install AutoSpectralRcpp. This is fully accessible from R and
integrates with AutoSpectral. But, when it gets to the slow bits in the
unmixing, it switches over to calculating in C++, so it can be 10-100x
faster. Do this after upgrading the BLAS.

[AutoSpectralRcpp](https://github.com/DrCytometer/AutoSpectralRcpp)

You will need [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
to compile this.

You can install AutoSpectralRcpp like so:

``` r
devtools::install_github("DrCytometer/AutoSpectralRcpp")
```

Third, turn on parallel processing in AutoSpectral. *Update*: This is
now supported via a `parallel` backend. Importantly, this moves away
from `future` and `future_lapply`, which were aborting sometimes on
Windows. More usefully, this is a bit faster in many cases, and should
support forking via `mclapply` on Mac and Linux systems, which will be
much faster. Perhaps most usefully, this appears to allow
parallelization of the native R per-cell fluorophore optimization in
`unmix.autospectral`, so you should see a nice speed up there.

The parallel processing in AutoSpectralRcpp operates via OpenMP and
works well. It is always activated, but the number of threads can be
configured.

To activate parallel processing, check the function arguments for a
`parallel` option and set it to `TRUE`. Additionally, there is control
over the number of threads used, which should be directly in the
function call via a `threads` argument. If you don’t know how many
threads to use, check the recommendation for your machine after running
`get.autospectral.param()`:

``` r
asp$max.worker.process.n
```

Multithreading will always default to this number if you do not
explicitly set a number of threads and set `parallel=TRUE`. This is one
less than `parallelly::availableCores()`, so it is designed to allow you
to keep working on minor stuff while AutoSpectral chugs along in the
background.

For unmixing larger data sets, you will do well to use a machine with
more CPUs. Suggestions for faster processing are welcome. Some modest
improvements are in the works.

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
  - Patch to `define.flow.control()` affecting universal negative
    definitions and impacting on `clean.controls()`.
  - Calculation of the unmixing matrix (Moore-Penrose pseudoinverse)
    will now be done using singular value decomposition `svd()` for
    numerical stability for all approaches. Up to now, it has been done
    with normal equations via `solve()`. This should be better in edge
    cases. In most cases, the only difference will be floating point
    error. Calculation time is equivalent because almost all of the
    computational effort is on projecting the raw data into the unmixed
    space via the unmixing matrix, not calculating the unmixing matrix.
  - New functions to `save.unmixing.matrix` and `calculate.weights`
  - Patches to `define.flow.control` that were causing redundant gates
    to be created.
  - Code legibility formatting.
