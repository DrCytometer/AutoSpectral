
<!-- README.md is generated from README.Rmd. Please edit that file -->

# AutoSpectral

<!-- badges: start -->

<!-- badges: end -->

AutoSpectral is AutoSpill updated for the spectral flow era.

It is also pretty complex and newly released, so there will be bugs.
Sorry. Please read the [help pages and
articles](https://drcytometer.github.io/AutoSpectral/) before submitting
bug reports. There’s a lot of info there.

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
spread.

At the moment, the following cytometers are supported:

- Cytek Aurora/Northern Lights (“aurora”)
- Sony ID7000 (“id7000”)
- BD FACSDiscoverS8 (“s8”)
- BD FACSDiscoverA8 (“a8”)
- Agilent NovoCyte Opteon (“opteon”)
- Beckman Coulter CytoFLEX mosaic (“mosaic”)
- ThermoFisher Attune Xenith (“xenith”)

Support for the Mosaic and Xenith is new in version 0.8.2.

There will likely be some unresolved issues with plotting data from
certain configurations of the ID7000 and Aurora. Edit: this should now
be fixed.

Please note that FCS 3.2 files from the S8 and A8 cytometers are not
fully supported in flowCore. You may receive warnings, but things should
still work.

If you want to use data from another cytometer and are wiling to provide
files for establishing the workflow, contact the author/maintainer.

This work has received funding from the KU Leuven C1 program, the
European Union’s Horizon 2020 research and innovation programme under
grant agreement No 874707 (EXIMIOUS), Wellcome Investigator Award,
222442/A/21/Z and UKRI Proactive Vaccinology Award, MR/Y004450/1
(IMMPROVE).

AutoSpectral is provided under an AGPL3 licence.

## Installation

You can install the development version of AutoSpectral from
[GitHub](https://github.com/) with:

``` r

# Optional: Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("flowWorkspace", "flowCore", "PeacoQC", "Biobase"))

# You'll need devtools or remotes to install from GitHub.
# install.packages("devtools")
devtools::install_github("DrCytometer/AutoSpectral")
```

## Example

This is a basic example of the workflow, using samples from the ID7000.
This only illustrates weighted least squares unmixing. For per-cell
autofluorescence extraction or per-cell fluorophore optimization, see
the articles on those topics.

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

Second, install AutoSpectralRcpp. This is fully accessible from R and
integrates with AutoSpectral. But, when it gets to the slow bits in the
unmixing, it switches over to calculating in C++, so it can be 10-100x
faster.

[AutoSpectralRcpp](https://github.com/DrCytometer/AutoSpectralRcpp)

You can install the development version of AutoSpectralRcpp like so:

``` r
devtools::install_github("DrCytometer/AutoSpectralRcpp")
```

Third, turn on parallel processing in AutoSpectral. At the moment, this
is not fully optimized in AutoSpectral. The parallel processing in
AutoSpectralRcpp operates via OpenMP and works well. It is always
activated, but the number of threads can be configured.

To activate parallel processing, either set the `parallel` flag to
`TRUE` in `asp` or check the function arguments for a `parallel` option.
Additionally, there is control over the number of threads used, either
directly in the function call or via `asp$max.worker.process.n`. Sorry,
this will be cleaned up soon.

``` r
asp$parallel <- TRUE
asp$max.worker.process.n
```

For unmixing larger data sets, you will do well to use a machine with
more CPUs. Suggestions for faster processing are welcome. Some modest
improvements are in the works.
