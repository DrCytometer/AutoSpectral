# 04 ID7000 WLS Example

## A simple ID7000 workflow

Below is a basic example of the workflow, using samples from the ID7000.
This only illustrates weighted least squares unmixing. For per-cell
autofluorescence extraction or per-cell fluorophore optimization, see
the articles on those topics. The point of this workflow is to show how
the control clean-up and automated determination of fluorophore spectra
can mimic or improve upon the unmixing obtained on the Sony using their
WLSM unmixing. The workflow here just uses weighted least squares (WLS)
for unmixing.

Please see the [Full
Workflow](https://drcytometer.github.io/AutoSpectral/articles/01_Full_AutoSpectral_Workflow.html)

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
univ.neg.spectra <- get.fluorophore.spectra(
   flow.control,
   asp,
   use.clean.expr = TRUE,
   title = "Universal Negative Cells" 
)

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
unmix.folder(
   sample.dir,
   no.af.spectra,
   asp,
   flow.control,
   method = "WLS" 
)
```
