# 10 Re-running AutoSpectral

## Reloading AutoSpectral to quickly return to your unmixing

In this article, we’ll cover some shortcuts to allow you to return to
AutoSpectral once you’ve already run an experiment. You may wish to
return to an experiment to unmix it again. If you want to create new
spectra, perhaps trying a new cleaning method or changing the gating,
start as you would normally. If however, all you want to do is re-unmix
or change the parameter names on the created FCS files, you don’t need
to go through the whole process of loading the single-stained control
files, gating, cleaning and extracting the spectra. In this article,
we’re going to look at the short, quick workflow to allow you to return
to a previous AutoSpectral run and go straight to unmixing. Note that
this workflow can be adapted to use spectral matrices obtained using
other methods, such as FlowJo or SpectroFlo, if you prefer. Some details
on that at the end.

We’re going to start, as always, by loading AutoSpectral and setting the
parameters for our cytometer. We also need the folder containing the
single-stained controls and the control file spreadsheet. You should
already have all of this.

``` r
library(AutoSpectral)
asp <- get.autospectral.param(cytometer = "id7000")
control.dir <- "path_to_my_single_stained_controls"
control.def.file <- "fcs_control_file.csv"
```

Then, we create a shell of the `flow.control` list, mostly to provide us
with the channels we’ll use for spectral unmixing.

``` r
flow.control <- reload.flow.control(control.dir, control.def.file, asp)
```

We also need to load in our spectra. If we’ve run AutoSpectral already,
these will be in `./table_spectra` and will, by default, be called
something like “Initial_autospectral_spectra.csv”. We need to read this
into R. For this, there’s a very simple helper function that reads the
named CSV file and sets the structure back as it was when written. If
your file isn’t in that folder, use the `spectra.dir` argument to define
the location.

``` r
spectra.file <- "Initial_autospectral_spectra.csv"
spectral.matrix <- read.spectra(spectra.file)
```

There’s an option in `read.spectra` to remove the autofluorescence `AF`
parameter while reading in. For this, set `remove.af` to `TRUE`. Note
that you can change which fluorophore(s) are excluded by naming them to
the `af.param` argument.

Now, we have all the information we need to unmix our files. We just
need to point AutoSpectral to the folder containing the FCS files we
want to unmix and tell it how we want to perform the unmixing. Here,
we’re using weighted least squares by calling `WLS`.

``` r
sample.dir <- "./Raw samples"
unmix.folder( sample.dir, spectral.matrix, asp, flow.control, method = "WLS" )
```

The unmixed FCS files go to the “AutoSpectral_unmixed” folder in your
current working directory, or, if you prefer, to the folder of your
choice if you specify this to unmix.folder as an argument to
`output.dir`.

------------------------------------------------------------------------

## Loading non-AutoSpectral spectral matrices

### FlowJo Spectra

You can also load a spectral matrix you’ve prepared using something
else. For instance, we can spectra that we isolate in FlowJo or directly
pull them from the SpectroFlo Experiment file. Why would you want to do
this? Well, perhaps you’re struggling with AutoSpectral and know you can
get decent spectra from your beads in FlowJo or know that you’ve already
gotten clean spectra on your cytometer’s software. Maybe you’d just like
to try out different unmixing approaches, perhaps checking the per-cell
autofluorescence extraction. Let’s look at how to do this.

You’ll need to run the steps above up to the point of loading in the
spectra. In other words, you need to create `asp`, `control.dir`,
`control.def.file` and `flow.control`. Hopefully this can be simplified
in the future.

Now, to load in a matrix from FlowJo, first create it using the Spectral
Population Viewer tool. Use negative and positive gating, as appropriate
for your sample.
<https://docs.flowjo.com/flowjo/advanced-features/spectral-populations/>

![FlowJo Spectral Viewer](figures/FlowJo_Spectral_Viewer.jpg)

FlowJo Spectral Viewer

To get the spectral (spillover) data, first set the Spectral Viewer Plot
to normalized values by going to Options/Normalized/Only Positive
Values. Then, right click on the plot and “Copy Content”. Paste this
into your favorite spreadsheet program (e.g., Excel), which will
probably give you a mess like this:

![FlowJo Matrix Formatting](figures/FlowJo_matrix_formatting.png)

FlowJo Matrix Formatting

![FlowJo Matrix Formatting](figures/Reformatted_FlowJo_Matrix.jpg)

FlowJo Matrix Formatting

To reformat it, click on the clipboard icon and select “Use Text Import
Wizard”. Select “delimited”, hit “Next”, then select “Comma” and finally
“Finish”.

![Text Import Wizard](figures/Text_import_wizard.jpg)

Text Import Wizard

![FlowJo Matrix Reformatted](figures/Reformatted_FlowJo_Matrix.jpg)

FlowJo Matrix Reformatted

The spectra are normalized to 100 in FlowJo, whereas in AutoSpectral and
in SpectroFlo everything is normalized to 1. You’ll need to re-normalize
to 1 (i.e., divide by 100) if you want to use this for unmixing. If you
start combining data from different sources with different normalization
schemes, re-normalize before unmixing or the normalization will carry
over into the scaling of the unmixed data, which you don’t want. This
renormalization happens automatically if you use the read.spectra
function below.

You should probably rename the stuff in column A to be the names of the
fluorophores, as this is what you’ll be importing into R.

Save the spreadsheet as a CSV file. Then you can read it into R as
follows:

``` r
spectra.file <- "FlowJo_Spectra.csv"
flowjo.spectra <- read.spectra(spectra.file, spectra.dir = "~/AutoSpectral_data/AlternativeMatrices/")
flowjo.spectra
```

### SpectroFlo Spectra

Pulling the spillover matrix directly from your SpectroFlo experiment is
also possible, provided you’ve gone through the unmixing process in
SpectroFlo. This should work whether you’ve used controls as part of a
Reference Group or used stored controls from the library. In other
words, we can obtain the spillover matrix that was used for the unmixing
in SpectroFlo.

List the name and location of your experiment file. The function will
extract the spillover matrix both to an object in R as well as writing a
copy to disk as a CSV file. Define where you’d like the CSV file to end
up as `output.dir`.

``` r
spectroflo.expt <- "~/AutoSpectral_data/AlternativeMatrices/20240620 Spectral Symposium-poor cell unmixed.Expt"
output.dir <- "./table_spectra"
spectral.matrix <- read.spectroflo.expt(spectroflo.expt, output.dir)
spectral.matrix
```

Note: this sometimes gives odd spectra. Check by plotting
`spectral.trace` or `spectral.heatmap`. I suspect this happens if you
re-import the experiment into an non-cytometer-connected version of
SpectroFlo, but I need to look into it a bit more.

Note 2: By default, the autofluorescence spectra stored in SpectroFlo
are not in a normalized format. I do not know why this is. When you read
them in using
[`read.spectroflo.expt()`](https://drcytometer.github.io/AutoSpectral/reference/read.spectroflo.expt.md),
they will be L-infinity (peak) normalized.

### BD (Chorus) Spectra

For BD instruments, the spectra are stored correctly in the `SPILL`
keyword in the unmixed flow files. So, these should be in any FCS file
from the A8 or S8 that is not a single-stained control and where you
have run the unmixing in Chorus.

``` r
s8.fcs <- "~/AutoSpectral_data/S8_data/BP0323502_1.fcs"
chorus.spectra <- read.bd.spectra(s8.fcs)
chorus.spectra
```

A plea to instrument manufacturers: Please, please start doing what BD
is doing. Including the fluorophore spectra used for the unmixing
calculation in the FCS file provides really important information for QC
and reproducibility.
