# AutoSpectral

**Version 1.0.0**

Welcome to the Development branch of AutoSpectral.

This will always be a work-in-progress, and will not necessarily be
stable at any given point. At the moment, most of the changes taking
place relate to speeding things up or fixing any reported bugs.

To install the development version, run the code below in your
R/RStudio. This will overwrite any prior installation of AutoSpectral.
To restore the stable version of AutoSpectral, simply re-run the
`install_github()` command without the @dev bit.

``` r
# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("flowWorkspace", "flowCore", "PeacoQC", "FlowSOM"))

# You'll need devtools or remotes to install from GitHub.
# install.packages("devtools")
# install.packages("remotes")
remotes::install_github("DrCytometer/AutoSpectral@dev")
```

## Introduction

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
your unmixing. AutoSpectral does not negate the requirement for good
controls or good panel design. It will work better with better controls.

AutoSpectral aims to provide reproducible unmixing, meaning that anyone
should be able to obtain reliably good (probably better) unmixing from
the same set of controls. The aim is to remove the human “fiddling”
part, which is slow and not always so scientific.

Plus, you can extract each cell’s individual autofluorescent background
in a manner specific to that cell, producing better unmixing with less
spread. Per-cell fluorescence spectral optimization can reduce unmixing
errors in some cases.

The central ideas are as follows:

- Encode Roederer’s Rules for controls as much as possible.
- Deal with cellular autofluorescence interfering with identifying
  fluorophore signatures.
- Autofluorescence is variable on a cell-to-cell basis, both in
  magnitude and in type. We can figure out, more or less, what this
  should be for each cell.
- Fluorophore emissions are variable, and this variability manifests on
  the level of the cell. Again, we can more or less figure this out and
  deal with it.

At the moment, the following cytometers are supported:

- Cytek Aurora (“aurora”)
- Cytek Northern Lights (“auroraNL”)
- Sony ID7000 (“id7000”)
- BD FACSDiscoverS8 (“s8”)
- BD FACSDiscoverA8 (“a8”)
- BD FACSymphony A5 SE (“a5se”) (probably still some work to be done for
  this one)
- Agilent NovoCyte Opteon (“opteon”)
- Beckman Coulter CytoFLEX mosaic (“mosaic”)
- ThermoFisher Attune Xenith (“xenith”)

AutoSpectral will probably not solve all your unmixing issues; hopefully
it will help. Here are some other unmixing tools that may also help:

- [flowUnmix](https://github.com/hally166/flowUnmix)
- [PanelBuilder](https://github.com/exaexa/panelbuilder)
- [Ozette Resolve](https://www.ozette.com/)

What AutoSpectral cannot do:

- Magic. Blind spectral unmixing of 50-color panels is not currently
  implemented.
- Fix bad controls. If you have two colors in your “single-color”
  sample, it probably won’t work (see control cleaning, though). If you
  have no signal, it won’t work. If your controls were run on a
  different machine, with different instrument settings or using a
  different antibody/fluorophore, it won’t work well. Note that this is
  likely to affect ID7000 samples, where Sony suggests running the
  controls on one set of voltages and the samples on another. I haven’t
  put any time into figuring out how that works yet.
- Fix poor panel design. If you have unmixing spread from using BV785
  and BV786 in combination, AutoSpectral won’t fix that. If you’re using
  BB515 and FITC, it’ll probably help quite a bit, though.
- Fix non-specific staining. If you’re struggling with this, check out
  the blocking articles on the Colibri blog.
- Fix instrument errors. If your laser wasn’t working, this won’t help.
- Fix all unmixing errors. AutoSpectral currently has some ability to
  fix unmixing errors, provided the controls are a representation of the
  sample. You can see some examples of errors being fixed in the
  pre-print, but you can also see some examples of errors not being
  fixed in the deliberately poorly designed 42-color Aurora data
  example. This is something I’m still working on.

## Installation

[![Stable](https://img.shields.io/badge/stable-master-blue)](https://github.com/DrCytometer/AutoSpectral)
[![Dev](https://img.shields.io/badge/dev-branch-orange)](https://github.com/DrCytometer/AutoSpectral/tree/dev)

### Latest Stable Release

**Version 1.0.0**

Version 1.0.0 is intended to greatly speed up the process of unmixing.
This is done by pre-screening the variation in autofluorescence and
fluorophore variation for each cell, identifying likely “best”
candidates rather than using a brute force approach. If you are using
AutoSpectralRcpp, you will need to update that for compatibility.

To install the latest, hopefully stable version, install using
`devtools` or `remotes`. You will need to install the Bioconductor
packages separately. As of version 0.9.1, AutoSpectral relies on
`FlowSOM` rather than `EmbedSOM`, so you will need to install `FlowSOM`
via Bioconductor.

``` r
# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("flowWorkspace", "flowCore", "PeacoQC", "FlowSOM"))

# You'll need devtools or remotes to install from GitHub.
# install.packages("devtools")
# install.packages("remotes")
remotes::install_github("DrCytometer/AutoSpectral")
```

As of version 0.8.7, there is a Shiny helper tool to assist you in
setting up your AutoSpectral control files. This is an interactive html
app that opens in RStudio. Hopefully this makes things easier. It is
new, so again, imperfect. To try it, visit
[AutoSpectralHelper](https://github.com/DrCytometer/AutoSpectralHelper).
If you update AutoSpectral, I recommend downloading a new version of the
app to ensure compatibility.

To install a specific release, e.g., a previous one, use the version
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
[Development](https://drcytometer.github.io/AutoSpectral/articles/14_Development.html)
for suggestions of where help is needed most.

You can install the development version of AutoSpectral from
[GitHub](https://github.com/) with:

``` r
remotes::install_github("DrCytometer/AutoSpectral@dev")
```

## Bug fixes and known issues

AutoSpectral is pretty complex and newly released, so there will be
bugs. Sorry. Thanks to all of you providing feedback. Please submit any
and all issues either using the Issues page or via email at
colibri-cytometry at gmail.

At this point, many of the issues arising are specific to either the
data set or how it is being defined in the “control file”.
Troubleshooting these often requires a copy of the raw data, at least a
partial set of the unstained and single-stained controls.

All of the functions available in AutoSpectral can be viewed
[here](https://drcytometer.github.io/AutoSpectral/reference/index.html).

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
Workflow](https://drcytometer.github.io/AutoSpectral/articles/01_Full_AutoSpectral_Workflow.html).

All functions available in AutoSpectral can be viewed
[here](https://drcytometer.github.io/AutoSpectral/reference/index.html).

Resolved issues, bug patches and improvements will be announced via the
NEWS and also tracked in the [Updates and
Issues](https://drcytometer.github.io/AutoSpectral/articles/13_Updates_And_Issues.html)
article.

### Known shortcomings

- Gating. The automated gating is not great. See the [help
  page](https://drcytometer.github.io/AutoSpectral/articles/06_Gating.html)
  for tips. I’m working on an alternative.
- Please note that FCS 3.2 files from the S8 and A8 cytometers are not
  fully supported in flowCore. You may receive warnings, but things
  should still work.
- More stuff in progress will appear in the [Development
  article](https://drcytometer.github.io/AutoSpectral/articles/14_Development.html).
- This is my first R package.

If you want to use data from another cytometer and are wiling to provide
files for establishing the workflow, contact the author/maintainer. See
existing information, which may also assist you in setting up your
control file, in the [cytometer
database](https://docs.google.com/spreadsheets/d/1wj7QPkgpsuPNeVKyt-WWdBu5R48aZTgEbH8-_bpKeBY/edit?usp=sharing).

AutoSpectral relies on a database of information containing fluorophore
emission details. If your fluorophore is not detected automatically by
[`create.control.file()`](https://drcytometer.github.io/AutoSpectral/reference/create.control.file.md)
and you want to add it, visit the Google sheet for the [fluorophore
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

## Installation and Runtime

Installation via GitHub should take only a minute or so. It takes less
than that on a Dell i7 core laptop, but that might be because the
dependencies are already installed.

Occasionally, the help gets corrupted. Just re-install if that happens.
If you know why this happens, let me know.

Installation of `AutoSpectralRcpp` will take a couple of minutes because
the code needs to compile. You will also first have to install Rtools to
have a compiler, and that will take longer, probably 10 minutes or so.

For more details and benchmarking of specific functions with an example
40-color cell control data set, see the article on
[Speed](https://drcytometer.github.io/AutoSpectral/articles/12_Speed_It_Up.html).

## Go Faster

See the article on
“[Speed](https://drcytometer.github.io/AutoSpectral/articles/12_Speed_It_Up.html)
for details on how to improve AutoSpectral’s performance on your system.

## Updates and news

See the article [Updates and
Issues](https://drcytometer.github.io/AutoSpectral/articles/13_Updates_And_Issues.html)
for more on this, or read the NEWS with the latest release.
