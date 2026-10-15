
<!-- README.md is generated from README.Rmd. Please edit that file -->

# AutoSpectral

<!-- badges: start -->

<!-- badges: end -->

## Introduction

AutoSpectral is AutoSpill (Roca et al., 2021) updated for the spectral
flow era.

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

- Encode Roederer’s Rules for controls (Roederer, 2002; Ferrer-Font et
  al., 2021) as much as possible.
- Deal with cellular autofluorescence interfering with identifying
  fluorophore signatures.
- Autofluorescence is variable on a cell-to-cell basis, both in
  magnitude and in type. We can figure out, more or less, what this
  should be for each cell. Manually defining multiple autofluorescence
  signatures is the established approach to this problem (Niewold et
  al., 2020; Wanner et al., 2022; Jameson et al., 2022; Roet et al.,
  2024; Pilkington, 2024); AutoSpectral automates the choice of
  signature for each cell.
- Fluorophore emissions are variable, and this variability manifests on
  the level of the cell. Again, we can more or less figure this out and
  deal with it. Spectral collinearity between fluorophores amplifies the
  resulting errors (Novo, 2022; Kmet & Novo, 2025; Mage et al., 2026).
- By formalizing the unmixing workflow, we make it reproducible and
  available to users of any experience level (although this is currently
  limited by the R interface).

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
- ChallenBio CytoStellar (“cytostellar”)

If you have data from a different instrument and are willing to share
files, contact the author/maintainer.

AutoSpectral will probably not solve all your unmixing issues; hopefully
it will help. Here are some other tools that may also help:

- [CompensAID](https://bioc.r-universe.dev/articles/CompensAID/CompensAID.html)
  (Olsman et al., 2026)
- [spectracle](https://github.com/nlaniewski/spectracle/)
- [flowstate](https://github.com/nlaniewski/flowstate)
- [spectreasy](https://github.com/pkheisig/spectreasy)
- [Ozette Resolve](https://www.ozette.com/)
- [flowUnmix](https://github.com/hally166/flowUnmix)
- [PanelBuilder](https://github.com/exaexa/panelbuilder)

What AutoSpectral cannot do:

- Magic. Blind spectral unmixing of 50-color panels is not currently
  implemented.
- Fix bad controls. If you have two colors in your “single-color”
  sample, it probably won’t work (see control cleaning, though). If you
  have no signal, it won’t work. If your controls were run on a
  different machine, that won’t work. If your controls were acquired
  with different instrument settings or using a different
  antibody/fluorophore, it may not work well. Note that this is likely
  to affect ID7000 samples, where Sony suggests running the controls on
  one set of voltages and the samples on another. I haven’t put any time
  into figuring out how that works yet. It may well be that the unmixing
  with AutoSpectral works fine in this situation.
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
- Adjust for running samples on different voltages/gains. Sony
  recommends running multi-colour samples on potentially very different
  voltage settings from the single-colour controls. This means the
  spectral profiles of your fluorophores will not be accurate except on
  in Sony’s software, where you have access to all the QC, voltage
  settings and PMT linearity curves needed to make the adjustment. None
  of this information is in the FCS files, so AutoSpectral can’t really
  do the same adjustment process. This is broadly true for other
  instruments as well–if you run your controls once a month and run
  samples later in the month without any new controls, your unmixing may
  not be very good.

## Installation

[![Stable](https://img.shields.io/badge/stable-master-blue)](https://github.com/DrCytometer/AutoSpectral)
[![Dev](https://img.shields.io/badge/dev-branch-orange)](https://github.com/DrCytometer/AutoSpectral/tree/dev)

### Latest Stable Release

**Version 1.8.7**

Version 1.6+ provides a new automated approach to extracting spectral
profiles. Most of the development here has been by Nathan Laniewksi.
This is fast, highly reproducible, clean and avoids all the issues with
gating. Try it out with `get.spectra.automated()`.

Version 1.6+ also provides a new unmixing pipeline for per-cell
fluorophore optimization to eliminate discontinuities that appeared in
some data sets. To use this approach, call `unmix.fcs()` or
`unmix.folder()` with `pipeline = "joint"`. Otherwise the unmixing will
default to the original approach.

This has been tested on Windows 11 and on macOS Sequoia 15.2 (Apple
Silicon).

Version 1.0+ is intended to greatly speed up the process of unmixing.
This is done by pre-screening the variation in autofluorescence and
fluorophores for each cell, identifying likely “best” candidates rather
than using a brute force approach. If you are using
[AutoSpectralRcpp](https://github.com/DrCytometer/AutoSpectralRcpp), you
will need to update that for compatibility.

To install the latest, hopefully stable version, install using `pak` or
`remotes`. As of version 1.7.0, AutoSpectral relies primarily on
self-organizing maps (SOM; Kohonen, 1995) via an in-house port of the
batched parallel SOM from `EmbedSOM` ([PubMed
32518625](https://pubmed.ncbi.nlm.nih.gov/32518625/)) via
`AutoSpectralRcpp`. If you do not install `AutoSpectralRcpp`, you will
need to install `FlowSOM` (Van Gassen et al., 2015) via Bioconductor or
using `pak`.

``` r
# You'll need pak or remotes to install from GitHub.
install.packages("pak")
pak::pak("DrCytometer/AutoSpectral")

# recommended
pak::pak("DrCytometer/AutoSpectralRcpp")

# required if you don't install AutoSpectralRcpp
pak:pak("FlowSOM")
```

Installation should take less than 5 minutes.

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
pak::pak("DrCytometer/AutoSpectral@v1.5.6")
```

### Dev branch

New features and patches will be implemented first in `dev` branch. So,
at any given point, this may not be working well and is not advised for
general use.

AutoSpectral is open source. If you are interested in contributing,
please visit
[Development](https://drcytometer.github.io/AutoSpectral/articles/16_Development.html)
for suggestions of where help is needed most.

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

Resolved issues, bug patches and improvements will be announced via the
NEWS, in
[Discussions](https://github.com/DrCytometer/AutoSpectral/discussions/24),
and also tracked in the [Updates and
Issues](https://drcytometer.github.io/AutoSpectral/articles/15_Updates_And_Issues.html)
article.

### Known shortcomings

- Not all unmixing errors will be fixed. Still working on this.
- AutoSpectral may not work well for unmixing of samples run on
  different PMT/APD voltage settings than what were used for the
  controls.
- Data from the A5SE may be more problematic–I’m looking into this.
- This would be better with a graphical interface.
- This is my first R package.
- R is not exactly a fast or modern programming language.
- The documentation is probably not amazingly clear. You can help with
  this by highlighting points that are unclear, or, even better,
  providing example workflows and explanations for other users.

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

## Installation and Runtime

Installation via GitHub should take only a minute or so. It takes less
than that on a Dell i7 core laptop, but that might be because the
dependencies are already installed.

Occasionally, the help gets corrupted. Just re-install if that happens.
If you know why this happens, let me know.

Installation of `AutoSpectralRcpp` will take a couple of minutes because
the code needs to compile. A C++ compiler is required; what you need
depends on your platform:

- **Windows**: Install
  [Rtools](https://cran.r-project.org/bin/windows/Rtools/). This
  provides the compiler and takes around 10 minutes to download and
  install. OpenMP is included and enabled automatically.
- **macOS**: Install the Xcode Command Line Tools
  (`xcode-select --install` in Terminal). OpenMP is not included in
  Apple’s compiler, but `AutoSpectralRcpp` will detect whether
  [Homebrew](https://brew.sh) and `libomp` are available and enable it
  automatically if so. If not, it will build without OpenMP
  (single-threaded C++ only) and tell you what to run. To enable OpenMP
  manually: `brew install libomp`, then reinstall the package.
- **Linux**: GCC is typically already present. OpenMP is enabled
  automatically.

For more details and benchmarking of specific functions with an example
42-color cell control data set, see the article on
[Speed](https://drcytometer.github.io/AutoSpectral/articles/14_Speed_It_Up.html).

## Go Faster

See the article on
[Speed](https://drcytometer.github.io/AutoSpectral/articles/14_Speed_It_Up.html)
for details on how to improve AutoSpectral’s performance on your system.

## Updates and news

See the article [Updates and
Issues](https://drcytometer.github.io/AutoSpectral/articles/15_Updates_And_Issues.html)
for more on this, or read the NEWS with the latest release.

## Citation

If you use AutoSpectral, please cite the manuscript:

Burton OT, Laniewski N, Rogers SS, Lopez Munoz-de-Leon A, Khan S,
Arunthavalingam J, Smith E, Bücken L, De Vuyst L, Humblet-Baron S,
Cerveira J, Scheible K, Dooley J, Liston A (2025). AutoSpectral improves
spectral flow cytometry accuracy through optimised spectral unmixing and
autofluorescence-matching at the cellular level. *bioRxiv* [doi:
10.1101/2025.10.27.684855](https://doi.org/10.1101/2025.10.27.684855)

The automated gating and robust linear modelling used to extract spectra
from single-stained controls were developed in AutoSpill (Roca et al.,
2021), so please cite that paper too.

## References

- Ferrer-Font L, Small SJ, Lewer B, Pilkington KR, Johnston LK, Park LM,
  Lannigan J, Jaimes MC, Price KM (2021). Panel optimization for
  high-dimensional immunophenotyping assays using full-spectrum flow
  cytometry. *Curr Protoc* 1(9):e222.
  [doi:10.1002/cpz1.222](https://doi.org/10.1002/cpz1.222)
- Jameson VJ, Luke T, Yan Y, Hind A, Evrard M, Man K, Mackay LK, Kallies
  A, Villadangos JA, McWilliam HEG, Perez-Gonzalez A (2022). Unlocking
  autofluorescence in the era of full spectrum analysis: Implications
  for immunophenotype discovery projects. *Cytometry A* 101(11):922-941.
  [doi:10.1002/cyto.a.24555](https://doi.org/10.1002/cyto.a.24555)
- Kmet R, Novo D (2025). Reducing spreading: removing the impact of
  irrelevant dyes improves unmixed flow cytometry data. *Cytometry A*.
  [doi:10.1002/cyto.a.24957](https://doi.org/10.1002/cyto.a.24957)
- Kohonen T (1995). *Self-organizing maps*. Springer.
- Mage PL, Konecny AJ, Mair F (2026). Measurement and prediction of
  unmixing-dependent spreading due to collinearity in spectral flow
  cytometry. *Cytometry A* 109(6):436-455.
  [doi:10.1002/cyto.a.70044](https://doi.org/10.1002/cyto.a.70044)
- Niewold P, Ashhurst TM, Smith AL, King NJC (2020). Evaluating spectral
  cytometry for immune profiling in viral disease. *Cytometry A*
  97(11):1165-1179.
  [doi:10.1002/cyto.a.24211](https://doi.org/10.1002/cyto.a.24211)
- Novo D (2022). A comparison of spectral unmixing to conventional
  compensation for the calculation of fluorochrome abundances from flow
  cytometric data. *Cytometry A* 101(11):885-891.
  [doi:10.1002/cyto.a.24669](https://doi.org/10.1002/cyto.a.24669)
- Olsman R, Bonte S, Hofmans M, Buysse M, Van der Borght K, Saeys Y, van
  der Velden VHJ, Van Gassen S (2026). CompensAID: An automated
  detection tool for reference errors. *Cytometry A* 109(2):124-134.
  [doi:10.1002/cyto.a.70016](https://doi.org/10.1002/cyto.a.70016)
- Pilkington KR (2024). Autofluorescence: From burden to benefit.
  *Cytometry A* 105(8):563-567.
  [doi:10.1002/cyto.a.24885](https://doi.org/10.1002/cyto.a.24885)
- Roca CP, Burton OT, Gergelits V, Prezzemolo T, Whyte CE, Halpert R,
  Kreft L, Collier J, Botzki A, Spidlen J, Humblet-Baron S, Liston A
  (2021). AutoSpill is a principled framework that simplifies the
  analysis of multichromatic flow cytometry data. *Nat Commun*
  12(1):2890.
  [doi:10.1038/s41467-021-23126-8](https://doi.org/10.1038/s41467-021-23126-8)
- Roederer M (2002). Compensation in flow cytometry. *Curr Protoc Cytom*
  Chapter 1:Unit 1.14.
  [doi:10.1002/0471142956.cy0114s22](https://doi.org/10.1002/0471142956.cy0114s22)
- Roet JEG, Mikula AM, de Kok M, Chadick CH, Garcia Vallejo JJ, Roest
  HP, van der Laan LJW, de Winde CM, Mebius RE (2024). Unbiased method
  for spectral analysis of cells with great diversity of
  autofluorescence spectra. *Cytometry A* 105(8):595-606.
  [doi:10.1002/cyto.a.24856](https://doi.org/10.1002/cyto.a.24856)
- Van Gassen S, Callebaut B, Van Helden MJ, Lambrecht BN, Demeester P,
  Dhaene T, Saeys Y (2015). FlowSOM: Using self-organizing maps for
  visualization and interpretation of cytometry data. *Cytometry A*
  87(7):636-645.
  [doi:10.1002/cyto.a.22625](https://doi.org/10.1002/cyto.a.22625)
- Wanner N, Barnhart J, Apostolakis N, Zlojutro V, Asosingh K (2022).
  Using the Autofluorescence Finder on the Sony ID7000 Spectral Cell
  Analyzer to identify and unmix multiple highly autofluorescent murine
  lung populations. *Front Bioeng Biotechnol* 10:827987.
  [doi:10.3389/fbioe.2022.827987](https://doi.org/10.3389/fbioe.2022.827987)
- Nguyen R, Perfetto S, Mahnke YD, Chattopadhyay P, Roederer M (2013).
  Quantifying spillover spreading for comparing instrument performance
  and aiding in multicolor panel design. *Cytometry A* 83\*3(:306-15).
  [doi: 10.1002/cyto.a.22251](https://doi.org/10.1002/cyto.a.22251)
- Cai X, Garcia-Garcia S, Kuhnen L, Gianniou M, Garcia Vallejo JJ
  (2026). Unmixing Spread Estimation Based on Residual Model in Spectral
  Flow Cytometry. *bioRxiv* (pre-print). [doi:
  10.64898/2026.01.27.701929](https://doi.org/10.64898/2026.01.27.701929)
