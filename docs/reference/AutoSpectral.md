# AutoSpectral: Tools for Unmixing Spectral Flow Cytometry Data

As both a refinement and advancement of the unmixing process of full
spectrum cytometry data, `AutoSpectral` provides a suite of
functions/tools that work in concert to provide the user with optimal
unmixing results.

Use `AutoSpectral` to:

- isolate and refine clean spectral signatures

- mitigate autofluorescence contamination in complex controls

- unmix fcs files using standard algorithms or the hallmark
  `AutoSpectral` approach

## `AutoSpectral` (v 1.5.7) Workflow

To maximize the success of `AutoSpectral`, the workflow is organized
into the following logical step-wise processes:

1.  [create.control.file](https://drcytometer.github.io/AutoSpectral/reference/create.control.file.md):
    Generate a .csv file that contains descriptive information about
    your single stained controls

    - [Control
      File](https://drcytometer.github.io/AutoSpectral/articles/01_Full_AutoSpectral_Workflow.html#creating-the-control-file)
      (see article)

2.  [define.flow.control](https://drcytometer.github.io/AutoSpectral/reference/define.flow.control.md):
    uses the result of
    [create.control.file](https://drcytometer.github.io/AutoSpectral/reference/create.control.file.md)
    to convert FCS files and their associated metadata into an optimized
    data structure for downstream workflows

    - [Flow
      Control](https://drcytometer.github.io/AutoSpectral/articles/01_Full_AutoSpectral_Workflow.html#loading-the-data)
      (see article)

3.  [clean.controls](https://drcytometer.github.io/AutoSpectral/reference/clean.controls.md):
    removes autofluorescent events from cell-based controls via
    PCA-based gating, selects top-expressing positive events, matches
    universal negatives by scatter, and downsamples for speed

4.  [get.fluorophore.spectra](https://drcytometer.github.io/AutoSpectral/reference/get.fluorophore.spectra.md):
    extracts normalized `[0, 1]` fluorophore reference spectra from the
    cleaned controls and performs cosine-similarity QC against the
    spectral reference library

5.  *(Optional)*
    [get.af.spectra](https://drcytometer.github.io/AutoSpectral/reference/get.af.spectra.md):
    extracts autofluorescence reference spectra for use with the
    AutoSpectral unmixing method

6.  *(Optional)*
    [get.spectral.variants](https://drcytometer.github.io/AutoSpectral/reference/get.spectral.variants.md):
    computes per-fluorophore spectral variants for per-cell fluorophore
    optimization during AutoSpectral unmixing

7.  [unmix.fcs](https://drcytometer.github.io/AutoSpectral/reference/unmix.fcs.md)
    /
    [unmix.folder](https://drcytometer.github.io/AutoSpectral/reference/unmix.folder.md):
    unmixes experiment FCS files using the extracted spectra. Supported
    methods are `"OLS"`, `"WLS"`, `"Poisson"`, and `"AutoSpectral"`

## Automated alternative (steps 2–4)

[get.spectra.automated](https://drcytometer.github.io/AutoSpectral/reference/get.spectra.automated.md)
replaces the `define.flow.control` → `clean.controls` →
`get.fluorophore.spectra` pipeline with a single function call,
requiring no scatter gating or interactive input.

## See also

Useful links:

- <https://github.com/DrCytometer/AutoSpectral>

- <https://drcytometer.github.io/AutoSpectral/>

- Report bugs at <https://github.com/DrCytometer/AutoSpectral/issues>

## Author

**Maintainer**: Oliver Burton <olivertburton@gmail.com>
([ORCID](https://orcid.org/0000-0003-3884-7373))

Authors:

- Oliver Burton <olivertburton@gmail.com>
  ([ORCID](https://orcid.org/0000-0003-3884-7373))

- Adrian Liston <al989@cam.ac.uk>
  ([ORCID](https://orcid.org/0000-0002-6272-4085))

Other contributors:

- Nathan Laniewski <Nathan_Laniewski@urmc.rochester.edu> \[contributor\]
