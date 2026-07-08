# Get Fluorophore Spectra - Automated Workflow

Single-function replacement for the three-step `define.flow.control` -\>
`clean.controls` -\> `get.fluorophore.spectra` pipeline. Extracts
normalised `[0, 1]` fluorophore reference spectra from a directory of
single-stained control FCS files without requiring scatter gating or
interactive input.

The algorithm is based on the `spectracle` approach, adapted to use
AutoSpectral's internal FCS reader, matrix-based data layout, and the
control-file metadata already required by the rest of the AutoSpectral
pipeline:

1.  Reads all single-stained controls with
    [`readFCS()`](https://drcytometer.github.io/AutoSpectral/reference/readFCS.md).

2.  Removes saturating events. \#' 3. For each fluorophore, determines
    the AF reference from the paired universal-negative file specified
    in `universal.negative` column of the control table, or - if that
    column is empty - uses internal-negative mode.

3.  **External-negative:** projection-based AF orthogonalisation to
    identify the empirical peak detector; cosine-similarity filter to
    select the least AF-like events; KNN scatter-matched AF subtraction;
    column median summary. **Internal-negative:** top-100 events by
    expected peak channel minus the mean of the bottom-10%; row mean
    summary; no cosine filter or kNN step.

4.  Applies QC: if the normalised signal cosine similarity against the
    spectral reference library is below `cosine.threshold`, the spectrum
    is refined using the legacy gating and cleaning approach.

5.  Returns the spectra matrix in AutoSpectral format (fluorophores x
    detectors).

The legacy workflow (`define.flow.control` + `clean.controls` +
`get.fluorophore.spectra`) remains available and is unchanged.

## Usage

``` r
get.spectra.automated(
  control.dir,
  control.def.file,
  asp,
  n.candidates = 1000L,
  n.spectral = 200L,
  k.neighbors = 2L,
  singlet.quantiles = c(0.85, 0.975),
  cosine.threshold = 0.9,
  peak.signal.threshold = 0.5,
  legacy.refinement = TRUE,
  top.expressing.override = NULL,
  return.af = FALSE,
  figures = TRUE,
  plot.cosine.filter = TRUE,
  plot.scatter.match = TRUE,
  verbose = TRUE
)
```

## Arguments

- control.dir:

  Character. Path to the directory containing the single-stained control
  FCS files.

- control.def.file:

  Character. Path to (or filename of) the control definition CSV. Must
  already exist and pass
  [`check.control.file()`](https://drcytometer.github.io/AutoSpectral/reference/check.control.file.md).

- asp:

  The AutoSpectral parameter list from
  [`get.autospectral.param()`](https://drcytometer.github.io/AutoSpectral/reference/get.autospectral.param.md).

- n.candidates:

  Integer, default `1000`. Number of top-expressing candidate events
  selected per fluorophore before cosine-similarity filtering. Ignored
  in internal-negative mode, where the top 5%% of events by peak channel
  are used directly.

- n.spectral:

  Integer, default `200`. Number of spectral events retained after
  filtering for low AF cosine similarity.

- k.neighbors:

  Integer, default `2`. Number of nearest neighbours in scatter space
  used for per-event AF subtraction.

- singlet.quantiles:

  Numeric, default `c( 0.85, 0.975 )`. Range of quantiles to use for
  doublet discrimination.

- cosine.threshold:

  Numeric, default `0.9`. Minimum cosine similarity against the spectral
  reference library to accept the automated spectrum; values below this
  threshold trigger legacy pipeline refinement.

- peak.signal.threshold:

  Numeric, default `0.5`. Minimum normalised signal in the expected peak
  detector. Used only for informational QC.

- legacy.refinement:

  Logical, default `TRUE`. Whether to run legacy pipeline on controls
  where the signature does not match the reference well.

- top.expressing.override:

  Named numeric vector or `NULL` (default). Override the event count for
  specific samples. Names should match the FCS filename (without
  extension) or the fluorophore name from the control file.

- return.af:

  Logical, default `FALSE`. Adds a single mean autofluorescence spectrum
  to the output to allow unmixing with OLS or WLS using autofluorescence
  extraction. Requires naming the unstained cell control sample as `AF`
  for the fluorophore in the control file.

- figures:

  Logical, default `TRUE`. Produce spectral trace, heatmap,
  cosine-similarity, and reference-library QC plots.

- plot.cosine.filter:

  Logical, default `TRUE`. When `figures = TRUE`, produce a multi-panel
  PDF showing the per-fluorophore cosine-similarity filter traces
  (binned by cosine similarity to AF, coloured from least to most
  AF-like). Set to `FALSE` to skip.

- plot.scatter.match:

  Logical, default `TRUE`. When `figures = TRUE`, produce a multi-panel
  PDF showing the KNN scatter-matched AF subtraction for each
  fluorophore (unstained background, selected spectral events, and their
  matched AF events). Set to `FALSE` to skip.

- verbose:

  Logical, default `TRUE`. Print progress messages.

## Value

A numeric matrix with fluorophores in rows and spectral detector
channels in columns, values normalised to `[0, 1]` (L-infinity norm,
peak = 1). Compatible with all downstream AutoSpectral functions.

## See also

[`get.fluorophore.spectra()`](https://drcytometer.github.io/AutoSpectral/reference/get.fluorophore.spectra.md)
for the legacy workflow.
[`spectral.reference.plot()`](https://drcytometer.github.io/AutoSpectral/reference/spectral.reference.plot.md)
for the QC report produced when `figures = TRUE`.
