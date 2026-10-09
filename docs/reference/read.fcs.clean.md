# Read And Clean An FCS File

Reads the scatter and spectral channels of an FCS file and removes
events that saturate any spectral detector, events at the scatter maxima
(`asp$scatter.data.max.x`/`asp$scatter.data.max.y`, when set) and,
optionally, doublets. Doublets are removed with a two-pass Area/Height
scatter-ratio gate: the FSC ratio cutoff first, then the SSC ratio
cutoff computed only on events that passed the FSC cutoff. This is the
cleaning applied to every control file by
[`get.spectra.automated()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectra.automated.md),
exported so that other packages can pre-process related samples
identically.

## Usage

``` r
read.fcs.clean(
  path,
  label,
  spectral.channels,
  scatter.channels,
  sat.value,
  singlet.quantiles,
  remove.doublets = TRUE,
  asp,
  verbose = TRUE
)
```

## Arguments

- path:

  File name and path of the FCS file.

- label:

  Character, name used in the progress message.

- spectral.channels:

  Character vector of spectral detector names.

- scatter.channels:

  Character vector of scatter channel names (Area). Matching Height
  channels are read, if present, for doublet removal.

- sat.value:

  Numeric, the detector saturation value. Events with any spectral
  channel at or above this value are removed. Typically
  `asp$expr.data.max`.

- singlet.quantiles:

  Numeric, length 2. Quantile cutoffs for the FSC and SSC Area/Height
  ratios, e.g. `c( 0.85, 0.975 )`.

- remove.doublets:

  Logical, default `TRUE`. Skipped automatically if Height channels are
  not present in the file.

- asp:

  The AutoSpectral parameter list. Prepare using
  `get.autospectral.param`.

- verbose:

  Logical, default `TRUE`.

## Value

Numeric matrix of retained events, with the scatter channels and
spectral channels present in the file (Height channels dropped).
