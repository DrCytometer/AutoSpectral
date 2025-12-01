# Define Flow Control

A complex function designed to convert the single-stained control FCS
files and the metadata `control.def.file` into a data structure for
AutoSpectral. It reads the metadata file, locates the corresponding FCS
files, determines the gates needed based on variables such as beads or
cells, `large.gate`, and `viability.gate`, and then creates gates for
each combination. It imports and gates the data in the FCS files and
assigns various factors to track the data.

## Usage

``` r
define.flow.control(
  control.dir,
  control.def.file,
  asp,
  gate = TRUE,
  parallel = FALSE,
  verbose = TRUE
)
```

## Arguments

- control.dir:

  File path to the single-stained control FCS files.

- control.def.file:

  CSV file defining the single-color control file names, fluorophores
  they represent, marker names, peak channels, and gating requirements.

- asp:

  The AutoSpectral parameter list defined using
  `get.autospectral.param`.

- gate:

  Logical, default is `TRUE`, in which case, automated gating will be
  performed. If `FALSE`, the FCS files will be imported without
  automatically generated gates applied. That is, all data in the files
  will be used. This is intended to allow the user to pre-gate the files
  in commercial software.

- parallel:

  Logical, default is `FALSE`, in which case parallel processing will
  not be used. Parallel processing will likely be faster when many small
  files are read in. If the data is larger, parallel processing may not
  accelerate the process much.

- verbose:

  Logical, default is `TRUE`. Set to `FALSE` to suppress messages.

## Value

A list (`flow.control`) with the following components:

- `filename`: Names of the single-color control files.

- `fluorophore`: Corresponding fluorophores used in the experiment.

- `antigen`: Corresponding markers used in the experiment.

- `control.type`: Type of control used (beads or cells).

- `universal.negative`: Corresponding universal negative for each
  control.

- `viability`: Logical factor; whether a control represents a viability
  marker.

- `large.gate`: Logical factor; large gate setting.

- `autof.channel.idx`: Index of the autofluorescence marker channel.

- `event.number.width`: Width of the event number.

- `expr.data.max`: Maximum expression data value.

- `expr.data.max.ceil`: Ceiling of the maximum expression data value.

- `expr.data.min`: Minimum expression data value.

- `channel`: Preliminary peak channels for the fluorophores.

- `channel.n`: Number of channels.

- `spectral.channel`: Spectral channel information.

- `spectral.channel.n`: Number of spectral channels.

- `sample`: Sample names (fluorophores).

- `scatter.and.channel`: FSC, SSC, and peak channel information.

- `scatter.and.channel.label`: Labels for scatter and channel.

- `scatter.and.channel.spectral`: FSC, SSC, and spectral channels.

- `scatter.parameter`: Scatter parameters used for gating.

- `event`: Event factor.

- `event.n`: Number of events.

- `event.sample`: Sample information for events. Links events to
  samples.

- `event.type`: Type of events.

- `expr.data`: Expression data used for extracting spectra.
