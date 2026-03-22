# Define Flow Control

A complex function designed to convert the single-stained control FCS
files and the metadata `control.def.file` into a data structure for
AutoSpectral. It reads the metadata file, locates the corresponding FCS
files, determines the gates needed based on variables such as beads or
cells, `large.gate`, and `viability.gate`, and then creates gates for
each combination. It imports and gates the data in the FCS files and
assigns various factors to track the data. Parallel processing
`parallel=TRUE` will likely speed up the run considerably on Mac and
Linux systems supporting forking, but will likely not help much on
Windows unless \>10 cores are available.

## Usage

``` r
define.flow.control(
  control.dir,
  control.def.file,
  asp,
  gate = TRUE,
  gating.system = c("landmarks", "density"),
  gate.list = NULL,
  parallel = FALSE,
  verbose = TRUE,
  threads = NULL,
  color.palette = NULL
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

- gating.system:

  Character string selecting the automated gating system to employ in
  defining the initial scatter gates for identifying cells in the FCS
  files. Options are `landmarks` and `density`. The `density` option
  uses the original gating system from AutoSpill, picking cell
  populations based on dense regions on FSC/SSC. The default `landmarks`
  system picks out the cell regions using the brightest events in the
  peak channel for the single-stained control(s). This approach is
  generally more robust. A good way to use this is to utilize
  `landmarks` in combination with specifying which single-stained
  control files should be used for defining the gates. For example,
  abundant, bright cell markers such as CD4 will reliably identify the
  lymphocyte region, as CD14 will identify the monocyte region, etc. For
  more instructions, see the help pages on GitHub.

- gate.list:

  Optional named list of gates. To use this, pre-define the gates using
  [`define.gate.landmarks()`](https://drcytometer.github.io/AutoSpectral/reference/define.gate.landmarks.md)
  and/or
  [`define.gate.density()`](https://drcytometer.github.io/AutoSpectral/reference/define.gate.density.md),
  ensure that the names of the gates correspond to the names in the
  `control.def.file`, and ensure that the `gate.name` column has been
  filled in for the `control.def.file`. Default `NULL` will revert to
  creating new gates.

- parallel:

  Logical, default is `FALSE`, in which case parallel processing will
  not be used. Parallel processing will likely be faster when many small
  files are read in. If the data is larger, parallel processing may not
  accelerate the process much.

- verbose:

  Logical, default is `TRUE`. Set to `FALSE` to suppress messages.

- threads:

  Numeric, number of threads to use for parallel processing. Default is
  `NULL` which will revert to `asp$worker.process.n` if `parallel=TRUE`.

- color.palette:

  Optional character string defining the viridis color palette to be
  used for the fluorophore traces. Use `rainbow` to be similar to FlowJo
  or SpectroFlo. Other options are the viridis color options: `magma`,
  `inferno`, `plasma`, `viridis`, `cividis`, `rocket`, `mako` and
  `turbo`.

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

## References

Roca, Carlos P et al. "AutoSpill is a principled framework that
simplifies the analysis of multichromatic flow cytometry data" *Nature
Communications* 12 (2890) 2021.

## See also

- [`tune.gate()`](https://drcytometer.github.io/AutoSpectral/reference/tune.gate.md)

- [`define.gate.landmarks()`](https://drcytometer.github.io/AutoSpectral/reference/define.gate.landmarks.md)

- [`define.gate.density()`](https://drcytometer.github.io/AutoSpectral/reference/define.gate.density.md)

- [`do.gate()`](https://drcytometer.github.io/AutoSpectral/reference/do.gate.md)

- [`gate.define.plot()`](https://drcytometer.github.io/AutoSpectral/reference/gate.define.plot.md)
