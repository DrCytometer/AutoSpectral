# Define Gate by Landmarks

Performs gating on scatter parameters using events selected on the basis
of high fluorescence intensity in the expected peak channel in the
single-stained control samples. This approach allows for identification
of cells fairly reliably by "landmarks", especially when using
well-defined, abundant populations such a T cells (e.g., CD3, CD4, CD8
single-stained samples), monocytes (CD14-stained sample) or neutrophils
(CD66, CD33).

## Usage

``` r
define.gate.landmarks(
  control.file,
  control.dir,
  asp,
  gate.name,
  gating.params = NULL,
  n.cells = 2000,
  percentile = 70,
  grid.n = 100,
  bandwidth.factor = 1,
  fsc.channel = NULL,
  ssc.channel = NULL,
  fsc.lims = NULL,
  ssc.lims = NULL,
  output.dir = "./figure_gate",
  filename = "landmark_gate_definition_",
  color.palette = "plasma",
  boundary.color = "black",
  points.to.plot = 1e+05,
  width = 5,
  height = 5,
  verbose = TRUE,
  control.table = NULL,
  check = TRUE
)
```

## Arguments

- control.file:

  File path and name for the CSV file defining the single- color control
  file names, fluorophores they represent, marker names, peak channels,
  and gating requirements.

- control.dir:

  File path to the single-stained control FCS files.

- asp:

  The AutoSpectral parameter list defined using
  `get.autospectral.param`.

- gate.name:

  Character, name for the gate. Useful for distinguishing gates when you
  have multiple types. Must match one string (name) in the `gate.name`
  column of the `control.file`.

- gating.params:

  Previously saved gating parameters. Load in the .rds file using
  `readRDS` and pass the result here if you wish to replicate a previous
  run. This is essentially just an updated version of `asp` containing
  any modifications based on information passed as arguments to
  `define.gate.landmarks`.

- n.cells:

  The number of cells to use for defining the gate boundary. This many
  cells will be selected from the peak channel (brightest first) in the
  single-color controls. For example, if you set `200` and marked files
  such as `CD3-PE.fcs` and `CD19-FITC.fcs` as `gate.define=TRUE` in the
  control file, the brightest 200 events in the YG1 channel from the
  CD3-PE file and the brightest 200 events in the B1 channel for the
  CD19-FITC file would be used to define the gate.

- percentile:

  Numeric 1 - 100, default `70`. The percentile cutoff for density in
  the scatter to use for defining the gate. For example, a value of `50`
  would take the 50% of cells closest to the density peak, more or less.
  Smaller numbers will define a tighter gate.

- grid.n:

  Numeric, default `100`. The binning grid for the kernel density
  estimation. If `n.cells` is very low, you may wish to lower this
  number to compress the search space.

- bandwidth.factor:

  Numeric, default `1`. A multiplier for the bandwidth for the kernel
  density estimation. Larger numbers will smooth the density, reducing
  discrimination between peaks in the density (such as between live and
  dead cells).

- fsc.channel:

  Channel to use for Forward Scatter. Default `NULL` will use the
  `asp$default.scatter.parameter[1]`, which is appropriate for your
  machine.

- ssc.channel:

  Channel to use for Side Scatter. Default `NULL` will use the
  `asp$default.scatter.parameter[2]`, which is appropriate for your
  machine. For machines with multiple side scatter measurements, you can
  change this.

- fsc.lims:

  Numeric vector. Limits for searching and plotting the FSC. The default
  `NULL` uses c`c(asp$scatter.data.min.x, asp$scatter.data.max.x)`.

- ssc.lims:

  Numeric vector. Limits for searching and plotting the SSC. The default
  `NULL` uses c`c(asp$scatter.data.min.y, asp$scatter.data.max.y)`.

- output.dir:

  File path where you want to save the results. Default is
  `./figure_gate`.

- filename:

  Character, name for the output files. Default is `gate_definition`.

- color.palette:

  Optional character string defining the viridis color palette to be
  used for the fluorophore traces. Default is `plasma`. Use `rainbow` to
  be similar to FlowJo or SpectroFlo. Other options are the viridis
  color options: `magma`, `inferno`, `plasma`, `viridis`, `cividis`,
  `rocket`, `mako` and `turbo`.

- boundary.color:

  Color for the gate boundary line on the plot. Default is `black`.

- points.to.plot:

  Numeric, default `1e5`. Maximum number of points to show on the plot.
  More points will take longer, but you really shouldn't have even close
  to this number when defining the gate with landmarks.

- width:

  Numeric, default `5`. Width of the saved plot.

- height:

  Numeric, default `5`. Height of the saved plot.

- verbose:

  Logical, default is `TRUE`. Set to `FALSE` to suppress messages.

- control.table:

  Dataframe or table of the control file, read in via
  [`define.flow.control()`](https://drcytometer.github.io/AutoSpectral/reference/define.flow.control.md)
  and cleaned up by that function.

- check:

  Logical, default is `TRUE`. Set to `FALSE` to skip consistency checks
  on gate defining columns in the control table: `large.gate`,
  `is.viability` and `sample.type`.

## Value

A set of points describing the gate boundary.

## References

Laniewski, Nathan. *flowstate*.
<https://github.com/nlaniewski/flowstate>

## See also

- [`tune.gate()`](https://drcytometer.github.io/AutoSpectral/reference/tune.gate.md)

- [`define.gate.density()`](https://drcytometer.github.io/AutoSpectral/reference/define.gate.density.md)

- [`gate.define.plot()`](https://drcytometer.github.io/AutoSpectral/reference/gate.define.plot.md)
