# Define Gate by Density

Performs gating on scatter parameters based on a density search.

The gating proceeds in three steps:

- Defines bounds by data trimming

- Defines a region around the target maximum found within the bounds

- Defines a gate around the target maximum, only within that region

The method uses numerical search of maxima over estimated densities and
Voronoi tessellations to improve density estimation around maxima.

## Usage

``` r
define.gate.density(
  control.file,
  control.dir,
  asp,
  gate.name,
  gating.params = NULL,
  n.cells = NULL,
  grid.n = NULL,
  bandwidth.factor = NULL,
  target.pop = NULL,
  neighbors = NULL,
  fsc.channel = NULL,
  ssc.channel = NULL,
  fsc.lims = NULL,
  ssc.lims = NULL,
  fsc.search.min = NULL,
  fsc.search.max = NULL,
  ssc.search.min = NULL,
  ssc.search.max = NULL,
  output.dir = "./figure_gate",
  filename = "density_gate_definition_",
  color.palette = "plasma",
  boundary.color = "black",
  points.to.plot = 1e+05,
  width = 5,
  height = 5,
  verbose = TRUE,
  control.table = NULL
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
  cells will be selected from the specified FCS files. The default,
  `NULL`, will revert to `asp$gate.downsample.n`. Around 1e5 is
  recommended.

- grid.n:

  Numeric, default `100`. The binning grid for the kernel density
  estimation. If `n.cells` is very low, you may wish to lower this
  number to compress the search space.

- bandwidth.factor:

  Numeric, default `NULL` will revert to
  `asp$gate.bound.density.bw.factor` (usually 1). A multiplier for the
  bandwidth for the kernel density estimation. Larger numbers will
  smooth the density, reducing discrimination between peaks in the
  density (such as between live and dead cells).

- target.pop:

  Numeric 0-n, default `NULL` will revert to
  `asp$gate.bound.density.max.target` (usually 1). In the output plots,
  the numbers correspond to density peaks. Whichever number is set here
  will select the population 1 less than that number as the target for
  defining the gate region.

- neighbors:

  Number of neighboring events to consider when finding local density
  maxima in the kernel density estimation. Default `NULL` will revert to
  `asp$gate.bound.density.neigh.size`, which is usually 3.

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

- fsc.search.min:

  The minimum bound for the search space for dense event populations for
  defining the gate on the x-axis. The default `NULL` reverts to
  `asp$default.gate.param$region.factor.x.low`, which is usually 0.05 or
  5% of the FSC axis.

- fsc.search.max:

  The maximum bound for the search space for dense event populations for
  defining the gate on the x-axis. The default `NULL` reverts to
  `asp$default.gate.param$region.factor.x.high`, which is usually 0.8 or
  80% of the FSC axis.

- ssc.search.min:

  The minimum bound for the search space for dense event populations for
  defining the gate on the y-axis. The default `NULL` reverts to
  `asp$default.gate.param$region.factor.y.low`, which is usually 0.05 or
  5% of the SSC axis.

- ssc.search.max:

  The maximum bound for the search space for dense event populations for
  defining the gate on the y-axis. The default `NULL` reverts to
  `asp$default.gate.param$region.factor.y.high`, which is usually 0.8 or
  80% of the SSC axis.#'

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

## Value

A set of points describing the gate boundary.

## References

Roca, Carlos P et al. "AutoSpill is a principled framework that
simplifies the analysis of multichromatic flow cytometry data" *Nature
Communications* 12 (2890) 2021.

## See also

- [`tune.gate()`](https://drcytometer.github.io/AutoSpectral/reference/tune.gate.md)

- [`define.gate.landmarks()`](https://drcytometer.github.io/AutoSpectral/reference/define.gate.landmarks.md)

- [`do.gate()`](https://drcytometer.github.io/AutoSpectral/reference/do.gate.md)

- [`gate.define.plot()`](https://drcytometer.github.io/AutoSpectral/reference/gate.define.plot.md)
