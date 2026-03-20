# Define Gate by Density

A short description...

## Usage

``` r
tune.gate(
  control.file,
  control.dir,
  asp,
  n.cells = c(100, 500, 2000),
  percentiles = c(30, 50, 70),
  grid.n = 100,
  bandwidth.factor = 1,
  fsc.channel = NULL,
  ssc.channel = NULL,
  fsc.lims = NULL,
  ssc.lims = NULL,
  output.dir = "./figure_gate_tuning",
  gate.name = "cell_gate",
  filename = "gate_tuning",
  color.palette = "mako",
  boundary.color = "red",
  points.to.plot = 1e+05,
  width = 9,
  height = 10
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

- n.cells:

  The number of cells to use for defining the gate boundary. The default
  is `c(100, 500, 2000)` and will test each of those values. This many
  cells will be selected from the peak channel (brightest first) in the
  single- color controls. For example, if you set `200` and marked files
  such as `CD3-PE.fcs` and `CD19-FITC.fcs` as `gate.define=TRUE` in the
  control file, the brightest 200 events in the YG1 channel from the
  CD3-PE file and the brightest 200 events in the B1 channel for the
  CD19-FITC file would be used to define the gate.

- percentiles:

  Numeric 1 - 100, default `c(30, 50, 70)`. The percentile cutoffs to
  test for density in the scatter to use for defining the gate. For
  example, a value of `50` would take the 50% of cells closest to the
  density peak, more or less. Smaller numbers will define a tighter
  gate.

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

  Numeric vector. Limits for plotting the FSC. The default `NULL` uses
  c`c(asp$scatter.data.min.x, asp$scatter.data.max.x)`.

- ssc.lims:

  Numeric vector. Limits for plotting the SSC. The default `NULL` uses
  c`c(asp$scatter.data.min.y, asp$scatter.data.max.y)`.

- output.dir:

  File path where you want to save the results. Default is
  `./figure_gate_tuning`.

- gate.name:

  Character, name for the gate. Useful for distinguishing gates when you
  have multiple types. Default is `cell_gate`.

- filename:

  Character, name for the output files. Default is `gate_tuning`.

- color.palette:

  Optional character string defining the viridis color palette to be
  used for the fluorophore traces. Default is `plasma`. Use `rainbow` to
  be similar to FlowJo or SpectroFlo. Other options are the viridis
  color options: `magma`, `inferno`, `plasma`, `viridis`, `cividis`,
  `rocket`, `mako` and `turbo`.

- boundary.color:

  Color for the gate boundary line on the plot. Default is
  `darkgoldenrod1`.

- points.to.plot:

  Numeric, default `1e5`. Maximum number of points to show on the plot.
  More points will take longer, but you really shouldn't have even close
  to this number when defining the gate with landmarks.

- width:

  Numeric, default `4`. Width of the saved plot.

- height:

  Numeric, default `4`. Height of the saved plot.

## Value

A combined plot of all gates. Saves plots of the gates defined using
combinations of the specified parameters.

## References

Laniewski, Nathan. *flowstate*.
<https://github.com/nlaniewski/flowstate>

## See also

- [`define.gate.landmarks()`](https://drcytometer.github.io/AutoSpectral/reference/define.gate.landmarks.md)

- [`define.gate.density()`](https://drcytometer.github.io/AutoSpectral/reference/define.gate.density.md)

- [`gate.sample.plot()`](https://drcytometer.github.io/AutoSpectral/reference/gate.sample.plot.md)
