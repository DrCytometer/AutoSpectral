# Get Autofluorescence Spectra

Extracts autofluorescence spectra from an unstained samples. Intended
for use with `unmix.autospectral`. Uses FlowSOM (EmbedSOM) clustering
for rapid identification of cells with similar AF profiles.

## Usage

``` r
get.af.spectra(
  unstained.sample,
  asp,
  spectra,
  som.dim = 10,
  figures = TRUE,
  plot.dir = NULL,
  table.dir = NULL,
  title = "Autofluorescence spectra",
  verbose = TRUE,
  refine = FALSE,
  problem.quantile = 0.99,
  remove.contaminants = TRUE,
  parallel = TRUE,
  threads = if (parallel) 0 else 1,
  heatmap.color.palette = "viridis",
  spectral.trace.color.palette = NULL,
  af.fill.color = "red",
  af.line.color = "black"
)
```

## Arguments

- unstained.sample:

  Path and file name for a unstained sample FCS file. The sample type
  and processing (protocol) method should match the fully stained
  samples to which the AF will be applied, ideally.

- asp:

  The AutoSpectral parameter list.

- spectra:

  Spectral signatures of fluorophores, normalized between 0 and 1, with
  fluorophores in rows and detectors in columns.

- som.dim:

  Number of x and y dimensions for the SOM. Default is `10`.

- figures:

  Logical, whether to plot the spectral traces and heatmap for the AF
  signatures. Default is `TRUE`.

- plot.dir:

  Directory (folder) where the plots will be saved. Default is `NULL`,
  which inherits from `asp$figure.af.dir`.

- table.dir:

  Directory (folder) where the spectra csv file will be saved. Default
  is `NULL`, which inherits from `asp$table.af.dir`.

- title:

  Title for the output spectral plots and csv file. Default is
  `Autofluorescence spectra`.

- verbose:

  Logical, controls messaging. Default is `TRUE`.

- refine:

  Logical, default is `FALSE`. Controls whether to perform a second
  round of autofluorescence measurement on "problem cells", which are
  those with the highest spillover, as defined by `problem.quantile`.
  When `FALSE`, behavior is identical to versions of AutoSpectral prior
  to 1.0.0. If you are working with samples containing complex
  autofluorescence, e.g., tissues or tumors, using `refine=TRUE` will
  improve autofluorescence extraction in the unmixing at the cost of an
  increase in unmixing time. The increase in time will depend on the
  method used to assign autofluorescence spectra per cell (residual
  based assignment is very fast) and whether you have installed
  `AutoSpectralRcpp`, which will speed up assignment and unmixing.

- problem.quantile:

  Numeric, default `0.99`. The quantile for determining which cells will
  be considered "problematic" after unmixing with per-cell AF
  extraction. Cells in the `problem.quantile` or above with respect to
  total signal in the fluorophore (non-AF) channels after per-cell AF
  extraction will be used to determine additional autofluorescence
  spectra, using a second round of clustering and modulation of the
  previously selected autofluorescence spectra. A value of `0.99` means
  the top 1% of cells, those farthest from zero, will be selected for
  further investigation.

- remove.contaminants:

  Logical, default is `TRUE`. A QC check is performed to exclude any
  autofluorescence spectra that are nearly identical to the fluorophore
  signatures in `spectra`. This helps deal with low-level contamination
  of unstained samples by single-stained control samples, which happens
  sometimes. To include these AF spectra, which can mess up unmixing if
  they are really fluorophore spectra, set `FALSE`.

- parallel:

  Logical, default is `TRUE`, which enables parallel processing for
  per-cell AF identification. Used when `refine=TRUE`.

- threads:

  Numeric, defaults to a single thread for sequential processing
  (`parallel=FALSE`) or all available cores if `parallel=TRUE`.Used when
  `refine=TRUE`.

- heatmap.color.palette:

  Optional character string defining the viridis color palette to be
  used for the fluorophore traces. Default is `viridis`. Options are the
  viridis color options: `magma`, `inferno`, `plasma`, `viridis`,
  `cividis`, `rocket`, `mako` and `turbo`.

- spectral.trace.color.palette:

  Optional character string defining the color palette to be used for
  the AF traces. Default is `NULL`, in which case default R Brewer
  colors will be assigned automatically. Options are the viridis color
  options: `magma`, `inferno`, `plasma`, `viridis`, `cividis`, `rocket`,
  `mako` and `turbo`.

- af.fill.color:

  Color for the shaded region indicating the range of variation in the
  autofluorescence. Feeds to `fill` in `geom_ribbon`. Default is "red".

- af.line.color:

  Color for the line representing the median autofluorescence spectrum.
  Default is "black".

## Value

A matrix of autofluorescence spectra.

## References

Van Gassen S et al. (2015). "FlowSOM: Using self-organizing maps for
visualization and interpretation of cytometry data." *Cytometry Part A*,
87(7), 636-645.
[doi:10.1002/cyto.a.22625](https://doi.org/10.1002/cyto.a.22625) Wehrens
R, Kruisselbrink J (2018). “Flexible Self-Organizing Maps in kohonen
3.0.” *Journal of Statistical Software*, *87*(7), 1-18.
[doi:10.18637/jss.v087.i07](https://doi.org/10.18637/jss.v087.i07)
