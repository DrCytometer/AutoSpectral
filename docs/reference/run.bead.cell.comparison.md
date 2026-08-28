# Compare single-color control spectra across particle types

Runs the full AutoSpectral spectral signature extraction pipeline
(legacy or automated) independently on each particle-type folder in
`particle.dirs` – typically one folder of cell-based single-stained
controls and one or more folders of bead-based single-stained controls –
then compares each non-reference particle type against the reference
type using cosine similarity, absolute spectral distance, and MAD-based
variability metrics. Produces per-particle-type diagnostic figures and a
consolidated results object.

## Usage

``` r
run.bead.cell.comparison(
  particle.dirs,
  asp,
  reference.type = "Cells",
  result.dir = "./results/bead_cell_comparison",
  control.filename = "fcs_control_file",
  n.cells.brightness = 500L,
  fluor.df = NULL,
  figures = TRUE,
  verbose = TRUE,
  legacy.pipeline = FALSE,
  legacy.gating.system = "density",
  legacy.gate = TRUE,
  n.candidates = 1000L,
  n.spectral = 200L,
  k.neighbors = 2L,
  allow.duplicate.controls = TRUE,
  mismatch.abs.alpha = 0.5,
  cluster.peak.threshold = 0.1,
  cluster.min.n = 3L,
  cluster.stat.threshold = 1.96,
  cluster.n.perm = 2000L,
  min.variant.rows = 3L,
  exclude.fluorophores = NULL
)
```

## Arguments

- particle.dirs:

  Named list of character scalars. Each name is a particle-type label
  (e.g. `"Cells"`, `"UltraComp Beads"`) and each value is the path to
  that particle type's folder of single-stained FCS files plus an
  unstained/AF reference file. Must contain an entry named
  `reference.type`.

- asp:

  An `AutoSpectralParam` object supplying cytometer, gating, and
  pipeline parameters. Per-folder figure and table directories are
  overridden internally for each particle type.

- reference.type:

  Character scalar. The name in `particle.dirs` to treat as the
  reference (typically cell-based) population that all other particle
  types are compared against. Defaults to `"Cells"`.

- result.dir:

  Character scalar. Root output directory for cached extractions,
  comparison figures, and the final results `.rds` file. Defaults to
  `"./results/bead_cell_comparison"`.

- control.filename:

  Character scalar. Base filename (without extension) used for the
  auto-generated control file in each particle-type folder. Defaults to
  `"fcs_control_file"`.

- n.cells.brightness:

  Integer scalar. Number of top (background- subtracted) events used per
  fluorophore when computing automated brightness via
  [`get.brightness.automated()`](https://drcytometer.github.io/AutoSpectral/reference/get.brightness.automated.md).
  Defaults to `500L`.

- fluor.df:

  Optional data frame with columns `Fluorophore` and `Class`, used to
  annotate mismatch/variability plots by dye class. If `NULL` (the
  default), a fallback is built from the package's
  `fluorophore_database.csv`, using `excitation.laser` as the class.

- figures:

  Logical scalar. Whether to generate diagnostic figures during spectral
  extraction. Defaults to `TRUE`.

- verbose:

  Logical scalar. Whether to print progress messages. Defaults to
  `TRUE`.

- legacy.pipeline:

  Logical scalar. If `TRUE`, spectral signatures are extracted with the
  legacy gated pipeline
  ([`define.flow.control()`](https://drcytometer.github.io/AutoSpectral/reference/define.flow.control.md)
  /
  [`clean.controls()`](https://drcytometer.github.io/AutoSpectral/reference/clean.controls.md)
  /
  [`get.fluorophore.spectra()`](https://drcytometer.github.io/AutoSpectral/reference/get.fluorophore.spectra.md));
  if `FALSE`, the automated (non-gated) pipeline
  ([`get.spectra.automated()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectra.automated.md))
  is used. Defaults to `FALSE`.

- legacy.gating.system:

  Character scalar. Gating system passed to
  [`define.flow.control()`](https://drcytometer.github.io/AutoSpectral/reference/define.flow.control.md)
  when `legacy.pipeline = TRUE`. Defaults to `"density"`.

- legacy.gate:

  Logical scalar. Whether to perform gating in
  [`define.flow.control()`](https://drcytometer.github.io/AutoSpectral/reference/define.flow.control.md)
  when `legacy.pipeline = TRUE`. Defaults to `TRUE`.

- n.candidates:

  Integer, default `1000`. Number of top-expressing candidate events
  selected per fluorophore before cosine-similarity filtering. Ignored
  in internal-negative mode, where the top 5%% of events by peak channel
  are used directly.

- n.spectral:

  Integer, default `50`. Number of spectral events retained after
  filtering for low AF cosine similarity.

- k.neighbors:

  Integer, default `2`. Number of nearest neighbours in scatter space
  used for per-event AF subtraction.

- allow.duplicate.controls:

  Logical, default `FALSE`. Set `TRUE` to permit multiple single-stained
  controls for the same fluorophore (diagnostic/QC use only). Each is
  tracked internally under a unique `sample` identifier and carried
  through to `marker.spectra`'s rownames; the true fluorophore identity
  is preserved as a `"fluorophore"` attribute for reference-library
  matching and for
  [`check.spectra.duplicates()`](https://drcytometer.github.io/AutoSpectral/reference/check.spectra.duplicates.md).

- mismatch.abs.alpha:

  Numeric scalar in `[0, 1]`. Line transparency used for the
  absolute-value spectral mismatch plot described below. Defaults to
  `0.5`.

- cluster.peak.threshold:

  Numeric scalar in `[0, 1]`. Passed to
  [`assess.mismatch.clusters()`](https://drcytometer.github.io/AutoSpectral/reference/assess.mismatch.clusters.md)
  as `peak.threshold`: a fluorophore/detector pair is only used in the
  cluster test where that fluorophore has on-peak signal. Defaults to
  `0.1`.

- cluster.min.n:

  Integer scalar. Passed to
  [`assess.mismatch.clusters()`](https://drcytometer.github.io/AutoSpectral/reference/assess.mismatch.clusters.md)
  as `min.n`: the minimum number of on-peak fluorophores required for a
  detector to receive a test statistic. Defaults to `3L`.

- cluster.stat.threshold:

  Numeric scalar. Passed to
  [`assess.mismatch.clusters()`](https://drcytometer.github.io/AutoSpectral/reference/assess.mismatch.clusters.md)
  as `cluster.threshold`: the per-detector \|t\| cutoff used to form
  candidate clusters. Defaults to `1.96`.

- cluster.n.perm:

  Integer scalar. Passed to
  [`assess.mismatch.clusters()`](https://drcytometer.github.io/AutoSpectral/reference/assess.mismatch.clusters.md)
  as `n.perm`: the number of sign-flip permutations used to build the
  null distribution of the maximum cluster mass. Defaults to `2000L`.

- min.variant.rows:

  Integer scalar. A fluorophore is automatically excluded from a given
  particle type's comparison against `reference.type` if its spectral
  variant matrix has fewer than this many rows (or is absent entirely)
  for `reference.type` or for that particle type – evaluated
  independently per comparison, so a fluorophore failing this check for
  one bead type does not remove it from comparisons against other bead
  types. Defaults to `3L` (i.e. more than two rows required).

- exclude.fluorophores:

  Character vector, or `NULL` (the default). Fluorophore names to
  exclude from every particle type's downstream comparison/plots, in
  addition to the automatic exclusion from `min.variant.rows`. Names not
  present in the data are ignored.

## Value

A list with components:

- Extraction:

  Named list (by particle type) of extraction results, each containing
  `spectra`, `brightness`, and `variants`.

- Comparison:

  Named list (by non-reference particle type) of `Cosine`,
  `Variability`, `VariabilityMAD`, and `Distance` results versus
  `reference.type`.

- Stats:

  Named list (by non-reference particle type) of
  [`mismatch.plot()`](https://drcytometer.github.io/AutoSpectral/reference/mismatch.plot.md)
  outputs.

- MAD:

  Named list (by non-reference particle type) of per-detector MAD
  spectra of the mismatch distance matrix.

- Clusters:

  Named list (by non-reference particle type) of the `clusters` data
  frame from
  [`assess.mismatch.clusters()`](https://drcytometer.github.io/AutoSpectral/reference/assess.mismatch.clusters.md).

The same list is also written to
`<result.dir>/bead_cell_comparison_results.rds`.

## Details

Processing proceeds in four stages:

1.  **Per-particle-type extraction.** For each entry in `particle.dirs`,
    a control file is created if absent, validated for unresolved
    fluorophore matches and a resolved AF row, and then run through the
    selected extraction pipeline to yield spectra, automated brightness,
    and spectral variants. Results are cached to
    `<result.dir>/<particle.type>/<particle.type>_extraction_<pipeline>.rds`
    and reused on subsequent calls; switching `legacy.pipeline` between
    runs uses a distinct cache file rather than reusing a stale
    extraction from the other pipeline. Particle types that fail
    validation or processing are skipped with a warning and excluded
    from downstream comparison; if `reference.type` itself fails, the
    function stops.

2.  **Fluorophore exclusion.** For each non-reference particle type, any
    fluorophore whose spectral variant matrix has `< min.variant.rows`
    rows (or is missing) for **that** particle type or for
    `reference.type` is dropped from **that particle type's** comparison
    only, along with any fluorophore named in `exclude.fluorophores`.
    Exclusion is evaluated independently per comparison, so a bead type
    that fails QC outright (e.g. from low acquisition counts) does not
    remove fluorophores from other, unrelated comparisons. A comparison
    that fails outright itself (e.g. no fluorophore survives exclusion
    for that pair) is skipped with a warning rather than aborting every
    other comparison.

3.  **Mismatch assessment.** Each non-reference particle type's spectral
    variants are first renormalized from the package's native L-infinity
    convention (peak detector fixed to 1.0) to unit L2 (Euclidean) norm
    via
    [`l2.normalize.spectra()`](https://drcytometer.github.io/AutoSpectral/reference/l2.normalize.spectra.md),
    so that distance-based comparisons are not artificially biased
    toward zero variance at each spectrum's peak detector. The
    renormalized variants are then compared against the reference type's
    renormalized variants via
    [`assess.mismatch()`](https://drcytometer.github.io/AutoSpectral/reference/assess.mismatch.md)
    (cosine similarity, unaffected by the choice of per-row
    normalization),
    [`assess.variability()`](https://drcytometer.github.io/AutoSpectral/reference/assess.variability.md),
    [`assess.variability.mad()`](https://drcytometer.github.io/AutoSpectral/reference/assess.variability.mad.md),
    and
    [`bead.cell.dist()`](https://drcytometer.github.io/AutoSpectral/reference/bead.cell.dist.md)
    (per-detector signed distance). `Extraction` in the return value
    retains the original L-infinity-normalized spectra/variants as
    produced by the extraction pipeline; only the comparison stages
    below use the L2-renormalized versions.

4.  **Error-metric plots.** `plot.mismatch()` is called per
    non-reference particle type, with shared axis limits computed across
    all particle types for comparability.

5.  **Spectral location of mismatch.** Per-detector MAD of the distance
    matrix is computed and visualized with
    [`spectral.variant.plot.dens()`](https://drcytometer.github.io/AutoSpectral/reference/spectral.variant.plot.dens.md)
    two ways: a signed view (raw `reference - other` differences,
    showing whether mismatch is systematically one-signed at a given
    detector) and a magnitude-only view (absolute differences, showing
    where mismatch concentrates in the spectrum irrespective of sign).
    Both use the same per-detector MAD of the signed distances as the
    reference line.

6.  **Cluster-based mismatch testing.** For each non-reference particle
    type,
    [`assess.mismatch.clusters()`](https://drcytometer.github.io/AutoSpectral/reference/assess.mismatch.clusters.md)
    tests whether mismatch is concentrated in a specific, contiguous
    region of the spectrum rather than spread evenly across detectors,
    using a sign-flip permutation test restricted to on-peak
    fluorophore/detector pairs. Detected clusters (with permutation
    p-values) are written to `<result.dir>/table_mismatch_clusters/`.
