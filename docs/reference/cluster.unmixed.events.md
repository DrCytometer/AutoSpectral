# Cluster Unmixed Events into Micro-Clusters

Groups events into small, homogeneous micro-clusters based on their
currently-unmixed abundance profile, and computes a per-cluster centroid
abundance vector and matched raw-detector-space centroid vector. This is
Phase 1 of the signature-error-correction workflow: it stabilises
cluster-level autofluorescence/background variation and prevents a few
dominant populations from swamping the downstream regression, while
still sampling the abundance range of every fluorophore.

## Usage

``` r
cluster.unmixed.events(
  unmixed,
  raw.data,
  asp,
  method = c("som", "stratify", "leiden", "kmeans"),
  cluster.id = NULL,
  n.clusters = 150,
  som.dim = 20,
  k.neighbors = 15,
  resolution = 1,
  unmixed.thresholds = NULL,
  n.levels = 10L,
  stratify.exclude = NULL,
  centroid.fun = c("median", "mean"),
  min.cluster.size = 5,
  cluster.cofactor = 500,
  threads = NULL,
  verbose = TRUE
)
```

## Arguments

- unmixed:

  Numeric matrix (events x fluorophores), OLS/WLS unmixed abundance for
  the sample under investigation. Typically the output of
  [`unmix.ols.fast()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.ols.fast.md)
  on `raw.data` with the current `spectra`.

- raw.data:

  Numeric matrix (events x detectors), the raw detector-space data that
  produced `unmixed`. Column order and names must match the detectors of
  `spectra`.

- asp:

  The AutoSpectral parameter list from
  [`get.autospectral.param()`](https://drcytometer.github.io/AutoSpectral/reference/get.autospectral.param.md).
  Used for `asp$bird.seed`.

- method:

  Character, `"som"` (default), `"leiden"`, or `"kmeans"`. `"som"` uses
  the package's own batch SOM
  ([`get.som.codes()`](https://drcytometer.github.io/AutoSpectral/reference/get.som.codes.md),
  the `AutoSpectralRcpp`-accelerated engine already used by
  [`get.af.spectra()`](https://drcytometer.github.io/AutoSpectral/reference/get.af.spectra.md)
  and
  [`get.fluor.variants()`](https://drcytometer.github.io/AutoSpectral/reference/get.fluor.variants.md),
  not `FlowSOM`/`EmbedSOM` directly), with events assigned to their
  nearest node via
  [`FNN::knnx.index()`](https://rdrr.io/pkg/FNN/man/knn.index.html).
  `"leiden"` builds a k-nearest-neighbour graph
  ([`FNN::get.knn()`](https://rdrr.io/pkg/FNN/man/get.knn.html)) and
  clusters it with
  [`igraph::cluster_leiden()`](https://r.igraph.org/reference/cluster_leiden.html) -
  `igraph` is not currently a package dependency, so this method needs
  it installed separately for now. `"kmeans"` gives an exact cluster
  count and is simplest to reason about, but slower at high event counts
  and high k.

- cluster.id:

  Optional integer vector, length `nrow(unmixed)`, a previously-computed
  cluster assignment (e.g. `$cluster.id` from an earlier call). When
  supplied, clustering is skipped entirely and this assignment is used
  directly to compute centroids - intended for reuse across correction
  iterations in
  [`correct.unmixing.signatures()`](https://drcytometer.github.io/AutoSpectral/reference/correct.unmixing.signatures.md),
  so cluster identity stays fixed while only `unmixed` (and therefore
  `x.clust`) changes as `spectra` is corrected. Default `NULL` (cluster
  fresh).

- n.clusters:

  Target cluster count for `method = "kmeans"`. Default `150`. Not used
  by `"som"` (see `som.dim`) or `"leiden"` (see `resolution`).

- som.dim:

  Integer, SOM grid side length for `method = "som"` (`som.dim^2`
  nodes). Default `20`.

- k.neighbors:

  Integer, neighbours per event for the k-NN graph in
  `method = "leiden"`. Default `15`.

- resolution:

  Numeric, resolution parameter passed to
  [`igraph::cluster_leiden()`](https://r.igraph.org/reference/cluster_leiden.html)
  for `method = "leiden"`. Higher values give more, smaller clusters.
  Default `1`. Unlike `n.clusters`/`som.dim`, this does not target an
  exact cluster count directly - see Exploration points.

- unmixed.thresholds:

  Named numeric vector of per-fluorophore positivity thresholds in
  unmixed space. Required for `method = "stratify"`, ignored otherwise.

- n.levels:

  Integer, number of abundance bins per fluorophore for
  `method = "stratify"`. Default `10`.

- stratify.exclude:

  Optional character vector of fluorophore names that are not eligible
  to be an event's dominant fluorophore under `method = "stratify"`.
  Autofluorescence belongs here: it is present in every event and would
  otherwise dominate every assignment. Default `NULL`.

- centroid.fun:

  Character, `"median"` (default) or `"mean"`. Function used to
  summarise each cluster's abundance and raw-data vectors. Median is
  more robust to within-cluster outliers.

- min.cluster.size:

  Integer. Clusters with fewer than this many events are dropped before
  the regression step, since their centroids are unstable. Default `5`.

- cluster.cofactor:

  Numeric, asinh cofactor applied to `unmixed` before clustering only.
  Abundances span several decades, so linear-scale distances are
  dominated by the brightest population; clustering on
  `asinh(unmixed / cluster.cofactor)` distributes nodes across the whole
  dynamic range. Centroids are always computed in linear units. Default
  `500`.

- threads:

  Integer, OpenMP threads for `method = "som"`. Default `NULL` (all
  available cores, matching
  [`get.som.codes()`](https://drcytometer.github.io/AutoSpectral/reference/get.som.codes.md)'s
  own default). Ignored by other methods.

- verbose:

  Logical, controls messaging. Default `TRUE`.

## Value

A named list:

- `cluster.id`:

  Integer vector, length `nrow(unmixed)`, the cluster assignment for
  every event (`NA` for events in clusters dropped by
  `min.cluster.size`).

- `x.clust`:

  Numeric matrix (clusters x fluorophores), the per-cluster centroid
  abundance vector. Row names are cluster IDs.

- `y.clust`:

  Numeric matrix (clusters x detectors), the per-cluster centroid
  raw-data vector, matched row-for-row to `x.clust`.

- `cluster.size`:

  Integer vector, length `nrow(x.clust)`, the number of events
  contributing to each cluster centroid.

## Exploration points

- SOM vs Leiden vs kmeans vs stratify: SOM grid cells are not guaranteed
  to be evenly populated, which is especially visible on data with
  strong block structure (e.g. concatenated single-stained controls,
  where most of the grid may be near-empty for any one fluorophore's
  positive population); Leiden should adapt cluster count to the actual
  density structure instead of a fixed grid, at the cost of not
  targeting an exact cluster count directly; kmeans gives the most
  uniform cluster sizes but assumes roughly spherical clusters in
  abundance space; stratify guarantees every fluorophore the same number
  of populated abundance levels regardless of rarity, at the cost of
  needing `unmixed.thresholds` up front. Worth comparing all four on the
  same real data.

- `resolution` (Leiden) doesn't target an exact `n.clusters` the way
  `som.dim`/`n.clusters` do for the other methods. A wrapper that
  adjusts `resolution` until `length(unique(cluster.id))` lands near a
  target could be added if a specific count is needed - not implemented
  here.

- The Leiden graph is currently unweighted (plain k-NN edges); an
  inverse-distance-weighted graph is a common refinement worth trying.

- median vs mean centroids: median is the codebase default elsewhere,
  but for skewed clusters, mean might track a systematic error pattern
  more faithfully. Worth comparing both.
