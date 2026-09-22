# Estimate Residual Spillover From a Known-Negative Mask

Thin, in-memory wrapper over the batched pair estimator, for callers
that already know which events are target-negative and do not need the
estimator to infer it.

## Usage

``` r
estimate.residual.spillover(
  unmixed,
  source,
  targets,
  negative.mask,
  threshold.source,
  threshold.target = NULL,
  spread.var = NULL,
  neg.var = NULL,
  ...
)
```

## Arguments

- unmixed:

  Numeric matrix, cells x fluorophores, already-unmixed abundances. Must
  contain a column named `source` and one column per entry of `targets`.

- source:

  Character scalar, the name of the source fluorophore's column in
  `unmixed`.

- targets:

  Character vector, the names of the target fluorophores' columns in
  `unmixed` to estimate residual spillover into.

- negative.mask:

  Logical vector, length `nrow(unmixed)`, `TRUE` for events already
  known to be target-negative. Used directly in place of the batched
  estimator's own negative-event inference.

- threshold.source:

  Numeric scalar or vector of length `nrow(unmixed)`, the source
  fluorophore's own per-event positivity boundary.

- threshold.target:

  Optional numeric matrix, cells x fluorophores, containing at least the
  columns in `targets`, giving each target's own per-event positivity
  boundary. When `NULL`, every event is treated as below threshold
  before `negative.mask` is applied.

- spread.var:

  Optional numeric vector, length `length(targets)`, the source's
  contribution to each target's spillover-spread variance. Defaults to
  zero for every target when `NULL`.

- neg.var:

  Optional numeric vector, length `length(targets)`, each target's
  negative-population variance. When `NULL`, computed as
  `stats::mad()^2` on each target's column of `unmixed`.

- ...:

  Additional arguments passed through to `.fix.envelope.slope.batch()`.
