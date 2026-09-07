# CONTEXT: Spectral signature error correction

Status as of the close of the gating-and-port session. This supersedes the
previous version of this document; the theory sections carried over are
condensed rather than restated in full.

---

## 1. Objective

Given stained sample data and a reference spectra matrix that is subtly
wrong for that sample (bead spectra applied to cells, spectra from another
day, lot, or instrument state), recover the correct per-fluorophore spectra
without making any row worse.

Test substrate: paired bead and cell single-stained control sets for the
same 16-colour Aurora panel, spectra extracted independently for each
particle type, tests run on concatenated data with ground truth available
in both directions. Production target: fully stained cell samples with a
same-day unstained control.

Outcome of this session: the method is implemented as a rewritten
`correct.unmixing.signatures()` and validated on the benchmark. See §7 for
production status.

---

## 2. Theory

### 2.1 Identifiability (carried over, condensed)

From mixed data alone, `spectra` is identifiable only up to left
multiplication: the row space is determined, not the rows. The OLS residual
is exactly orthogonal to every row, so a full-panel residual correction can
recover only the out-of-span error component; on the real bead/cell pair
the parallel (invisible) fraction is 0.84. The full-panel correction path
is closed (§2.2 of the previous document): one step reaches its fixed
point, which is not the truth.

Restricting the design to the fluorophores actually present makes the
problem solvable: in the rank-1 limit the only unidentifiable direction is
the row's own scale, removed by renormalisation. The update sign is settled
algebra: regressing the residual on abundance estimates the error directly,
and the update is `spectra + beta`.

### 2.2 The restricted path has its own ceiling (new)

The restricted residual is orthogonal to every *active* row, so the
component of a dye's true error lying in the span of its nuisance rows is
invisible to the restricted path too. Measured with ground truth
(diagnostics Stage I): near zero on cells for every dye, but substantial on
beads for the tandems - BV711 64% invisible under the flat co-activity
mask, 50% under the spread-scaled mask (which also shrank the nuisance set
from 4 dyes to 1). The spread-scaled mask therefore buys real
recoverability and is the default. The component along the dye's own row is
pure scale and is absorbed by renormalisation; it can be ignored.

### 2.3 Descent geometry of the residual line search (new)

With the restricted design and a perturbation of row `j` along `beta`, the
envelope theorem gives the derivative of the re-fit residual objective at
step zero as

```
J'(0) = -2 [ (sum x_j) (alpha . beta) + (sum x_j^2) ||beta||^2 ]
```

Measured (Stage G): every dye descends, including the ones the line search
vetoed. The `t.hat = 0` failures were therefore curvature, not sign: the
objective is a narrow valley and the old step grid's smallest step (0.125)
overshot it. Extending the grid down to 0.03125 resolved BV711 on cells.

### 2.4 The no-intercept slope is the right estimator (new)

If the residual is `r_c = x_c E(x_c)` with an abundance-dependent error
direction (a variant mixture), the no-intercept slope estimates

```
sum( x^2 E(x) ) / sum( x^2 )
```

the abundance-weighted mean shift - exactly the optimal target for a single
reference row unmixing a heterogeneous population. A free intercept fit
through curved data steals part of the signal. Empirically the no-intercept
direction was closer to the true error on ~14 of 16 rows on both
substrates (BV711 cells: cosine 0.72 versus 0.28). The production fit is
anchored instead by background subtraction plus the zero-abundance anchor
bin.

### 2.5 Curvature and variant mixtures, confirmed (new)

Stage H: the per-bin residual direction is not stable across abundance for
the tandems (BV711 cells minimum bin-to-bin cosine -0.72, i.e. the
direction reverses), and the bin-to-bin drift lies almost entirely in the
dye's own spectral-variant span (drift.in.delta 0.87 for BV711, 0.97 for
PerCP-eFluor 710, cells). The variant-mixture hypothesis for tandem
under-recovery is confirmed.

### 2.6 Variant-subspace projection: defensive only (new)

Constraining the correction to the span of the dye's variant deltas was
tested (Stage J) and rejected as the estimator: the true bead-versus-cell
error is only 45-75% inside the variant span for the dyes that recover
well, so the projection taxes genuine corrections heavily (PE recovery
0.97 unconstrained to 0.14 constrained) while the within-control
tandem-degradation axis explains curvature *within* a control, not the
cross-substrate shift itself. The projection's demonstrated value is
damage limitation (cells catastrophes of -8 to -17 collapsed to -0.4 to
+0.1), superseded in that role by the gate stack of §3.

Also rejected in the same stage: significance gating with naive
homoscedastic standard errors (z of 8-900 on everything, wrong directions
included - residuals are dominated by structured variant scatter, not
white noise), and a trust region from within-control delta norms (the cap
fired on the best corrections; within-control variant magnitudes
underestimate cross-substrate shifts).

### 2.7 Jacobian (carried over)

The corrected two-term formula in `compute.signature.jacobian()` is exact
(finite-difference error 2.8e-6). It is now diagnostic-only; nothing in the
production path uses it.

---

## 3. Gates

### 3.1 Failed gates (extended)

| Gate | Failure |
|---|---|
| r-squared on the bin regression | Inverted: collapses toward zero when a dye converges. |
| Split-half slope agreement | A systematically mis-specified model reproduces itself. |
| Held-out residual reduction (as sole gate) | Blind to in-span error; rewards repurposed rows. |
| AF PC projection | No-op after scatter-matched subtraction. |
| Step/ALS concordance | ALS is the iterated one-step; they agree at cosine 0.999 while both diverge. |
| Naive z-significance | Structured residuals make everything "significant". |
| Delta-norm trust region | Caps the best corrections; wrong scale for cross-substrate shifts. |

### 3.2 The working stack

All applied per fluorophore per iteration; any failure leaves the row
untouched.

1. `x.span > min.span * |threshold|` - enough abundance range to identify
   a slope.
2. `explained > 0.5` - at its brightest bin, the dye's own term accounts
   for most of the background-subtracted signal. Correctly excludes the
   too-dim controls (PerCP, NovaFluor Blue 610-30S).
3. Held-out step search: `t.hat > 0` and `gain > min.gain`, with the step
   grid extended down to 0.03125 for narrow-valley (tandem) dyes.
4. `rel.step <= 0.08` - a correction larger than ~8% of the row it
   corrects is rejected outright, never scaled down.
5. `span.drift <= 1.10` - a row rotating toward whatever else the sample
   contains inflates its own apparent abundance.
6. `bg.align > -0.9` (new) - the background-confound gate. A genuine
   spectral error gives the event-level regression an intercept and a
   slope with independent physical origins. A common-mode background
   residual whose magnitude tracks brightness is one physical direction
   split in two by the regression, so intercept and slope come out
   anti-collinear. Separation on this data: divergers at -0.93 to -0.999
   on both substrates, genuine corrections no worse than -0.66 (cells) /
   -0.56 (beads). Known miss: Spark Blue 550 on beads (-0.73). Known
   near-boundary pass that should pass: PerCP-eFluor 710 cells (-0.849).
   The gate is a cosine - unitless, panel- and instrument-independent -
   and doubles as the backstop that detects residual background confounds
   whatever subtraction mode produced them.

---

## 4. Main-population gating

Density-relative scatter gate (keep events above a fraction, default 0.1,
of the modal 2D scatter density; no per-run landmarks). Applied to the
sample and to the unstained control. It removes the noise/debris fraction
that was contaminating three things at once: the dominance populations
(noise events go dominant for whichever dye their AF direction leaks
into), the global-mean background estimate (the likely origin of the bead
background confounds), and the evaluation medians (the likely origin of
Spark Blue 550's misleadingly negative bead score - by skew plots the
correction is good).

Caveat: on beads the gate also removes multiplets, which are bright and
informative. Acceptable because the production target is fully stained
cells, not beads; for bead-based testing set the gate off. Toggles:
`diag.gate.main` in the harness, `gate.main` in the package function.

---

## 5. What was built this session

Harness (`run_signature_correction_diagnostics_1.R`, Stage D fit):
multivariate slope (all active columns, preventing a co-varying nuisance
dye's error being booked against the dominant dye), spread-scaled
co-activity mask, extended step grid, `bg.align` gate, dominance
assignment returned so Stages E/F evaluate the fitted populations, optional
main-population gating of sample and unstained.

New diagnostic scripts: `run_signature_correction_diagnostics_2.R`
(Stage G descent audit, Stage H ALS/concordance/curvature, Stage I
restricted identifiability ceiling) and
`run_signature_correction_diagnostics_3.R` (Stage J variant-subspace
constrained estimator, own- and cross-basis).

Skew script: fitted-population evaluation, degenerate-median guard,
attribution controls auto-selected by disagreement.

Package: `correct_unmixing_signatures.R` rewritten in full. The exported
`correct.unmixing.signatures()` now implements the dominance path with the
complete gate stack, optional density gating, spread-scaled co-activity,
`bg.mode` in {global.mean, scatter.knn, none}, and optional ground-truth
recovery diagnostics. The cluster-residual full-panel loop and its
arguments (`cluster.method`, `rlm.method`, `learning.rate`, `direction`,
`collinear.damping`, `smooth.window`, `recluster.every`) are removed - a
deliberate API break, since that path is structurally closed.

No longer called from production and pending deprecation notes:
`fit.signature.error.model()`, `update.signature.matrix()`,
`calculate.cluster.residuals()`, `restrict.cluster.unmixing()`,
`compute.signature.jacobian()`, `diagnose.single.positive.clusters()`.
`cluster.unmixed.events()` remains in use by the diagnostics.

---

## 6. Results (gated benchmark run)

Cells: accepted set now includes BV711 (new, via the finer step grid).
Median leakage 1.186 -> 1.151 against a correct-spectra floor of 1.059.
Beads: median leakage 0.133 -> 0.101 against a floor of 0.094. By skew
plots, PerCP-eFluor 710 remains uncorrected on both substrates, BV711
uncorrected on beads (out of production scope), Spark Blue 550 on beads
mostly corrected with a small residual.

Metric caveats: gated-run medians are not comparable to earlier ungated
runs (the evaluation population changed). Several accepted cells dyes show
higher MAD-scaled abundance error despite lower leakage (FITC, Spark Blue
550, PE) - the suspected cause is benign: a correction that moves the
row's peak channel rescales all abundances via L-infinity renormalisation,
and the median-absolute-difference metric counts that constant scale shift
as error. Stage E needs a scale-matched error column (fit and remove a
per-dye scale against the truth-unmixed values) before those numbers are
used to judge anything.

---

## 7. Production status and workflow

Status: production candidate. Functionally complete and validated on the
bead/cell benchmark. Before calling it production-ready:

1. Parity check - run the package function on the benchmark with
   `true.spectra` supplied and confirm the `recovery` table matches the
   harness Stage D.
2. Fully stained sample test - dominance under genuine co-expression is
   the untested regime. On single-stain concatenates, co-activity is
   mostly spillover; on real samples it is real biology, and the nuisance
   machinery has not yet been exercised against that.
3. AF-heavy validation - this dataset (PBMCs and beads) barely exercises
   the background machinery.

Workflow position - runs after variant extraction, before unmixing:

1. `get.spectra.automated()` -> reference spectra from the controls.
2. `get.spectral.variants()` -> `thresholds`, `spillover.spread`,
   variant matrices.
3. `correct.unmixing.signatures( raw.data = fully stained raw,
   scatter = its scatter, spectra, unmixed.thresholds =
   variants$thresholds, spillover.spread = variants$spillover.spread,
   unstained + unstained.scatter = same-day unstained,
   bg.mode = "scatter.knn" for cells, gate.main = TRUE )`.
4. Feed the corrected `$spectra` to `unmix.fcs()` /
   `unmix.autospectral.joint()`. The variant matrices are absolute
   spectra and remain valid alongside the corrected reference; for full
   consistency of thresholds and deltas, variants can be re-extracted
   against the corrected spectra, but this is optional.

Input provenance rule: when `gate.main` is on, `unmixed.thresholds` must
come from a gated unstained control (the harness now does this); thresholds
from ungated unstained data will be slightly off for noisy files.

---

## 8. Open problems and next steps

### 8.1 PerCP-eFluor 710 - the remaining priority

Uncorrected on both substrates; on cells, its residual leakage is
dominated by the persisting BV711 -> PerCP-eFluor 710 channel even after
BV711's own row corrected. Two candidate mechanisms, with different fixes:
its variant mixture is too heterogeneous for a single-row correction
(fix: per-variant correction - fit the residual on the dye's delta basis
per abundance bin, or emit the corrected row as a candidate variant into
the joint per-cell pipeline and let per-cell selection arbitrate, with the
adoption fraction among bright dominant cells as the promotion gate - the
strongest available "never worse" deployment); or the remaining error is
in-span, which belongs to the `fix.my.unmix()` stage below.

### 8.2 `fix.my.unmix()` integration - the in-span complement

Estimates the F x F in-span matrix from population negativity on a fully
stained sample: exactly the component the residual path cannot see, with
complementary data requirements (negativity events versus dominance
populations). Runs *after* row-shape correction, because its unmixed-space
slopes conflate shape error with spillover error. Required fixes, from
code review of the current file:

- `ratioed.coefficients <- marker.spillover * similarity.matrix` biases
  the fixed point through both the update and the convergence test; the
  fixed point becomes `marker.spillover * similarity = I`. Trust belongs
  in a damping weight on the update step (an SE-based weight from
  `fit.robust.linear.model()`), not on the coefficient estimate.
- `compensation.matrix <- solve( ratioed.coefficients )` inside the loop
  is dead code.
- The `max.iter` argument is ignored; the loop tests `asp$rs.iter.max`.
- `downsample` lacks guards for `FALSE` and for counts exceeding the
  event count.
- `plot.unmix.fix()` still expects the removed `spread.estimate`
  argument.
- Known one-directional bias: gating to below-threshold events truncates
  the response from above, shrinking slopes toward zero (the safe
  direction); the spread-scaled boundary reduces it.

### 8.3 AF deconvolution as a `bg.mode`

Per-event unmixing against `rbind( af.spectra, spectra )` using the
`get.af.spectra()` library, subtracting the fitted AF component, with a
per-dye hotspot guard (a dye whose spectrum sits close to the AF span
falls back to scatter-knn rather than letting the projection eat real
signal). Preferred long-term because it does not depend on fluorophore
abundances and handles low-rank panels. Not implemented: it cannot be
validated on this dataset (AF is a non-issue here) and needs an AF-heavy
set. The `bg.align` gate stays as the backstop either way.

### 8.4 Carried-over preconditions

- AF rows in both benchmark tables are still untrustworthy (previous
  §7.2); re-extract before quantitative benchmarking of anything
  intercept-related.
- Cells PerCP and NovaFluor Blue 610-30S controls are too dim to
  characterise; the `explained` gate correctly excludes them. Reacquire.
- `laser.groups` auto-derivation for any future smoothing across detector
  array boundaries.

### 8.5 Housekeeping

- Stage E scale-matched abundance error column (§6).
- Deprecation notes on the retired phase functions; NEWS entry for the
  `correct.unmixing.signatures()` rewrite and API break.
