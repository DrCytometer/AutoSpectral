# CONTEXT: fix.my.unmix / correct.unmixing.signatures — consolidated status

Supersedes, for narrative purposes, `CONTEXT_fix_my_unmix.md`,
`CONTEXT_fix_my_unmix_session_summary.md`, `CONTEXT_A10_A11_session_summary.md`,
`CONTEXT_fast_diagnostics_session_summary.md`,
`CONTEXT_round3_review_and_convergence_freeze.md`,
`CONTEXT_envelope_slope_openmp_session_summary.md` and
`CONTEXT_synthetic_null_session_summary.md`. Updated after a follow-on session
that closed out the batched-pair-estimator correctness gap this document's
previous version flagged as blocking (section 4): the bulk-population cap was
ported into `fix_envelope_truncated_batch_rcpp.cpp`, validated by a dedicated
A/B test built for the purpose, and the production diffs applied and
confirmed loaded. The subsequent real/synthetic null-fit re-run surfaced two
further findings, both now resolved by direct code inspection rather than
inference: neither null-fit script can actually exercise the bulk cap
(both size their input at exactly `max.truncated.events`), and the
tandem-enrichment Fisher's test in the prior session's report used the wrong
column (a stale hardcoded candidate list rather than
`fluorophore_database.csv`'s `class` field) — corrected here. Also new this
session: `correct.spectra.glasso()`, a complete, previously-unstarted
candidate third estimator (item 11, section 5), now implemented and awaiting
its first validation run.

**Updated after a follow-on session that ran `run_correct_spectra_glasso_validation.R`
for the first time.** A scatter-column-naming bug in the script itself was
found and fixed first (`.sim.scatter()` always names its output `c("FSC",
"SSC")`; the script now renames to `flow.control$scatter.parameter` after
calling it). With the script running, the validation surfaced a real,
still-unresolved phase-one problem shared by both estimators (neither
recovers an injected spillover coefficient at anything close to its true
size), a real and fixed phase-two bug in `extract_raw_signature.R` that was
silently rejecting good `fix.my.unmix()` candidates, and a second,
independent phase-one problem specific to `fix.my.unmix()`'s own pairwise
envelope estimator. None of this was resolved by inference or by reading
source alone — it took several rounds of diagnostics, most of which
disconfirmed the working hypothesis at the time, to get here. See section 9,
which is the primary record of that session and should be read before
trusting either estimator's phase-one spillover-detection step on
library-error-style corruption.

**Updated after a follow-on session that ran Stage J of
`run_signature_correction_diagnostics_aurora.R` for the first time,
extended the suite with four new stages (K–O), and ran all of them to
completion.** Stage J closes out the `max.step = 0.08` ceiling question
flagged three sessions ago (section 1): it is the binding constraint for
exactly one fluorophore in this benchmark, and a `0.15` default is now
recommended. Stages K–O chase a different, harder question the Stage J
run surfaced along the way — several fluorophores that don't correct at
any `max.step`, and two on Beads that correct confidently but in the
*wrong direction* — through five specific, sequentially-tested hypotheses
(held-out split noise, two real bugs in the diagnostic reimplementation
itself, background-anchor leverage, aggregate/tail contamination, and
autofluorescence data hygiene). Four of the five are now ruled out with
direct evidence; the fifth (AF) is confirmed as a real, independent data
problem worth fixing but is *also* ruled out as the explanation for the
specific failure it was suspected of causing. Root cause for the Beads
direction-flip is not yet identified. See section 10, which is the primary
record of this session.

**Updated after a follow-on session that used a Python reproduction to find
the root cause of section 9.4's "complete failure to recover either
injected pair," fixed it in `correct.spectra.glasso()`, confirmed the fix
at the pairwise-coefficient level, and then hit a second, separate, still-
unresolved problem before the corrected coefficients ever reached the
function's actual output.** Section 9.4's numbers were reading the wrong
cell of a matrix that was never wrong: the injected corruption's signal
sits in the transposed pairing from the one that section checked, for
reasons specific to how OLS unmixing changes basis, not a defect in either
estimator's fitting code. Restructuring `active.set` construction around
this fixed the pairwise recovery check (confirmed in R, matching
`fix.my.unmix()`'s own numbers closely). It did not fix `deg.glasso`:
phase two now rejects both corrected candidates for a different reason
(`resid.rel` over threshold, in a code path the larger active set makes
more likely to trigger), and the diagnostic built to explain why was never
read this session. See section 11, which is the primary record of this
session and supersedes section 9.4's "not identified" for root cause
specifically, while leaving the function's practical status unchanged
(still not recommended for production).

**Updated after a follow-on session that used the same Python-reproduction
workflow on `correct.unmixing.signatures()` itself, found a general
mechanism behind several of section 10's never-correcting fluorophores,
fixed it, and validated the fix in Python against real ground-truth
spectra (not yet re-confirmed in R).** A dye's own emission footprint —
how many of the panel's detectors it actually emits into — turns out to
predict correctability better than the size of its cross-substrate error
alone: the held-out step search's objective sums over every detector, so
a real shape error confined to a handful of UV channels is diluted by
noise from the ~55 channels a narrow dye never touches. Restricting the
slope fit and the held-out objective to a dye's own footprint
(`footprint.frac`, new parameter), combined with supplying a real
`spillover.spread` for Beads (previously absent; now supplied by Oliver
as `bead_spillover_spread.csv`), raised the simulated accepted count from
6/16 to 11-12/16 per direction with no fluorophore moving away from
truth. Both mechanisms plausibly generalise to `fix.my.unmix()` and
`correct.spectra.glasso()`, which share `extract_raw_signature.R`'s
structurally identical `resid.rel` gate and already guard against a
missing `spillover.spread` with their own warning — neither checked this
session. See section 12, which is the primary record of this session.

**Updated after a follow-on session that ran `run_footprint_validation.R`
for the first time (confirming section 12 on the real production function
and real gated data), retracted an intervening `spillover.spread`-decoupling
proposal, and tested true AF deconvolution (`get.af.basis()` /
`deconvolve.af.background()`) for the first time as a candidate `bg.mode`
for `correct.unmixing.signatures()`.** Section 12's `footprint.frac`
recommendation is now confirmed on real data (Cells 10/16 → 11/16, Beads
8/16 → 10/16); a proposal from an intervening, otherwise-unrecorded session
to decouple `spillover.spread` from nuisance-set gating is retracted, since
its original, unmodified behaviour is a clear net positive on the same real
data. A second, mechanistically distinct AF-removal approach — the SVD/
joint-OLS mechanism `fix.my.unmix()` already uses by default as
`bg.mode = "af.deconv"` — was tested for the first time as an option for
`correct.unmixing.signatures()` specifically: net negative when used as a
first background-removal step ahead of any spectral correction, essentially
inert (one small exception) when layered on top of already phase-1-corrected
spectra, matching the role it actually plays in `fix.my.unmix()`. Also
corrects a terminology error from an intervening session: what that session
called "AF deconvolution" (a SOM-based candidate-AF-library plus per-cell
discrete assignment) is "per-cell AF extraction," a different mechanism,
now tracked separately. See section 13, which is the primary record of this
session.

**Updated by a consolidation pass that folds in four sessions' worth of
`fix.my.unmix()`-side work not previously recorded in this document:**
`CONTEXT_hypernegative_penalty_and_cluster_audit.md`,
`CONTEXT_hypernegative_penalty_implementation.md`,
`CONTEXT_gating_downsampling_and_af_basis_fixes.md`, and
`CONTEXT_pair_restricted_residual_validation.md`. The first three are a
single continuous thread — a soft over-correction penalty that didn't
engage strongly enough, a redesign into a hard-reject fraction-based gate,
and, chasing why the redesign still didn't visibly fix the motivating
case, two much larger and unrelated bugs in scatter gating and
downsampling that were silently excluding the exact population the
correction needed most — now consolidated as section 14 below, with all
four documents' file-level changes cross-checked directly against current
source rather than taken on the documents' word. Every specific
`fix.my.unmix.R` and `correct_unmixing_signatures.R` change those four
sessions describe as applied is confirmed present in the copies in this
project as of this pass, **except** `run_hypernegative_cluster_audit.R`
and `run_pair_restricted_residual_validation.R` themselves, both of which
still contain the pre-fix logic their own originating documents describe
as superseded — see section 14.6 and section 15's opening note. The fourth
document is new, unresolved analytical work (a detector-space residual
check for phase two, restricted to pairs with both suspect rows removed
from the design) rather than a bug-fix thread, and is kept live as section
15 rather than folded into 14's narrative. Nothing in sections 1–13 below
is touched by this pass except where a stale "not yet applied" claim is
corrected against current source (section 1, section 3, section 5, section
7); the historical narrative in sections 9–13 is left exactly as written.

Benchmark in current use: FACSDiscover, 33 fluorophores + AF, paired Cells /
Beads single-stained concatenates (`a8.*` objects), cross-substrate ground
truth in both directions. The older 16-colour Aurora bead/cell pair
(`c.results`, `run_signature_correction_diagnostics_1/2.R`, continued this
session as `run_signature_correction_diagnostics_aurora.R` with Stages K–O
appended) remains the benchmark for `correct.unmixing.signatures()`'s own
row-shape diagnostics and has not been ported to the Discover panel. This
session's Python reproduction (section 13) used the same 16-colour Aurora
bead/cell benchmark, against newly-provided FSC/SSC main-population-gated
concatenated files (`small_Concatenated_cells.fcs`,
`small_Concatenated_BD_beads.fcs`). Separately, sections 14–15 use a
different benchmark entirely: a real, gated seven-colour lung/spleen panel
(BUV395, BUV805, BV421, PE, PE-Cy7, APC, eFluor 780 + AF) with genuine
biological co-expression, not concatenated single-stained controls — see
section 14's opening note and section 15.5 for why that distinction matters
to how far those sections' findings can be trusted.

---

## 1. Headline status

Three estimators now, three maturity levels:

- **`correct.unmixing.signatures()`** (row-shape, out-of-span-blind residual
  correction) — production candidate, functionally complete, validated on the
  Aurora benchmark. The `max.step = 0.08` ceiling flagged in
  `CONTEXT_signature_correction2.md` **was tested (Stage J, two sessions
  ago) and only partially confirmed**: it is the binding constraint for
  exactly one fluorophore in this benchmark
  (`PerCP-eFluor 710`, both substrates), not the broader group of
  never-correcting fluorophores that motivated re-running Stage J in the
  first place. `max.step = 0.15` is recommended as the new default —
  recovers `PerCP-eFluor 710` cleanly, plateaus beyond that, and one
  fluorophore (`NovaFluor Blue 610-30S`, already known to be too dim to
  trust) shows a real regression at `max.step ≥ 0.25`, arguing against going
  higher. The never-correcting group Stage J's sweep couldn't rescue is
  **now substantially resolved (section 12), and the fix is now confirmed
  on the real production function against real data (section 13.1)**: a
  Python reproduction found that a dye's own emission footprint (how many
  detectors it actually emits into), not the size of its cross-substrate
  error, is what predicts whether the held-out gain gate can detect a real
  correction — the objective sums over all 64 detectors, diluting
  narrow-emission dyes' signal with noise from channels they never touch.
  Restricting the fit and the held-out objective to a dye's own footprint
  (`footprint.frac`, new parameter, default `0.02`), combined with supplying
  a real `spillover.spread` for Beads (previously absent; now supplied),
  raised the accepted count from 10/16 to 11/16 (Cells) and 8/16 to 10/16
  (Beads) on the real, gated benchmark via `run_footprint_validation.R`,
  with one small, isolated regression (`PE`, Beads) against several larger
  gains. **This consolidation pass confirms all of it is now live in
  production `correct_unmixing_signatures.R`, not merely recommended**:
  direct reading shows `max.step = 0.15`, `footprint.frac = 0.02`,
  `footprint.min.channels = 3L`, and the `n.split.trials`/`min.split.frac`
  parameters (default `1L`/`0.6`, backward-compatible, not yet promoted to
  a non-`1` default) all present as of this pass — see section 3, which
  previously and incorrectly listed these as unapplied. **A separate,
  intervening proposal to also decouple
  `spillover.spread` from nuisance-set gating is retracted (section
  13.1)** — real gated data shows the original, unmodified behaviour is a
  net positive, not a bug. **`bg.mode = "af.deconv"` (the SVD/joint-OLS AF
  deconvolution mechanism `fix.my.unmix()` already uses by default) was
  tested for the first time as an option for this function's own
  background-removal step and is not recommended there (section 13.4-13.5)**
  — net negative as a first-phase step, essentially inert as a second-phase
  step; the existing two-phase split with `fix.my.unmix()` already places
  the mechanism correctly.
- **`fix.my.unmix()`** — production function's three patches
  (`gate.on.bias` default `FALSE`, `leakage.margin` + seed-matching, the
  batched pair estimator wired into phase one) are now all confirmed working
  as intended. The batched estimator's bulk-cap gap (previously this
  document's top-priority blocking item) is **resolved and validated** — see
  section 4. Two further things surfaced while confirming this, both now
  understood and neither blocking: the standard null-fit scripts cannot
  actually exercise the bulk cap at their current settings (section 4,
  "Why the null-fit re-run couldn't confirm this directly"), and a real,
  separate gap between this session's `n.fitted` counts and the numbers this
  document previously recorded is open at low priority (section 5, item 3).
  **New this session (section 9): a fourth patch, to `extract_raw_signature.R`'s
  `explained.total` calculation, is applied and confirmed working — but a
  second, separate, unresolved problem in phase one (the pairwise envelope
  slope estimator itself) was found underneath it. Production `fix.my.unmix()`
  currently produces good final corrected spectra on the one benchmark this
  session tested, but only because phase two's full-panel design routes
  around phase one's broken coefficient, not because phase one works.**
- **`fix.my.unmix()`'s over-correction ("hypernegative") problem — four
  sessions of work, now resolved and deployed (section 14).** A real,
  visually- and cluster-confirmed over-correction on a lung sample
  (BUV805→PE-Cy7) turned out to have two independent causes, neither of
  them the phase-one estimator itself: `fix.my.unmix()`'s internal scatter
  gate and its default stratified downsample were both silently excluding
  the large, high-SSC BUV805+/PE-Cy7+ cells the correction needed to see,
  from both the stained fit and the AF reference. With `scatter.gate =
  FALSE` and `downsample = FALSE` (new parameters), the over-correction is
  "perfectly solved" by Oliver's own direct comparison against an
  independent OLS unmix. Alongside this, a genuine estimator gap was also
  fixed: the old `max.hypernegative.delta` soft damper is replaced by a
  hard-reject fraction-based gate (`min.hypernegative.frac`/
  `max.hypernegative.frac`/`min.hypernegative.events`), per Oliver's
  design mandate that a negative source-target correlation is unambiguous
  evidence of an unmixing error and should not merely be discouraged. All
  of this is confirmed present in current `fix_my_unmix.R` by direct
  reading this pass. **Not fully closed**: phase two's own signature
  re-measurement still runs on the (possibly downsampled) `stained.fit`,
  not the new `stained.fit.full`, so the corrected spectra themselves —
  not just the reported `unmixed.final` — could still be under-informed by
  the same thinned tail in a fluorophore-rich panel; and the two auxiliary
  scripts this work produced, `run_hypernegative_cluster_audit.R` and (see
  next bullet) `run_pair_restricted_residual_validation.R`, both still
  contain logic their own originating documents describe as superseded —
  see section 14.6.
- **A new, unresolved analytical thread: pair-restricted residual
  validation of `fix.my.unmix()`'s phase two (section 15).** A standalone
  script regresses a detector-space residual — with *both* members of a
  suspect pair removed from the design, not just one, which is shown to be
  the minimum restriction under which an in-span error becomes visible at
  all — on a seven-colour lung sample, as an independent check on phase
  two's own signature-acceptance gates. Flagged one function
  (`eFluor 780`) that is best-in-panel on every existing acceptance gate
  yet shows a consistent, split-half-significant *worsening* against four
  of five informative partners under this test, plus a smaller, real
  mutual-coupling effect between the panel's two most collinear rows
  (`BUV395`/`BUV805`). **Explicitly flagged by Oliver as read too far into
  a substrate that cannot yet support the interpretation**: the test's
  population definitions come from full-panel unmixing on a real,
  biologically co-expressing sample, not from concatenated single-stained
  controls, so an observed effect cannot yet be distinguished from real
  co-expression. Not wired into production in any form; a rerun on
  concatenated single-stains for the same panel is the identified
  prerequisite before treating any of this thread's specific findings as
  confirmed.
- **`correct.spectra.glasso()`** — run for the first time two sessions ago
  (section 9), where it appeared to completely fail to recover either of
  two injected spillover coefficients and the root cause was not
  identified. **That root cause is now identified and fixed (section 11):
  section 9.4's numbers were reading the wrong cell of a matrix that was
  never wrong.** The pairwise recovery check itself is now confirmed
  working — `recovered.glasso` for both injected pairs lands within ~10% of
  the true coefficient, matching `fix.my.unmix()`'s own numbers on the same
  data closely enough to confirm this is a shared, general identification
  property and not specific to either estimator. **But this did not fix the
  function's actual output**: corrected `spectra` is bit-identical to the
  corrupted starting spectra for both injected rows in every run this
  session, including a run with 5x the events. Phase two rejects both
  candidates for `reason = "fit"` (`resid.rel` over threshold), routed
  through a disadvantaged code path (`joint = FALSE`) that the corrected,
  larger active set makes more likely to trigger, not less. Whether the
  elevated `resid.rel` reflects real coupling outside the lasso-selected
  active set or ordinary bin noise in a still-modest per-target population
  is not determined — a diagnostic built to answer exactly that question is
  in the code but its output was never captured this session. **Still not
  recommended for production or further benchmarking; `fix.my.unmix()`
  remains the production estimator.**
- **Autofluorescence data hygiene — new finding this session (section 10.5),
  independent of any of the three estimators above.** `bd.spectra["AF",]`
  and `cell.spectra["AF",]` in the Aurora benchmark are bit-identical across
  every detector — the same AF spectrum is being used for both particle
  types, which is physically implausible for at least one of them.
  Confirmed to matter in practice: `Spark Violet 538` is strongly collinear
  with this shared AF row (`cos = 0.774`). Regenerating both particle
  types' AF independently via `get.af.spectra()` is recommended as a
  standing data-quality fix regardless of the item below, since correct,
  substrate-specific AF is a precondition for several other things in this
  pipeline (per-cell AF assignment, the hotspot-vs-AF check already used in
  Stage D) even where it doesn't explain a specific bug.
- **A reproducible, wrong-direction correction on Beads — found last
  session (section 10.4), deprioritised this session on Oliver's
  operational judgement ("hasn't been a problem in the past").** Two
  fluorophores (`BUV805`, `Spark Violet 538`) get a held-out-validated,
  non-noisy correction from `correct.unmixing.signatures()`'s dominance-fit
  approach on Beads that moves *away* from ground truth. Four candidate
  mechanisms were ruled out last session; this session's Python
  reproduction — idealised single-spectrum-per-dye, IID Gaussian/Poisson/
  detector noise, both directions, population sizes up to 10x this
  benchmark's own — could not reproduce the flip at all, for either dye, at
  any population size. That rules out a fifth mechanism (generic
  finite-sample bias in the point estimator, independent of real
  structured/correlated residual noise) without identifying the true cause,
  which most likely requires something an idealised IID model cannot
  produce. **Not chased further this session per Oliver's steer; not
  resolved.** See section 12.3.
- **Autofluorescence-removal mechanisms — two now tested for
  `correct.unmixing.signatures()`, neither recommended for that function's
  own background-removal step (section 13).** Per-cell AF extraction
  (`get.af.spectra()`/`assign.af.fluorophores()`, previously mislabelled
  "AF deconvolution") measurably regresses several already-well-corrected
  Cells fluorophores. True AF deconvolution (`get.af.basis()`/
  `deconvolve.af.background()`, the SVD/joint-OLS mechanism, tested for the
  first time this session) is net negative used as a first background-removal
  step and essentially inert used as a second-phase step on top of
  already-corrected spectra — consistent with, and no argument against,
  `fix.my.unmix()`'s own placement of this mechanism as a second-phase
  default. Neither finding is a criticism of either mechanism in its own
  intended role; both are specifically about whether either belongs inside
  `correct.unmixing.signatures()`'s own pipeline, and the answer for both is
  no.
- The larger architecture gap from the previous version of this document
  still stands: the outer restart loop, restricted-design rescue arm, angle
  clamping, plateau retirement, brightest-N refinement, per-cell AF, and
  continuous null-trust scaling all remain diagnostic-only, not merged into
  the package function. Now that section 4 is closed out, this is the single
  largest remaining item.

---

## 2. Files: what to retire, what to keep

### Safe to fold into this document and remove
- `CONTEXT_fix_my_unmix.md`
- `CONTEXT_fix_my_unmix_session_summary.md`
- `CONTEXT_A10_A11_session_summary.md`
- `CONTEXT_fast_diagnostics_session_summary.md`
- `CONTEXT_round3_review_and_convergence_freeze.md`
- `CONTEXT_envelope_slope_openmp_session_summary.md`
- `FixMyUnmix_Progress.Rmd`
- `run_fix_my_unmix_diagnostics_A8.R`, `A9.R`

### Keep — still-live technical content
- `CONTEXT_information_sources.md` (ideas 1 and 5 still unscored against the
  corrected estimator — see section 6)
- `CONTEXT_signature_correction2.md` (identifiability theory and gate
  derivation still current; its §8.2 bug list is resolved and can be struck)
- `run_fix_my_unmix_diagnostics_A11.R` (Stages Z/W/G/G2 still the only tool
  for those questions)
- `correct_spectra_glasso.R`, `run_correct_spectra_glasso_validation.R` — new
  last session; the candidate third estimator and its recovery test. See
  section 5, item 11.
- `run_signature_correction_diagnostics_aurora.R` — the Aurora bead/cell
  benchmark script, Stages A–J from prior sessions plus Stages K–O (held-out
  split stability, repeated-split validation, fitted-direction check,
  anchor/leverage check, AF-substitution test — see section 10). Stages
  K/M/N/O are self-contained (do not call `.diag.dominance.fit()`) and
  reimplement the relevant fragment of production logic directly; Stage L
  calls `.diag.dominance.fit()` with the `n.split.trials`/`min.split.frac`
  parameters (section 10.3).
- `bead_spillover_spread.csv` — **new this session, supplied by Oliver.**
  Same 16x16 layout as `cell_spillover_spread.csv`. Previously absent; its
  absence measurably degraded nuisance-set quality for several Beads
  fluorophores in this session's Python reproduction (section 12.2). Wire
  into `run_signature_correction_diagnostics_aurora.R` and
  `run_footprint_validation.R` alongside the cell one going forward.
- `/mnt/user-data/outputs/sigcorr_python/` (this session) — a Python
  reproduction of `correct.unmixing.signatures()`'s dominance-fit
  algorithm (`corr_sig_core.py`, footprint-restricted;
  `corr_sig_core_baseline.py`, the pre-session unrestricted port), a pooled
  single-stained-control simulator (`sim_pooled.py`), an empirical
  spillover-spread estimator for use when a real one isn't available
  (`est_spillover_spread.py`, not needed for Beads now that a real file
  exists, but still useful for spot-checking one), and driver/diagnostic
  scripts (`run_baseline.py`, `run_final_comparison.py`,
  `diag_bidirectional.py`). See section 12.6 for scope and reuse notes.
  Companion R script `run_footprint_validation.R` runs the two changes this
  environment found against the real production function on the real
  benchmark — see section 7, item 1. **Run this session against real gated
  data; see section 13.1.**
- `small_Concatenated_cells.fcs`, `small_Concatenated_BD_beads.fcs` — **new
  this session, supplied by Oliver.** FSC/SSC main-population-gated versions
  of the existing `Concatenated_cells.fcs`/`Concatenated_BD_beads.fcs`,
  replacing this session's Python reproduction's own uniform downsample.
  `code/load_real_data.py` (Python reproduction) now loads these directly.
- `corr_sig_core.py` (Python reproduction, continued this session) — added
  `get_af_basis()`, `cosine_similarity()`, `calculate_hotspot_matrix()`,
  `deconvolve_af_background()` (true AF deconvolution, `bg_mode =
  "af_deconv"`, section 13.4); renamed the earlier per-cell-library port
  from `bg_mode = "af_deconvolution"` to `"af_library"` to avoid confusion
  between the two mechanisms (section 13.3); fixed a port-fidelity gap in
  `fit_slope()` that had drifted out of sync with production's footprint
  masking (section 13.2).

### Downgrade to "historical, do not rely on"
- `CONTEXT_parameter_sweep_results.md` and `run_fix_my_unmix_parameter_sweep_3.R`
  — both predate the A10/A11 dominance-mask fix; treat cell-substrate verdicts
  as unverified.

### Safe to fold into this document and remove (added this pass)
- `CONTEXT_hypernegative_penalty_and_cluster_audit.md`
- `CONTEXT_hypernegative_penalty_implementation.md`
- `CONTEXT_gating_downsampling_and_af_basis_fixes.md`

  All three are a single continuous thread, fully consolidated into section
  14 below with every file-level claim re-verified against current source.
  Safe to retire as standalone documents; nothing in them is live that
  isn't now in section 14.

### Keep — still-live technical content (added this pass)
- `CONTEXT_pair_restricted_residual_validation.md` — new, unresolved
  analytical thread (section 15), explicitly not yet confirmed against the
  right substrate (Oliver's own caveat, section 15.5). Keep until the
  concatenated-single-stain rerun it calls for has been done.
- `run_hypernegative_cluster_audit.R` — the cluster-centroid audit script
  built alongside the hypernegative-penalty work (section 14.3). **Currently
  contains logic its own originating document describes as superseded**:
  direct reading this pass shows `neg.threshold <- -fit$threshold.matrix.final`
  (the old mirrored boundary) still in the script, not
  `fit$neg.threshold.matrix.final` (the real, directly-measured one) as
  `CONTEXT_gating_downsampling_and_af_basis_fixes.md` describes as already
  applied. See section 14.6 — re-apply before trusting a re-run.
- `run_pair_restricted_residual_validation.R` — the phase-two residual
  check (section 15). **Also contains logic its own originating document
  describes as already fixed**: the current verdict block still collapses
  to "improved"/"worsened" via an `ifelse` chain that checks `delta.j` or
  `delta.k` with `OR`, with no `"mixed"` case reachable, where
  `CONTEXT_pair_restricted_residual_validation.md` describes a fix already
  applied that computes each member's status independently and reports
  `"mixed"` when they disagree — see section 15's opening note. The
  session's own reported results (8 improved / 4 worsened / 1 mixed / 8
  unchanged) are not reproducible from the copy of the script in this
  project as it currently stands.

---

## 3. Production deployment status of the three patches

1. **`gate.on.bias`, default `FALSE`.** Deployed. No adverse effect observed
   in the subsequent null-fit re-runs — this toggle only affects phase two's
   `"bias"` rejection path and is orthogonal to the issue in section 4.
2. **`leakage.margin` (0.05) + seed-matched before/after calls.** Deployed.
   Working as intended: in the real-null re-run, `BV421` cleared the leakage
   gate at a 3.3% held-out-leakage increase that the old unmargined gate
   would have blocked. Effect size on this run was modest (this benchmark
   doesn't reproduce the larger false-positive rate the original
   FACSDiscover fast-diagnostics run found) but directionally confirmed.
3. **Batched pair estimator (`.fix.envelope.slope.batch()` /
   `fix_envelope_truncated_batch_rcpp()`) wired into phase one.** Deployed,
   and now **fixed and validated**. The fix: `fix_envelope_truncated_batch_rcpp.cpp`'s
   per-target loop now buckets every target-negative event into bright
   (source-positive, always kept) and bulk (subsampled down to
   `max_truncated_events` once the total exceeds it, one `std::mt19937`
   stream per target seeded from R), reproducing `select.negative()`'s cap
   instead of always fitting the full negative-selected population.
   `threshold_source`, `seed`, and `max_truncated_events` were added to the
   Rcpp signature and threaded through `.fix.envelope.slope.batch()` and
   `fix.my.unmix()`'s own `max.truncated.events` parameter, which previously
   reached the scalar path only. Confirmed loaded in production:
   `"threshold_source" %in% names(formals(AutoSpectralRcpp::fix_envelope_truncated_batch_rcpp))`
   returns `TRUE`. Validated by a dedicated A/B test
   (`test_batch_bulk_cap_ab.R`, real `a8.cell.spectra` panel, 80,000
   synthetic cells sized so per-pair negative populations run 55,000–65,000 —
   comfortably above the 20,000 cap, unlike either null-fit script, see
   section 4) comparing scalar, the unpatched batch estimator, and the
   patched one head to head:

   | | vs. scalar, max \|diff\| | vs. scalar, median \|diff\| | one-sweep elapsed |
   |---|---|---|---|
   | batch, unpatched | 0.00479 | 5.91e-5 | 0.689s |
   | batch, patched | **0.00044 (11x tighter)** | **2.14e-5 (2.8x tighter)** | **0.425s** |
   | scalar (reference) | — | — | 1.264s |

   `span` (from the binned-envelope block, untouched by this fix) matched
   the unpatched batch estimator's output exactly once the A/B test harness's
   own simplified wrapper was corrected to include that block — see the
   session's working notes for how a first A/B run's span mismatch was
   traced to the test harness, not the fix, before being re-run and
   resolved. The patched estimator is not just correct but faster than what
   was shipping, since capping the bulk means the Huber fit works on ~20,000
   events instead of the ~60,000+ the unpatched version was fitting.

**A fourth patch, previously recorded here as "not yet applied," is
confirmed deployed as of this consolidation pass: `correct.unmixing.signatures()`'s
`max.step` default and the `n.split.trials`/`min.split.frac` parameters**
(section 10). Direct reading of the current `correct_unmixing_signatures.R`
shows `max.step = 0.15` as the default and `n.split.trials`/`min.split.frac`
present as parameters (default `1L`/`0.6`, backward-compatible — `1L`
reproduces the original fixed, un-reseeded split exactly). `n.split.trials`
has not been promoted above its backward-compatible default of `1`, so the
split-robustification behaviour exists but isn't engaged by default —
raising it remains a live option, not something still to be wired in.

**A fifth patch is now deployed and confirmed on real data:
`footprint.frac`/`footprint.min.channels` and the `spillover.spread`-missing
warning** (section 12, confirmed via `run_footprint_validation.R` in
section 13.1). Current `correct_unmixing_signatures.R` (confirmed by direct
reading this session) already carries `footprint.frac = 0.02` as its
default. An intervening proposal to additionally decouple
`spillover.spread` from nuisance-set gating was tested against real gated
data this session and retracted (section 13.1) — it should not be applied.

**A sixth, seventh, and eighth patch, all to `fix.my.unmix()`, are deployed
and confirmed present as of this pass — see section 14 for the full
narrative:**

6. **Hypernegative correction penalty**, redesigned from a soft delta-based
   damper into a hard-reject fraction-based gate: `min.hypernegative.frac`
   (default `0.01`), `max.hypernegative.frac` (default `0.05`),
   `min.hypernegative.events` (default `200L`), scoped to each source's own
   positive population. `coefficient.decay` (default `0.05`, scaled
   continuously by `(1 - trust)`) and `keep.history`/`spillover.history`
   (opt-in per-iteration diagnostic snapshot) are deployed alongside it.
   `thresholds.final`, `neg.thresholds.final`, `threshold.matrix.final`,
   `neg.threshold.matrix.final`, `unmixed.final`, `residual.final`, and
   `dominant.final` are all exposed in the return value. Confirmed present
   by direct reading of `fix_my_unmix.R`. **Not independently re-validated
   against the full 33-fluorophore panel** — see section 14's open items.
7. **`scatter.gate` (default `TRUE`) and `landmark.quantile` (default
   `NULL`).** `scatter.gate = FALSE` skips scatter gating entirely for both
   the unstained and stained samples; `landmark.quantile` is a lighter
   rescue threshold computed once from the unstained file and reused
   unchanged on the stained one. Confirmed present.
8. **`stained.fit.full`**, a pre-downsample copy of the projected stained
   sample, used so `unmixed.final`/`residual.final` reflect every gated
   event regardless of what the (possibly downsampled) phase-one loop
   actually fit on. Confirmed present, **but confirmed only partially
   applied**: phase two's own signature re-measurement (the code that
   decides dominance and feeds `extract.raw.signature()`) still reads from
   the downsampled `stained.fit`, not `stained.fit.full` — see section 14.5,
   item 2 for the open half of this fix.

---

## 4. RESOLVED — the batched pair estimator did not reproduce the scalar
   estimator's bulk-population cap

### Symptom (as originally observed)
Re-running the real/synthetic null-fit comparison with the patched
`fix_my_unmix.R` (batch estimator wired in, bulk cap still missing at that
point) produced a large, substrate-asymmetric shift in how many pairs the
phase-one estimator fitted at all:

| | pre-patch (scalar) | post-patch, batched, uncapped |
|---|---|---|
| synthetic-null `n.fitted` | 213 | **59** |
| real-null `n.fitted` | 242 | 237 |
| synthetic-null `span` (summary quantiles) | ~17.5M–30.4M | ~7.0M–8.6M |

Real barely moved. Synthetic collapsed to roughly a quarter of its previous
fitted-pair count, with span shrinking 2–4x across the board.

### Mechanism, confirmed by reading source directly
`.fix.envelope.slope()` (scalar) and `fix_envelope_truncated_batch_rcpp.cpp`
(batched, pre-fix) both ran the same `max.mask.passes`-bounded iterative
refinement loop with identical `settled`/`valid` logic — that part of the
port was faithful. The difference was what each fit *at every pass*:

- **Scalar**: `select.negative()` caps the near-origin bulk at
  `max.truncated.events` (20,000) on every pass, keeping all bright events
  plus a bounded random sample of the bulk.
- **Batched, pre-fix**: fit the Huber slope on the *entire* negative-selected
  population, every pass, uncapped.

The batched estimator's own (pre-fix) documentation justified dropping the
cap with a classical-leverage argument ("the bulk sits at the origin and
carries no leverage"), true for OLS-style leverage on the fitted
coefficient but not for the Huber M-estimator's iterative MAD-based scale
re-estimation: when the full bulk outnumbers the informative bright tail by
orders of magnitude, the residual scale — and therefore which points get
down-weighted as outliers — is set almost entirely by the bulk, diluting the
tail's relative influence and biasing the fitted slope toward zero. A slope
shrunk toward zero under-corrects target-negativity for bright source
events, directly shrinking observed `span` and failing more pairs' span/
coverage/min-rise gates — the exact pattern seen above.

**Why the original correctness validation missed this**:
`test_fix_envelope_slope_vectorized.R` validated batched-vs-scalar parity on
a 20-fluorophore synthetic test with subsampling *disabled on both sides* —
a deliberate methodological choice to isolate the core math from RNG-driven
sampling differences. That choice meant the validation never exercised the
code path where the two implementations actually differed in practice. Not a
flawed test — a test that answered a narrower question than "is this safe to
deploy on data with large per-pair populations." The test script also
called `.fix.envelope.slope.batch()` with a `max.truncated.events` argument
that didn't exist on the function at the time — it had been written against
the intended, post-fix signature and would not even have run; it should now
run cleanly against the current production code and is worth executing as
an independent cross-check (section 7).

### Resolution
Ported `select.negative()`'s bright/bulk partition and cap into
`fix_envelope_truncated_batch_rcpp.cpp`'s per-pass loop; see section 3, item
3 for the fix mechanics and validation numbers. Production diffs applied to
`fix_envelope_truncated_batch_rcpp.cpp` and `fix_my_unmix_vectorized.R`
(Rcpp signature, `.fix.envelope.slope.batch()`'s call site and formals,
`fix_my_unmix.R`'s phase-one call site) and confirmed loaded.

### Why the null-fit re-run couldn't confirm this directly
Both null-fit scripts size their phase-one input at exactly 20,000 total
events — `run_fix_my_unmix_synthetic_null.R`'s `sn.n.cells <- 20000L`, and
`run_fix_my_unmix_real_null.R`'s `rn.fix.args$downsample <- 20000L` — the
same value as `max.truncated.events`'s default. Since a pair's
negative-selected population can never exceed the file it's drawn from, and
the cap only fires once that population *exceeds* `max.truncated.events`,
**neither script can trigger the bug in the first place, whether or not the
fix is applied.** `ni_total <= n <= 20000 <= max.truncated.events` holds
unconditionally at these settings. This document's original acceptance
criterion for this section — "a synthetic-null re-run's `n.fitted` and
`span` distribution should return to something close to the pre-batch
numbers (213 fitted, span ~17.5M–30M)" — assumed the null-fit re-run was the
right instrument to confirm the fix with; it structurally is not, at either
script's current settings. This was never checked against the actual script
parameters when originally written, only asserted as "very likely." The fix
is validated instead by the dedicated A/B test in section 3 (built
specifically to exceed the cap), which is authoritative for this question.
If a null-fit-framework confirmation is still wanted, `sn.n.cells` and
`rn.fix.args$downsample` would need raising well above 20,000 in both
scripts — optional, not required, see section 7.

The synthetic-null `n.fitted` this session's re-run actually produced (76,
against the 213 originally recorded) is real but is not this issue — see
section 5, item 3.

### Consequence for `sn_pair_comparison.csv`'s `excess` column
The original concern here — that `bias.synthetic` was inflated by the
uncapped estimator's toward-zero bias, making `excess = bias.real -
bias.synthetic` unreliable in magnitude — is moot for the current null-fit
scripts specifically, since (per above) they never exercised the uncapped
path either before or after the fix. The `excess` ranking's shape (same top
pairs, same signs, reproduced across runs) was and remains informative; the
magnitudes were never actually contaminated by this specific bug for these
two scripts, though they may still be affected by whatever is producing the
`n.fitted` gap in section 5, item 3, which is a live question.

---

## 5. Open / unresolved questions, by priority

1. **A reproducible, wrong-direction correction on Beads for `BUV805` and
   `Spark Violet 538` (section 10.4). Deprioritised this session on
   Oliver's operational judgement ("hasn't been a problem in the past") —
   kept as item 1 for reference continuity, not as a priority ranking; see
   section 7 for the current priority order.**
   `correct.unmixing.signatures()`'s dominance-fit
   approach accepts a held-out-validated correction for both — not noisy,
   not a leverage artifact, not explained by AF collinearity — that moves
   the row further from ground truth (`recovered = -0.10`, `-1.57`).
   Independently corroborated for `Spark Violet 538` by the existing Stage H
   ALS estimator (`cos.als.E = -0.861`, essentially matching this session's
   own `cos.slope.true = -0.824`). Five mechanisms tested last session, four
   ruled out; root cause open. **A sixth mechanism (generic finite-sample
   bias in the point estimator, independent of real structured/correlated
   residual noise) was ruled out this session (section 12.3)**: an
   idealised single-spectrum-per-dye, IID-noise Python reproduction could
   not reproduce the flip for either dye, at any population size up to 10x
   this benchmark's own. See section 10 for the full trail through the
   first five mechanisms.
2. **Autofluorescence data hygiene — new this session, real and
   independent of item 1 (section 10.5).** `bd.spectra["AF",]` and
   `cell.spectra["AF",]` are bit-identical across every detector in the
   Aurora benchmark. `Spark Violet 538` is strongly collinear with this
   shared row (`cos = 0.774` against the shared AF, `0.487` against Beads'
   own unstained-derived AF — reduced but still substantial). Regenerating
   both substrates' AF independently via `get.af.spectra()` is recommended
   regardless of item 1, since this session's specific test (swapping in a
   Beads-derived AF direction) did **not** move `cos.slope.true` for any of
   the three fluorophores tested — so this is confirmed as a data problem
   worth fixing on its own merits, not confirmed as item 1's explanation.
3. **`correct.unmixing.signatures()`'s `max.step` ceiling — Stage J run this
   session, mostly resolved (section 10.1). Confirmed applied to production
   as of this pass — see section 3.** `max.step = 0.15` recovers
   `PerCP-eFluor 710` on both substrates and plateaus beyond that.
   `max.step ≥ 0.25` shows a real regression on `NovaFluor Blue 610-30S`
   (already excluded from this benchmark as too dim), arguing against going
   higher than `0.15`.
4. **Held-out split-search robustification (`n.split.trials`,
   `min.split.frac`) — implemented and validated this session, narrow
   practical impact so far (section 10.2–10.3). Confirmed present in
   production as of this pass, still at its backward-compatible default —
   see section 3.** Confirmed to correctly distinguish genuine
   split-to-split noise (one dye, `BUV805` on Cells, 17/40 passing splits)
   from the two dominant failure modes it does *not* help with:
   fluorophores with exactly-zero achievable held-out gain at every split
   (four Cells fluorophores), and fluorophores with a tightly reproducible
   gain that sits just under `min.gain` regardless of split (several Beads
   fluorophores). At this session's `min.split.frac = 0.6`, it did not flip
   any accept/reject outcome in this benchmark; it is backward-compatible at
   `n.split.trials = 1` (still the shipped default) and safe to leave
   available. Worth keeping as a diagnostic lever (e.g. testing a lower
   `min.split.frac` specifically for `BUV805`) rather than as a general
   fix for the stuck-fluorophore problem, which turned out to have other
   causes.
5. **`n.fitted` gap between this session's null-fit re-run and the numbers
   this document previously recorded — real, separate from section 4, not
   yet understood.** This session's fixed-estimator re-run: synthetic-null
   `n.fitted` 76, real-null `n.fitted` 192, against 213 / 242 recorded
   earlier in this document's history. Confirmed **not** attributable to the
   bulk-cap issue (section 4 — neither script can exercise that path at
   either point in time). The convergence-loop-discard bug is confirmed
   fixed in the current `fix_my_unmix.R` (the spillover update at lines
   ~761-771 now runs before the convergence break check). Most likely
   explanation: the compounding effect of `gate.on.bias`/`leakage.margin`
   (deployed since the 213/242 baseline was recorded) tightening or shifting
   the phase-one acceptance boundary — not yet confirmed. Suggested check:
   re-run the synthetic-null with `gate.on.bias`/`leakage.margin` reverted to
   their pre-patch settings, same `sn.seed`, and see how much of the gap
   closes. Low priority; doesn't block anything else in this document.
6. **Tandem-degradation hypothesis: still underpowered, previous session's
   Fisher's test corrected.** `run_fix_my_unmix_synthetic_null.R` computes
   two different tandem flags and this document's prior report used the
   wrong one. `sn.deg$is.tandem` (line 446,
   `grepl("tandem", class, ignore.case = TRUE)` against
   `fluorophore_database.csv`'s `class` column) is correct and is what the
   script's own `fisher.test()` call (line 457) uses. `sn.real.deg$is.tandem.hyp`
   (line 487) is dead code — a hardcoded 9-name list left over from the
   original, already-discredited ad hoc hypothesis (the one this document
   already flagged for wrongly including `RB670`/`RB705`/`RR688`). Lines
   466-494 of that script (the whole `sn.real.deg` block) duplicate the
   `deg.real` computation via a different path and should be deleted — see
   section 7. Recomputed correctly, on the larger accepted set this
   session's leakage-margin patch and working estimator together produced
   (17 of 33 fluorophores, up from 9): 11 of 17 accepted are tandem-class
   (64.7%) against the 57.6% panel base rate; Fisher's exact test odds
   ratio 1.83, p = 0.49, approximate 95% CI on the odds ratio [0.45, 7.41].
   Still no significant enrichment, direction now mildly consistent with
   (not opposed to) the original hypothesis, interval still wide. Confirm
   against the script's own `fisher.test()` output directly rather than a
   hand reconstruction once the dead code above is removed and it's re-run.
7. **Convergence-freeze mechanism (`fx.converge.ratio`) still unimplemented**
   anywhere, including the diagnostic harness. Fully specified, not yet
   built — see the fix-side next step in section 7.
8. **Reproducibility/non-determinism** under the OpenMP-only architecture
   is still untested (single-core vs. multi-core diff). Somewhat narrowed by
   this session's fix: the batched estimator's bulk sampling is now seeded
   per target from R's own RNG stream rather than depending on thread
   scheduling, so at least that source of non-determinism across
   `n.threads` settings is closed. The broader single-core-vs-multi-core
   check is still unrun.
9. **Whether past production `fix.my.unmix()` runs converged in one
   iteration** is still unchecked against any saved real-world result.
10. **Beads' reliance on the univariate fallback path** remains unexplained
    — possibly related to item 1 above (both are Beads-specific,
    small-population effects), not yet connected.
11. **BV650/AF-coupling theoretical risk in the restricted-rescue arm**
    remains watched via `applied.helped`, not resolved either way.
12. **Idea 5 (shear) and idea 1 (hypernegativity) re-scoring** against the
    corrected estimator — still never run (A11 Stages Z/W).
13. **Per-cell AF untested on Beads.** Per-cell AF *extraction* (not
    deconvolution — see section 13.3) has now been tested on Cells
    (section 2 of this session's report, folded into section 13.3): net
    negative there. Still untested on Beads specifically as its own
    question (as opposed to the AF-deconvolution tests in section 13.4,
    which did cover both substrates).
14. **Graphical lasso / neighborhood selection — implemented two sessions
    ago, run for the first time last session, fails its own recovery test
    for a reason not yet identified. See section 9 for the full trail.**
    `correct.spectra.glasso()`
    (`correct_spectra_glasso.R`, ~1200 lines) is a complete, documented,
    exported function performing the same two-phase task as
    `fix.my.unmix()` with a different phase-one estimator. Mechanism: for
    each target fluorophore, its compensated abundance over its own
    target-negative event set is regressed jointly on every other
    fluorophore's abundance at once, with an L1 penalty (coordinate-descent
    lasso path, base R, `.glasso.lasso.path()`) that drives an uninvolved
    fluorophore's coefficient to exactly zero — neighbourhood selection
    (Meinshausen & Bühlmann 2006) applied row by row to the spillover
    matrix instead of `fix.my.unmix()`'s one-source-at-a-time envelope fit
    plus dominance masking. `lambda` is chosen per target by a two-way
    split-half holdout, with the sparser of two disagreeing halves
    preferred within `margin.frac` of the minimum error — the same
    safe-direction bias `fix.my.unmix()`'s target-negative truncation uses.
    Phase two carries only each target's lasso-selected co-fluorophores
    (`active.set`) into `extract.raw.signature()`, rather than the whole
    panel through a ridge penalty. Shares `fix.my.unmix()`'s
    read/gate/background-basis/downsample architecture and most of its
    phase-two acceptance-gate vocabulary (`leakage.margin`, `max.hotspot`,
    `max.angle`, `min.explained`/`max.explained`, `null.fit` bias
    subtraction) directly — parameter names and defaults match where the
    concept is shared. `max.truncated.events` (default 20,000) caps the
    lasso fit's own per-target working set the same way `fix.my.unmix()`
    caps the Huber fit's, for the same cost reason — worth keeping in mind
    given section 4's history if this estimator is ever benchmarked against
    a data source that can exceed that population size. A validation script
    (`run_correct_spectra_glasso_validation.R`) is written and has been run:
    a controlled recovery test, methodologically analogous to
    `run_signature_correction_recovery_test.R` — inject a known spectral
    corruption (`target_j <- normalize(true_j + eps * true_k)`) into a copy
    of a real reference spectra table, simulate a fully-stained sample from
    the *true, uncorrupted* spectra via `sim.flow.data()`, and check whether
    `correct.spectra.glasso()` started from the corrupted spectra recovers
    both the injected coefficient and the row shape it came from. Result:
    both estimators fail to recover the injected coefficients; see section 9
    for the full trail. Understanding *why* is the natural next milestone
    for this estimator, but is currently lower priority than section 5,
    item 1 above, since that failure mode is closer to production
    (`correct.unmixing.signatures()` is the current production candidate;
    glasso is not).
15. **`correct.unmixing.signatures()`'s never-correcting fluorophore
    group — found and fixed two sessions ago, validated in Python, and now
    confirmed in R this session (section 12, section 13.1). Resolved and
    applied.** Emission footprint width (detectors a dye actually emits
    into), not cross-substrate error size, predicts whether the held-out
    gain gate can detect a real correction — confirmed by direct
    computation against the real Aurora spectra CSVs (Spearman ρ = 0.66
    between footprint width and measured gain plateau; the four
    narrowest-footprint dyes in the panel are exactly the four with
    near-zero gain regardless of population size). `footprint.frac`
    (new parameter, restricts the slope fit and held-out objective to a
    dye's own footprint) combined with supplying Beads a real
    `spillover.spread` (previously absent, now supplied) raised the
    accepted count from 10/16 to 11/16 (Cells) and 8/16 to 10/16 (Beads) on
    the real production function against real gated data
    (`run_footprint_validation.R`, section 13.1).
16. **AF deconvolution proper (`get.af.basis()`/`deconvolve.af.background()`)
    as a `correct.unmixing.signatures()` background-removal option — new
    this session, tested, not recommended (section 13.4-13.5).** Net
    negative as a first background-removal step ahead of any spectral
    correction; essentially inert (one small, unexplained exception,
    `eFluor 450`) as a second-phase step on top of already-corrected
    spectra. Consistent with the package's existing placement of this
    mechanism as `fix.my.unmix()`'s second-phase default; no change
    recommended to that architecture. The `eFluor 450` exception is a
    minor open thread, not currently prioritised (section 7).
17. **Phase two's signature re-measurement still reads the downsampled
    `stained.fit`, not `stained.fit.full` (section 14.5, item 2).** Only
    `unmixed.final`/`residual.final` were fixed to reflect every gated
    event; `extract.raw.signature()`'s own input, and therefore the
    corrected spectra themselves, could still be under-informed by the
    same thinned tail on a fluorophore-rich panel. Real compute-cost
    increase to weigh against a fix; not yet tested in isolation with
    downsampling left on.
18. **`run_hypernegative_cluster_audit.R` and
    `run_pair_restricted_residual_validation.R` both still contain logic
    their own originating documents describe as already superseded** — the
    audit script's negative boundary and the residual script's verdict
    logic respectively (section 2, section 14.6, section 15's opening
    note). Neither script's prior reported numbers should be treated as
    reproducible from the copies currently in this project until
    reconciled.
19. **The eosinophil AF-representation gap (section 14.4).** A small
    population reads hypernegative in PE-Cy7/PE against BUV805 under the
    shared, low-rank `af.deconv` basis but is positioned correctly under
    per-cell AF assignment — a case where the shared basis under-represents
    a spectrally distinct population's own background, not a spillover-
    correction error. Natural next candidate for `bg.mode = "per.cell"`,
    alongside the alveolar macrophage population that motivated it
    originally.
20. **`.fix.stratified.sample()`'s quota math on fluorophore-rich panels
    (section 14.5, item 4).** `floor.n = n.fluorophores *
    downsample.min.stratum` easily exceeds the positive budget on a panel
    with many fluorophores, forcing every stratum's retention far below
    its nominal floor by uniform scaling — confirmed sufficient on its own
    to produce a visibly biased-looking result with no boundary logic
    involved. `stained.fit.full` works around this for `unmixed.final`; the
    quota math itself is unchanged and could still bite phase two (item 17
    above) or any other consumer of the downsampled fit.
21. **`som.dim` sensitivity in `cluster.unmixed.events()`-based auditing
    (section 14.3).** 20 vs. 30 gave different cluster boundaries for what
    appears to be the same lung population; not investigated.
22. **Pair-restricted residual validation's specific findings
    (`eFluor 780`, `BUV395`/`BUV805` mutual coupling, `APC`) are not yet
    confirmed against the right substrate (section 15).** Oliver's own
    assessment is that this session's interpretation reads more into the
    results than a multi-colour sample (real biological co-expression) can
    currently support. A rerun against concatenated single-stained controls
    for the same seven-colour panel is the identified prerequisite, and is
    also entangled with a separate, not-yet-reconciled session's rejection
    of the null-fit approach generally (section 15.5) — the calibration
    plan proposed alongside the `eFluor 780` finding should be treated as
    open to revision, not a settled next step.

---

## 6. Information sources to reassess

- **Idea 2 (structural row completion)** — the tandem-donor-axis version
  originally failed per A9/A10 verdicts. The Fisher-test result (section 5,
  item 6, now corrected) is a second, independent line of evidence pointing
  the same general direction at the whole-panel level: still no clearly
  detectable relationship between a dye's tandem status and whether/how much
  its signature gets corrected, though the interval is wide enough not to
  rule out a moderate effect. Doesn't rule out a donor-axis effect at the
  individual-pair level (which idea 2 actually targets), but it's one more
  reason not to prioritize reviving this idea without new evidence.
- **Idea 5 (shear) and idea 1 (hypernegativity)** — unchanged recommendation:
  check overlap against the restricted-design rescue's successes before
  re-running Stage Z/W.
- **Ideas 3 and 4 (spread/VIF, SSI)** — cleanly refuted, remain closed.

---

## 7. Next steps for testing, in order

1. **Closed.** `run_footprint_validation.R` has been run (section 13.1),
   and this consolidation pass confirms the two diffs this item originally
   called for — `max.step` default `0.08 → 0.15` (section 5, item 3 /
   section 10.1) and the `n.split.trials`/`min.split.frac` parameters
   (section 5, item 4 / section 10.2–10.3) — are both present in current
   `correct_unmixing_signatures.R` alongside `footprint.frac`. Nothing
   further required for this item.
2. **Check whether `extract_raw_signature.R`'s `resid.rel` gate has the
   same footprint-dilution property `correct.unmixing.signatures()`'s
   `gain` objective had (section 12.5) — directly relevant to section
   11.3's still-open question.** `resid.rel` (line
   278-279 of that file) is computed the same way the pre-session `gain`
   objective was: a full-panel residual norm over a full-panel signal
   norm, with no restriction to the target fluorophore's own emission
   footprint. This function is phase two for *both* `fix.my.unmix()` and
   `correct.spectra.glasso()` (section 9.1, section 11.3), so if the same
   mechanism applies here, a fix would help both estimators at once — and
   would give a third candidate mechanism (alongside section 11.3's bin
   noise and real-unmodelled-coupling) for why `BUV805`/`BV480` are
   rejected there at `resid.rel = 0.0494`/`0.0558`. Suggested check: for
   the two rejected Discover-panel rows, compute each row's own emission
   footprint width the way this session did for the Aurora panel (section
   12.2's method), and check whether `resid.rel` computed over only that
   footprint clears `0.03`. Not yet tested — the Python reproduction built
   for section 12/13 covers `correct.unmixing.signatures()` only, not
   `extract_raw_signature.R` (section 12.6).
3. **Regenerate `bd.spectra["AF",]` and `cell.spectra["AF",]` independently**
   via `get.af.spectra()`, one call per substrate against that substrate's
   own unstained FCS file, with `use.unmixed = FALSE` (per that function's
   own documented warning about collinear panels — directly applicable
   here) and `deduplicate = TRUE`. Compare the result's row-1 (population
   mean) AF against both the current shared row and the earlier
   `.diag.af.pcs()`-derived stand-in before committing back to the CSVs.
   Independent of item 4 below, but do it anyway (section 5, item 2).
4. **Chase the Beads direction-flip (section 5, item 1 / section 10.4) —
   deprioritised this session on Oliver's operational judgement, kept here
   for completeness rather than as an active recommendation.** If picked
   back up: (a) check whether `BUV395`'s partial `cos.no.anchor` shift in
   Stage N (0.352 → 0.562, the one variant in that stage that moved more
   than noise) is worth a dedicated look; (b) check whether the same
   pattern appears on the FACSDiscover benchmark's Beads context; (c)
   consider whether Beads' much smaller per-pair populations make the
   10-11-bin abundance discretization itself unreliable for these two dyes
   specifically. This session's idealised-noise Python reproduction could
   not reproduce the flip (section 12.3), which rules out one more
   mechanism but doesn't point at what to try next among (a)-(c).
5. **Closed.** The dead `sn.real.deg`/`is.tandem.hyp` block is confirmed
   absent from the current `run_fix_my_unmix_synthetic_null.R` (467 lines
   total, no match for either name) — the housekeeping this item called
   for is done. Re-running and recording the script's own `fisher.test()`
   output directly, superseding the hand-reconstructed number in section 5
   item 6, is still worth doing but is no longer a code-cleanup item.
6. **Find the mechanism behind `correct.spectra.glasso()`/`fix.my.unmix()`
   phase one's failure to recover an injected spillover coefficient in
   OLS-abundance space (section 9, "Open — root cause not identified").**
   Shared by both estimators, not explained by regularization,
   source-panel redundancy, or panel-oblique sensitivity (all three tested
   and disconfirmed). Suggested starting point: derive, rather than
   empirically probe, the expected relationship between an `eps`-sized row
   perturbation and the resulting slope in target-negative OLS-abundance
   space for a known design, then check the two injected pairs against
   that derivation instead of against intuition.
7. **Convergence-freeze test.** Implement `fx.converge.ratio` in
   `run_percell_4_6.R` (fully specified, not yet built), restore both
   substrates to the target loop, run, and check `applied.helped` on
   `converged == TRUE` rows.
8. **Single-core determinism check.** Two identical runs of the fast
   harness, diff `fast_ledger.csv` and `fast_spillover_matrix.csv`. Note the
   batched pair estimator's own thread-count-invariance is now confirmed by
   construction (section 5, item 8) — this check is about the rest of the
   pipeline.
9. **Promote the alternating loop into production `fix.my.unmix()`** —
   sequence after item 7 above, since it depends on the convergence-freeze
   mechanism.
10. **Per-cell AF on Beads**, **Beads univariate-reliability investigation** —
    lower priority, no dependency on the above, possibly connected to
    section 5, item 10.
11. **Optional: chase the `n.fitted` gap** (section 5, item 5) — revert
    `gate.on.bias`/`leakage.margin` on a synthetic-null re-run, same seed,
    and see how much of the 213-vs-76 gap closes. Not required for anything
    else in this document.
12. **Optional: raise `sn.n.cells`/`downsample` above 20,000** in both
    null-fit scripts if a null-fit-framework confirmation of the section 4
    fix (on top of the A/B test already done) is wanted for the standard
    diagnostic harness going forward.
13. **Optional, low priority: chase `eFluor 450`'s small, reproducible
    wrong-direction move under AF deconvolution (section 13.4).** Not
    explained by hotspot collinearity; not currently affecting any
    production recommendation, since `af.deconv` is not recommended inside
    `correct.unmixing.signatures()` regardless. Worth a note if `af.deconv`
    is ever considered for a role closer to that function's own pipeline.
14. **Optional: port `get.spectral.variants()`'s Beads' `gate.main`
    2D-density pre-gate into the Python reproduction** — the leading
    candidate for the residual Beads accepted-count gap between this
    session's Python reproduction and the real R run (8/16 vs. 10/16;
    section 13.2). Not required for any pending R diff, since the real R
    numbers are already the authoritative ones.
15. **Reconcile `run_hypernegative_cluster_audit.R` and
    `run_pair_restricted_residual_validation.R` against their own
    originating documents (section 5, item 18 / section 14.6 / section 15's
    opening note).** The audit script's negative boundary
    (`-fit$threshold.matrix.final` → `fit$neg.threshold.matrix.final`) and
    the residual script's verdict logic (independent `status.j`/`status.k`
    with a reachable `"mixed"` outcome) are both described as already fixed
    in their respective documents but are not present in the copies
    currently in this project. Do this before re-running either for
    anything that will be trusted — this is a prerequisite for items 16 and
    19 below, not an independent item.
16. **Re-run `run_hypernegative_cluster_audit.R` on lung and spleen with the
    hard-reject hypernegative gate live** (section 14.3's original next
    step, still not done), once item 15 above is resolved — specifically
    checking whether cluster 240 (lung, BUV805/PE-Cy7) stops crossing the
    boundary, and whether any genuinely co-expressed pair elsewhere in the
    33-fluorophore Discover panel gets wrongly killed by the new hard cap.
    This is the real test of whether a hard reject was the right design
    choice — sound reasoning, not yet validated against a full panel.
17. **Decide whether to extend the `stained.fit.full` fix to phase two**
    (section 5, item 17 / section 14.5, item 2) — test with downsampling
    left on (for speed) but the fix applied, confirming the corrected
    spectra themselves, not just `unmixed.final`, match a full-data
    reproduction. The confirmed win this session used `downsample = FALSE`
    directly, which bypasses this specific code path.
18. **Revisit `.fix.stratified.sample()`'s quota math** (section 5, item 20)
    for panels with many fluorophores, now that the severity of the
    floor-vs-budget scaling is quantified — either raise defaults or
    protect each stratum's own extreme tail explicitly, mirroring the
    scatter gate's landmark rescue, rather than relying on downstream
    full-data recomputation alone.
19. **Rerun `run_pair_restricted_residual_validation.R` against
    concatenated single-stained controls for the seven-colour lung panel**
    (section 5, item 22 / section 15.5, section 15.7), once item 15 above
    is resolved, and compare directly against the multi-colour results
    already recorded in section 15.3 — the identified prerequisite before
    treating `eFluor 780`, the `BUV395`/`BUV805` coupling, or the `APC`
    finding as more than suggestive.
20. **Investigate the eosinophil AF-representation gap** (section 5, item
    19 / section 14.4) as a candidate for `bg.mode = "per.cell"`, alongside
    the alveolar macrophage population that originally motivated that
    option.

---

## 8. Verified current-state facts (spot-checked against source, not assumed)

- `fix_my_unmix.R`'s convergence-discard bug: fixed, confirmed present in
  the current copy (spillover update at lines ~761-771 runs before the
  convergence break check).
- All old dead-code bugs from `CONTEXT_signature_correction2.md` §8.2:
  absent from current source.
- `gate.on.bias`, `leakage.margin` (with seed-matched `.fix.leakage()`
  calls), `n.threads`, and the batched pair estimator are all present in the
  current production `fix_my_unmix.R`.
- **The batched pair estimator's bulk-cap gap (section 4) is fixed and
  validated.** `fix_envelope_truncated_batch_rcpp.cpp` now reproduces
  `select.negative()`'s bright/bulk partition and cap; confirmed loaded via
  `"threshold_source" %in% names(formals(AutoSpectralRcpp::fix_envelope_truncated_batch_rcpp))`
  returning `TRUE`, and validated by a dedicated 80,000-cell A/B test (see
  section 3, item 3, for the numbers).
- **Neither `run_fix_my_unmix_synthetic_null.R` nor `run_fix_my_unmix_real_null.R`
  can exercise the bulk-cap fix at their current settings** — both size
  phase-one's input at exactly 20,000 events, equal to
  `max.truncated.events`'s default, so the cap can never fire regardless of
  whether the fix is applied. Confirmed by reading `sn.n.cells <- 20000L`
  and `rn.fix.args$downsample <- 20000L` directly.
- Real-null re-run this session (fixed estimator): `n.fitted` 192. Synthetic:
  `n.fitted` 76. Both differ from the 242/213 recorded earlier in this
  document's history; confirmed unrelated to section 4 (see above); cause
  not yet established — section 5, item 5.
- `fluorophore_database.csv`'s `class` column is the authoritative
  tandem/non-tandem source for this panel: 19 of 33 fluorophores are
  tandem-class (`SIRIGEN tandem` or `Protein tandem`); `RB670`, `RB705`, and
  `RR688` are class `Real`, not tandem.
- `run_fix_my_unmix_synthetic_null.R` computes tandem status two different
  ways: `sn.deg$is.tandem` (line 446, correct, `class`-derived) and
  `sn.real.deg$is.tandem.hyp` (line 487, dead code, a hardcoded 9-name list
  from the original discredited hypothesis). This document's prior report
  of the tandem-enrichment result used the wrong one; corrected in section
  5, item 6.
- Tandem-enrichment test on this session's real-null accepted set (n=17,
  using the correct `class`-derived flag): 11/19 tandem accepted (57.9%
  rate) vs. 6/14 non-tandem accepted (42.9% rate); Fisher's exact test
  p = 0.49; accepted-tandem fraction 0.647 vs. panel base rate 0.576;
  approximate 95% CI on odds ratio [0.45, 7.41].
- `correct.spectra.glasso()` (`correct_spectra_glasso.R`) is a complete,
  exported, documented function — not a draft. Its validation script
  (`run_correct_spectra_glasso_validation.R`) has been run, at
  two different `complexity` settings (`1/3` and `0.1`). Both runs: neither
  injected pair (`BUV805 <- BUV395`, `eps = 0.12`; `BV480 <- BV510`,
  `eps = 0.08`) is recovered — `gv.glasso.fit$spillover[source, target]`
  comes back within ±0.005 of 0 both times, and `deg.glasso` (the corrected
  row's angle to ground truth) is bit-identical to `deg.start` (i.e. no
  correction at all) both times. See section 9 for the full diagnostic
  trail and what's been ruled out.
- The validation script itself had a real bug, now fixed: `.sim.scatter()`
  always names its output columns `c("FSC", "SSC")` regardless of the
  cytometer, so `run_correct_spectra_glasso_validation.R` now renames
  `gv.fs.scatter`/`gv.un.scatter`'s columns to
  `flow.control$scatter.parameter` immediately after generating them, and a
  `length(flow.control$scatter.parameter) < 2L` guard was added alongside
  the existing column-match checks. This is a script fix, not a package fix.
- `run_percell_4_6.R`'s target loop still runs `c("Cells")` only.
- **`bd.spectra["AF",]` and `cell.spectra["AF",]` are bit-identical across
  every detector in the Aurora benchmark** (confirmed by direct
  subtraction: all 52 UV1-A…R8-A values are exactly 0). Same AF spectrum in
  use for both particle types.
- **`Spark Violet 538`'s spectral row is strongly collinear with that
  shared AF row**: `cos = 0.774`. Against a Beads-specific AF direction
  derived from Beads' own unstained control (`.diag.af.pcs()`, top singular
  vector), collinearity drops to `cos = 0.487` — reduced but still
  substantial. `BUV805` and `BUV395` are not meaningfully collinear with
  either AF version (`cos ≤ 0.34` in all four combinations tested).
- **Swapping a Beads-specific AF direction into the dominance-assignment
  unmix does not resolve the Beads direction-flip.** Dominance population
  sizes for `Spark Violet 538`/`BUV805`/`BUV395` shift by at most ~2%
  (2374→2372, 1571→1575, 2297→2248), and the fitted-slope-vs-truth cosine
  is essentially unchanged for all three (`-0.824→-0.834`,
  `0.350→0.349`, `0.352→0.338`). This rules out AF-collinearity-via-
  dominance-assignment as the explanation for the direction-flip, though
  the underlying AF data-hygiene issue remains real (section 10.5).
- **`correct.unmixing.signatures()`'s `max.step` ceiling only binds one
  fluorophore in the Stage J sweep of this benchmark** (`PerCP-eFluor 710`,
  both substrates): `rel.step` for the other never-correcting fluorophores
  sits either well under the current `0.08` default or is rejected for a
  different reason entirely (exactly-zero held-out gain, or a
  reproducible-but-wrong-direction accepted step). `max.step ≥ 0.25`
  produces a confirmed regression on `NovaFluor Blue 610-30S`
  (`recovered = -1.06` to `-1.09`).
- **`.diag.dominance.fit()` now accepts `n.split.trials` and
  `min.split.frac`** (added this session, defaults `1L`/`0.6`, confirmed
  backward-compatible: `n.split.trials = 1` reproduces the original fixed,
  un-reseeded split exactly). Not yet ported to production
  `correct_unmixing_signatures.R`.
- **Current production `correct_unmixing_signatures.R` (re-read in full
  this session) already has `footprint.frac = 0.02` and
  `footprint.min.channels = 3L` as defaults, and its `fit.slope()` closure
  restricts the `lm.fit()` call itself to `footprint.mask`-selected
  columns** (`slope <- rep(0, ncol(r.bin)); slope[footprint.mask] <-
  stats::coef(fit.m)[2,]`) — channels outside a dye's own footprint carry
  an exact zero slope, not merely a de-emphasised one in a diagnostic
  statistic. `bg.mode` is `c("global.mean", "scatter.knn", "none")`; no
  `"af.deconv"` option exists in this function. `max.step` is still `0.15`
  in this reading (already raised from `0.08` at some point since section
  5, item 3 was written — confirmed present, not confirmed when it was
  applied).
- **`get.af.basis()` and `deconvolve.af.background()`, in
  `deconvolve_af_background.R`, read in full this session.** `get.af.basis()`
  builds a components-x-detectors basis from the top right singular vectors
  of a raw, uncentred unstained matrix; `n.pc = "auto"` (default) retains
  every component whose singular value exceeds the largest a
  column-independently-permuted version of the same matrix produces, floored
  by an electronic read-noise term (`read.var = 125^2` default). Sign
  convention: each row's largest-magnitude detector reads positive.
  `deconvolve.af.background()` fits `[af.basis; spectra]` jointly by OLS for
  every event, clips only the AF-component coefficients non-negative, and
  subtracts the reconstructed background; an optional `target` argument
  drops AF components whose hotspot coupling to that fluorophore
  (`calculate.hotspot.matrix()`) exceeds `max.hotspot` (default `5`), always
  keeping the leading (mean-background) component. This is the mechanism
  `fix.my.unmix()` uses by default as `bg.mode = "af.deconv"`.
- **`get_spectral_variants.R`, read in full this session, uses the same SVD
  idea independently and for a different purpose.** `af.pcs.list` — one
  per unique universal-negative file listed in the control table, `svd(dat,
  nu=0, nv=4)` with a fixed `nv=4` (not `get.af.basis()`'s own
  `n.pc = "auto"`) — feeds `af.collinear.threshold`'s later
  fluorophore/AF-component collinearity check inside `get.fluor.variants()`,
  not `get.spectral.variants()`'s own dominance-assignment or background
  removal, which use `get.af.spectra()`'s discrete library instead
  (`af.spectra <- get.af.spectra(...)`, passed to `unmix.autospectral()`).
  `get.spectral.variants()` itself therefore uses per-cell AF extraction
  (`get.af.spectra()`), not AF deconvolution (`get.af.basis()`), for its own
  dominance/background step — `af.pcs.list` is a separate, narrower use of
  the SVD-basis idea within the same file.
- **This session's Python port of `correct.unmixing.signatures()`
  (`corr_sig_core.py`) had a port-fidelity gap in `fit_slope()`, now
  fixed (section 13.2).** It previously fit the per-detector regression on
  every detector, restricting only the diagnostic R-squared to
  `footprint_mask` — a deliberate design choice recorded in that file's own
  comments at the time, which turned out to diverge from current
  production's actual behaviour (previous bullet). After the fix, the
  Python reproduction's Cells accepted count under `footprint.frac = 0.02`
  + `spillover.spread` matches the real R run's exactly (11/16); Beads
  remains 2 fluorophores short (8/16 vs. 10/16), with every individual
  fluorophore's direction and approximate magnitude agreeing — see section
  13.2 for the residual-gap discussion.

---

## 9. Prior session: `correct.spectra.glasso()` validation run, the
   `extract_raw_signature.R` fix, and a second, unresolved phase-one problem

This section is the primary record of the session before this one. Read it
before trusting either estimator's phase-one spillover-detection step on
library-error-style corruption.

The `run_correct_spectra_glasso_validation.R` script initially failed with a
column-mismatch error before ever reaching `correct.spectra.glasso()` or
`fix.my.unmix()`: `.sim.scatter()` always names its output columns
`c("FSC", "SSC")` regardless of the cytometer, colliding with
`flow.control$scatter.and.channel.spectral`. Fixed by renaming
`gv.fs.scatter`/`gv.un.scatter`'s columns to `flow.control$scatter.parameter`
right after generating them, plus a `length(flow.control$scatter.parameter)
< 2L` guard. Script-only fix; no package change. See section 8 for the
confirmed-applied note.

### 9.1 The `extract_raw_signature.R` bug — found, fixed, confirmed working

**Symptom.** With the script running, both `correct.spectra.glasso()` and
`fix.my.unmix()` returned `spectra` bit-identical to the corrupted starting
spectra for every row in the panel, not just the two injected ones —
`deg.glasso`/`deg.fix` matched `deg.start` to ~15 significant digits. Both
functions' own `signature.log` showed a candidate row had in fact been
computed and then rejected (`accepted = FALSE`).

**Mechanism.** `extract.raw.signature()` has two branches depending on
whether the population supports a joint fit (`fit.joint`, gated on
`n.levels.use >= n.floor` where `n.floor = length(active) + 3`). In the
non-joint fallback branch, nuisance fluorophores are subtracted from the
raw signal at their *current* spectra before the target is fit
univariately:

```r
y.res <- y.bin
if ( length( nuisance ) > 0 )
  y.res <- y.bin - x.bin[ , nuisance, drop = FALSE ] %*%
    spectra[ nuisance, , drop = FALSE ]
```

But `explained.total`'s denominator (`top.norm`) is taken from the
*un-subtracted* `y.bin`, while its numerator (`fitted`) comes from a model
fit to the *subtracted* `y.res`. Whenever the population has real
co-expression of other panel members and `active` spans many fluorophores,
this numerator/denominator mismatch mechanically biases `explained.total`
low, independent of whether the target's own slope was estimated correctly.
`fix.my.unmix()` always passes `active = fluorophores` (the whole panel), so
it hits this branch hardest; `correct.spectra.glasso()` passes a curated,
usually small `active.set`, so the bias is present but much smaller. This
tracked exactly across the four `signature.log` rows checked that session:
glasso/`BV480` (joint fit, no mismatch) scored `explained.total = 1.004`;
glasso/`BUV805` (non-joint, 6 nuisance fluorophores subtracted) scored
`0.955`; `fix.my.unmix()`/`BUV805` and `/BV480` (non-joint, 32 nuisance
fluorophores subtracted) scored `0.677`/`0.690` — well below the default
`min.explained = 0.8` gate, causing rejection.

**Fix, applied to `extract_raw_signature.R`.** In the non-joint branch,
reconstruct what the whole candidate model (target term plus the nuisance
terms it subtracted) predicts, and use that reconstruction — not the
target-only `fitted` — as `explained.total`'s numerator:

```r
nuisance.fitted <- if ( length( nuisance ) > 0 )
  x.bin[ , nuisance, drop = FALSE ] %*% spectra[ nuisance, , drop = FALSE ] else
  matrix( 0, nrow = nrow( y.bin ), ncol = ncol( y.bin ),
         dimnames = dimnames( y.bin ) )
```

added inside the non-joint branch (both the intercept and no-intercept
sub-branches keep `y.res <- y.bin - nuisance.fitted` unchanged), with
`fitted.total <- fitted + nuisance.fitted` computed at the end of that
branch, and `explained.total`'s numerator switched to
`if (fit.joint) fitted else fitted.total`. `resid.rel` is untouched (it's
correctly scoped to `y.res` already).

**Status: applied and confirmed working.** Patched, sourced, and re-run
that session with `min.bin.events` reverted to its default `50` (see 9.2 for
why that reversion mattered). Result: `deg.fix` for `BUV805` dropped from
`3.097°` (= `deg.start`, i.e. no correction) to `0.212°`; `BV480` dropped
from `2.679°` to `0.245°` — roughly 92-93% of the injected angular error
removed. This is a genuine, working fix, independent of everything else in
this section.

### 9.2 A dead end tried along the way: `min.bin.events` — do not repeat

Before the `extract_raw_signature.R` fix was applied, `min.bin.events` was
temporarily dropped to `10` (from the script's default `50`) as a
zero-code-change way to force `fit.joint = TRUE` for the two corrupted
targets, to confirm the `explained.total` diagnosis without editing source.
**This is a useful one-off diagnostic and should not be used as a
production or default setting.** With `active = fluorophores` (33 in this
panel) and only 419-460 events in a target's dominant population,
`n.levels.use` at `min.bin.events = 10` becomes ~41-46 — just barely above
`n.floor = 36` — leaving only 5-10 residual degrees of freedom for a joint
regression on bin means that are themselves noisy (only ~10 events per
bin). That's a nearly critically-determined regression on noisy inputs: it
introduced spurious `deg.fix` of 0.4-0.9° on rows that started at
`deg.start = 0` (already correct), i.e. it made previously-correct rows
worse, while the two rows that actually needed correcting were *still*
unchanged. `min.bin.events` should stay at `50` in production and in this
validation script; the real fix is 9.1, not a bin-count change.

### 9.3 A second, separate problem in `fix.my.unmix()` phase one

Even with the 9.1 patch working and `deg.fix` looking good, the *pairwise*
phase-one coefficient for both injected pairs is still wrong —
`gv.fix.fit$spillover["BUV395", "BUV805"]` and
`gv.fix.fit$spillover["BV510", "BV480"]` never approached `0.12`/`0.08` in
any run that session, including the run where `deg.fix` recovered well.
`fix.my.unmix()`'s final row correction only looks good because phase two
(`extract.raw.signature()`, called with `active = fluorophores` regardless
of what phase one found) re-measures the row directly from raw
detector-space data and doesn't depend on phase one's spillover matrix
being accurate for this particular correction — phase one's matrix only
decides which events count as target-negative elsewhere in the pipeline.
**This means phase one's own pairwise envelope estimator has a real,
uninvestigated problem, currently masked in this benchmark by phase two's
full-panel design.** It would not be masked for `correct.spectra.glasso()`,
whose phase two only sees phase one's `active.set` (section 9.4), and it may
not be masked in production scenarios where phase two's candidate is
rejected for other reasons and the spillover matrix itself is relied on.

`fix.my.unmix()`'s own `coefficient.log` for the two pairs (`complexity =
1/3` run, `min.bin.events = 50`):

| source | channel | slope.truncated | coverage | trust |
|---|---|---|---|---|
| BUV395 | BUV805 | -0.00076 | 0.218 | 1 (accepted) |
| BV510  | BV480  | -0.00460 | 0.026 | 0 (rejected) |

Two distinct failure modes, not one:

- **`BV480 <- BV510`, `trust = 0`**: rejected by the `coverage` gate
  (`0.026 < min.negative.frac = 0.10`). `coverage` is specifically the
  fraction of the source's brightest third of source-positive events that
  remain target-negative after the fitted slope is removed — the gate meant
  to distinguish real co-expression from spillover. Whether this rejection
  is *correct* (this pair really does look like genuine co-expression to the
  estimator) or a *symptom* of the same underlying estimation problem as the
  next bullet is not resolved — see 9.4, since this pair's raw OLS-abundance
  correlation (`-0.485`, section 9.4) is real and stable but the wrong sign
  and far too small relative to `eps`, which is consistent with either
  explanation.
- **`BUV805 <- BUV395`, `trust = 1`**: passed the coverage gate outright and
  was accepted, but the fitted slope (`-0.00076` against a true `0.12`) is
  wrong by roughly two orders of magnitude and the wrong sign. This is not a
  gate-calibration question — the estimator was given the pair, accepted
  its own answer, and the answer is wrong.

**Open, not yet investigated**: whether `.fix.envelope.slope()`'s Huber
fit and masking-pass loop is itself flawed, or whether (per 9.4) there is
simply very little real signal for it to find in this specific pair's
OLS-abundance space, in which case near-zero is close to the correct
answer given its inputs and the bug is upstream. Section 7, item 6.

### 9.4 `correct.spectra.glasso()`: complete failure to recover either
   injected pair — root cause not identified

**The headline result, unchanged across two very different settings.**
`gv.pair.recovery$recovered.glasso` (read directly from
`gv.glasso.fit$spillover[source, target]`, the same accumulated matrix used
throughout the outer iteration, not an intermediate value) and
`gv.spectra.recovery$deg.glasso` at `complexity = 1/3`:

| target | source | eps | recovered.glasso | deg.start | deg.glasso |
|---|---|---|---|---|---|
| BUV805 | BUV395 | 0.12 | 0 | 3.097° | 3.097° |
| BV480  | BV510  | 0.08 | -0.00301 | 2.679° | 2.679° |

And again at `complexity = 0.1` (an order of magnitude sparser random
co-expression):

| target | source | eps | recovered.glasso | deg.start | deg.glasso |
|---|---|---|---|---|---|
| BUV805 | BUV395 | 0.12 | -0.00042 | 3.097° | 3.097° |
| BV480  | BV510  | 0.08 | -0.00427 | 2.679° | 2.679° |

`deg.glasso` is bit-identical to `deg.start` in both runs — no correction
whatsoever. `convergence.log` shows the outer loop converging in a single
iteration; `lambda.log` shows `n.selected` in the 8-22 range (out of 32
candidates) for every one of the 33 targets in the `complexity = 0.1` run,
including the 31 with no injected error, which is itself worth further
scrutiny but is not, per the diagnostics below, the explanation for these
two pairs' failure.

**Hypotheses tested and disconfirmed, in the order tried — recorded so a
future session doesn't repeat them:**

1. **Wrong corruption basis.** `deg.start` (3.097°/2.679°) tracks `eps`
   sensibly, confirming the injection itself (`target_j <-
   normalize(true_j + eps * true_k)`) landed as intended. Not the cause.
2. **Lasso over-regularization hiding a real coefficient.** Disconfirmed by
   tracing the full coordinate-descent path (`.glasso.lasso.path()`)
   directly, at every one of 40 `lambda` values from `lambda.max` down to
   `lambda.max * 1e-3` — essentially unregularized at the low end. For
   `BUV805 <- BUV395`, `beta.source` never exceeds `-0.00029` at *any* point
   on the path. For `BV480 <- BV510`, `beta.source` is stable across nearly
   the whole path (`-0.0003` to `-0.0021`) rather than being suppressed at
   high `lambda` and released at low `lambda` — i.e. it isn't a
   regularization artifact, it's genuinely what the unregularized data
   supports. Confirmed independently by plain OLS with zero regularization
   at all: `cor = -0.049`, `slope = -0.00016` for `BUV805 <- BUV395`; `cor =
   -0.485`, `slope = -0.00198` for `BV480 <- BV510` — both far below `eps`
   in magnitude, one of them wrong-signed, before the lasso is involved at
   all.
3. **Source fluorophore redundant with the rest of the panel** (i.e. its
   true contribution gets absorbed elsewhere during OLS unmixing before it
   ever reaches the target's negative population). Tested by regressing
   each source's own spectrum on the rest of the panel (excluding target and
   source): `BUV395` R² = 0.299, `BV510` R² = 0.497 against the rest of the
   panel — mid-range, not the near-1 this hypothesis would need.
   Disconfirmed.
4. **Low panel-oblique sensitivity between this specific target/source
   pair** (i.e. `calculate.hotspot.matrix()`'s `sqrt(VIF)`-scale entry for
   this pair is unusually small relative to the panel, so even a real `eps`
   perturbation barely reaches abundance space for this pair specifically).
   Tested directly: panel-wide off-diagonal hotspot median `0.362` (IQR
   `[0.184, 0.667]`); `hotspot["BUV805", "BUV395"] = 1.526`;
   `hotspot["BV480", "BV510"] = 2.437`. **Both pairs are well above the
   panel's typical sensitivity, not below it** — the opposite of what this
   hypothesis needs. Disconfirmed, decisively.

**Where this leaves things: root cause not identified.** Four specific,
independently testable mechanisms were checked and none explain why an
`eps`-sized single-row perturbation produces a target-negative-population
slope roughly 40-100x smaller than `eps` (and wrong-signed for both pairs)
in `unmixed.comp` — for `BV480 <- BV510`, a pair with above-median hotspot
sensitivity, an R² against the rest of the panel that rules out simple
redundancy, and a raw correlation (`-0.485`) too strong and too stable
across the entire regularization path to be noise. Something about how an
`eps * source` perturbation to a single row's L∞-renormalized composition
propagates through OLS/pseudo-inverse unmixing into the compensated
abundance space that phase one reads from is not behaving the way the
`eps`-sized-error framing of this test assumes, and the mechanism is not
yet understood. This is not a lasso bug, not a source-redundancy issue, and
not a low-sensitivity-pair issue. See section 7, item 6 for the suggested
next step.

**Root cause identified two sessions later — see section 11.1.** The
`beta` this section checked (`response = target, regressor = source`)
correctly estimates a quantity that really is close to zero for this
corruption model; that is not a failure to recover, it is the right answer
to that pairing. The recoverable `eps`-sized signal is in the transposed
pairing (`response = source, regressor = target`), which the per-target
loop was already computing correctly for every fluorophore in the panel —
it was simply never read. No lasso, redundancy, or sensitivity mechanism
needed explaining after all; the estimator's own fitted matrix was right
the whole time.

### 9.5 Diagnostic code from that session, for reuse

Kept here rather than in a throwaway script since a future session will
likely need at least one of these again. All three assume
`run_correct_spectra_glasso_validation.R` has been run through section 5
(`gv.glasso.fit`, `flow.control`, `asp`, `gv.spectra.wrong`, `gv.fs.path`,
`gv.un.path`, `gv.af.name` all exist) and that `devtools::load_all()` was
used so `correct_spectra_glasso.R`'s non-exported `.glasso.*` helpers are
callable directly.

**Raw (unregularized) correlation and OLS slope for a target/source pair**,
recomputed from scratch via the same read/gate/project steps
`correct.spectra.glasso()`'s `read.gated()`/`project()` use internally —
useful any time the question is "does a real relationship exist before any
regularization or masking is applied":

```r
gv.diag.read.gated <- function( file.name, gate.polygon = NULL, label ) {

  expr.data <- readFCS( file.name, columns = flow.control$scatter.and.channel.spectral )
  gate.data <- expr.data[ , flow.control$scatter.parameter ]

  if ( is.null( gate.polygon ) )
    gate.polygon <- do.gate(
      gate.data, viability.gate = FALSE, large.gate = TRUE,
      samp = label,
      scatter.and.channel.label = flow.control$scatter.and.channel.label,
      control.type = "cells", asp )

  keep <- which( sp::point.in.polygon( gate.data[ , 1 ], gate.data[ , 2 ],
                                       gate.polygon$x, gate.polygon$y ) != 0 )

  list( data = expr.data[ keep, flow.control$spectral.channel, drop = FALSE ],
        gate = gate.polygon )
}
```

**Full lasso path for one target/source pair** (not just the two-way
holdout's chosen point), via `.glasso.lasso.path()` and
`.glasso.neighborhood.row()` directly — useful for telling "regularization
hid it" apart from "it isn't there at any regularization level."

**Panel-wide hotspot check** for a target/source pair against the panel's
own off-diagonal distribution, via `calculate.hotspot.matrix()` on the
actual (corrupted) design being unmixed with — useful for testing the
panel-oblique-sensitivity hypothesis on any future pair before relying on
it as an injection choice.

Full code for all three, as actually run that session, is in that session's
chat log rather than reproduced in full here; ask if a standalone script
consolidating them is wanted.

---

## 10. This session: Stage J run, held-out split diagnostics (Stages K–L),
    a Beads-specific direction-flip (Stages M–N), and an AF hygiene finding
    that didn't explain it (Stage O)

This section is the primary record of this session. The starting question
was narrow — "does raising `max.step` fix the fluorophores that don't
correct in Stage D/J's Aurora bead/cell benchmark?" — and the honest answer
turned out to be "only one of them," which opened up the harder question
this section spends most of its length on. As with section 9, several
plausible-looking mechanisms were tested and ruled out before landing on
the two genuinely open findings (10.4, 10.5's disconfirmed connection to
10.4).

### 10.1 Stage J: `max.step` sweep, run for the first time

Swept `max.step ∈ {0.08, 0.15, 0.25, 0.50}` against `.diag.dominance.fit()`
on both Cells and Beads.

**Result: `max.step` is the binding constraint for exactly one
fluorophore**, `PerCP-eFluor 710` (both substrates) — `rel.step` sits just
above the `0.08` default (`0.084` Cells, similar order on Beads) and clears
at `0.15`. `recovered` jumps from `0` to `0.43` (Cells) / `0.50` (Beads,
needing two accepted steps rather than one) and plateaus there; `0.25` and
`0.50` add nothing further for this dye.

**Regression risk at higher `max.step`**: `NovaFluor Blue 610-30S` (already
excluded from the accept-criteria discussion as too dim to trust) shows
`recovered` going from `0` to `-1.06`/`-1.09` once `max.step ≥ 0.25` lets a
large, wrong step through unmodified — exactly the failure mode this
document's own Stage J interpretation note anticipated. This argues for
`max.step = 0.15`, not something more aggressive.

**The remaining never-correcting fluorophores are not a `max.step`
story at all.** Five on Cells (`Spark UV 387`, `BUV395`, `BUV805`,
`eFluor 450`, `FITC`) and a substrate-shifted but analogous set on Beads
show `n.steps = 0` identically across the entire grid, with `rel.step`
for several of them sitting *under* the current `0.08` default already —
raising the ceiling cannot help a row whose held-out step search never
proposes a nonzero step in the first place. This became the rest of the
session's investigation.

**Recommendation: `max.step` default `0.08 → 0.15`** in both
`.diag.dominance.fit()` and production `correct_unmixing_signatures.R`.
Diff written, not yet applied to production (section 3).

### 10.2 Stages K/M/N/O: two real bugs in this session's own diagnostic code

Before any conclusion could be trusted, two mistakes in this session's
first-draft diagnostic reimplementations were found and fixed — both worth
recording as a methodology note for future sessions, since neither was a
production or existing-diagnostic-function bug:

1. **Flat thresholds instead of spillover-spread-widened thresholds.**
   Stage D/J's actual calls to `.diag.dominance.fit()` (and production
   `correct.unmixing.signatures()`) always pass `spillover.spread`, so
   co-activity (`above`, and therefore the nuisance/`active` set) is always
   computed via `get.spread.thresholds()`, not a flat `>` comparison. This
   session's first Stage K draft used a flat comparison. Confirmed to
   matter: three Beads fluorophores (`BV711`, `RB705`, `PerCP-eFluor 710`)
   have `n.nuisance = 1` in Stage D's real log, which the flat-threshold
   version could not reproduce. Fixed in all subsequent stages (K, M, N, O)
   by computing `above` via `get.spread.thresholds()` whenever
   `ctx$ss.mat` is available, matching production exactly.
2. **Unseeded background-reference resampling.** `.diag.bg.subtract()`
   subsamples its unstained reference pool to `max.reference` (default
   50,000) whenever the pool exceeds that — true for Cells' ~102,000 gated
   unstained events — via an un-reseeded `sample()`. This session's first
   Stage K draft called it once, fresh, right after its own `set.seed()`,
   producing a different reference draw than whatever Stage D's original
   run happened to use, then reused that single draw identically across
   all 40 trials for every Cells fluorophore. This is indistinguishable
   from split-noise in a naive read of the output and was the direct cause
   of `PerCP-eFluor 710` showing 0/40 in the first draft despite Stage D
   accepting it outright (`t.hat = 1.0`). Fixed by raising
   `max.reference` for this diagnostic's call so no subsampling occurs.

**General lesson, worth restating for any future standalone diagnostic
reimplementation of a fragment of production logic**: mirror every
gate-relevant parameter (not just the ones that seem obviously relevant to
the question being asked) and every seeding/RNG dependency exactly, or a
diagnostic bug can produce a result that looks exactly like the phenomenon
under investigation. Both bugs above were only caught by cross-checking
`n.nuisance` and specific dye outcomes against Stage D's own already-printed
log, not by the new stage's output looking implausible on its own.

### 10.3 Stages K (corrected) and L: held-out split noise explains one
    fluorophore, not the pattern

With both bugs fixed, Stage K re-drew each fluorophore's 50/50 held-out
split 40 times and recorded the resulting `t`/`gain` distribution.

**The diagnostic tell**: under genuine split-to-split noise, unanimous
results (0/40 or 40/40 passing) across 40 independent draws are
vanishingly unlikely (`P ≈ 2 × 0.5^40` at a 50/50 true rate). Seeing that
pattern on 6/7 Cells fluorophores and 7/8 Beads fluorophores tested meant
most of them are close to *deterministic* under resampling, not noisy —
the opposite of what the original "single unlucky split" hypothesis
predicted.

**Cells**: four fluorophores (`Spark UV 387`, `BUV395`, `eFluor 450`,
`FITC`) show `gain = 0.000000` at *every one* of 40 independent splits —
not usually zero, exactly zero — meaning no candidate step in `step.grid`
ever beats the baseline for these dyes, regardless of split. `BUV805`
(17/40 passing, real gain up to `0.0089` in the upper quartile) is the one
genuine split-noise case in this benchmark. `PerCP-eFluor 710` (40/40,
confirming the earlier bug) and `BV711` (40/40, positive control) behave
as expected once the bugs above were fixed.

**Beads**: four fluorophores (`Spark UV 387`, `BUV395`, `eFluor 450`,
`RB705`) show a *tightly reproducible* gain that sits just under
`min.gain = 0.002` at essentially every split (e.g. `eFluor 450`:
`0.00187`–`0.00193` across all 40 draws) — narrowly and consistently
rejected by the gain gate, not by `t.hat = 0`, and not helped by
resampling since there's no real split-to-split variability to average
over. `BV711` (`n.nuisance = 1` once correctly computed) shows genuinely
zero gain once its real nuisance term is included — a different dye's
"stuck" reason than the flat-threshold bug's first draft suggested. `PE`
(positive control) and, separately, `BUV805`/`Spark Violet 538` behave
as described in 10.4 below.

**`n.split.trials`/`min.split.frac` added to `.diag.dominance.fit()`**
this session (defaults `1L`/`0.6`, confirmed backward-compatible) and
validated via Stage L (`n.split.trials = 15`, `min.split.frac = 0.6`)
against the same ground truth Stage D/J use. Result matches the Stage K
analysis exactly: the four zero-gain Cells fluorophores and the four
sub-threshold-gain Beads fluorophores are unchanged (`recovered = 0` for
all of them, identical to the single-split baseline); `BUV805` (Cells)
stays rejected because 17/40 (42.5%) doesn't clear the `0.6` threshold —
lowering `min.split.frac` specifically for this one dye is a plausible,
low-priority follow-up, not a general fix.

**Conclusion: split noise is real but narrow.** It explains one
fluorophore (`BUV805`, Cells) and is a validated, safe-to-leave-available
tool, but is not the explanation for the broader never-correcting group,
which splits into "genuinely no detectable signal" (Cells) and "signal
consistently just under threshold" (Beads) — neither of which repeated
splitting can address, since there's nothing split-dependent to average
over in either case.

### 10.4 Stages M/N: a reproducible, wrong-direction correction on Beads —
    unresolved

Two Beads fluorophores don't fit either pattern above: `BUV805` and
`Spark Violet 538` show 40/40 passing splits with *tight, comfortably
above-threshold* gain (`~0.0101`, `~0.0094`) — genuinely, reproducibly
accepted — and yet the accepted correction moves the row **away** from
ground truth: `recovered = -0.10` (`BUV805`) and `-1.57`
(`Spark Violet 538`), unchanged by `n.split.trials = 15` (Stage L).

**Stage M (fitted-direction check)**: computed the cosine between the
fitted slope's direction and the true error direction (`s.true - s.start`)
directly.

| fluorophore | cos.slope.true | cos.slope.meanbright | best.other.match | cos.slope.best.other |
|---|---|---|---|---|
| BUV805 | 0.350 | 0.023 | Spark UV 387 | -0.237 |
| Spark Violet 538 | **-0.824** | 0.091 | PerCP | 0.117 |
| BUV395 (contrast) | 0.354 | 0.084 | Spark Violet 538 | 0.135 |

`Spark Violet 538`'s fitted correction is not just imprecise, it's close
to a clean reflection of the true direction. `cos.slope.meanbright` near
zero for all three rules out a generic brightness/scale confound;
`cos.slope.best.other` low and physically implausible (`PerCP` for a
violet-excited dye) rules out an unrecognised co-varying partner not
caught by `nuisance.frac`. `BUV805`/`BUV395` show a milder version of the
same pattern — positive but far from 1, meaning most of the fitted slope's
magnitude is orthogonal to the true error direction.

**Independent corroboration, already in this document**: Stage H (a
completely different estimator — ALS, not the binned regression Stage M
uses) reported, for the same dye on the same substrate,
`cos.als.E = -0.86135` — within 0.04 of Stage M's `cos.slope.true = -0.824`,
from unrelated code. Two independently-implemented methods agreeing this
closely is strong evidence the anti-alignment is a real property of this
dye/substrate pair, not a bug in either estimator. Worth noting: Stage H's
own interpretation note predicted a *different* signature for a genuine
problem (low `cos.dir.min` with high `drift.in.delta`, diagnosing tandem
variant-mixture drift). `Spark Violet 538`'s `cos.dir.min = 0.968` is high,
not low — so whatever this is, it is a *different* failure mode than the
one that note was written to catch.

Also notable: `Spark Violet 538` corrects cleanly on **Cells** in this same
Stage J run (`recovered = 0.538`, correct direction) — the failure is
specific to the Beads population, not a broadly bad reference row.

**Stage N (anchor-leverage and aggregate-contamination check)**: refit
each slope three ways — baseline, without the single background-anchor
bin, and after trimming the brightest 2% of events (a plausible
doublet/aggregate proxy, since beads are more prone to this than cells).

| fluorophore | cos.baseline | cos.no.anchor | cos.trim.top2pct |
|---|---|---|---|
| BUV805 | 0.350 | 0.401 | 0.326 |
| Spark Violet 538 | -0.824 | -0.809 | -0.825 |
| BUV395 | 0.352 | 0.562 | 0.307 |

Neither variant moves `Spark Violet 538` or `BUV805` meaningfully — the
anti-alignment is a property of the bulk of the assigned population, not a
few high-leverage points or a handful of aggregate events. **Both
hypotheses ruled out for these two dyes.** `BUV395`'s `cos.no.anchor`
moved further than the others (0.352 → 0.562) — the one data point across
Stages M/N/O that didn't stay flat, and not yet followed up (section 7,
item 4).

### 10.5 Stage O: AF data hygiene confirmed real, disconfirmed as the
    explanation

Prompted by the observation that `Spark Violet 538` is visibly collinear
with autofluorescence, and that `bd.spectra["AF",] - cell.spectra["AF",]`
is exactly zero across every detector — the same AF spectrum is in use for
both particle types, which cannot be correct for both.

**Confirmed real and quantified**: `cos(AF, Spark Violet 538) = 0.774`
against the shared (wrong-for-at-least-one-substrate) AF row. Against a
Beads-specific AF direction derived from Beads' own unstained control
(`.diag.af.pcs()`, the top singular vector of the raw uncentred unstained
matrix — the same lightweight technique already used elsewhere in this
script's Stage D preamble), collinearity drops to `0.487` — reduced, but
still substantial. `BUV805`/`BUV395` are not meaningfully collinear with
either AF version (`cos ≤ 0.34` throughout).

**Mechanism reasoned through before testing**: `af.name` is always excluded
from `panel`, so AF never enters the per-dye restricted regression's
`active` design directly. The only channel through which a wrong AF row
could corrupt this pipeline is the full-panel unmix used for dominance
assignment and thresholding (`unmixed <- unmix.ols.fast(raw.data,
spectra)`), which does include the AF row. If AF and `Spark Violet 538`
are collinear and AF is wrong for Beads, that unmix's collinear subspace
is exactly where OLS is least stable — a plausible route to systematically
misassigning which events count as `Spark Violet 538`-positive vs.
background, which would explain a population-wide (not leverage-driven)
contamination consistent with Stage N.

**Tested directly — disconfirmed as the explanation for 10.4.** Swapped
Beads' own AF direction into `ctx$spectra["AF",]`, recomputed thresholds,
recomputed dominance assignment, and reran the Stage M direction check
with only that one row changed:

| fluorophore | n.old (shared AF) | n.new (own AF) | cos.slope.old | cos.slope.new |
|---|---|---|---|---|
| Spark Violet 538 | 2374 | 2372 | -0.824 | **-0.834** |
| BUV805 | 1571 | 1575 | 0.350 | 0.349 |
| BUV395 | 2297 | 2248 | 0.352 | 0.338 |

Dominance population sizes shift by at most ~2%, and `cos.slope.new` is
essentially unchanged for all three — if anything, marginally worse for
`Spark Violet 538`. The corrected-AF full-panel unmix does not produce a
meaningfully different dominance assignment or fitted direction. **The
AF-collinearity-via-dominance-assignment route, reasoned through above, is
ruled out as the cause of the direction-flip.**

**What this leaves**: AF data hygiene remains a real, independent
data-quality problem — using an unmeasured, substrate-mismatched AF is
wrong regardless of whether it explains this particular symptom, and it's
a precondition for other parts of this pipeline (per-cell AF assignment,
the AF-vs-dye hotspot check already used in Stage D) even where it isn't
the cause of a specific bug. Regenerating both substrates' AF properly via
`get.af.spectra()` (with `use.unmixed = FALSE`, per that function's own
documented warning about collinear panels — directly applicable here) is
recommended as a standing fix (section 7, item 3), independent of section
10.4's still-open question.

### 10.6 Where this leaves things

**Resolved**: `max.step` — one fluorophore genuinely capped, `0.15`
recommended, diff written and not yet applied.

**Resolved, narrow**: held-out split noise — real for one Cells
fluorophore (`BUV805`), validated tool (`n.split.trials`/`min.split.frac`)
added to the diagnostic function, not a general explanation for the
never-correcting group.

**Confirmed real, independently worth fixing**: AF data hygiene — shared,
unmeasured AF row between substrates, real collinearity with `Spark Violet
538` specifically.

**Open, highest priority for next session**: the Beads
`BUV805`/`Spark Violet 538` direction-flip. Four mechanisms tested and
ruled out with direct evidence (split noise, background-anchor leverage,
aggregate/tail contamination, AF-collinearity-via-dominance-assignment).
Corroborated independently by an unrelated estimator (Stage H's ALS,
`cos.als.E`) for `Spark Violet 538`. Not yet explained by anything tested
this session. See section 7, item 4 for suggested next steps — in
particular, `BUV395`'s partial movement under Stage N's no-anchor variant
is the one loose thread across Stages M/N/O worth following up before
generating new hypotheses from scratch.


---

## 11. This session: root cause of section 9.4 found and fixed; a second,
    unresolved problem stops it from reaching the function's output

Follow-on to section 9. Section 9.4 left `correct.spectra.glasso()`'s
complete failure to recover either injected pair unexplained after four
disconfirmed hypotheses. This session found the actual cause using a
Python reproduction of just the phase-one estimator (not ported: phase
two's `extract.raw.signature()` gating, which is where the session's
second, unresolved problem turned out to live — see 11.4).

### 11.1 Root cause: the recovery check, and `active.set`, were both
    reading the wrong cell of a matrix that was never wrong

Corrupting fluorophore `j`'s reference row with `eps` of fluorophore `k`'s
shape (`spectra.wrong[j,] <- normalize(spectra[j,] + eps*spectra[k,])`)
does not produce a positive coefficient of `k` on `j` — the pairing
section 9.4 checked, and the pairing both estimators' phase-one loops are
built to detect for ordinary spillover. It produces a *negative*
coefficient of `j` on `k`: unmixing is a change of basis, and if `S.wrong`
still spans the same row space as `S.true` (a linear combination of two
existing rows still is one), then noise-free signal from real `j`-only
cells decomposes in the corrupted basis as `+1/c` on `j`'s own channel and
`-eps` on `k`'s — because `S.wrong[j]` now contains `eps` of `k`'s shape,
and that borrowed component has to be subtracted back out to represent
pure `j` correctly. Verified both algebraically (noise-free, exact to five
decimal places: `M[source,target] = 0`, `M[target,source] = -eps`) and
against noisy simulated data regressing the reversed pairing (`slope =
-0.106` against true `-0.12`, `-0.077` against true `-0.08`).

This is not symmetric confounding the way positive-coefficient spillover
is with genuine co-expression — a marker-negative population cannot read
below zero from real emission, so a negative coefficient is unambiguous
(the same hypernegativity argument in `CONTEXT_information_sources.md`
section 4, here derived mechanistically for library-error-style
corruption specifically rather than argued generally).

Critically, **the per-target lasso loop, run over every fluorophore in the
panel as intended, was already estimating this correctly** — the `-eps`
value shows up in the fitted matrix at `spillover[target, source]` (row =
the corrupted row, column = the donor), filled in automatically while the
*donor* is being fit as its own target. Nothing about the lasso, the mask
passes, the gating, or the outer convergence loop was wrong. What was
wrong: `active.set[[j]]` was built only from `j`'s own row-fit (a negative
coefficient found there was folded in as if it meant genuine spillover
into `j`), never from other fluorophores' fits that named `j` as a
negative-coefficient source — so the true donor of a corrupted row's own
error never made it into that row's own signature re-measurement. The
validation script's recovery check had the identical bug: it read
`spillover[source, target]`, which is correctly ~0, and reported that as
failure to recover.

### 11.2 Fix applied and confirmed at the pairwise level

`correct_spectra_glasso.R`: `active.set` construction restructured. A
fluorophore's own row-fit now contributes only its *positive* survivors
(`own.set`, `b > 0` instead of `b != 0` — a negative own-fit coefficient
was never evidence about that fluorophore's own row). A new pass
(`reversed.donors`) scans every fit for negative coefficients and routes
each one to the *named source's* active set, not the fit's own target.
`active.set` is their union. Returned `reversed.donors` alongside it as a
named-list diagnostic. `run_correct_spectra_glasso_validation.R`: recovery
check reads `spillover[target, source]` (expects `-eps`) instead of
`spillover[source, target]`, reports the old cell alongside for contrast,
and adds a direct `donor.in.active.set` boolean.

Confirmed in R, `complexity = 0.1`, `n.cells = 100000`:

| target | source | eps | recovered.glasso | other.cell.glasso | donor.in.active.set | recovered.fix | other.cell.fix |
|---|---|---|---|---|---|---|---|
| BUV805 | BUV395 | 0.12 | -0.1089 | -0.0004 | TRUE | -0.1073 | -0.0004 |
| BV480  | BV510  | 0.08 | -0.0807 | -0.0045 | TRUE | -0.0795 | 0.0000  |

Both estimators recover within ~10% of the true coefficient once the
correct cell is read, and agree with each other closely — strong evidence
this is a general property of the identification strategy both share, not
specific to `correct.spectra.glasso()`'s lasso. **This also very likely
explains section 9.3's finding** that `fix.my.unmix()`'s own pairwise
envelope estimator produced a coefficient "wrong by two orders of
magnitude" for the identical `BUV805 <- BUV395` pair — 9.3 checked
`spillover["BUV395","BUV805"]`, the same cell this section now shows is
genuinely close to zero. Not directly re-confirmed at 9.3's own
`complexity = 1/3` setting, so stated as likely rather than closed; the
check is one line if anyone wants it settled outright.

### 11.3 Unresolved: the fix does not reach `spectra`

`gv.spectra.recovery$deg.glasso` is bit-identical to `deg.start` for
`BUV805` (`3.097°`) and `BV480` (`2.679°`) in every run this session,
including after raising `gv.n.cells` from `20000` to `100000` (a change
that, by itself, tells us nothing about whether the underlying gate value
moved — a rejected row leaves `spectra.new` untouched, so its angle to
ground truth is definitionally identical to `deg.start` regardless of how
close the rejected fit actually came). `fix.my.unmix()`, same corrupted
spectra, same synthetic data: `deg.fix` for the same two rows dropped to
`0.146°`/`0.150°`.

`gv.glasso.fit$signature.log` for both rows (`n.cells = 20000` run):
`reason = "fit"`, `resid.rel = 0.0494` (`BUV805`) / `0.0558` (`BV480`)
against `max.resid = 0.03`. Both `joint = FALSE`: `BUV805`'s corrected
active set has 9 members (`n.floor = 12`), but its dominant population
(433 events at `n.cells = 20000`) caps `n.levels.use` at `floor(433/50) =
8`, below `n.floor`, forcing the non-joint sequential-subtraction branch
instead of the joint ridge fit. `BV480`: 13 active members, `n.floor = 16`,
same shortfall.

Two candidate mechanisms, not yet distinguished:

1. **Bin noise.** ~54 events/bin at 8 bins, with real detector noise, may
   simply be enough to push `resid.rel` past a tight `0.03` bar on its
   own, independent of any missing signal.
2. **Real, unmodelled coupling.** `fix.my.unmix()` always subtracts *all*
   32 other fluorophores at their current spectra before fitting;
   `correct.spectra.glasso()`'s sparse active set subtracts only the 8-13
   the lasso selected. Real weak cross-talk from fluorophores correctly
   left out of that selection would inflate `resid.rel` for the sparse
   estimator specifically and not for the always-everything one.

A diagnostic (`resid.rel.full`) was added to `correct_spectra_glasso.R`'s
phase-two loop specifically to distinguish these — an extra,
logging-only `extract.raw.signature()` call with `active = fluorophores`
on rows rejected for `"fit"`, never affecting acceptance. **Confirmed
present in the current file. Its value was never reported this session**,
so mechanism 1 vs. 2 remains open. Raising `n.cells` 5x did not move
`deg.glasso` (see caveat above about why that specific observation is not
informative either way), so the population-size story is at best
unconfirmed, not supported.

### 11.4 Process note: scope of the Python reproduction

The Python reproduction built this session (synthetic data generator,
`unmix.ols.fast` port, and a faithful port of `.glasso.lasso.path()` /
`.glasso.select.lambda()` / `.glasso.neighborhood.row()`) covered phase
one only, and it did what it was for: found 11.1's mechanism from first
principles, and the R port of the resulting fix was confirmed working
against that same mechanism on the first R run (11.2). It did not cover
`extract_raw_signature.R` — the joint/non-joint branching, `n.floor`
logic, or gate calibration that 11.3's problem lives in — so when that
second problem surfaced, there was no fast-iteration environment for it,
and diagnosing it fell back to a slower loop of hypothesizing from code
and asking for R output. If this is picked back up, porting
`extract.raw.signature()`'s two branches (or at minimum a standalone
repro of the non-joint sequential-subtraction path, since that's
specifically where both rejected rows are landing) before touching R
again would avoid repeating that.

### 11.5 Where this leaves things

**Resolved**: section 9.4's "root cause not identified" — it was a
matrix-reading direction bug, not an estimation problem, and it is fixed
in both the production function and the validation script. Likely also
resolves section 9.3 (not independently re-confirmed at that section's own
settings).

**Open, next step already built and waiting**: read `resid.rel.full` vs.
`resid.rel` in `gv.glasso.fit$signature.log` for `BUV805`/`BV480` from the
next run. If `resid.rel.full` clears `0.03`, the problem is active-set
design (sparse-but-correct is leaving real coupling unmodelled) and the
fix is a design decision — whether phase two's nuisance removal should
stay sparse-selected or go back to ridge-inclusive-of-everything the way
`fix.my.unmix()` does it, independent of phase one's identification.
If it doesn't clear `0.03` either, the problem is closer to `max.resid`
being uncalibrated for populations this size, or the populations
themselves being smaller than realistic production controls would supply.

**Function status unchanged**: `correct.spectra.glasso()` is still not
recommended for production or further benchmarking. The bug fixed this
session was real and is now closed, but it was never the bug standing
between this function and usable output.

---

## 12. This session: a Python-based troubleshooting workflow applied to
    `correct.unmixing.signatures()`, a general mechanism found behind
    several of section 10's never-correcting fluorophores, fixed and
    validated in Python (not yet in R)

Follow-on to section 10, and methodologically a sibling of section 11 (same
Python-reproduction-first approach, applied to a different function).
Scope: `correct.unmixing.signatures()`'s dominance-fit algorithm only —
`fix.my.unmix()` and `correct.spectra.glasso()` were not touched this
session, though section 12.5 discusses what does and doesn't transfer.

### 12.1 The workflow, for reuse

The pattern that worked, in order:

1. **Read the production R function in full before writing any Python**,
   not just its docstring — `correct_unmixing_signatures.R`'s dominance
   assignment, background subtraction, per-fluorophore binned regression,
   held-out step search, and all five acceptance gates were ported
   line-for-line (`corr_sig_core.py`), with variable names tracking the R
   source (`spectra_new`, `x_bin`, `r_bin`, `t_hat`, `bg_align`) so a diff
   against the R function stays direct. Checked against the diagnostic
   script's own copy (`.diag.dominance.fit()` in
   `run_signature_correction_diagnostics_aurora.R`) for consistency, since
   that copy is what section 10's findings were generated from.
2. **Scope-check what the target function actually calls before assuming
   auxiliary files are needed.** `cluster_unmixed_events.R` and
   `extract_raw_signature.R` were initially flagged as required inputs but
   turned out not to be called by `correct.unmixing.signatures()` at all —
   that function has its own inline binning logic. Confirmed by reading
   the source, not by assumption, before building anything.
3. **Build the simulator to match the real benchmark's actual data
   structure, not a generic one.** The real Aurora benchmark concatenates
   single-stained control files (one fluorophore per file) plus one
   unstained file — not a fully-stained multi-positive sample. The
   simulator built this session (`sim_pooled.py`) reproduces that
   structure specifically (per-dye titration blocks plus an unstained
   block), reusing the noise pipeline already validated in the prior
   glasso session's `sim_flow_data.py` rather than rebuilding it.
4. **Validate the simulation reproduces known qualitative patterns before
   trusting it to test anything new.** Before touching the algorithm, the
   unmodified baseline reproduction was checked against section 10's own
   findings: the fluorophores it left uncorrected matched the real
   benchmark's never-correcting group closely enough (by cross-substrate
   angle, computed directly from `Cells_spectra_corrections.csv` /
   `BD_spectra_corrections.csv`) to trust the environment for further
   testing.
5. **Test the specific mechanistic hypothesis directly against real
   project data, not synthetic proxies for it.** Oliver's wavelength/
   footprint hypothesis (this section) was checked by computing each real
   fluorophore's actual emission footprint from the real spectra CSVs and
   correlating it against measured gain, before writing any fix.
6. **Sweep the candidate fix's own parameter and look for where it stops
   being safe, not just where it starts working.** `footprint.frac` was
   swept 0–0.05; 0.02 was chosen specifically because 0.05 showed a
   measurable regression on `PE` (an already-good correction) while 0.02
   did not, on any fluorophore, in either direction.
7. **Hand back working code plus a same-shaped R script that exercises the
   real production function on the real benchmark**, not just the Python
   result — the R side is the only side that can confirm anything, given
   this session's noise model is idealised (12.6). **Confirmed effective
   this session: Oliver ran `run_footprint_validation.R` and it matched
   section 12's Python prediction (section 13.1).**

### 12.2 Finding: emission footprint width, not error size, predicts
    correctability

Oliver's hypothesis (longer emission wavelength → more signature
uncertainty generally, except on red/APD detectors; `Spark UV 387`
specifically restricted to ~3 UV detectors, overlapping `BUV395`, leaving
"minimal extra information") was tested directly against the real Aurora
panel rather than argued from first principles. For each fluorophore,
footprint was counted as the number of detectors where
`Cells_spectra_corrections.csv`'s row exceeds 2% of its own peak:

| fluorophore | footprint (detectors) | gain plateau (n→∞) |
|---|---|---|
| Spark UV 387 | 7 | 0.0013 |
| BUV395 | 9 | 0.00009 |
| eFluor 450 | 9 | 0.00001 |
| FITC | 11 | 0.00016 |
| BUV805 | 15 | 0.0017 |
| Spark Violet 538 | 16 | 0.0020 |

("gain plateau" = the held-out `gain` statistic swept from `n.per.dye =
4,000` to `40,000` events per single-stain block; it converges to a stable
non-zero value for some fluorophores and to ~0 for others, distinguishing
a real-but-small effect from no effect at all — see 12.6 for why this
matters more than raw sample size.) Spearman correlation between footprint
width and gain plateau across the full 16-dye panel: **ρ = 0.66**. More
tellingly, the four narrowest-footprint dyes in the table above are
*exactly* the four with near-zero plateau, regardless of population size —
not a coincidence at this panel's specific angle sizes, since `Spark UV
387` (1.31° cross-substrate error) plateaus at a real, detectable 0.0013
while `Spark Violet 538` (0.96°, a *smaller* angle) plateaus at a larger
0.0020. Error size alone does not order these; footprint does, better.

Mechanism: the held-out step search's objective (`gain`, and the `fit.slope`
regression feeding it) sums over every detector in the panel. A real shape
error confined to a handful of UV channels moves that sum by only as much
as those few channels change, while the other ~55-57 channels — which the
narrow-footprint dye never meaningfully emits into — contribute pure noise
to both the numerator and denominator, diluting the signal-to-noise ratio
of the objective specifically for narrow-emission dyes. This is a property
of the objective's denominator, not a `Spark UV 387`-vs-`BUV395`
collinearity problem specifically — the fix that follows doesn't reference
`BUV395` at all and still resolves `Spark UV 387` partially (at
`footprint.frac = 0.05`; not at `0.02`, see below).

Also checked and not the explanation on its own: `nuisance.frac`-based
co-activity flagging. None of the four near-zero-plateau dyes were being
assigned nuisance partners in the baseline run — the dilution happens in
the objective's own detector sum, independent of the active-set
composition.

### 12.3 Finding: `spillover.spread` availability materially changes
    Beads' nuisance-set quality and detectability

`bead_spillover_spread.csv` did not exist in the project when this session
started; only `cell_spillover_spread.csv` did. Before Oliver supplied a
real one, an empirical estimator was built (`est_spillover_spread.py`,
kept for future use — e.g. spot-checking a real one, or a substrate that
genuinely lacks one) from the simulated single-stain pool itself, and
tested for effect: on `PerCP-eFluor 710` (Beads), moving from flat
thresholds (`spillover.spread = NULL`, `n.nuisance = 2`) to a
spread-scaled threshold (`n.nuisance = 0`) moved `gain` from `0` to
`0.048` and flipped the fluorophore from rejected to accepted. This
reproduces and quantifies, on a concrete example, what
`CONTEXT_signature_correction2.md` §2.2 already argued in general ("the
spread-scaled mask... buys real recoverability"). Once Oliver supplied the
real `bead_spillover_spread.csv`, all further runs this session used it in
place of the estimator.

Combined with 12.2's `footprint.frac = 0.02`, both directions of the
Aurora benchmark, idealised noise, `n.per.dye = 4,000`, seed 42:

| | accepted, baseline (flat thresholds, no footprint restriction) | accepted, both changes |
|---|---|---|
| Cells (starting spectra = Beads) | 6/16 | 11/16 |
| Beads (starting spectra = Cells) | 6/16 | 12/16 |

No fluorophore moved away from ground truth in either direction at
`footprint.frac = 0.02` (checked directly: minimum `recovered` value
across both tables is `0.0`, never negative). Both `PerCP` and `NovaFluor
Blue 610-30S` — flagged by Oliver at the start of this session as
plausibly uncorrectable — picked up partial, positive correction as a side
effect of these changes on Cells (`recovered = 0.23`, `0.70`). Treat that
as a bonus, not a guarantee: the simulator draws every dye across the same
bright titration range, so it does not reproduce these two dyes' real
physical dimness, which the real benchmark's `explained` gate is
specifically there to catch. **This idealised-simulation prediction is now
confirmed, in direction and rough magnitude, on the real production
function against real gated data — see section 13.1.**

**Item 4 from the previous message (the reproducible wrong-direction
correction on Beads, section 10.4/section 5 item 1): tested and not
reproduced.** An idealised single-spectrum-per-dye, IID Gaussian/Poisson/
detector-noise simulation could not reproduce the `BUV805`/`Spark Violet
538` flip for either dye, at population sizes up to 10x this benchmark's
own. A bidirectional-gain safety check was also prototyped (evaluate the
held-out objective at both `+t` and mirrored `-t` of the fitted slope;
reject if both directions show comparable gain, on the reasoning that a
genuine one-signed spectral error should not benefit equally from being
reversed) and found to never trigger a false rejection on any currently-
working correction in this simulation — cheap and safe to add, but its
ability to catch the real bug is unproven, since the bug itself could not
be reproduced to test detection against. Per Oliver's steer, not pursued
further this session; see section 5, item 1 and section 7, item 4 for
status.

### 12.4 R diffs written, not yet applied

Five anchored edits to `correct_unmixing_signatures.R`, given to Oliver as
before/after blocks in this session's chat (not applied to the source by
Claude, per this project's standing convention):

1. Roxygen documentation for two new parameters, `footprint.frac` (default
   `0.02`) and `footprint.min.channels` (default `3L`).
2. The two parameters added to the function signature.
3. A `verbose`-gated `message()` when `spillover.spread` is `NULL`, noting
   the cost demonstrated in 12.3 — bringing this function in line with
   `fix.my.unmix()` and `correct.spectra.glasso()`, which already warn on
   this (see 12.5).
4. Computation of `footprint.mask` per fluorophore per iteration
   (`spectra.new[j,] > footprint.frac * max(spectra.new[j,])`, falling
   back to every detector if fewer than `footprint.min.channels` survive).
5. `fit.slope()`'s regression and `residual.gain()`'s objective both
   restricted to `footprint.mask`. Deliberately *not* restricted:
   `explained`, `x.span`, `bg.align`, `span.drift` — these gates exist
   specifically to catch background/scale confounds across the whole
   panel, which narrowing would undermine.

**Status: applied to production and confirmed on real data this session
(section 13.1).** `run_footprint_validation.R` (built two sessions ago)
has now been run against the real production function and the real
benchmark, confirming the direction and rough magnitude of section 12.3's
Python prediction (Cells 10/16 → 11/16, Beads 8/16 → 10/16). See section
13.1 for the full numbers and section 8 for the confirmed current-defaults
reading of `correct_unmixing_signatures.R`.

### 12.5 Applicability to `fix.my.unmix()` and `correct.spectra.glasso()`

Not tested this session — both points below are checked against the
current source, not against a run — but both transfer more directly than
a generic "might also apply" caveat would suggest:

- **`spillover.spread`: already guarded, more thoroughly than
  `correct.unmixing.signatures()` was before this session.** Both
  `fix_my_unmix.R` and `correct_spectra_glasso.R` take
  `variants$spillover.spread` and already emit `warning("No
  'spillover.spread' in 'variants'; positivity ...`" when it's `NULL`
  (`fix_my_unmix.R` line 385-386, `correct_spectra_glasso.R` line
  387-388). This session's finding (12.3) is corroborating evidence that
  an earlier decision to warn on this was well-founded, not a new problem
  for those two functions. `correct.unmixing.signatures()` was the
  outlier in not warning; diff 3 above brings it into line.
- **`extract_raw_signature.R`'s `resid.rel` gate has the same structural
  shape as the `gain` objective 12.2 found the problem in — untested, but
  a strong candidate.** `resid.rel <- sqrt(sum((y.res -
  fitted)^2))/max(sqrt(sum(y.res^2)), eps)` (lines 278-279) is a full-panel
  residual norm over a full-panel signal norm, with no restriction to the
  target fluorophore's own emission footprint — the same shape as the
  pre-session `gain` computation this session restricted. This function is
  phase two for *both* `fix.my.unmix()` and `correct.spectra.glasso()`
  (section 9.1, section 11.3), so a real fix here would help both
  estimators from one change, not one. It is also directly relevant to
  section 11.3's still-open question: `BUV805` and `BV480` are rejected
  there at `resid.rel = 0.0494`/`0.0558` against `max.resid = 0.03`, and a
  footprint-dilution effect would be a third candidate mechanism alongside
  11.3's two (bin noise; real unmodelled coupling from the sparser active
  set). Not checked this session whether either rejected row is in fact
  narrow-footprint on the Discover panel, or whether restricting
  `resid.rel` to each target's own footprint would move either past
  `0.03` — see section 7, item 2 for the suggested check. `explained` and
  `explained.total` in the same file already use a `top.bin`/`top.norm`
  full-panel calculation deliberately, the same choice this session made
  for `correct.unmixing.signatures()`'s own `explained` — that gate is not
  a candidate for the same fix.
- **The footprint-vs-error-size finding itself (12.2) is not specific to
  the dominance-fit algorithm's mechanics** — it is a general
  consequence of computing any full-panel residual metric for a
  narrow-emission dye on a wide detector array, and should be expected
  wherever such a metric gates acceptance, which both other estimators'
  phase two does via `extract_raw_signature.R`.

### 12.6 Process note: scope and limits of this session's Python
    reproduction

Covered: `correct.unmixing.signatures()`'s dominance-fit algorithm in
full, ported line-for-line. Not covered, and worth flagging before trusting
any specific number from this session too far:

- **Noise model is IID and single-spectrum-per-dye.** No spectral variant/
  tandem heterogeneity, no scatter-correlated AF, a single `readout_sd`/
  `counts_per_unit` applied uniformly across all 64 detectors rather than
  the real per-detector values `estimate_noise_model.R` would fit (Oliver's
  own note this session: red/APD detectors show no measured noise
  increase despite longer wavelength, plausibly due to detector
  physics — not represented in this session's flat noise preset). This is
  the most likely reason the wrong-direction Beads bug (12.3) could not be
  reproduced: this project's own prior finding (section 3.1, via a
  different investigation) is that real residuals here are structured
  ("structured variant scatter"), not IID — exactly the property an
  idealised simulation cannot produce by construction. Any number from
  this session's simulation should be read as "what an idealised version
  of this mechanism does," not as a prediction of the real benchmark's
  exact effect size.
- **`bg.mode = "scatter.knn"` not implemented** — both directions used
  `"global.mean"` this session, even though the real benchmark uses
  `"scatter.knn"` for Cells (confirmed from
  `run_signature_correction_diagnostics_aurora.R`'s own Stage D call).
  Defensible for this session's specific questions (both 12.2's and
  12.3's mechanisms are demonstrated independent of background-subtraction
  mode), but not exercised. **`bg_mode = "scatter_knn"` is now implemented
  in the Python reproduction (`signature_knn_subtract()`) and was used
  throughout section 13's Cells testing.**
- **Every dye simulated across the same bright titration range**,
  regardless of real physical brightness — noted above (12.3) as the
  likely reason `PerCP`/`NovaFluor Blue 610-30S` look more correctable in
  this session's simulation than they likely are in reality.
- **`extract_raw_signature.R` and `cluster_unmixed_events.R` not ported** —
  turned out not to be needed for this function specifically (12.1, point
  2), but this means 12.5's `resid.rel` hypothesis is exactly that, a
  hypothesis, not something this session's environment could test.

If this is picked back up: porting `extract_raw_signature.R`'s `resid.rel`
calculation (or at minimum a standalone repro of it) before touching R
again would let 12.5's hypothesis be tested the same way 12.2's was, and
would let both open questions in 12.5 be answered from the same
environment rather than by inference.

### 12.7 Where this leaves things

**Resolved, validated in Python, and now confirmed in R against real
data**: the mechanism behind several of section 10's never-correcting
fluorophores (footprint dilution, 12.2) and the cost of a missing Beads
`spillover.spread` (12.3). Diffs written (12.4) and confirmed via
`run_footprint_validation.R` (section 13.1).

**Ruled out, not resolved**: one more candidate mechanism (generic
finite-sample estimator bias, independent of structured noise) for the
Beads wrong-direction bug, which remains open but deprioritised per
Oliver's steer (section 5, item 1).

**Tested and answered this session (13)**: `extract_raw_signature.R`'s
`resid.rel` hypothesis from this section remains untested (still open,
section 7 item 2), but a different, related question this session opened
— whether AF-removal mechanisms belong inside
`correct.unmixing.signatures()`'s own pipeline — is now answered for two
specific mechanisms (per-cell AF extraction, AF deconvolution proper): no
for both, in that function's own background-removal role.

**Function status unchanged pending further R confirmation of the
remaining open items**: `correct.unmixing.signatures()` remains the
production candidate; the `footprint.frac`/`spillover.spread` changes are
now confirmed and applied; `max.step` and `n.split.trials`/`min.split.frac`
remain recommended but not yet applied (section 7, item 1).

---

## 13. This session: `run_footprint_validation.R` run for the first time
    (confirms section 12); true AF deconvolution (`get.af.basis()` /
    `deconvolve.af.background()`) tested for the first time as a candidate
    `bg.mode` for `correct.unmixing.signatures()`

Follow-on to section 12. Two things happened this session, both against a
new, better benchmark input: Oliver produced FSC/SSC main-population-gated
versions of the concatenated single-stained control files
(`small_Concatenated_cells.fcs`, 622,522 events; `small_Concatenated_BD_beads.fcs`,
81,260 events), replacing the ungated files section 12's own Python
reproduction had been using via a uniform downsample. All numbers in this
section use the gated files throughout, both in R and in this session's
Python reproduction.

### 13.1 `run_footprint_validation.R` run for the first time — section 12's
    recommendation confirmed on real data, with `spillover.spread`'s own
    contribution isolated and independently reconfirmed

Section 12.4 left two things un-run: `footprint.frac`/`footprint.min.channels`
against the real production function on the real benchmark, and (from an
intervening session not otherwise recorded in this document) a proposed
further change decoupling `spillover.spread` from nuisance-set gating,
tested only in an ungated, heavily downsampled Python reproduction and
provisionally called a second bug ("Bug 2") in that session's own report.

Oliver ran `run_footprint_validation.R` himself, with only `footprint.frac`
applied (not the proposed `spillover.spread` decoupling), against the new
gated `small_` files:

| substrate | accepted, old defaults → new (`footprint.frac = 0.02` + `spillover.spread` supplied) | notable moves |
|---|---|---|
| Cells | 10/16 → 11/16 | APC 0.498 → 0.920, PE 0.076 → 0.717, FITC newly accepted at 0.322 |
| Beads | 8/16 → 10/16 | RB744 newly accepted at 0.679, Super Bright 600 newly accepted at 0.250, APC 0.561 → 0.741, PerCP-eFluor 710 0.015 → 0.483; PE the one mild regression, 0.809 → 0.742 |

**This confirms section 12's recommendation directly on the real
production function**, on real gated data, for the first time. It also
settles the intervening "Bug 2" proposal: the "new" column already
includes `spillover.spread`'s original, unmodified nuisance-gating
behaviour (never decoupled in this run), and it is a clear net positive —
several large genuine gains against one small give-back (`PE` on Beads).
The Python reproduction that motivated the "Bug 2" proposal was built
against an ungated, more heavily downsampled dataset; re-run against the
same gated `small_` files this session (Python port, `spillover.spread`'s
original behaviour restored, `bg_mode = "scatter_knn"` for Cells /
`"global_mean"` for Beads), it now agrees with the real R run's overall
direction and, on Cells specifically, with its exact accepted count
(11/16) once a separate port-fidelity gap was fixed (13.2). **The
`spillover.spread`-decoupling change is retracted. It should not be
applied.**

### 13.2 Port-fidelity fix: `corr_sig_core.py`'s `fit_slope()` had drifted
    out of sync with production `correct_unmixing_signatures.R`

Before this session's re-run, the Python reproduction's `fit_slope()`
fitted the per-detector regression on every detector and restricted only
the diagnostic R-squared to a dye's own footprint mask. Current production
`correct_unmixing_signatures.R` (confirmed by direct reading) restricts the
`lm.fit()` call itself to the footprint-masked columns, so channels outside
a dye's own footprint carry an exact zero slope and are never touched by a
correction step:

```r
fit.m <- stats::lm.fit( x = cbind( 1, x.bin ),
                        y = r.bin[ , footprint.mask, drop = FALSE ] )
slope <- rep( 0, ncol( r.bin ) )
slope[ footprint.mask ] <- stats::coef( fit.m )[ 2, ]
```

Fixed in the Python port to match exactly (fit restricted to the masked
columns, slope zero-filled elsewhere). This is a Python-port-only fix —
production R already had this right — but it mattered for trusting this
session's own diagnostic numbers: after the fix, the Python reproduction's
Cells accepted count (11/16) matches the real R run's exactly, resolving a
previously-unexplained single-fluorophore discrepancy (APC) noted in an
earlier session's report. Beads still runs 2 fluorophores short of the
real R count (8/16 vs. 10/16) with every individual fluorophore's
direction and approximate magnitude agreeing; the `gate.main` 2D-density
pre-gate, still not ported, is the leading candidate for that residual gap
and is more likely to matter for Beads' smaller populations than for
Cells'.

### 13.3 A terminology correction, load-bearing for the rest of this
    section: "AF deconvolution" vs. "per-cell AF extraction"

An earlier session (not otherwise recorded in this document) built and
tested `get.af.spectra()`'s SOM-based candidate-AF-library machinery plus
`assign.af.fluorophores()`/`assign.af.residuals()`'s per-cell discrete
assignment, and reported it under the name "AF deconvolution." That is
**per-cell AF extraction**, a different mechanism from AF deconvolution:
each event is assigned to ONE of a small dictionary of candidate AF
spectra and has that one spectrum's fitted abundance subtracted. Tested
against the real gated Cells data, this measurably regresses several
already-well-corrected fluorophores (most notably `Spark Violet 538`,
0.954 → 0.242 recovered), traced to `assign.af.fluorophores()`'s
"minimise apparent fluorophore signal" heuristic mis-assigning a small
minority (4%) of one dye's own strongly-positive population to the wrong
AF candidate — a heuristic built for sparse real-tissue positivity, which
a single-stained control's large single-dye-positive populations violate
by design. Confirmed not specific to that one assignment rule:
`assign.af.residuals()` (which collapses almost the entire file to the
population-mean AF candidate) produces nearly the same regression. Neither
result is retracted; both stand as documentation of that specific,
correctly-named mechanism's behaviour on this benchmark. Renamed in the
delivered code from `bg_mode = "af_deconvolution"` to
`bg_mode = "af_library"` to avoid confusion with 13.4 below.

**AF deconvolution proper** is a mechanistically distinct thing:
`get.af.basis()` (SVD of the raw, uncentred unstained matrix into a small
basis of autofluorescence principal components) plus
`deconvolve.af.background()` (that basis fit *jointly* with the
fluorophore panel by OLS for every event, only the AF-component
coefficients clipped non-negative before subtraction) — both in
`deconvolve_af_background.R`, and the same mechanism `fix.my.unmix()`
already uses by default as `bg.mode = "af.deconv"`. `get_spectral_variants.R`
uses the same SVD idea independently (`af.pcs.list`, one basis per
universal-negative file, `nv = 4` fixed rather than `get.af.basis()`'s own
`n.pc = "auto"`) for a different downstream purpose (per-file AF-PC
directions feeding `af.collinear.threshold`, not background removal).
Neither function has previously been discussed in this document, and
neither is currently an option in `correct.unmixing.signatures()`
(`bg.mode` there is `"global.mean"` / `"scatter.knn"` / `"none"`
only) — this session is the first to test it in that role.

### 13.4 AF deconvolution proper, tested for the first time: net negative
    as a first-phase background-removal step, inert-but-safe as a
    second-phase step

Ported `get.af.basis()` and `deconvolve.af.background()` to Python exactly
(component-count selection including the permutation-based `n.pc = "auto"`
heuristic, sign convention, hotspot-based component freezing), wired in as
`bg_mode = "af_deconv"`, and tested two ways against the real gated
Cells/Beads data, on top of the same `footprint.frac = 0.02` +
`spillover.spread` configuration as 13.1.

**As a replacement for the existing background-removal step** (the same
role `global_mean`/`scatter_knn` play now, i.e. before any spectral
correction runs): net negative on Cells (8 of 16 fluorophores move
meaningfully away from the baseline, including several — `FITC`,
`Spark Violet 538`, `APC`, `Spark Blue 550`, `PE`, `Super Bright 600` —
that were already well-corrected; `eFluor 450` flips to an accepted,
wrong-direction correction; `PerCP-eFluor 710` loses its one accepted
correction entirely), roughly neutral on Beads (one mild regression,
several mild improvements, no large moves).

**As a second-phase step layered on top of already phase-1-corrected
spectra** — the role `af.deconv` actually plays in `fix.my.unmix()`,
which only runs after a first pass has already brought the panel spectra
close to correct: almost entirely inert. 10 of 11 previously-accepted
Cells fluorophores see no further change (correctly rejected, no held-out
gain left to find); all 8 accepted Beads fluorophores likewise unchanged.
One exception, reproducible across every configuration tested
(`n_pc = "auto"`, `n_pc = 4`, before or after phase 1): `eFluor 450`
consistently picks up a small, real, wrong-direction accepted correction
(0.62° → 0.9-1.0° from ground truth) that phase 1 never touches either way
(it never accepts a correction for this dye under any tested
configuration). Not explained by AF-basis/fluorophore hotspot collinearity
— the actual coupling, computed the same way `deconvolve.af.background()`
does, tops out at 3.7 (Cells) / 3.3 (Beads) against a `max.hotspot = 5`
default, for every configuration including `eFluor 450` specifically.

**Mechanism, common to both AF-removal approaches tested in this document
(this session's and the earlier session's per-cell extraction, 13.3),
despite reaching it by different routes**: both compute a per-event
background estimate that is itself a function of the very (deliberately
wrong, being corrected) panel spectra. `deconvolve.af.background()`'s
estimate is `af_coef @ af_basis`, where `af_coef` comes from a single
joint OLS solve against `[af_basis; spectra]` — when a starting row is
badly wrong (the deliberate setup of this whole cross-substrate
validation), whatever error that row can't explain has somewhere to go in
that joint solve, and some of it lands in the AF-component coefficients,
contaminating the background estimate with abundance-correlated content
before the correction loop ever runs. Confirmed this isn't simply "losing
scatter-based per-event specificity": a `global_mean` baseline on Cells
(also spectra-light, no scatter matching either) holds up about as well as
`scatter_knn` does, while `af_deconv` still regresses badly under the
identical comparison — `spectra`'s only role in `global_mean`/`scatter_knn`
is deciding which events count as background, a far weaker dependency than
being solved against directly. This is consistent with, and explains,
`fix.my.unmix()`'s own choice to use `af.deconv` only as a *second*-phase
mechanism, after a first pass has already gotten the panel spectra close:
testing it as a *first* background-removal step, ahead of any spectral
correction, inverts that dependency and runs into exactly the fragility
the two-phase split already guards against.

**A separate, secondary finding, worth recording independent of the
above**: `get.af.basis()`'s `n.pc = "auto"` selection retains only 1
component (the mean-background direction) on both Cells and Beads at this
benchmark's real values — not for lack of real multi-component structure
(components 2-4's singular values run 600,000-1,100,000 for Cells, well
above a plain electronic-read-noise floor of ~29,000), but because the
permutation-derived noise threshold itself comes out enormous (~2.6
million for Cells) on raw, untransformed spectral flow data: raw detector
values span orders of magnitude with a heavy right skew, and permuting
each detector column independently preserves that skew, which alone is
enough to give the permuted matrix a large second singular value with no
real cross-detector correlation involved. This looks like a real caveat of
the `n.pc = "auto"` heuristic specifically on raw (rather than
transformed) detector data. Forcing `n.pc = 4` (matching
`get_spectral_variants.R`'s own fixed `nv = 4`) does not change the
qualitative result in either test above.

### 13.5 Recommendation

**Do not add `bg.mode = "af.deconv"` to `correct.unmixing.signatures()`'s
own background-removal step.** The evidence here argues against it in
exactly the role it would play there — a first background-removal pass,
ahead of any spectral correction. This is not a criticism of the mechanism
itself: the package's existing two-phase architecture
(`correct.unmixing.signatures()` first, `fix.my.unmix()`'s `af.deconv`
second) already places it where the evidence says it belongs, and nothing
here argues for changing that split. **No R diff is written for this
finding** — the recommendation is to leave the current architecture as it
is, not to change it.

**The one open, real, reproducible finding worth a future session's
attention**: `eFluor 450`'s small but consistent wrong-direction move
under `af_deconv`, present in every configuration tested and not explained
by hotspot collinearity. Low priority relative to everything else in this
document (this dye never corrects under the current production pipeline
either), but worth a note in case `af.deconv` is ever considered for a use
closer to `correct.unmixing.signatures()`'s own pipeline in the future.

### 13.6 Where this leaves things

**Resolved, confirmed on real data**: section 12's `footprint.frac`/
`spillover.spread` recommendation (13.1) — `run_footprint_validation.R`
has now been run, for the first time, against the real production
function on real gated data, and the result matches section 12's Python
prediction. This closes section 7, item 1 of this document's prior
version.

**Retracted**: the intervening "Bug 2" proposal (decoupling
`spillover.spread` from nuisance-set gating) — real gated data shows the
original, unmodified behaviour as a net positive, not a bug. Do not apply.

**New, tested for the first time, not recommended for production**: AF
deconvolution proper (`get.af.basis()`/`deconvolve.af.background()`) as a
`correct.unmixing.signatures()` background-removal option. Net negative
before spectral correction, inert-but-safe (one small exception) after —
consistent with, and no argument against, the package's existing
placement of this mechanism as `fix.my.unmix()`'s second-phase default.

**Terminology now corrected going forward**: "AF deconvolution" refers
only to the `get.af.basis()`/`deconvolve.af.background()` SVD/joint-OLS
mechanism (this section). The SOM-library/discrete-per-cell-assignment
mechanism from `get.af.spectra()`/`assign.af.fluorophores()`/
`assign.af.residuals()`, previously reported under that name in error, is
"per-cell AF extraction" (13.3) and is tracked separately.

---

## 14. Hypernegative correction penalty, cluster-based audit, and gating/downsampling fixes to `fix.my.unmix()`

Consolidates three documents — `CONTEXT_hypernegative_penalty_and_cluster_audit.md`,
`CONTEXT_hypernegative_penalty_implementation.md`, and
`CONTEXT_gating_downsampling_and_af_basis_fixes.md` — covering four
sessions of continuous work, all on a real seven-colour lung/spleen panel
(BUV395, BUV805, BV421, PE, PE-Cy7, APC, eFluor 780 + AF), not the
FACSDiscover or Aurora benchmarks used elsewhere in this document. Every
file-level claim below has been re-verified against the current copies of
`fix_my_unmix.R`, `create_biplot.R`, and `run_hypernegative_cluster_audit.R`
in this project as of this consolidation pass, with any discrepancy
flagged explicitly rather than silently trusted (see 14.6).

### 14.1 The starting problem

`fix.my.unmix()` was observed to over-correct on multi-colour lung
samples, producing hypernegative events — the running test case throughout
was BUV805 vs. PE-Cy7 in lung, where Oliver confirmed the two dyes are
genuinely co-expressed biologically (minimally so in spleen), which
complicated distinguishing a real modelling error from real biology for
most of this thread.

### 14.2 The penalty mechanism: three iterations

1. **First attempt**: `max.hypernegative.delta`, a soft trust-damper
   mirroring the existing positive spread threshold about zero, penalising
   a candidate in proportion to how much it *increased* the hypernegative
   fraction relative to before. Engaged only weakly on the real test case
   (`trust = 0` at the final iteration) and did not resolve the
   over-correction.
2. **Retracted hypothesis**: that the coefficient was accepted once early
   and never revisited. `coefficient.decay` was added on this theory,
   gated on `trust == 0` exactly — barely moved the number even from a
   from-scratch run, ruling this out. Corrected understanding: `trust == 0`
   is a much narrower condition than "weakly supported" (the inverse-
   variance weight is essentially never exactly zero for a pair that
   clears the hard acceptance gates), so `coefficient.decay` was
   reformulated to scale continuously by `(1 - trust)` on every
   off-diagonal cell, every iteration — this version works as designed and
   is not itself in question, though it will not meaningfully erode a
   coefficient the estimator keeps independently re-finding, which is
   correct behaviour, not a bug.
3. **Per-iteration history** (`keep.history`/`spillover.history`) showed
   `trust` consistently 0.16–0.26 across all 15 iterations for
   BUV805→PE-Cy7, both the fresh candidate and the running total converging
   smoothly — this fully retracted the "orphaned coefficient" theory. The
   coefficient reads as genuinely, repeatedly, independently re-identified,
   which reframed the question as whether ~0.03 is real spectral overlap,
   real co-expression, or both, rather than an estimator artefact.

**Final design, per Oliver's explicit mandate**: a negative correlation
between a source's abundance and a target channel's corrected value is
treated as **unambiguous** evidence of an incorrect coefficient — unlike
coverage/disagreement, which this thread repeatedly showed can look
identical to real co-expression, there is no biological mechanism that
produces a genuine over-correction past the negative boundary. It should
therefore be a hard reject at a cap, not a proportional penalty. This
replaced `max.hypernegative.delta` with `min.hypernegative.frac` (default
`0.01`) / `max.hypernegative.frac` (default `0.05`) / a hard-fires-below
`min.hypernegative.events` (default `200L`) gate, scoped to each source's
own positive population rather than all events (the main fix over the
original design — an all-events denominator dilutes a real over-correction
whenever the source's own positive population is a small share of the
total sample), and judged on the *absolute* post-candidate fraction rather
than the delta over baseline, deliberately not rewarding a coefficient for
merely improving on a bad starting point. The penalty is linear between
floor and cap, reaching exactly zero (hard rejection) at the cap. All of
this is confirmed present in current `fix_my_unmix.R`
(`min.hypernegative.frac`/`max.hypernegative.frac`/
`min.hypernegative.events` as parameters; `hypernegative.base`/
`hypernegative.after`/`hypernegative.new` in `coefficient.log`).

### 14.3 The cluster-based audit

`cluster.unmixed.events()` (pre-existing, previously unused in production)
was wired into a new workflow: `fix.my.unmix()` now exposes
`unmixed.final`, `residual.final`, `threshold.matrix.final`,
`neg.threshold.matrix.final`, `thresholds.final`, `neg.thresholds.final`,
and `dominant.final` in its return value — confirmed present by direct
reading — computed once right after the final compensation is fixed, so
they match exactly what phase two itself works from. A new standalone
script, `run_hypernegative_cluster_audit.R`, consumes these: Section 1
clusters each tissue's `unmixed.final`/`residual.final`
(`min.cluster.size = 20`, well above the library's own default of 5, since
here the centroid itself is the evidence rather than a regression-input
stabiliser); Section 2 flags clusters whose channel centroid crosses the
negative boundary only *after* the accepted coefficient is applied, not
before; Section 3 reports, per cluster, how much of a target's dominant
population is also positive for a suspected partner.

Run on lung + spleen, this found one real, small over-correction (lung
cluster 240, 22 of 19,631 events, a genuine PE-Cy7-negative population
pushed from −2342 to −4626 against a boundary of −4127), and a composition
asymmetry that ran opposite to naive expectation — spleen's much smaller
BUV805-dominant population is proportionally *more* PE-Cy7-co-positive
than lung's, despite spleen having less overall tissue co-expression. This
is treated as an open biological question about what defines that
population in each tissue, not evidence of a modelling artefact.

A separate, real bug was found and fixed in `create.biplot()` along the
way: its first-draft spillover-spread reference curves were flat
`geom_hline`/`geom_vline` lines, independent of the other axis's signal —
this did not reflect how `get.spread.thresholds()` actually behaves, and
the bug was confirmed scoped entirely to the new plotting code, not the
core estimator. Corrected to `geom_line` curves (intercept the flat
threshold, slope `spread.kappa · √(spillover.spread)`, computed in raw
space and biexponentially transformed for display). Using the corrected
curves, Oliver visually confirmed BUV805+ events in lung crossing the
negative boundary under the corrected spectra — direct visual confirmation
of the same problem the cluster audit had already quantified.

**Known, unresolved limitation**: cluster assignment is sensitive to
`som.dim` — 20 vs. 30 gave different cluster boundaries for what appears
to be the same lung population (cluster 240 splits into three smaller
clusters, 1/61/31, at `som.dim = 30`). Not investigated further.

### 14.4 Root cause of the motivating case: gating and downsampling, not the estimator

Implementing the new hard-reject penalty (14.2) did not visibly fix the
BUV805/PE-Cy7 lung case on its own, which motivated tracing why. Two
independent, much larger bugs were found upstream of the estimator
entirely:

**Scatter gating.** Both the unstained and fully-stained samples are gated
through the same polygon, fit once from the unstained file and reused
verbatim on the stained one. `large.gate = TRUE` widens this by
extrapolating the strict, density-defined population's convex hull outward
past a quantile — an extrapolation from the main population's shape, not a
guarantee of reaching a genuinely separate large-cell population.
BUV805+/PE-Cy7+ cells in this panel are large, high-SSC events, and were
excluded from *both* files by this same mechanism — from the stained fit
(so the pair-fit never saw them) and from the unstained fit (so the AF
basis never saw their autofluorescence either). `scatter.gate` (master
on/off, default `TRUE`) and `landmark.quantile` (a lighter-weight
scatter-only rescue threshold, computed once from the unstained file and
reused unchanged on the stained file so both stay matched, default `NULL`)
were added, threaded through the shared `read.gated()` closure. Confirmed
present in current `fix_my_unmix.R`.

**Downsampling.** Even with `scatter.gate = FALSE`, the correction still
didn't reproduce a direct OLS unmix of the full gated file — the
brightest BUV805 events and brighter PE-Cy7 events were simply absent from
`unmixed.final`. Root cause: `unmixed.final` is computed from
`stained.fit$abundance`, which by default has already been overwritten by
a stratified downsample (`downsample = 20000`) before the phase-one loop
ever runs. No deliberate hypernegative-based exclusion exists anywhere in
`.fix.stratified.sample()` — every event is either in the background
bucket or in some fluorophore's positive stratum, and sampling within a
stratum is uniform, not brightness-based — but the *quota* math is severe
on a panel with many fluorophores: `floor.n = n.fluorophores *
downsample.min.stratum` easily exceeds the positive budget, and every
stratum's retained count gets scaled down by the same ratio, which alone
is enough to visibly thin a rare, heavy right tail at a low retention rate
with no boundary logic involved at all. `stained.fit.full`, a copy of the
pre-downsample projection, was added specifically so `unmixed.final`/
`residual.final` reflect every gated event regardless of what the
(necessarily downsampled, for speed) iterative phase-one loop fits on.
Confirmed present in current `fix_my_unmix.R`.

With **both** `scatter.gate = FALSE` and `downsample = FALSE` together (the
strongest, simplest test), the BUV805/PE-Cy7 over-correction was
"perfectly solved," confirmed by Oliver's own direct comparison against an
independent OLS unmix. This specific test did not exercise the
`stained.fit.full` decoupling in isolation, since with `downsample = FALSE`
the stratified-sample block never runs and `stained.fit.full` is identical
to `stained.fit` anyway — that diff's value is specifically for keeping
downsampling on for speed while still getting full-data fidelity in the
reported matrices, a usage mode not yet independently tested (section 5,
item 17).

**The one remaining discrepancy**, found while manually reproducing
`unmixed.final` via `unmix.ols()` plus a fresh `get.af.basis()` call, is
attributed to AF, not to the spillover correction: a small eosinophil
population reads hypernegative in PE-Cy7/PE against BUV805 under the
shared, low-rank `af.deconv` basis, but is positioned correctly under
per-cell AF assignment — confirming the corrected spectra themselves are
right, and that this is an AF *representation* gap (the shared basis
under-representing a spectrally distinct population's own background), not
a spillover-correction error (section 5, item 19).

### 14.5 Key learnings

- A negative correlation between a source's abundance and a channel's
  corrected value is unambiguous evidence of an incorrect coefficient, so
  the penalty judges the absolute resulting fraction, not the change from
  baseline, and does not reward a coefficient for merely improving on a
  bad starting point — judging by direction-of-change rather than
  absolute end-state was exactly the loophole the penalty was built to
  close.
- A downsample's stratification floor is a per-stratum *count*, not a
  fraction of that stratum's true size; on a panel with many fluorophores
  the shared budget forces every stratum's actual retention far below the
  nominal floor once quotas are scaled down to fit, and this alone can
  produce a visibly biased-looking result with no boundary-based logic
  involved anywhere (section 5, item 20).
- A single scatter gate polygon, fit once from an unstained control and
  reused on the stained sample, keeps the two consistent with each other —
  but a population genuinely outside that gate is excluded from *both*
  simultaneously, corrupting the AF reference the same way it corrupts the
  correction's own input data. Fixing one without the other reintroduces
  the mismatch in a new direction.
- `get.af.basis()`'s `n.pc = "auto"` component-count selection is not
  reproducible across separate calls on the same data (its permutation-based
  noise-edge estimate and its `max.events` subsampling are both unseeded);
  always reuse a returned `af.basis` rather than recomputing one when exact
  reproduction matters.
- `unmixed.final` = joint OLS abundance (panel + AF-basis design) times
  `compensation`, not the joint OLS fit alone — the phase-one residual
  spillover correction is a separate multiplication applied after
  `project()`, easy to omit when manually reproducing the pipeline.
- `unmix.ols()` reorders/validates columns by name against its `spectra`
  argument and warns if it has to; `unmix.ols.fast()` (what
  `fix.my.unmix()` calls internally) assumes positional alignment silently
  — a real silent-failure mode worth remembering for other manual
  comparisons, though not the cause of anything in this thread.

Open items from this thread are tracked in section 5 (items 17, 19, 20,
21) and section 7 (items 16–18, 20).

### 14.6 Discrepancy: the audit script does not match its own documentation

`CONTEXT_gating_downsampling_and_af_basis_fixes.md` states that
`run_hypernegative_cluster_audit.R`'s Section 2 was changed to use
`fit$neg.threshold.matrix.final` directly, replacing the old
`-fit$threshold.matrix.final` mirrored boundary. **Direct reading of the
current script this pass shows the old line still in place**:

```r
neg.threshold <- -fit$threshold.matrix.final
```

with no reference to `neg.threshold.matrix.final` anywhere in the file.
`create_biplot.R`, by contrast, does correctly use the real,
directly-measured `variants$neg.thresholds` rather than a mirrored value —
so the fix was applied at the panel-level/`create.biplot()` layer but not
carried into this specific audit script, or the script in this project
predates that edit. Either way, any cluster-audit re-run described in the
gating-fixes document's "open questions" (re-run on lung/spleen with the
new hard-reject penalty live, checking cluster 240) used the *old* mirrored
boundary, not the corrected one, and should be treated as needing a re-run
rather than as confirming anything about the new penalty's interaction
with the real boundary. See below for the diff to bring the script back in
line with its own documentation.

---

## 15. Pair-restricted residual validation of `fix.my.unmix()` phase two

New, unresolved analytical thread from `CONTEXT_pair_restricted_residual_validation.md`,
built on top of the gating/downsampling fixes in section 14 (the test
substrate was run with `scatter.gate = FALSE, downsample = FALSE`, section
14.4's confirmed-working configuration). **Not yet confirmed against the
right substrate — see 15.5 — and the script's own verdict logic does not
currently match what this section's reported results describe, per the
note below.** Treat everything in 15.3 as suggestive, not confirmed, until
both are resolved.

### 15.1 The question

`fix.my.unmix()`'s two phases correct in-span spillover error and
row-shape error respectively; `correct.unmixing.signatures()` corrects
out-of-span row-shape error from detector-space residuals under a design
restricted to the dominant fluorophore plus co-active nuisance dyes. This
thread asked whether a similar detector-space residual check could confirm
or refute what `fix.my.unmix()`'s own phase-two signature step does, given
that phase's own acceptance gates are otherwise the only check on whether
a candidate correction is real.

The full-panel OLS residual cannot do this, for an exact reason. Let `S`
be the current spectra, `W = (SS')⁻¹S` the unmixing matrix, `P = S'W` the
projector onto `row(S)`, and `r = (I-P)y` the OLS residual. If a row's
error lies entirely inside `row(S*)`, then `row(S) = row(S*)`, `P` is
unchanged, and `r` is identically unchanged — not approximately, exactly.
Since in-span error is precisely what `fix.my.unmix()` exists to correct,
a correctly-operating phase-two step is invisible to the full-design
residual by construction; under isotropic noise, regressing the full
residual on any abundance is guaranteed to return zero regardless of how
wrong `S` is. In lay terms: if you ask a model "how much does the data
disagree with everything I already believe," and the model already
believes the (wrong) correction, the disagreement can't see the thing it
already absorbed — you have to take the belief itself back out before
asking the question.

Leave-one-out restriction is insufficient for the same reason
`correct.unmixing.signatures()`'s `nuisance.frac` deliberately keeps
co-active partners in its own design: if row `j` is contaminated with
`s_j = s_j* + ε·s_k*`, removing only `j` from the design leaves `s_k*` in
it, and the contamination term vanishes because `s_k*` is still spanned.
**Both members of a suspect pair must be removed together** for the error
to become visible; this is the design `run_pair_restricted_residual_validation.R`
implements.

### 15.2 What the script does

For each of the top `prv.n.pairs` most collinear pairs (ranked by
starting-spectra cosine similarity): builds a restricted design (panel
minus the pair, plus the AF background basis, always retained so
background structure is projected out rather than mistaken for pair
signal); pools events single-positive for one member, single-positive for
the other, and double-negative for both (capped at 2x the combined
single-positive count); regresses the restricted residual jointly on both
abundances with an intercept and compares the fitted response per row
against what that row's candidate spectrum predicts, both before and after
`fix.my.unmix()`'s correction; reports a `separability` statistic (cosine
between the two rows' model-predicted responses) and a `regressor.corr`
statistic (correlation between the two abundances over the scored
events), scoring a pair as `identifiable` only if both clear threshold;
puts an error bar on each angle via repeated 50/50 split-half
re-estimation, calling a change "improved" or "worsened" only if it exceeds
that row's own split-half SD; and can build the rest of the restricted
design from either the starting spectra or the corrected spectra, which
isolates whether an effect is intrinsic to the pair under test or
contingent on what the rest of the panel's own corrections did.

### 15.3 Results on the seven-colour lung sample

21 pairs scored, 0 unidentifiable. Panel-level weighted score decreased
overall (a net improvement). Per-pair: 8 improved, 4 worsened, 1 mixed, 8
unchanged using the corrected-spectra design; 9/4/1/7 using the
starting-spectra design — the two designs agreeing on nearly every verdict
is itself informative.

- **BUV395 and BUV805** both showed a consistent, several-SD improvement
  against 5 of their 6 pairing partners. Tested against each other
  directly, the result came back mixed in both designs: BUV805 still
  improves cleanly, but BUV395 shows a small, reproducible worsening. This
  is the single largest pairwise cosine in the panel and consistent with a
  genuine, small mutual coupling between the two most collinear rows
  rather than a contradiction of the 5-partner result.
- **APC**'s candidate was rejected on `offset` with a strong background-
  confound signature (`bg.align = -0.656`) and the panel's worst VIF; the
  pair-restricted test agreed qualitatively, but this is exactly the case
  flagged in 15.5 as needing the multi-colour caveat, since APC's own
  reference population is not guaranteed pure single-positive.
- **BV421 and PE-Cy7** were rejected at phase two with the smallest event
  counts in the panel and the largest proposed rotations; **PE** was
  rejected specifically because its candidate would have measurably
  increased leakage into other channels, a cleaner and more specific
  mechanism than the other three.
- **eFluor 780** is the standout finding and the one most likely to
  generalise: best-behaved row in the signature log by every
  self-referential metric, yet a consistent worsening against 4 of its 5
  informative partners, reproducing almost identically whether the
  restricted design uses starting or corrected spectra for the other rows
  — which rules out "leakage from still-uncorrected partners" as the
  explanation. The two measurements are asking different questions:
  `fix.my.unmix()`'s own gates check whether a candidate fits its own
  population and doesn't increase leakage into other channels; they do not
  check whether removing the row from the design changes what other
  channels' bright populations say about it.

### 15.4 Full result tables

`pair_restricted_residual_after.csv`/`_before.csv` (per-pair) and
`pair_restricted_summary_after.csv`/`_before.csv` (panel-level).

### 15.5 Two caveats, flagged by Oliver, that limit how far to trust 15.3

**This is a real multi-colour sample, not concatenated single-stained
controls.** The single-positive/double-negative populations are defined by
threshold-crossing in the full-panel unmixing, not by acquisition-time
single staining. An event crossing one threshold while reading negative
for its partner is not guaranteed free of real, independent biological
signal — a confound the test cannot distinguish from a true spectral
error, and one that concatenated single-stained controls would remove by
construction. Oliver's assessment is that this thread's interpretation,
particularly the `eFluor 780` and `APC` findings, reads more into the
results than this substrate can currently support, and that a rerun on
concatenated single-stains is the correct next step before treating
section 15.3 as more than suggestive. That rerun has not been done
(section 5, item 22 / section 7, item 19).

**The null-fit approach has been separately and substantially rejected as
not useful, in a different session not summarised in this document.** This
directly affects the calibration plan floated for the `eFluor 780` finding
(running this script against a concatenated single-stain "real null" to
establish its own noise floor), since that plan is a variant of the same
null-fit logic. Until a reconciliation happens, that calibration step
should be treated as open to revision, and `fix.my.unmix()`'s own
`null.fit` parameter should not be assumed to be the right vehicle for it.

### 15.6 Implementation options, not decided

**(a) Check-and-warn**, a standalone function reporting pair-restricted
angles as a diagnostic table without altering acceptance behaviour —
reversible, lower risk, the natural near-term choice given 15.5.

**(b) A new gate inside phase two's acceptance stack** — real
per-iteration cost, and risks rejecting genuinely good corrections
(`eFluor 780`'s candidate is the clearest example: best-in-panel on every
existing gate, flagged only by this one) on a metric whose false-positive
rate on real multi-colour data is exactly what 15.5 leaves open. **Not
recommended until (a) has been run on concatenated single-stains and the
null-fit question is resolved.**

### 15.7 Discrepancy: the script's verdict logic does not match its own documentation

`CONTEXT_pair_restricted_residual_validation.md` states that a verdict-
logic bug — calling a pair "improved" whenever *either* member cleared the
improvement threshold, even when the *other* member simultaneously
worsened in the same row — was found and fixed mid-session, with the fix
computing `status.j`/`status.k` independently and reporting `"mixed"` when
they disagree. **Direct reading of the current script this pass shows the
pre-fix logic still in place**:

```r
prv.result$verdict <- with( prv.result, ifelse(
  !identifiable, "unidentifiable", ifelse(
    is.finite( delta.j ) & delta.j < -deg.j.sd |
      is.finite( delta.k ) & delta.k < -deg.k.sd, "improved", ifelse(
        is.finite( delta.j ) & delta.j > deg.j.sd |
          is.finite( delta.k ) & delta.k > deg.k.sd, "worsened",
        "unchanged" ) ) ) )
```

This chain has no path to a `"mixed"` result at all — the OR in the first
branch means any pair where `delta.j` improves and `delta.k` worsens (or
vice versa) is simply reported as `"improved"`, exactly the bug the
document describes as fixed. The BUV395/BUV805 "mixed" result reported in
15.3 is therefore **not reproducible from the current script** — either
the run that produced it used a different, since-reverted copy, or the
"mixed" label in the saved CSVs was assigned by some other means not
visible in this file. Either way, 15.3's specific per-pair verdicts should
be treated as unconfirmed until the script is fixed and re-run. See below
for the diff.