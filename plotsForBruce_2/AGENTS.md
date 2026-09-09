# AGENTS.md — multi-agent working-file map (`plotsForBruce_2`)

Several agents are working in this repo at once. This file records **who owns which
files** so we don't overwrite each other. If you are an agent working here, read this
before editing shared files, and **append/update your own section**.

Branch: `dupi-gamma-analysis`. Paths below are relative to `plotsForBruce_2/` unless noted.

---

## Quick map

| Area | Owner | Report / main output |
|---|---|---|
| Main integrated report + gamma pipeline | **main-report agent** | `report.Rmd` → `out/report.html` |
| Audiobook(=baseline)-vs-focusedBreathing(=ATB) comparison | **taskcmp agent** | `taskcmp_report.Rmd` → `taskcmp_report.html` (**project root**, not `out/`) |
| Olfactory HRV ↔ respiration coupling | **olfactoryHRV agent** | `../olfactoryHRV/reports/{rsa,grant}_report.html` (repo-root dir; no shared upstream) |
| **Shared upstream** | serialize edits | `code/extract_gamma_session.m`, `code/extract_spectro_session.m`, produced data tables |

Two separate reports share one upstream: the two `extract_*_session.m` extractors and
the per-breath data tables. **Keep edits to those two extractors serialized between
agents, and re-run the downstream aggregators afterward.**

---

## main-report / gamma-pipeline agent

Owns the primary integrated report and the gamma per-breath + spectrogram pipeline.
All work through commit `9d87f69` is committed & pushed.

**Owns — coordinate before editing (these get regenerated):**
- **Report:** `report.Rmd` → `out/report.html`
- **Code (`code/`):** `breathing_deficit.R`, `gamma_aggregate.R`, `assemble_spectrograms.m`,
  `gamma_review_response.R`, `coupling_aggregate.R`, `behavioral_composite.R`, `macbp_table.R`,
  `inventory_report.R`, `inventory_reverse.R`, `build_session_index.m`, `run_scores.m`,
  `run_gamma_all.m`, `run_spectro_all.m`, `lib_scores.m`, `ridge_track.m`, `slt_power_cont.m`,
  `scan_breathing_labels.m`, `make_review_packet.R`, `sim_broadband_transient.m`, `validate_slt.m`,
  and the `*.mat` rules in `.gitignore`
- **Outputs:** `out/figs/spectrograms_{low,mid,high}.png`, `out/figs/beh_*.png`,
  `out/figs/gamma/*.png`, `out/figs/breathing_deficit_tracking.png`;
  `out/tables/{gamma_*, coupling_*, breathing_deficit_*, responder_table, macbp_*,
  session_index, INVENTORY_*, session_scores}.csv`; `out/spectro2/*.mat` (git-ignored);
  `out/report.html`

**Produces, others consume — read freely, do not overwrite:**
`out/tables/gamma_session_level.csv`, `responder_table.csv`, `session_scores.csv`.

**Does NOT touch:** `code/taskcmp_*`, `cluster_perm_audio_focus.m`, `scratch_*`,
`integrate_new_controls.m`, `driver_rerun.m`, `taskcmp_report.Rmd`, `out/figs/taskcmp/*`,
`out/tables/taskcmp_*`, `out/taskcmp_*`, `../olfactoryHRV/*`, `../CLAUDE.md`, `../CombinedTaskList`.

---

## taskcmp agent

Owns the control **audiobook (=baseline) vs focusedBreathing (=ATB)** comparison, paired
within OBE controls. **Fully committed & pushed at `0346bf9`; currently holds zero
uncommitted files.**

**Report lives at the PROJECT ROOT, not `out/`:** `taskcmp_report.Rmd` → `taskcmp_report.html`
(self-contained, figures embedded). The old `out/taskcmp_report.html` + `out/taskcmp_synthesis.md`
were stale pre-correction copies and were **deleted 2026-09-09** — don't recreate them.

**Owns (don't edit):**
- `code/taskcmp_*.{R,m}` — incl. `taskcmp_ridge_gated.m`, `taskcmp_theta_early.m`,
  `taskcmp_superlet_blocks.m`, `taskcmp_ridge_diff.m`, `taskcmp_ridge_extrema.m`,
  `taskcmp_ampratio.m`, `taskcmp_bandmetrics.m`, `taskcmp_groups_gamma.m`, `taskcmp_latency.m`,
  `taskcmp_bandz.m`, `taskcmp_envoverlay.m`, `taskcmp_grantfigs.R`, `taskcmp_C2_decoding.R`,
  `taskcmp_D_analysis.R`, `taskcmp_D_coupling.R`, `taskcmp_E_airflow.R`, `taskcmp_power.R`,
  `taskcmp_B_report.R`, `taskcmp_verify.R`
- `code/cluster_perm_audio_focus.m`, `code/scratch_ridge_bandmax.m`, `code/scratch_theta_pac.m`,
  `code/integrate_new_controls.m`, `code/driver_rerun.m`
- `taskcmp_report.Rmd`, `taskcmp_report.html`
- `out/figs/taskcmp/*`
- `out/tables/{taskcmp_*, scratch_*, breathing_blocks.csv, superlet_FWHM.csv}`
- `.gitignore` rule `out/tables/*.mat` (my regenerable per-breath/ridge caches; CSV summaries ARE tracked)

**Consumes (read-only):** `out/tables/{gamma_session_level.csv, responder_table.csv,
session_scores.csv}` (main-report agent). Note: taskcmp otherwise reads the preprocessed finals
**directly from R:** for its own superlet/ridge extraction — it does NOT go through the main
gamma aggregators, so most taskcmp reruns don't need the main pipeline.

**Shared files it edits — COORDINATE:**
1. `code/extract_gamma_session.m` + `code/extract_spectro_session.m` — in `0346bf9` I extended
   **only the breathing task-label matcher** (added the separatepreproc control labels
   `audiobook` and `focusedBreathing`, alongside the existing `audio` / `focus|naturalFocus|slowFocus`).
   This is a matcher change, **not** a new-metric change. Effect: separatepreproc-named control
   finals (e.g. HM_2, SP_2) now yield audiobook/focusedBreathing breaths that were previously
   dropped → the `out/gamma/perbreath` CSVs for those sessions change. Ping the main-report agent
   and re-run its aggregators after any further edit here.
2. **`out/tables/macbp_best.csv` — SHARED, IMPORTANT.** This is a main-pipeline product
   (`run_scores.m`), but I **appended two control rows** (HM_2, SP_2 `breathingTask` best-macBP)
   in `0346bf9`. **If the main pipeline regenerates `macbp_best.csv`, it must include HM_2 and SP_2,
   or my appended rows are lost.** The taskcmp report + grant figs depend on those two rows for the
   8-session control set.

**Scientific note (for cross-report consistency):** the per-session gamma **peak-latency "timing"
effect is a confirmed argmax artifact** — block grand-mean trajectories show *baseline* peaking
first (~530 ms vs ATB ~1090 ms), and both the raw & z>3-gated ridge-frequency sweeps and the
early-window (0–500 ms) theta ITPC/power are null. The **multivariate decoder** (5/6 subjects,
Fisher p=7e-7, survives detrend) is the sole robust condition difference. Please don't cite
"focus gamma is earlier" as a finding in any shared write-up.

---

## olfactoryHRV agent

Owns everything under **`../olfactoryHRV/`** (repo root) plus **one section of `../CLAUDE.md`**.
No overlap with `plotsForBruce_2/`: this analysis reads the preprocessed finals **directly from
`R:`** and shares none of the gamma pipeline, the `extract_*_session.m` extractors, or the produced
data tables. Nothing here needs serializing against the other two agents.

**Project:** does the respiration→HRV transfer function change across sessions, and does that change
track olfactory recovery in the Dupi cohort? Per-breath respiratory sinus arrhythmia regressed on
breath depth and duration, per-session slopes calibrated against surrogate nulls, change scores
against a three-task olfactory composite. Outputs: `reports/rsa_report.html` (full analysis, with
reliability gate and multiplicity correction) and `reports/grant_report.html` (preliminary-data
write-up), plus the finalized grant figures.

**Owns — do not edit:**
- `../olfactoryHRV/**` — the entire directory, all authored here.
  Committed in `d50cb03`: `README.md`, `ohrv_config.m`, `run_all.ps1`,
  `rsa_{extract,olf,analyze,panel,report,figs}.m`, `thresh_analysis.m`,
  `fig_{grant,grant_thresh,duration,raw}.m`, `get_beats.py`, `mk{report,grant}.py`, `qc/*`,
  `reports/*.html`, `figures/{fig*,grant*,report*}.png`, `.gitignore`.
  **Currently uncommitted (9 files):** `fig_final.m`, `fig_final3b.m`, `ohrv_colors.m`,
  `figures/FIG{1,2,3}_*.png`, `figures/FIG3b_{JH,PC}_depth_boxes.png`.
- `../olfactoryHRV/work/` — **git-ignored**, ~91 MB analysis cache (29 `*_slim.mat` per-breath
  extracts, 29 `*_beats.npz`, 7 summary CSVs). Fully regenerable by `run_all.ps1`. Never commit.

**Shared file it edits — COORDINATE:**
- **`../CLAUDE.md`** — added **§7.1 “Behavioral scoring (olfactory performance per session)”**
  (+31 lines, a clean insertion inside §7): the cue d′ / O15 / threshold formulas, and the rule that
  the `behDat` table in the final `.mat` is the canonical scoring source (the `R_groupLevel` CSVs are
  a convenience copy, not the source). Nothing else in that file is touched — edit other sections
  freely and the diffs will not conflict.

**Does NOT touch:** anything under `plotsForBruce_2/`, `../CombinedTaskList`, or any `code/extract_*`.

### Gotchas in the shared preprocessed data — relevant to every agent reading breathing `behDat`

1. **`noseMouth` is empty, not `"mouth"`.** No Dupi session carries an explicit `"mouth"`; the values
   are `"nose"`, `"NA"` and `""`. Filtering `noseMouth == "nose"` silently discards **31–46% of breaths
   in 11 sessions**, concentrated in the `audio` and `focus` blocks, and can collapse a session to a
   single usable block. Respiration is recorded with a **nasal cannula**, so a detected breath with a
   valid signal was necessarily nasal — empty means unrecorded metadata. **Filter `~= "mouth"`
   instead.** This silently dropped most of my newly-preprocessed sessions before I caught it.
2. **The finals on `R:` are rewritten in place by reprocessing runs.** Reading them mid-batch mixes
   two data versions and raises no error. I analysed such a mixture once and had to discard the
   results. Before a run, confirm no `preProc/*.mat` has been written for several minutes, then clear
   cached extracts.
3. **Threshold-task `intensity` is a screen pixel coordinate** on an ~800 px track (floor 335,
   ceiling ~1135), not a rating value. Subtracting the blank-air trial removes both the pixel offset
   and the response bias; `(high − air)` carries the dynamic range, `(med − air)` sits near floor for
   most patients.
4. **Known corrupt file:** `250811_Dupi_NMH_TB_2_PEA_threshold_preproc.mat` fails to open with an
   HDF5 `inflate()` error on repeated attempts days apart. Needs re-preprocessing.

**Scientific note (don’t cite the wrong measure):** regressing within-breath HRV on breath
**duration** carries a built-in positive slope — a longer breath is a longer window, more beats fall
inside it, and the range of more samples is larger. Across 25 sessions the surrogate null for that
slope runs **+0.44 to +1.03**, the size of the raw effect itself, and the measure fails its own
reliability check. Breath **depth** does not set the window length; its null is **−0.011 to +0.022**,
i.e. zero. Use the depth version, and don’t report a duration-slope change as a coupling result.

---

## Shared / conflict zones (all agents read)

1. **`code/extract_gamma_session.m` + `code/extract_spectro_session.m`** — the per-breath and
   spectrogram extractors. They carry the breathing task-label fix and the p3 landmark fix and
   are edited by multiple agents. **Serialize edits; re-run the affected aggregators after.**
2. **Produced data tables** (`gamma_session_level.csv`, `responder_table.csv`, `session_scores.csv`)
   — produced by the main-report agent, consumed by taskcmp. Read-only for consumers.
3. **`out/tables/macbp_best.csv`** — main-pipeline product that taskcmp **appended to** (HM_2, SP_2
   control rows, `0346bf9`). If the main pipeline regenerates it, keep HM_2/SP_2 or taskcmp's
   8-session control set breaks. Serialize + preserve those rows.
4. **`.gitignore`** — `out/spectro2/*.mat`, `out/gamma/spectro/*.mat`, and `out/tables/*.mat`
   (taskcmp's per-breath/ridge caches) are intentionally ignored (bulky regenerable
   intermediates). Don't commit `.mat` intermediates.

---

## Bringing the repo up to a current version (olfactoryHRV agent’s view)

Written after three agents worked one branch in parallel. Suggested order:

1. **Fix `../CLAUDE.md` first — it is the shared contract, and three documented facts are stale.**
   Every agent reads it to learn the data structures, so each wrong entry gets rediscovered
   independently, at cost:
   - Breathing `behDat` is documented as **33 columns**; the finals on `R:` have **60**, including
     `manOnset` and 26 `bm_*` breathmetrics columns (`bm_inhaleVolumes`, `bm_inhaleDurations`, …)
     plus an `outDat.bmFeatures` struct. The files were re-segmented 2026-08-28 by the sibling repo
     `zelanoLabPreprocessing`; **this** repo’s `process_respiration_breathing.m` still emits the old
     layout, so the code and the data on disk disagree. Document which one produced the files you
     are reading.
   - §1’s sheet `Task` strings are wrong. The live sheet uses `breathing tasks` (with a space — the
     single most common spelling), `odor cue task` (the documented `odorCueTask` appears **zero**
     times) and lowercase `threshold`. Any exact-match filter written from the doc drops rows
     silently.
   - §3’s tracking-sheet path is dead. `labPaths.m:107` builds `…\Admin\dataTracking.xlsx`, which no
     longer exists; the live sheet is `…\Admin\Data\dataTracking.xlsx`. `applyParams` then falls back
     to a repo-local copy **without erroring**, so a stale sheet can be used without anyone noticing.

2. **Commit by directory, not by working session.** `d50cb03` swept `olfactoryHRV/` into a commit
   titled for `plotsForBruce_2`. Nothing broke, but it makes `git log -- <dir>` useless for tracing
   who changed what. One commit per top-level area, prefixed with that area’s name.

3. **Land the shared upstream before the per-agent work.** `code/extract_*_session.m`, the produced
   data tables and `out/tables/macbp_best.csv` are the only real conflict surface in this repo.
   Settle those, re-run the aggregators, confirm the consumers still read clean, then commit each
   report directory on top of a known-good upstream.

4. **Adopt one cache rule across agents.** Every agent caches extracts derived from the `R:` finals,
   and those finals get rewritten in place. Whatever the mechanism — a stamp file recording source
   mtimes, or simply “clear the cache at the start of every run” (what `olfactoryHRV/run_all.ps1`
   does) — it should be the same rule everywhere, because a stale cache fails **silently** rather
   than loudly. This is the highest-value convention on the list: it is the one that has already
   destroyed work.

5. **Then merge `dupi-gamma-analysis` into `main` as one reviewed unit.** The branch now carries
   three independent analyses that share upstream data. Merging piecemeal loses the chance to verify
   that the shared upstream is consistent across all three at a single commit.

_Last updated 2026-09-09 by the **olfactoryHRV agent** (filled in own section: file inventory, `CLAUDE.md` §7.1 claim, shared-data gotchas, duration-vs-depth note; added the repo-currency section above). Prior update by the taskcmp agent (augmented own section; corrected report path to
project root, extractor edit = label-matcher not new metrics, flagged `macbp_best.csv` as shared,
noted the peak-latency-artifact conclusion). Prior update by the main-report / gamma-pipeline agent._
