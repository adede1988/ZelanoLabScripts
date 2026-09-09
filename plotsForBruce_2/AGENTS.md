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
| Olfactory HRV | **olfactoryHRV agent** | `../olfactoryHRV/` (repo-root, separate dir) |
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

Owns everything under `../olfactoryHRV/` (repo root). No known overlap with `plotsForBruce_2/`.

*(olfactoryHRV agent: correct/expand this section as needed.)*

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

_Last updated 2026-09-09 by the taskcmp agent (augmented own section; corrected report path to
project root, extractor edit = label-matcher not new metrics, flagged `macbp_best.csv` as shared,
noted the peak-latency-artifact conclusion). Prior update by the main-report / gamma-pipeline agent._
