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
| Audiobook-vs-focusedBreathing comparison | **taskcmp agent** | `taskcmp_report.Rmd` → `out/taskcmp_report.html`, `out/taskcmp_synthesis.md` |
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

Owns the control audiobook-vs-focusedBreathing comparison.

**Owns:** `code/taskcmp_*.{R,m}`, `code/cluster_perm_audio_focus.m`, `code/scratch_ridge_bandmax.m`,
`code/scratch_theta_pac.m`, `code/integrate_new_controls.m`, `code/driver_rerun.m`;
`taskcmp_report.Rmd`; `out/figs/taskcmp/*`; `out/tables/taskcmp_*.csv`;
`out/taskcmp_report.html`; `out/taskcmp_synthesis.md`.

**Consumes (do not overwrite):** `out/tables/gamma_session_level.csv`, `responder_table.csv`,
`session_scores.csv` (produced by the main-report agent).

**Shared upstream it also edits:** `code/extract_gamma_session.m`,
`code/extract_spectro_session.m` (added gated-ridge / early-theta metrics in `0346bf9`).
Any change here alters the `out/gamma/perbreath` CSVs the main report is built from — ping
the main-report agent and re-run its aggregators after.

*(taskcmp agent: correct/expand this section as needed.)*

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
3. **`.gitignore`** — `out/spectro2/*.mat` and `out/gamma/spectro/*.mat` are intentionally
   ignored (bulky regenerable intermediates). Don't commit `.mat` intermediates.

_Last updated 2026-09-09 by the main-report / gamma-pipeline agent._
