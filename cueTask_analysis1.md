# cueTask_analysis0 — analysis plan (macBP gamma + event-locked spectrograms)

Plan only — **no implementation yet**. This catalogues five tasks that consume the
**cue-task preprocessed `.mat` files** and produce FOOOF summaries, per-subject and
group-mean spectrograms, and an R Markdown report. It is written against the data
structures documented in `CLAUDE.md` §6.

---

## 0. Scope, inputs, and shared conventions

### 0.1 Data source (what each file gives us)
- **Files:** every cue final `<root>\<id>\preProc\<id>_cueTaskPreproc.mat`, top-level var
  **`outDat`** (load via `fieldnames`; cue is always `outDat`). Session list from
  `applyParams('cueTask','main')` (≈35 sessions), filtered to those that actually exist on
  disk **and** have ≥1 `macBP` channel.
- **Signal:** `outDat.data` `[nChan × nSamp]` @ **`outDat.fs = 500 Hz`**; channels addressed
  by `outDat.labels`.
- **macBP channels:** `find(cellfun(@(x) contains(x,'macBP'), outDat.labels))` → `macBP1..5`
  (bipolar macro pairs; **not all files have all 5**; count = nMacro−1).
- **Events (all in 500 Hz samples):**
  - `outDat.TTL.trialStart` — trial-onset sample per trial (cue `TTL` table is `[nTrial×3]`
    `{trialStart, response, sniff}`).
  - `outDat.behDat.finalOnset` — phase-refined sniff onset per trial (cue = one sniff/trial,
    `wiTriali==1`); also `outDat.behDat.n` (trial number) and `type` (`hit/miss/cr/fa`).
- **Identity:** `outDat.sessID` (e.g. `250623_Dupi_NMH_KS_2`), `outDat.type`
  (`Dupi`/`OBE`/`EEG`). Parse `subID`/`sessNum` from the sessID like
  `splitToSingleChan_allTasks.m`: `bits=strsplit(sessID,'_')` → `subID=bits{4}`,
  `sessNum=str2double(bits{5})` (default 1 if absent).
- **Figure folders:** per-subject = recompute portably as
  `fullfile(labPaths().figPath, sessID, 'cueTask')` (≡ stored `outDat.figs`); create if
  absent. Group = `R:\Neurology\Zelano_Lab\Lab_Common\Adam\Dupi_processing\groupStatFigs`.

### 0.2 Toolboxes & coding rules (per request)
- Use **EEGLAB** (`C:\Users\Adam\Documents\eeglab2026.0.0`) and **FieldTrip**
  (`C:\Users\Adam\Documents\fieldtrip-20260518`) functions; write only glue/formatting code.
- **Pipe to FieldTrip *through* EEGLAB** (`eeglab2fieldtrip`) so formatting stays consistent.
- Prefer EEGLAB unless a FieldTrip function is specifically needed.
- **Resolved tool choices (from clarifying Q&A):**
  - **FOOOF → FieldTrip** `ft_freqanalysis` (`cfg.output='fooof_peaks'` / `'fooof_aperiodic'`;
    uses FieldTrip's bundled brainstorm `process_fooof`, no Python).
  - **Spectrograms → EEGLAB `newtimef`** (log-spaced freqs, dB baseline); group means by
    averaging the per-subject dB ERSP matrices.
- **Path hygiene (gotcha):** EEGLAB and FieldTrip share many function names. Plan to
  `ft_defaults` for the FOOOF step and run EEGLAB (`eeglab nogui`) for the spectrogram step,
  managing `addpath` order so neither shadows the other (and `eeglab2fieldtrip` lives in the
  EEGLAB `dipfit` plugin). Document the working path order in the implementation.

### 0.3 Time/sample conventions
- `1 sample = 2 ms` (fs=500). Helper: `ms→samp = round(ms/1000*fs)`.
- Spectrogram **display** window: **−1000 … +3000 ms**. **Epoch** (extraction) window is
  padded to **−1750 … +3750 ms** to avoid wavelet edge effects at 2 Hz (≈3 cycles = 1500 ms);
  display is cropped to −1000…+3000 ms after TF estimation.
- **Baseline** (both lockings): **−700 … −200 ms relative to `trialStart`**.
- Drop any event whose padded epoch would run off either end of the recording.

### 0.4 Groups (updated to FOUR groups)
Assign each cue session by `type` + `sessNum`:
| Group | Rule |
|---|---|
| `DupiS1` | `type=='Dupi'` & `sessNum==1` |
| `DupiS2` | `type=='Dupi'` & `sessNum==2` |
| `DupiS3` | `type=='Dupi'` & `sessNum==3` |
| `Control` | **all non-Dupi** cue sessions (OBEControl + EEG studies pooled) |
Dupi sessions with `sessNum>3` (if any) are reported in the data table but excluded from the
4-group means unless you say otherwise.

### 0.5 Trial inclusion (assumption to confirm)
Use **all trials** for spectrograms (not split by `hit/miss/cr/fa`). Flagged in §6.

---

## Task 1 — redo baseline normalization

**Goal:** z-score power time series for use in spectrograms using a z-score permutation method against the same baseline period as before

**Steps**
1. Load `outDat`; pick the `bestMac` channel (from Task 1's stored label); events =
   `outDat.TTL.trialStart` (valid, padded epochs only).
2. **Epoch via EEGLAB.** Build a 1-channel `EEG` (the bestMac signal), insert `trialStart`
   events, `pop_epoch` to **[−1.75, +3.75] s** → `EEG.data` `[1 × pnts × trials]`.
3. **TF via `newtimef`** on that channel:
   - `'freqs',[2 120], 'nfreqs',100, 'freqscale','log'` (log-spaced 2–120 Hz);
   - `'cycles',[3 0.8]` (Morlet, increasing cycles — tunable);
   - `'baseline',[-700 -200]` (ms, rel. trialStart), `'trialbase','full'`;
   - Unlike previous plots, this time use myChanZscore.m to z-score data against a baseline. Prior to submitting data to z-scoring, do a simple z-score at every time frequency point within the baseline period. Any trials with any z-scored values greater than 10 or less than -10 should be eliminated from use in the z-scoring
   - After z-scoring and before plotting, eliminate trials with z-score values greater than 15 at any time frequency point
   - make a table reporting how many trials were removed at each of these steps for all subjects in the final report. Facilitate this by adding noise trial elimination information into the .csv file that you already have been building for fooof output. This .csv will become a general place for all summary info across the analysis. 
   - **Persist `powbase`** (per-freq baseline power for this subject's bestMac) — Task 3 reuses it. A good strategy might be to extract power values from a long epoch that covers both the trialStart and the finalOnset locked plots, then z-score all of it at once.
4. **Plot ourselves** for axis control: `imagesc(times, 1:nFreq, ersp)` with **log-frequency
   y-axis labeled with raw Hz** (e.g. ticks at 2,4,8,16,32,64,120), x-axis cropped to
   **−1000…+3000 ms**, symmetric diverging colormap, colorbar = "dB vs baseline", vertical
   line at t=0. Title = `sessID` + bestMac label.
5. Save to subject cue figs, e.g. `<id>_<bestMac>_TFR_trialStart.png`; also save the numeric
   `ersp/times/freqs/powbase` to a per-subject `.mat` for Task 4 averaging. Overwrite figures and values from prior run.
6. Add mean respiration trace as an overlay on the figure. 


**Outputs:** per-subject trialStart-locked spectrogram PNG + numeric TFR `.mat` 


---

## Task 2 — `bestMac` spectrogram locked to `finalOnset` (same trialStart baseline)

**Goal:** identical spectrogram but time-locked to `behDat.finalOnset`, **normalized to the
same trialStart baseline** so the two figures share a color scale.

**Key method (cross-locked baseline).** `newtimef`'s built-in baseline would be relative to
`finalOnset`; we must instead reuse Task 2's **`powbase`** (the trialStart baseline). So:
1. Epoch the bestMac channel around `finalOnset` (**[−1.75, +3.75] s**, EEGLAB `pop_epoch`).
2. `newtimef` with the **same** freq/cycle settings but **`'baseline', NaN`** → returns the
   **raw dB power** TFR (`erspRaw`, no baseline removal), same `times`/`freqs` grid as Task 1.
3. **Apply the stored baseline:** Applying the stored baseline now could be done in a few ways. You could z-score using a single long epoch per trial that's big enough to cover both plot windows. That way you do it all in one shot, then carve it up for plotting. This would mean being careful with variable time offsets between trialStart and finalOnset in order to get the plots aligned correctly. Or you could tack on the baseline period to the finalOnset period and then do the z-scoring anew
4. Use the same noise removal rules as in Task 1
5. Plot exactly as Task 1 (same axes, same color limits) and save
   `<id>_<bestMac>_TFR_finalOnset.png` + numeric `.mat`.
6. Add mean respiration trace as an overlay on the figure. 

**Outputs:** per-subject finalOnset-locked spectrogram PNG + numeric TFR `.mat`. Overwrite figures and values from prior run

---

## Task 4 — Group-mean spectrograms (4 groups × 2 lockings)

**Goal:** mean spectrograms aggregating the per-subject TFRs from Tasks 1 & 2, per group.

**Steps**
1. Assign each subject to `DupiS1 / DupiS2 / DupiS3 / Control` (§0.4).
2. For each (group × locking): load the per-subject numeric TFR `.mat`s (identical
   `times`/`freqs` grids by construction), **average the z-score matrices across subjects**
   (track *n* per cell; equal-weight subjects). Optionally also save the across-subject SEM.
3. Plot with the **same axes/scale** as Tasks 1–2 (log-freq raw labels, −1000…+3000 ms,
   shared color limits across all group plots for comparability), title = group + locking + *n*.
4. Save **8 figures** (4 groups × {trialStart, finalOnset}) to
   `…\Dupi_processing\groupStatFigs`, e.g. `group_<Group>_TFR_<locking>.png`, plus the numeric
   group means as `.mat`. Overwrite prior run outputs. 
5. Add mean respiration trace as an overlay on the figures. 

**Outputs:** 8 group-mean spectrogram PNGs (+ numeric means) in `groupStatFigs`.

---

## Task 5 — R Markdown report (`.Rmd`)

**Goal:** update the R Markdown with the new figures and outputs

**Sections**
1. **Methods** — preprocessing recap (downsample 500 Hz, macBP bipolar), FOOOF via FieldTrip
   (settings, gamma rule + fallback), spectrograms via EEGLAB `newtimef` (freqs, cycles,
   baseline, cross-locked baseline for finalOnset), group definitions.
2. **Data availability table** — one row per session: `subID, sessID, group, #macBP,
   bestMac, bestMacMethod, nTrials` (from the CSV + a small MATLAB-exported `data_manifest.csv`).
3. **FOOOF findings table** — rendered from `cueTask_fooof_summary.csv` (e.g. `knitr::kable`/
   `DT`), highlighting `isBestMac` rows and gamma-detected channels.
4. **Single-subject examples** — for **each of the 4 groups**, **2 example subjects**, show
   **both** spectrograms (trialStart + finalOnset) ⇒ 4×2×2 = **16 embedded PNGs**
   (`knitr::include_graphics`).
5. **Group averages** — the **8** group-mean spectrograms.
6. (Optional) brief gamma-power group comparison from the FOOOF CSV.

**Inputs the Rmd needs from MATLAB:** `cueTask_fooof_summary.csv`, `data_manifest.csv`, and
the figure files at known paths (the Rmd takes a config block with the `groupStatFigs` path
and a list of example subject IDs per group).

**Outputs:** `cueTask_report.Rmd` (+ knitted HTML).

---

## 6. Assumptions:
1. **Session set:** all cue finals with ≥1 `macBP` channel (skip macro-less sessions). 
2. **Trials:** all trials pooled (not split by `hit/miss/cr/fa`).
4. **Group color scaling:** share one color-limit across all group plots (max comparability)
   vs per-plot autoscale. Default = shared.
5. **`newtimef` params:** `nfreqs=100`, `cycles=[3 0.8]`, epoch pad ±1.75 s
6. **CSV granularity:** one master CSV (all channels, all subjects).
7. **`bestMac` re-save:** overwrite the cue `.mat` in place adding `outDat.bestMac`
   (+`bestMacMethod`). bestMac should already be in the outDat structs so no overwriting is needed on this step
8. **DupiS>3 / EEG-study cue sessions:** reported in the manifest, excluded from the 4 group
   means. 

## 7. Deliverables inventory
- Updated cue `.mat` files: `outDat.bestMac`, `outDat.bestMacMethod`.
- Per subject: `<id>_macBP_fooof_periodic.png`, `<id>_<bestMac>_TFR_trialStart.png`,
  `<id>_<bestMac>_TFR_finalOnset.png` (+ numeric TFR `.mat`s).
- Group: 8 `group_<Group>_TFR_<locking>.png` (+ numeric means) in `groupStatFigs`.
- Tables: `cueTask_fooof_summary.csv`, `data_manifest.csv`.
- Report: `cueTask_report.Rmd` (+ knitted HTML).

## 8. Proposed code layout (for the implementation pass)
- `cueAnalysis/run_task1_fooof.m` · `run_task2_trialStart_TFR.m` ·
  `run_task3_finalOnset_TFR.m` · `run_task4_group_means.m` · `cueTask_report.Rmd`
- Shared helpers: `cue_session_table.m` (list+group+macBP inventory),
  `cue_make_EEG.m` (outDat→EEGLAB EEG for chosen channels), `tfr_newtimef_channel.m`
  (epoch+newtimef wrapper returning `ersp/times/freqs/powbase`), `plot_TFR.m`
  (shared log-freq dB plotter so Tasks 2–4 look identical).
