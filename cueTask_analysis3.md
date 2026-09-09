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

**Goal:** evaluate possible noise 

**Steps**
1. Load `outDat`; pick the `bestMac` channel (from Task 1's stored label); events =
   `outDat.TTL.trialStart` (valid, padded epochs only).
2. **Epoch via EEGLAB.** Build a 1-channel `EEG` (the bestMac signal), insert `trialStart`
   events, `pop_epoch` to **[−1.75, +5.75] s** → `EEG.data` `[1 × pnts × trials]`.
3. Make a plot showing all single trial epochs on top of each other separated vertically by 50 uV
4. Flag trials with sharp deflections by plotting them in red in the plot. Sharp deflections are defined as a greater than 80 uV max to min difference in a 10 ms period. In addition to visually flagging these trials, eliminate them from use throughout the rest of the analyses. This replaces the noise rejection rules of previous analyses. 
5. save to the per subject figs folder under the name "singleTrialRawMac.png"
6. Propagate the change in noise detection throughout all analyses in the pipeline.
7. Add 4 example plots to the html output, 2 with a high number of noise trials and 2 with a low number of noise trials


**Outputs:** per-subject png, updated noise detetion, update of all subsequent analyses


---
