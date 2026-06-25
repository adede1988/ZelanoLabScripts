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

## Task 1 — Quantifying oscillations in epochs

**Goal:** apply fooof to epochs between -500 ms and + 1500 ms relative to finalOnset to track oscillatory changes

**Steps**
1. Load `outDat`; pick the `bestMac` channel (from earlier work's stored label); events =
   `outDat.behDat.finalOnset` (valid, padded epochs only).
2. **Epoch via EEGLAB.** Build a 1-channel `EEG` (the bestMac signal), insert `finalOnset`
   events, `pop_epoch` to **[−1.75, +5.75] s** → `EEG.data` `[1 × pnts × trials]`.
3. Utilize Morlet Wavelets to extract power time series for frequencies spanning 5 Hz to 58 Hz in 100 linear spaced frequencies spanning the range. 
4. Take the mean power across trials and timepoints in 250 ms windows to make epoch specific power spectra spanning the range -500 ms to + 2000 ms. This should yield 10 power spectra
5. Apply fooof to each power spectrum, plot the periodic power spectra at the subject level with the different time epochs represented by different colored lines that progress through a logical progression from yellow ocher to purple. Save this per-subject figure into the subject figs folder under the name gammaTimeProgression.png. Also, save the largest detected peak from each epoch in the 30 to 58 Hz range into the overall .csv file that you've been building throughout this analysis.  
6. Construct a plot showing the number of participants with a detectable oscillatory peak in the gamma range in each epoch. Make this plot for each group. Add these plots to the .rmd summary
7. Subtract from each participant's epoch-detected oscillatory peaks the value of the channel's overall fooof peak that was detected earlier. Then plot the group average of these channel-centered values across the epochs. Show individual points as lighter connected dots so all data are shown in the one plot. 
8. Add 4 examples (one from each group) of the gammaTimeProgression.png plots to the .rmd. 
9. Add the group average channel-centered epoch gamma peak plot to the .rmd. Perform stats using the epoch peaks in the .csv to ask if the gamma peak frequency changes as a function of epoch. Do this analysis within each group independently using linear mixed effects modeling. Treat epochs with no fooof detected gamma peak as missing values in this analysis. 


**Outputs:** per-subject png, updated values to summary .csv file, group channel-centered gamma peak figure, statistical analysis of gamma peak changes over epochs, all added into .rmd output render


---
