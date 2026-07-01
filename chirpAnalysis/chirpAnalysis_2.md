# chirpAnalysis.md — OB sniff-locked gamma: one chirping oscillator, or two?

**Analysis-stage spec for Claude Code.** Sits alongside `CLAUDE.md` (the pipeline +
data-structure reference). This document defines a new downstream analysis that *consumes*
the preprocessed finals
It does not touch the preprocessing pipeline.

---

## 0. The question and why it is hard

In OB sniff-locked single-trial electrophysiology, the sniff drives a gamma power increase
that exhibits a **downchirp**: instantaneous frequency starts high and descends. The question:

> Is this **one oscillator modulating its frequency** down over the trial, or **two
> oscillators** — a higher-frequency one earlier and a lower-frequency one later — whose
> aggregate *looks* like a chirp but is actually two separate generators?

### 0.1 The identifiability problem (read before building anything)

A single descending sweep and two fixed-frequency components that overlap in time and
phase-couple can produce nearly identical short signal segments. That degeneracy is *why* the
aggregate reads as a chirp. So the bare "one vs two generators" question has a region where no
fit can separate the hypotheses. We therefore answer the sharper, decidable question:

> **During the frequency transition, is there a single spectral component whose center
> frequency moves, or two simultaneously-present components at different (slowly-varying)
> frequencies?**


## 1. Scope & current targets

- **Now:**  `cueTask`. **Eventually:** all four tasks (`breathingTask`, `cueTask`,
  `threshTask`, `O15`).
- Code must be **task-agnostic** in its signal processing; branch only where the data layout
  forces it (sniff selection, epoch guarding). Mirror how the existing pipeline keeps the
  shared signal path identical across tasks.

---

## 2. Integration contract (data, paths, idioms)

### 2.1 Input

- **the session finals** `<root>\<id>\preProc\<id>_<task>preproc.mat`. 

- **`fs = 500 Hz` for all final signals** (2 ms/sample). Gamma 25–58 Hz is well within Nyquist.

### 2.2 Channel selection

- The outDat struct loaded from the session preproc file should include a field called `bestMac` which indicates the string label for the target channel. The labels field indicates all channel labels in the same order as the chanels X time data matrix. 
  **Never hard-code a row index** — montages vary by session.


### 2.3 Epoching

- Build sniff-locked epochs around **`behDat.finalOnset`**. finalOnset indexes directly into the channels X time data matrix.
  Use `finalOnset`, not `sniffOnset`.
- **Window:** Use a window of -1s to +3s around the finalOnset index. That is, -500 timepoints to +1500 timepoints at 500 Hz sampling rate
- **Sniff selection (task-specific config, not signal-path branching) key information lives in behDat table:**
  - `cueTask`, `breathingTask`, and `threshTask`: one sniff/trial simply take `behDat.finalOnset` values and align to them.
  - `O15`: multiple sniffs/trial (`sniffLabel ∈ {start,free,confirm}`, `wiTriali` is variable). Make
    the selected sniff(s) a config field (default to the task-relevant sniff).
- **O15 long-window guard:** truncate each epoch at the **next** `finalOnset` so a long window
  doesn't run into the following sniff; set a `burstTruncated` flag where the window is cut
  short, so a truncated coexistence window cannot masquerade as a genuinely brief one.
- **breathingTask bad data:** for the breathingTask, use `behDat.goodBreath` as an indicator for which finalOnsets to use. Only use finalOnsets where goodBreath==1
- Again, at current, we are only focusing on the cueTask. We will loop back to other tasks later 

### 2.4 Trial QC — defer to the existing pipeline

Utilize the cue_noise_trials.m function found in the cueAnalysis folder of ZelanoLabScripts to do noise rejection. Calibrate K to maintain at least 80% of data across the whole dataset, but individual sessions may maintain less if they are noisey. Apply noise detection to the `bestMac` channel. 


### 2.5 Figures

Save under the session figure folder (`chanDat`/`outDat.figs`, i.e.
`…\Lab_Common\Adam\Dupi_processing\<id>\`) in a new subfolder **`singleTrialSpectrograms\`**.
Filename: `sub-<sessID>_ch-<chi>_trial-<NNN>.png`. It's okay to overwrite files from prior runs in this folder. 

### 2.6 Driver & code layout

- Write a new overall script called "chirpAnalysisV2.m" All other functions/scripts should be called from here to make it easy for me to follow the code after you've written it. Any of the helper functions that you do not use from the earlier analysis that are already in the chirpAnalysis folder, move them to a new subfolder called "oldAnalysis"

### 2.7 Compute

- Finals are **`-v7.3` (HDF5)** — partial-load with `matfile`/`h5info`; don't pull the whole
  `data` matrix into RAM if you can avoid it.
- **Run on the lab desktop over SSH** (`ssh labdesktop`; 8-core / ~96 GB; no GPU). `parfor`
  over trials within a channel. Write outputs to `E:\` (never `C:\`). Code synced via git. You'll need to update paths to prefer new GitHub folder on `E:\GitHub` on the lab machine. Profile **one** session locally before launching the batch.

---

## 3. External tools (clone + register paths)

Clone into the GitHub dir on each machine (`C:\Users\Adam\Documents\GitHub\` home /
`G:\My Drive\GitHub` lab), sync via git, and register paths in `labPaths.m` (a new block beside
the eeglab root) or a `setup_chirpAnalysis_paths.m` so both machines resolve identically. All
MATLAB except where noted.

| Tool | URL | Role |
|---|---|---|
| **Superlets** | `https://github.com/TransylvanianInstituteOfNeuroscience/Superlets` | single-trial TFRs; use the **FASLT** (fractional adaptive) variant for a smooth sweep |
| **MPACT** | `https://github.com/jiecui/mpact` | adaptive chirplet transform via matching pursuit + EM refinement (MPEM); noise-robust, built for biosignals — the chirplet engine, do **not** hand-roll MP |
| **Synchrosqueezed chirplet** | `https://github.com/ziyuchen7/Synchrosqueezed-chirplet-transforms` | Chen & Wu 2023; disentangles modes with curved/crossing instantaneous frequencies — the right tool for the converging-IF read over long bursts |
| **Frequency_ridge_tracking** | `https://github.com/DavidBondesson/Frequency_ridge_tracking` | lightweight dynamic-programming ridge tracker on the superlet TFR (recommended default ridge extractor) |
| **MODA** (+ **PyMODA**) | `https://github.com/luphysics/MODA` · | principled ridge fallback — pull Iatsenko `ecurve`/path-optimization from `scripts/`, ignore the GUI in `allguis/` |


> If the MATLAB Wavelet Toolbox is licensed, its built-in `tfridge` + `wsst`/`fsst` are a
> zero-clone option for ridge extraction and synchrosqueezing.

---

## 4. Analyses: 

- time frequency decomposition at the single trial level using superlets method: 
  - Using epoched data from the bestMac channel: Demean, detrend. Keep: broadband-gamma (zero-phase FIR **15–75 Hz**, `ft_preproc_bandpassfilter`
    with `'firws'`) and its `hilbert()` (analytic signal → phase + envelope). **No in-band notch**
    (58 Hz ceiling already clears line noise). Utilize epochs that are padded by 1.5 seconds on either end of the target epochs to avoid edge artifacts during processing. Cut down to the target epochs for after time frequency extraction is complete. 
  - FASLT, **20–70 Hz** (pad beyond 25–58 for context), base cycles ~3, adaptive order ~3–30,
    frequency step ≤1 Hz. Cache the superlet kernel set across trials within a channel. 
  - Across trials, utilize the time period -1000ms to -500ms as the baseline period. Use the function myChanZscore to z-score the time frequency representations. These z-scored TF representations will be used for ridge extraction and plotting below.
- Apply ridge extraction to the per trial time frequency representations: 
  - Run the DP ridge tracker (Frequency_ridge_tracking; fallback MODA `ecurve`) on the FASLT TFR
    within **25–58 Hz**. Output primary **f̂(t)** and **ridge-power(t)**. 
  - Iterative **peeling**: remove the primary ridge's TF support in the TF representation. To do this: At each time point along the ridge, find the the full width half maximum (FWHM) extent of the ridge along the frequency axis. Find a gaussian that fits the observed power along the frequency axis at the FWHM points and the central peak. Subtract this gaussian from the TFR at that time point. Do this across all time points. re-extract a second ridge from the remaining TFR
  - Write a new field into the outDat struct called 'ridgeInfo'
   - Write a field in ridgeInfo called 'primaryRidge' and record info for the primary ridge: Record a trials X time matrix of f̂(t) values called "f", record a trials X time matrix of ridge-power(t) values called 'p'
   - Write a field in ridgeInfo called 'secondaryRidge' and record info for the secondary ridge: Record a trials X time matrix of f̂(t) values called "f", record a trials X time matrix of ridge-power(t) values called 'p'
   - Consider if there are other key parameters to be saved in ridgeInfo. Add anything you think would be useful here to help summarize. 
   - For trials skipped because of noise, fill the space in the matrices that would have been occupied by these trials with NaN values so that the matrices maintain the same number of trials as the behDat table. 
  - plot a series of single trial TF plots to be saved to the subject specific figs path found in outDat (okay to write directly to R drive even if working from data copied to E on lab machine). On these plots, overlay both the primary and secondary ridges as different colored lines.
  - plot all primary ridge traces on a single plot over the average TF representation with semi translucent lines to show the across trial pattern. Save this figure to the subject specific figs path found in outDat (okay to write directly to R drive even if working from data copied to E on lab machine)
  - plot all secondary ridge traces on a single plot over the average TF representation with semi translucent lines to show the across trial pattern. Save this figure to the subject specific figs path found in outDat (okay to write directly to R drive even if working from data copied to E on lab machine)
  - Do not add the full TFR matrices to the outDat struct as this will explode the file size. However, you will need the z-scored values for the next analysis, so maintain either in RAM or as temporary files on E:.
  

- phase/power continuity across time and frequency
 - Make a new field in the outDat struct called "powerPhaseContinuity". It will be a struct with key descriptive stats from this analysis
 - Take the smoothed (use: smoothdata(ridgeInfo.primaryRidge.p(trial, :), 'gaussian', 50)) ridge power time series for each trial. Find the time point of maximum ridge power in the primary ridge that occurs in the window between 0 ms and +700 ms relative to finalOnset. Call this the gammaPeakTime. The associated frequency of the ridge at this time point is the gammaPeakFrequency. The associated power of the ridge (unsmoothed) at this time point is the gammaPeakPower. Store gammaPeakTime, gammaPeakFrequency, and gammaPeakPower as trials X 1 vectors in powerPhaseContinuity leaving values as Nan for noise trials 
 - using the TFR information from the previous analysis, take the z-scored power time series for each trial at the trial specific gammaPeakFrequency, smooth it using the same 50 sample gaussian as above, looking backwards and forwards in time assess when did the gamma power at the gammaPeakFrequency rise above z=2 and when did it fall back below z=2. Call these time points relative to finalOnset the peakBurstOnset and peakBurstOffset. Also calculate the peakBurstLength as peakBurstOffset - peakBurstOnset. If the smoothed gammaPeakPower never rises above z=2, then the onset, offset, and length values are not calculable. Create a flag variable called hasPeak which is 1 when a trial has a z>2 peak and 0 when it does not. Store peakBurstOnset, peakBurstOffset, peakBurstLength, and hasPeak as trials X 1 vectors in powerPhaseContinuity leaving values as Nan for noise trials. peakBurstOnset, peakBurstOffset, and peakBurstLength will also be Nan were hasPeak==0
 - Using the primary ridge frequency at that time point, narrowband filter the data around the ridge frequency at the detected peak time point with a bandwidth of 5 Hz (+- 2.5): `ft_preproc_bandpassfilter` with `'firws'`) and its `hilbert()` (analytic signal → phase + envelope)
 - Also create a sin wave that is the same length as the trial and is aligned such that its phase equals the phase of the narrowband signal at the gammaPeakTime. 
 - Looking forwards and backwards in time, at what time point is the circular difference in phase between the sin wave and the narrowband signal greater than pi/4 and at what time point is it greater than pi/2. Call these time points peakPhaseOnsetNarrow, peakPhaseOffsetNarrow, peakPhaseOnsetWide, and peakPhaseOffsetWide. Record these values into the powerPhaseContinuity struct as vectors with Nan values for noise trials. These are still calculable for trials with hasPeak==0, so calculate for these anyway. To make this more robust, calculate the average phase difference across 5 consecutive time steps (10ms) and look for where the thresholds are crossed in the average. In addition, save the full phase difference vector from each trial into a matrix that is trials X time called peakPhaseConsistency in powerPhaseContinuity with Nans for noise trials. 
 - Next, construct a series of narrowband signals using the same method as above and spanning the frequency range covered by the primary ridge's f values across the trial in steps of 1 Hz. Look forward and backward from gammaPeakTime. Start with the instantaneous phase from the narrowband signal at gammaPeakTime for gammaPeakFrequency. For each time point moving forward and backward from the gammaPeakTime, take the primary ridge frequency at that time point and use that frequency to calculate the estimated phase progression for that time step. Compare this calculated phase to the observed phase at the ridge frequency at that time point. Carry forward your estimated phase without recentering to observed values. Continue to compare the calculated phase to the observed at each time step. Keep in mind, the observed phase for comparison will be determined by the ridge frequency, so different narrowband signals will be used as the ridge moves in frequency. Find the time point where the estimated phase and the observed phase are greater than pi/4 apart and the point where they are pi/2 apart. Again, take averages of the phase difference across 5 time steps (10ms) to make the analysis more robust. Record the points where the thresholds are crossed as ridgePhaseOnsetNarrow, ridgePhaseOffsetNarrow, ridgePhaseOnsetWide, ridgePhaseOffsetWide. Record these values into the powerPhaseContinuity struct as vectors with Nan values for noise trials. These are still calculable for trials with hasPeak==0, so calculate for these anyway. In addition save the full phase difference vector for each trial into a matrix that is trials X time called ridgePhaseConsistency in powerPhaseContinutity with Nans for noise trials. 
 - For each subject, plot on a single set of axes the trial level peakPhaseConsistency and ridgePhaseConsistency traces along with thicker lines for the mean across trials
 - For each subject, plot a histogram of peakBurstLength
 - save figures to the subject specific figs paths found in outDat (okay to write directly to R drive even if working from data copied to E on lab machine)

- save the outDat structs with the new fields you've calculated, overwriting the previous preproc files. Copy from E: to R: to overwrite and update files on R: as well. 

---


