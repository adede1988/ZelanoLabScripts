# chirpAnalysis.md — OB sniff-locked gamma: one chirping oscillator, or two?

**Analysis-stage spec for Claude Code.** Sits alongside `CLAUDE.md` (the pipeline +
data-structure reference). This document defines a new downstream analysis that *consumes*
the preprocessed finals
It does not touch the preprocessing pipeline.

Read chirpAnalysis_2.md to refamiliarize with the prior analysis step.

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


- Regenerate the TFR plots, when saving limit the z axis to [-10 +20] on both single trial and mean plots
- Create group level mean TFR plots by averaging the subject level matrices together. Save these to "...\Lab_Common\Adam\Dupi_processing\groupStatFigs"

- Refine phase/power continuity plots:
  -  for each subject, plot on a single set of axes (not two separate panels, but literally a single set of axes) both the peakPhaseConsistency and ridgePhaseConsistency traces from individual trials along with thicker lines indicating the across trial means. Use different colors for peak versus ridge based consistency. Center the time axis at 0 = gammaPeakTime and plot for +- 1000 ms in both directions. Overwrite previous version of these plots in place
  - Create means across subjects for these plots. Plot subject level means as individual lighter lines and group level means as thicker lines. Make one plot per group. Save these to "...\Lab_Common\Adam\Dupi_processing\groupStatFigs"


- Ridge power continuity
 - In the previous analysis we calculated peakBurstOnset, peakBurstOffset, and peakBurstLength. These were all calculated at the gammaPeakFrequency. I'd like to get analogous values for ridgeBurstOnset, ridgeBurstOffset, and ridgeBurstLength. These values need to be added into the powerPhaseContinuity struct. To obtain them, take the primaryRidge power timeseries p, smooth it (use: smoothdata(ridgeInfo.primaryRidge.p(trial, :), 'gaussian', 50). Look forward and backward along the smoothed primaryRidge power timeseries starting from gammaPeakTime. Find when the ridge values fall below z = 2. These threshold crossings will determine ridgeBurstOnset and ridgeBurstOffset which you can use to calcualte ridgeBurstLength. 
 - overwrite the burstLength histogram with a new figure that displays both the ridgeBurstLength distribution and the peakBurstLength distribution
 - Make a group level version of this plot as a histogram of all burstLength and ridgeBurstLength values for each group. Save these to "...\Lab_Common\Adam\Dupi_processing\groupStatFigs"

- Perform peak aligned analysis
 - Create new epochs that span +- 1000 ms around the gammaPeakTime on each trial
 - execute the same TFR and primary ridge extraction that you performed in the previous task set
 - make a new field inside ridgeInfo called "peakLockedRidge", record the same information into it as the primary/secondary ridge information. 
 - Make a new average across trials time frequency plot that displays the frequency axis with 0 centered at gammaPeakFrequency (even if that's different for different trials) and then +- 10 Hz from there and time centered at gammaPeakTime with +- 1000ms from there. Save this to the subject figs folder. You may need to increase the headroom of the initial TFR decomposition so that you have the space to pull +- 10 Hz around the gammaPeakFrequency for every trial.
 - Make averaged plots for each group displaying the equivalent time frequency plot averaged over subjects. These plots should be constructed by averaging the subject level averages together. Save these to "...\Lab_Common\Adam\Dupi_processing\groupStatFigs"
 - remember to avoid using noise trials throughout. 

--- 

## 5. Output report: 

- Edit the chirpAnalysisReport.rmd document. This will be a wholesale rewrite. I want this to reflect the version 2 analysis that we're doing now. I don't want any of the analytics or methods from the old analysis. The only thing to keep is the data being used. Everything else likely needs to go. 
- Show single trial examples of the TFR plots with ridge overlays, get at least 5 trials each from BW and AD as well as 20 additional trials from various subjects. 
- Show the group level TFR plots 
- show 2 subject level examples (OBE should be BW and AD) from each group of the peakPhaseConsistency/ridgePhaseConsistency plots. Also show the group level of these plots. Include stats comparing peakPhaseConsistency versus ridgePhaseConsistency at the group level across time with permutation across time for multiple comparison correction. Do stats in R, save the per subject values needed for the stats from Matlab in a .csv to pass between. Do stats for each group independently
- Show examples from the same subjects for the peakBurstLength v. ridgeBurstLength histograms. Also show the group level plots. Again, do stats on these comparing the peak v. ridge values within each group. Again, do stats in R on values passed from Matlab via .csv file. 
- Show examples from the same subjects for the peak aligned TFR plots. Also show the group averaged plots. 
- knit the report when done  

- save the outDat structs with the new fields you've calculated, overwriting the previous preproc files. Copy from E: to R: to overwrite and update files on R: as well. Ensure that the repo is aligned across local and lab machines with all code changes committed and pushed on both machines and the .html report is included in the repo.

---


