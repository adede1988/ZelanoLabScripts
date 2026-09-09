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

- **Now expanding to:**  `cueTask`, `threshTask`, and `O15`. **Eventually:** all four tasks (`breathingTask`, `cueTask`,
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
    the selected sniff(s) a config field (default to the task-relevant sniff). Analyze data split by sniffLabel. 
- **O15 long-window guard:** truncate each epoch at the **next** `finalOnset` so a long window
  doesn't run into the following sniff; set a `burstTruncated` flag where the window is cut
  short, so a truncated coexistence window cannot masquerade as a genuinely brief one.
- **cueTask** Split analysis by trial type `behDat.type`. Specifically, analyze 'hit' and 'cr' trials separately.
- **threshTask** There are three types of trials: "air", "med", and "low". Analyze data split between these types
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


- This will be a lot of work to apply the work done in prior analysis rounds to more data and more separate categories of trials. 

- Specifically, generate all plots at the subject and group level for cueTask, O15, and threshTask data. 
- For the cueTask, treat outDat.behDat.type == 'hit' and outDat.behDat.type =='cr' as two separate categories of trials to be averaged separately for all measures. Disgard other trial types.
- For the O15 task, treat outDat.behDat.sniffLabel == 'start', outDat.behDat.sniffLabel == 'free', and outDat.behDat.sniffLabel == 'confirm' as three different categories of trial types to be averaged separately for all measures. 
- For the threshTask, treat outDat.behDat.type == 'air', outDat.behDat.type == 'low', and outDat.behDat.type == 'med' as three different categories of trial types to be averaged separtely for all measures. 

- Run all analyses that have been done previously but generate new subject level averages for different trial types and apply the analysis pipeline to the two new task sets. Update figure names to reflect task and trial type in the file names.


--- 

## 5. Output report: 

- Edit the chirpAnalysisReport.rmd document. Integrate the addition of different trial types and of the new tasks. 
- Carefully verify that statistical numeric outputs match values shown in figures. 
- knit the report when done  

- save the outDat structs with the new fields you've calculated, overwriting the previous preproc files. Copy from E: to R: to overwrite and update files on R: as well. Ensure that the repo is aligned across local and lab machines with all code changes committed and pushed on both machines and the .html report is included in the repo.

---


