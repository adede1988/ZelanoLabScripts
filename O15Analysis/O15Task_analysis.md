# O15Task analysis — design record (macBP gamma + sniff-locked spectrograms)

Mirrors the cue-task analysis (`cueAnalysis/`, `cueTask_analysis0-4.md`) for the **O15 odor-
identification task**, consuming the preprocessed `<id>_O15preproc.mat` finals. Entry point:
`run_O15Analysis_all(sessFilter, doKnit)`. Report: `O15Task_report.Rmd` → `O15Task_report.html`.

## What is different from the cue task (and why)

1. **Events are sniffs, split into three types.** O15 has no cue-onset TTL to lock to — every
   marker is a sniff. Each trial = a **start** sniff (odor onset), several **free** sniffs
   (exploration/identification), and a **confirm** sniff (button press). We split
   `behDat.sniffLabel ∈ {start, free, confirm}` (= `sniffType` 1/2/3; counts ≈ 15 / ~120 / 15)
   and lock each type on `behDat.finalOnset`. (`TTL` table `[15×20]` is informational only.)

2. **A single shared baseline before the first trial-start sniff.** Instead of cue's per-trial
   trialStart baseline, all three sniff types are z-scored against ONE window taken before the
   first start sniff: `[firstStart − 1 s − ≤30 s … firstStart − 1 s]`. One `newtimef` pass over
   that segment gives the baseline TF power; wavelet edge frames (±1500 ms) are trimmed → the
   baseline frame distribution per frequency (~180 frames).

3. **Plain z (comparable across sniff types).** `o15_ztfr_multi` reuses `myChanZscore` but passes
   the trial-mean power as a **single column** (`ntrials = 1`), so it returns a PLAIN z =
   `(power − baselineMean) / baselineSD` (baseline-SD units), not cue's bootstrap SEM-z. This is
   intentional: a SEM-z (∝ √n) would make the 120-trial *free* map ~2.8× hotter than the
   15-trial *start* map purely from trial count. One shared baseline + plain-z put the three
   maps on a single comparable scale. (Verified: start & free share clim ≈ 0.59 a.u.)

## Stages (each is also an independent driver)

| Stage | Driver | Output |
|---|---|---|
| 1 FOOOF / bestMac | `run_o15_fooof_all` / `run_o15_fooof_one` (reuse `cue_fooof_macBP`) | `bestMac` into each `.mat`; `O15Task_fooof_summary.csv`; `<id>_macBP_fooof_periodic.png` |
| 2 Baseline-z spectrograms | `run_o15_ztfr` (+ `o15_ztfr_multi`; relative noise via `cue_noise_trials` with `K=o15_noise_K()`; reuse `cue_plot_ztfr`/`cue_plot_singletrial`) | 3× `<id>_<bestMac>_O15TFR_<start\|free\|confirm>.png`, `singleTrialRawMac_start.png`, `<id>_O15_bestMac_TFR.mat`, noise/QC cols |
| — Noise K calibration | `o15_calibrate_noise` (sweep K over O15 sniffs) → set `o15_noise_K` | `_calib_sniffMaxZ_O15.csv` + console K-sweep table |
| 3 Time-resolved gamma | `run_o15_gamma_epochs` (computed **per sniff type**: start/free/confirm) | `O15Task_gammaEpochs.csv` (col `locking`), `gammaTimeProgression_<type>.png` |
| 4 Group means + manifest | `run_o15_task_group` + `o15_make_manifest` | `group_<G>_O15TFR_<locking>.png`, `O15Task_group_means.mat`, `O15Task_data_manifest.csv` |
| 5 Report | `O15Task_report.Rmd` | `O15Task_report.html` |

## Conventions

- **Fresh finals only.** `o15_session_table` flags `fresh` = final modified on/after 2026-06-24
  (this preprocessing run); drivers use `T(T.onDisk & T.fresh,:)`. Older finals are ignored.
- **Groups** (same as cue): `DupiS1/S2/S3` = Dupi & sessNum 1/2/3; `Control` = all non-Dupi;
  Dupi sessNum>3 → `DupiSx` (excluded from the 4-group means).
- **Sessions without `macBP` are excluded** (O15 EEG/macros are sometimes absent).
- **Noise (relative rule):** ported from the cue analysis via `cue_noise_trials` — over the whole
  recording, robust-z the sliding 10 ms max-min (`zd = (d−median)/(1.4826·MAD)`); a sniff (epoch
  −1.75…+5.75 s) is noisy if any 10 ms window has `zd > K`. `K` is channel-relative so it transfers
  across datasets; set in `o15_noise_K` (default 10), tuned for O15 by `o15_calibrate_noise`. Noisy
  sniffs are excluded from the averages. (Replaces the earlier absolute >80 µV / 10 ms rule, which
  penalized uniformly high-amplitude channels.)
- **Where:** per-subject figs `…\Dupi_processing\<id>\O15\`; group CSVs/figs/means
  `…\Dupi_processing\groupStatFigs\` (all O15 artifacts prefixed `O15Task_*` / `*_O15TFR_*` so
  they never clobber the cue outputs).

## Reused unchanged from `cueAnalysis/`
`cue_fooof_macBP`, `cue_noise_trials`, `cue_plot_ztfr`, `cue_plot_singletrial`, `cue_plot_fooof`
(+ repo `myChanZscore`, `labPaths`, `applyParams`).
