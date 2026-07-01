# threshTask analysis — task description / build plan

**Goal.** Produce a full `threshTask` report in the exact style of the `cueAnalysis/`
pipeline (FOOOF → bestMac → noise rejection → z-score spectrograms → group means →
time-resolved gamma → knitted HTML). The **only scientific change** from cue: instead of
splitting the event-locked spectrograms by **locking** (`trialStart` vs `finalOnset`), split
by **odor condition** — the `finalOnset`-locked response of trials where
`outDat.behDat.type == 'low'`, `== 'med'`, and `== 'air'`. Everything else (methods,
noise rule, FOOOF settings, z-scoring, plotting, grouping, gamma-epoch analysis, report
layout) is ported **unchanged in spirit** from cue.

This document is the executable spec. Build the code under `threshAnalysis/` mirroring
`cueAnalysis/`, run it on the lab desktop, copy results back to `R:`, knit.

---

## 0. Reference: what makes thresh different from cue

Read `CLAUDE.md` §6 (the preprocessed data structure) first. The thresh-specific facts that
drive this port:

| Item | cue | thresh |
|---|---|---|
| Final file | `<id>_cueTaskPreproc.mat` | **`<id>_PEA_threshold_preproc.mat`** |
| Top-level var | `outDat` | `outDat` |
| `applyParams` list call | `applyParams('cueTask','main')` | **`applyParams('threshTask','main')`** |
| `TTL` table | `[nTrial × 3]`: `trialStart, response, sniff` | **`[45 × 3]`: `start (=sniff−1000), trial, sniff`** |
| Per-trial anchor for noise/baseline | `TTL.trialStart` | **`TTL.start`** (the pre-sniff marker; = sniff − 1000 samples = 2 s pre-sniff @500 Hz) |
| behDat rows | 40 sniffs, 1/trial (`sniffLabel=="cued"`) | **45 sniffs, 1/trial** (`sniffLabel=="cued"`, `moreThan1==0`) |
| behDat condition column | `type` = hit/miss/cr/fa (SDT) | **`type` ∈ {`air`,`low`,`med`}** (odor intensity) — *this is the split variable* |
| Split variable for spectrograms | locking (trialStart / finalOnset) | **`behDat.type` (air / low / med), all finalOnset-locked** |
| Trials per split | ~all trials in each locking | ~15 trials each (45 / 3), minus noise |
| Analysis channels | `macBP1..5` | `macBP1..5` (thresh always runs `preprocess_macros`) |
| Groups | DupiS1/S2/S3/Control | **same** (Dupi sessions 1/2/3 by trailing `_N`; Control = all non-Dupi) |
| Session count | ~35 | **~29** (live count; confirm with the session table) |

**The single rename that touches several files:** every cue reference to `od.TTL.trialStart`
becomes `od.TTL.start` for thresh (noise screen anchor + shared baseline anchor). There is no
separate "response" locking in thresh; the sniff/`finalOnset` is the only response event, and
it is what we split three ways.

**Baseline (important).** Keep the cue mechanism byte-for-byte: each condition's `finalOnset`
epochs are z-scored against **those same trials' pre-`start` baseline** (the −700..−200 ms
window relative to `TTL.start`), tacked on per-trial before `myChanZscore` (exactly as cue
z-scores `finalOnset` against the `trialStart` baseline). So all three condition maps use the
identical baseline definition; they differ only in *which trials* feed them. Because
`start = sniff − 2 s`, the −700..−200 ms pre-`start` window sits ~2.2–2.7 s before the sniff —
a clean pre-stimulus baseline.

---

## 1. File-by-file port plan (`cueAnalysis/` → `threshAnalysis/`)

Copy each cue file, rename `cue_`→`thresh_` / `run_cue_`→`run_thresh_`, and apply the edits in
the right column. Files marked **identical** need only the rename + the `cue`→`thresh` string
swaps in paths/CSV names. Keep the heavy commenting style.

| cue file | → thresh file | edits required |
|---|---|---|
| `cue_init_paths.m` | `thresh_init_paths.m` | identical (EEGLAB + FieldTrip via `L.fieldtrip`; keep the **`PLUGINLIST` guard** so `eeglab nogui` runs once per session — see §6 Java-IO gotcha). |
| `cue_session_table.m` | `thresh_session_table.m` | swap `applyParams('cueTask',…)`→`'threshTask'`; final filename → `<id>_PEA_threshold_preproc.mat`; fuzzy fallback glob → `[id '*PEA_threshold*.mat']`. Group logic unchanged. |
| `cue_noise_trials.m` | `thresh_noise_trials.m` | **identical** (relative K robust-z sharp-deflection rule). Start K=10; re-calibrate (below). |
| `cue_calibrate_noise.m` | `thresh_calibrate_noise.m` | identical; anchor epochs on `TTL.start`. Re-run to pick K for thresh amplitudes; write `_calib_trialMaxZ.csv`; target >0% and ≤20% trial loss dataset-wide. |
| `cue_fooof_macBP.m` | `thresh_fooof_macBP.m` | identical FOOOF (2–58 Hz, knee, gamma 30–58). **Noise screen anchor:** change `ev = od.TTL.trialStart` → `ev = od.TTL.start`. bestMac selection logic (periodicPeak / flattenedFallback / allNoisy_leastBad, 30% rejRate exclusion) unchanged. |
| `run_cue_tasks123.m` | `run_thresh_tasks123.m` | identical driver (FOOOF all sessions, rewrite CSV). CSV → `threshTask_fooof_summary.csv`; writes `bestMac` back into each final; per-subject `<id>_macBP_fooof_periodic.png`. |
| `run_cue_fooof_one.m` | `run_thresh_fooof_one.m` | identical single/subset driver (append/replace rows). |
| `cue_plot_fooof.m` | `thresh_plot_fooof.m` | identical. |
| `cue_plot_singletrial.m` | `thresh_plot_singletrial.m` | identical. |
| `plot_session_singletrials.m` | `plot_session_singletrials_thresh.m` | identical (per-trial stacked raw, bestMac red where zd>K), anchor `TTL.start`. |
| `run_cue_singletrial.m` | `run_thresh_singletrial.m` | identical. |
| `cue_ztfr_pair.m` | **`thresh_ztfr_conds.m`** | **REWRITE — the core change.** See §2. |
| `run_cue_ztfr.m` | `run_thresh_ztfr.m` | rewrite the locking loop as a condition loop; see §2. |
| `cue_plot_ztfr.m` | `thresh_plot_ztfr.m` | identical (map + respiration overlay). Just called 3× (low/med/air). |
| `run_cue_task4_group.m` | `run_thresh_task4_group.m` | `lockings = {'trialStart','finalOnset'}` → `conds = {'low','med','air'}`; group PNGs `group_<G>_TFR_<cond>.png`; means `.mat`. 4 groups × 3 conds = 12 PNGs. |
| `run_cue_gamma_epochs.m` | `run_thresh_gamma_epochs.m` | split the finalOnset gamma progression **by condition** (see §3); add a `cond` column to `threshTask_gammaEpochs.csv`; per-subject `gammaTimeProgression_<cond>.png` (or one overlaid figure). |
| `cue_make_manifest.m` | `thresh_make_manifest.m` | derive `threshTask_data_manifest.csv` from the thresh FOOOF CSV (status ok / noMacBP / excluded(allNoise)). |
| `run_cueAnalysis_all.m` | `run_threshAnalysis_all.m` | master wrapper, same stage order; env override `THRESH_GROUPDIR`; knit `threshTask_report.Rmd`. |
| `cueTask_report.Rmd` | `threshTask_report.Rmd` | mirror sections; spectrogram + gamma sections become condition-split (see §4). |

New CSV names (all under `groupDir`): `threshTask_fooof_summary.csv`,
`threshTask_gammaEpochs.csv`, `threshTask_data_manifest.csv`. Per-subject figures under
`<figRoot>/<id>/threshTask/`. Group PNGs + CSVs under `groupDir` (shared
`…/Dupi_processing/groupStatFigs/` is safe: thresh CSVs are prefixed and thresh group PNGs use
`{low,med,air}` suffixes that never collide with cue's `{trialStart,finalOnset}`).

Add a per-job env override **`THRESH_GROUPDIR`** (mirror `CUE_GROUPDIR`) in every driver:
`groupDir = getenv('THRESH_GROUPDIR'); if isempty(groupDir), groupDir = fullfile(L.figPath,'groupStatFigs'); end`.

---

## 2. The core change — `thresh_ztfr_conds.m` + `run_thresh_ztfr.m`

**`thresh_ztfr_conds(sig, rsp, fs, startEv, foByCond)`** replaces `cue_ztfr_pair`.

- `startEv` : per-trial `TTL.start` sample indices (the baseline anchor), clean trials only.
- `foByCond`: struct `.low .med .air`, each a vector of `finalOnset` sample indices for that
  condition's **clean** trials, **paired** with the matching entries of `startEv` (so each
  condition trial can find its own pre-`start` baseline frames).

Method (identical engine to cue, three times):
1. For each condition c ∈ {low, med, air}: epoch that condition's `finalOnset` events and the
   same trials' `start` events with `newtimef` (log freqs 2–120 Hz, Morlet `[3 0.8]`,
   `baseline NaN` → raw power `|alltfX|²`).
2. Baseline = −700..−200 ms of the trials' `start` epochs. Tack the per-trial baseline frames
   onto the front of each `finalOnset` trial and `myChanZscore` per frequency (reuse the
   `zpipe` subfunction verbatim). Average trials → the condition's z-map.
3. Respiration: mean of that condition's `finalOnset`-locked `rsp` epochs (overlay).

Return `out.freqs`, `out.low/.med/.air` = `struct(map,times,resp,respT,nFinal)` (empty if
<3 clean trials), and `out.qc` (per-condition nTrials, nOOB, nFinal).

**`run_thresh_ztfr.m`** (mirror `run_cue_ztfr` structure):
1. Load final, get `bestMac` (skip if none), get rsp (`rspIDX`,`rspFlip`), `figDir = <figRoot>/<id>/threshTask`.
2. Noise on **all `TTL.start` epochs** (`thresh_noise_trials`), draw `singleTrialRawMac.png`,
   write `nTrials`/`nNoiseTrials` onto the bestMac CSV row immediately (survive a later skip).
3. Build paired vectors from `behDat`: for each row j, `k = bd.n(j)`; anchor `s = TTL.start(k)`,
   `f = bd.finalOnset(j)`, `cond = char(bd.type(j))`; drop if noisy(k) or non-finite. Bucket
   `f`,`s` into `foByCond.(cond)` / `startByCond.(cond)`.
4. `R = thresh_ztfr_conds(sig, rsp, fs, startEv, foByCond)` (pass paired start per condition).
5. Per-condition color scale: `a = prctile(abs([low;med;air] maps), 98)` shared across the
   three so they're comparable; `clim=[-a a]`. Plot 3 PNGs
   `<id>_<bestMac>_TFR_{low,med,air}.png` via `thresh_plot_ztfr`.
6. Save `<id>_thresh_bestMac_TFR.mat` with `tfrOut.low/.med/.air`, `freqs`, `clim`, identity,
   `nNoiseTrials`. Write `nFinal_low/med/air` (+ `nLow/nMed/nAir` clean-trial counts) to the CSV.

CSV QC columns (replace the cue locking columns): `nTrials, nNoiseTrials, nLow, nMed, nAir,
nFinal_low, nFinal_med, nFinal_air`.

---

## 3. Time-resolved gamma — `run_thresh_gamma_epochs.m`

Cue runs one finalOnset-locked gamma progression over all clean trials. Because thresh splits
finalOnset three ways, **run the progression per condition**:

- For each condition, take that condition's clean `finalOnset` trials; Morlet power 2–58 Hz
  (100 linear freqs); average into ten 250 ms windows (−500..+2000 ms); FOOOF each window
  (2–58 Hz, knee); record the largest 30–58 Hz peak per (condition × epoch).
- Output long CSV `threshTask_gammaEpochs.csv` with an added **`cond`** column:
  `sessID,subID,sessNum,type,group,cond,epoch,centerMs,gammaPeakFreq,gammaPeakPower,gammaDetected`.
- Per-subject figure: either three `gammaTimeProgression_<cond>.png` or one figure with the
  three conditions as line-style families. Keep the ochre→purple epoch ramp.
- **Low-N caveat:** ~12–15 trials/condition after noise. If a condition has <~8 clean trials,
  skip its progression for that session and log it (don't silently emit a noisy FOOOF). The
  Rmd's per-epoch "participants with a peak" count then naturally reflects available sessions.
- Report LME (§6) gains a condition factor: `gammaPeakFreq ~ epoch * cond + (1|sessID)` per
  group (or a simpler `~ epoch + cond` if the interaction is singular), testing whether gamma
  peak frequency/its epoch-trajectory differs by odor intensity.

*(Alternative if trial counts prove too thin: pool all trials for the gamma progression
—odor-agnostic, exactly like cue— and keep the condition split only in the spectrograms. Decide
after seeing the calibrated clean-trial counts; document whichever is chosen in the Rmd.)*

---

## 4. Report — `threshTask_report.Rmd`

Mirror the cue Rmd sections; change only the split dimension and helper globs:

- **Params** `groupDir`, `figRoot` (same defaults, R: paths). `conds <- c('low','med','air')`.
- Helpers: `subFig(sid,bestMac,cond) -> file.path(fr, sid, 'threshTask', sprintf('%s_%s_TFR_%s.png', sid, bestMac, cond))`; `groupFig(g,cond) -> group_<g>_TFR_<cond>.png`.
- **§1 Methods** — copy cue methods; replace the "two lockings share the trialStart baseline"
  paragraph with "three odor conditions (air/low/med), all finalOnset-locked, each z-scored
  against its own trials' pre-`start` baseline"; note `TTL.start` as the baseline anchor and the
  odor-intensity design (air = clean-air control, low/med = PEA intensities). Keep the FOOOF
  (2–58 Hz, knee), gamma band, and relative-noise-rule text verbatim.
- **§2 Data available** — from `threshTask_data_manifest.csv`.
- **§3 FOOOF findings** — `threshTask_fooof_summary.csv`; dynamic column select incl. `rejRate`,
  `selectionMethod`, `isBestMac` (reuse cue's `intersect(...)` idiom).
- **§3b Noise rejection** — relative rule table + example single-trial plots (2 high / 2 low noise).
- **§4 Single-subject example spectrograms** — one example per group, **3 panels** (low/med/air).
- **§5 Group-average spectrograms** — 4 groups × 3 conditions (12 panels), shared ±5 z scale.
- **§6 Time-resolved gamma** — example progressions, per-epoch peak counts, channel-centered
  peak frequency, and the per-group LME **with the condition factor** (§3).
- Knit self-contained: `output_options=list(self_contained=TRUE)`, `RSTUDIO_PANDOC` set.

---

## 5. Staging & compute workflow (E: ↔ R:) — REQUIRED

Runs on the lab desktop (`ssh labdesktop`; **GlobalProtect VPN must be connected**). See
`CLAUDE.md` "Remote compute" and the `lab-remote-compute-workflow` memory for the infra; this is
the thresh-specific recipe. **`E:` is the workspace; never do job I/O on `C:` (~4 GB free).**

### 5a. Stage inputs R: → E: (no credential; scp from home over the existing key-based SSH)
Home already has R: mounted at `/r/` (Git Bash). A WMI-detached lab job has **no** R: mount, so
finals must live on local `E:`. For each thresh session, copy the final into an E: mirror of the
lab-common tree:
```bash
# from HOME (R: visible as /r/). <root> ∈ {Dupi, OBEControl, AllStudyData/EEGbreathing}
scp "/r/Neurology/Zelano_Lab/Lab_Common/<root>/<id>/preProc/<id>_PEA_threshold_preproc.mat" \
    "labdesktop:E:/Lab_Common/<root>/<id>/preProc/<id>_PEA_threshold_preproc.mat"
```
(scp is ~4.2 MB/s over the VPN; a full thresh set is a one-time stage. Only re-stage sessions
whose finals changed.) Put a git-ignored **`E:\labPaths_local.m`** on the lab (labCommon=`E:\Lab_Common\`,
fieldtrip=`E:\fieldtrip-20260518`) and `addpath('E:\')` first in every job so `labPaths` resolves
to E: in the detached session. FieldTrip-with-FOOOF must exist at `E:\fieldtrip-20260518`
(release zip, has `external/brainstorm/process_fooof`).

### 5b. Run detached (survives VPN/harness drops)
Write a job script `E:\run_thresh_job.m` that: `addpath('E:\')`; adds the repo
(`G:\My Drive\GitHub\ZelanoLabScripts` + `threshAnalysis`); `setenv('THRESH_GROUPDIR','E:\Lab_Common\Adam\Dupi_processing\groupStatFigs')`;
`diary` to `E:\_thresh_job.log`; runs the stages; writes a sentinel `E:\_thresh_job_DONE.txt`
in a `try/catch`. Launch it **WMI-detached** from a LOCAL PowerShell (base64 `-EncodedCommand`,
nothing leaks):
```powershell
$script = @'
$matlab = "C:\Program Files\MATLAB\R2024b\bin\matlab.exe"
$r = ([wmiclass]"Win32_Process").Create('"' + $matlab + '" -sd "E:\" -batch "run(''E:/run_thresh_job.m'')"')
"PID $($r.ProcessId)"
'@
$enc=[Convert]::ToBase64String([Text.Encoding]::Unicode.GetBytes($script))
ssh labdesktop "powershell -NoProfile -NonInteractive -EncodedCommand $enc"
```
Poll the sentinel + tail the log with `ssh -n labdesktop "..."`. **Do not** rely on a harness
background job — the IDE harness tears those down; the detached MATLAB survives.

### 5c. ⚠️ WMI-detached MATLAB Java-IO gotcha (hit on 2026-06-30 — design around it)
In `-batch` detached mode stdout is a broken pipe with a finite kernel buffer. After ~1–2
sessions of output, `fprintf` throws `IOException: Error writing to output stream` (diary uses
C-stdio and is unaffected). Mitigations, in order:
1. **Keep the `PLUGINLIST` guard** in `thresh_init_paths` so `eeglab nogui`'s internet check runs
   **once** (repeated calls accelerate the buffer fill). Already in cue's `cue_init_paths`.
2. Redirect Java at job top (delays failure ~1 session):
   `System.setOut(java.io.PrintStream(java.io.FileOutputStream('E:\_javaout.txt')))`.
3. **Split the batch into one-session-per-MATLAB invocations** (most robust): loop the session
   list in the *launcher*, spawning a fresh `-batch` MATLAB per session for the FieldTrip stages
   (FOOOF/ztfr/gamma). Each process prints little enough to stay under the buffer.
4. Computations finish before the `fprintf` dies — **verify by output files** (CSV rows, PNGs,
   TFR mats), then resume from the next stage. `run_thresh_ztfr`/`_gamma` persist per-session, so
   they're resumable; the group step is cheap to re-run.
5. Avoid re-running the FieldTrip PSD path (`ft_freqanalysis`/`ft_checkdata`) inside a detached
   job when the FOOOF result already exists — patch the CSV/mat directly (as done for the cue
   bestMac fix).

### 5d. Copy results back E: → R: (scp from home; no credential)
After the sentinel and an output check:
```bash
# report + group CSVs/PNGs
scp     "labdesktop:E:/Lab_Common/Adam/Dupi_processing/threshTask_report.html" \
        "/r/Neurology/Zelano_Lab/Lab_Common/Adam/Dupi_processing/threshTask_report.html"
scp -r  "labdesktop:E:/Lab_Common/Adam/Dupi_processing/groupStatFigs/." \
        "/r/Neurology/Zelano_Lab/Lab_Common/Adam/Dupi_processing/groupStatFigs/"
# per-subject threshTask figures + TFR mats (per changed <id>)
scp -r  "labdesktop:E:/Lab_Common/Adam/Dupi_processing/<id>/threshTask/." \
        "/r/Neurology/Zelano_Lab/Lab_Common/Adam/Dupi_processing/<id>/threshTask/"
# finals with bestMac written back (only if FOOOF re-ran; large)
scp     "labdesktop:E:/Lab_Common/<root>/<id>/preProc/<id>_PEA_threshold_preproc.mat" \
        "/r/Neurology/Zelano_Lab/Lab_Common/<root>/<id>/preProc/<id>_PEA_threshold_preproc.mat"
```
Knit either on the lab (self-contained, E: params) or on **home** with R: params
(`params=list(groupDir='R:/…/groupStatFigs', figRoot='R:/…/Dupi_processing')`). Home has R:
mounted, so a home knit reads the just-copied figures directly.

*(In-place alternative: if working an interactive lab session instead of detached, you can run
against R: directly using the single-invocation `net use R: … & <job> & net use R: /delete`
credential pattern from `CLAUDE.md`. The detached scp-staging path above is preferred because it
needs no credential.)*

---

## 6. Validation / acceptance criteria

- `thresh_session_table(false)` lists ~29 sessions on disk with sane group counts; each has
  `macBP` (combo/no-macBP sessions excluded, as in cue).
- `threshTask_fooof_summary.csv` has one row per macBP channel, exactly one `isBestMac==1` per
  session, with `rejRate`, `selectionMethod`, and per-condition noise/final columns populated.
- For a spot-check session, the three `_TFR_{low,med,air}.png` exist, share a color scale, and
  each shows the correct `nFinal` in its title; the sum of `nFinal_{low,med,air}` ≈ clean trials.
- Group PNGs: 12 files `group_{DupiS1,DupiS2,DupiS3,Control}_TFR_{low,med,air}.png`.
- `threshTask_gammaEpochs.csv` has the `cond` column; per-subject gamma progression figure(s)
  render; the §6 LME runs per group with the condition factor.
- `threshTask_report.html` knits self-contained (~few MB) with all sections populated.
- Re-run parity: running one session via `run_thresh_fooof_one({id})` then
  `run_thresh_ztfr({id})` → `run_thresh_gamma_epochs({id})` → `run_thresh_task4_group([])`
  reproduces that session's rows/figures without disturbing others.

## 7. Execution checklist

1. Build `threshAnalysis/` files per §1 (start by copying cue files + rename + string swaps;
   then rewrite `thresh_ztfr_conds.m`, `run_thresh_ztfr.m`, `run_thresh_gamma_epochs.m`,
   `run_thresh_task4_group.m`, and `threshTask_report.Rmd` per §2–§4).
2. `git add threshAnalysis/ && git commit && git push` from home.
3. On lab: `git pull` (ignore Google-Drive `.git` corrupt-loose-object warnings).
4. Stage thresh finals R: → E: (§5a). Ensure `E:\labPaths_local.m` + `E:\fieldtrip-20260518`.
5. Calibrate K (`thresh_calibrate_noise`) → set K in `thresh_noise_trials`.
6. Run detached job (§5b), one-session-per-MATLAB for the FieldTrip stages (§5c#3): FOOOF →
   ztfr → gamma → group means + manifest → knit. Poll sentinel; verify output files.
7. Copy results E: → R: (§5d). Knit from home with R: params if not already self-contained.
8. Update the `cue-macbp-gamma-spectrogram-analysis` / `lab-remote-compute-workflow` memories
   with any new thresh-specific gotchas.
