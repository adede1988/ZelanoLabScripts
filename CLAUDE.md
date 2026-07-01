# CLAUDE.md — Respiration / EEG / iEEG preprocessing pipeline & preprocessed‑data reference

## What this project is

MATLAB preprocessing pipelines for intracranial‑EEG / scalp‑EEG + respiration (+ ECG)
data across **four tasks**: `breathingTask`, `cueTask`, `threshTask`, `O15`. Every
session is tracked in **`dataTracking.xlsx`**, which is the single source of truth for
session lists and per‑session parameters — no lists or params are hard‑coded in scripts.

This document is written for the **next stage of work**: editing the pipeline, adding a
new task to it, or — most importantly — **analysing the preprocessed `.mat` files**. The
pipeline itself (the big refactor + the "simplify / standardize" pass) is **complete and
validated**; the sections below describe what it produces and how to consume it.

> If you only want to *run* preprocessing, read `preprocessingReadme.md` (reference) and
> `tutorialPreprocessing.md` (step‑by‑step). This file is the architecture + **data‑structure**
> reference. `REFACTOR_NOTES.md`, `taskList.md`, `SimplifyStandardize*.md` are the historical
> change logs for how the current code came to be.

---

## 1. The four tasks

The same family of scripts preprocesses four experiment types. They differ only in raw
layout, photodiode/TTL structure, and behavioral file; the downstream signal processing is
shared.

| Task | sheet `Task` value(s) | `outDat.task` | What's unique |
|---|---|---|---|
| `breathingTask` | `breathingTasks`, `waveBreathing` | `breathingTask` (legacy files: `breathing`) | paced‑breathing blocks; **ECG/HRV**, **per‑breath metrics (`bmObj`)**, **paced target‑trace alignment**, **emotion ratings**. The richest task. |
| `cueTask` | `odorCueTask` | `cueTask` | odor cue / sniff / response TTLs; hit/miss/cr/fa behavior |
| `threshTask` | `Threshold` | `threshTask` | PEA intensity/pleasantness threshold; 45 single‑sniff trials |
| `O15` | `O15` | `O15` | loads genuinely raw data directly; photodiode TTLs parsed by `detect_ttls_O15` |

---

## 2. Data flow & file naming

```
RAW  (Neuralynx .mat + behavioral .mat/.csv)
  │   <task>_makeOutDat.m            (breathing / cue / thresh ONLY)
  │     • parse photodiode → TTLs, load raw behavior, stitch blocks
  ▼
<root>\<id>\preProc\<id>_<task>PreProc.mat        ← INTERMEDIATE  ("raw" outDat)
  │   <task>PreProc_main.m
  │     • applyParams(task,id) → P                 (per‑session params from the sheet)
  │     • assemble_outDat_all(S,P) → outDat,raw,TTL (load intermediate/raw → outDat)
  │     • shared pipeline + task‑specific onset/behavior
  ▼
<root>\<id>\preProc\<id>_<task>preproc.mat         ← FINAL  (the file you analyse)
```

**O15 has no `_makeOutDat`** — `O15PreProc_main.m` loads `raw_O15.mat` directly, and
`assemble_outDat_all`'s O15 branch runs `detect_ttls_O15` inside it.

### Exact filenames (load name → save name)

| Task | intermediate loaded by `assemble_outDat_all` | final saved by `_main` | top‑level var in `.mat` |
|---|---|---|---|
| breathing | `<id>_breathingPreProc.mat` | `<id>_breathingPreproc.mat` | **`chanDat`** (via `parSave`) |
| cue | `<id>_cueTaskPreProc.mat` | `<id>_cueTaskPreproc.mat` | `outDat` |
| thresh | `<id>_PEA_threshold_preproc.mat` | `<id>_PEA_threshold_preproc.mat` | `outDat` |
| O15 | `<id>\raw\raw_O15\raw_O15.mat` (`curDat`) | `<id>_O15preproc.mat` | `outDat` |

> **Windows / case‑insensitive FS quirk (intentional, do not "fix"):** for
> breathing/cue/thresh the intermediate (`…_PreProc.mat`) and the final (`…_preproc.mat`)
> are the **same file** — `_main` overwrites the intermediate in place. So a fully
> processed session has only the final on disk. Older breathing finals may store the
> variable as `out` or `chanDat`; loaders try `outDat → chanDat → out`.

**To load a final robustly:**
```matlab
s = load(finalPath);
fn = fieldnames(s);                 % {'outDat'} or {'chanDat'} or {'out'}
outDat = s.(fn{1});
```
A file is **fully processed** iff it has `outDat.moreThan1` (and, for breathing,
`outDat.bmObj` / `outDat.baseEmotion`). 

---

## 3. `dataTracking.xlsx` — the source of truth

Default path `R:\Neurology\Zelano_Lab\Lab_Common\Admin\dataTracking.xlsx`, `Sheet1`,
**header row = row 2, data from row 3**. Read through one function, **`applyParams.m`**:

```matlab
cfg = applyParams(task, 'makeOutDat'|'main')   % Mode A: session list for a loop
P   = applyParams(task, sessID)                % Mode B: one session's parameter struct
% task ∈ {'breathingTask','cueTask','threshTask','O15'}
```
A row is used only if its `Task` maps to one of the four pipelines **and** `Raw Data
Extracted` is non‑blank and not `INCOMPLETE`. Rows are deduped by `Subject ID`
(case‑insensitive, first wins), in sheet order. Live counts: ~55 breathing / 35 cue /
29 thresh / 38 O15 sessions.

### Parameter columns (read by `applyParams`, written by `writeParams`)

| Column | Tasks | Meaning / default |
|---|---|---|
| `datPre` | all | per‑session data root (blank → Type→root default; rebased from the canonical `R:\…\Lab_Common\` prefix onto this machine's `labCommon`) |
| `rspIDX` / `rspFlip` | all | which respiration channel (among labels containing `rsp`), and ±1 polarity. Defaults `1` / `1`. |
| `hasEEG` | all | run `preprocess_eeg` (default true; O15 false) |
| `spikeClean` | all | targeted‑ICA spike clean of macros (default true; O15 false) |
| `spikeThresh` / `spikeWin` | all | spike detector params (default `20` / `11`) |
| `macroRemove` | all | macro channels to drop before bipolar (`""`=none → `[]`) |
| `hasMacros` | breathing | run `preprocess_macros` (cue/thresh/O15 always run it) |
| `respThresh` / `cuedBackBuff` / `adjWin` | cue/thresh/O15 | sniff‑detection windows (default `500` / `150` / `500`) |
| `beatSpec` | breathing | ECG beat‑detection spec for `detectBeats` (default `1,0,gt,3`) |
| `ttlRemoveIdx` / `ttlNote` | O15 | aberrant photodiode TTL indices to drop / note |
| `isNewStd` | breathing/cue/thresh | selects the "new standard" ingestion branch (`newList`/`newSet`) |
| `paramSource` | target rows | `curated` (trusted, runs unattended) or `guess` (carried forward → forces interactive verification before any save) |
| `Data Preprocessed` | all | set to `X` by `writePreProcX` when a final is saved |

`P` (Mode B) always carries: `task, type ('Dupi'|'OBE'|'EEG'), fs_target=500, debug,
computeResp, rspIDX, rspFlip, hasEEG, spikeClean, spikeThresh, spikeWin, macroRemove,
paramSource` + task extras (breathing: `hasMacros, beatSpec, getBeats`; cue/thresh/O15:
`respThresh, cuedBackBuff, adjWin` + `ttlMap`; O15 also `pd`, `ttl`).

---

## 4. Running on any machine — `labPaths.m`

Every machine‑specific path comes from `labPaths.m`, which auto‑detects the machine by
`USERNAME` (cases: `adam` = home desktop, `dtf8829` = lab workstation) and derives
everything else (repo, eeglab, lab‑common roots, Admin sheet, behavioral dirs, figure dir,
breathing target‑trace dir). **New machine = add one `case`** (or a git‑ignored
`labPaths_local.m`); unknown machines error with a copy‑paste template. The four lab data
roots:

- Dupi → `…\Lab_Common\Dupi\`  · OBEControl → `…\Lab_Common\OBEControl\`
- EEGbreathing → `…\Lab_Common\AllStudyData\EEGbreathing\`  · figures → `…\Lab_Common\Adam\Dupi_processing\<id>\`

---

## 5. What the pipeline computes (stage by stage)

`assemble_outDat_all(S,P)` builds the initial `outDat` (see §6), then `_main` runs:

**Shared (identical across all four tasks — never edit per‑task):**

1. **`downsample_data(outDat, 500)`** — `resample` to **`fs_target = 500 Hz`**, then per
   channel: 4th‑order IIR **high‑pass 0.03 Hz**, 4th‑order IIR **low‑pass ≈Nyquist**, and
   **notch 60/120/180 Hz** (Q≈35). Sets `fs=500, origFS, downsampled=1`.
2. **`preprocess_eeg(outDat, EEGLOC, P)`** *(only if `P.hasEEG`)* — validates the **first
   32 channels** against `eegLocs_standard_coords.csv`, attaches coords, detects+interpolates
   noisy channels (`removeNoiseChansVolt`), removes blinks on the good channels (when >10
   survive, via `blinkRemoveWrapper`), and computes a **surface Laplacian**. Appends QC
   channels and EEG fields (see §6.3).
3. **`preprocess_macros(outDat, P)`** *(breathing only if `P.hasMacros`; others always)* —
   finds `macro` channels, optional `macroRemove`, **bipolar re‑references adjacent pairs →
   `macBP1..5`**; if `P.spikeClean`, splits the >10 Hz component and removes spike artifacts
   via targeted ICA; appends `macBP*` + `spikeCleanVec`.

**Respiration / onsets:**

- **cue / thresh / O15:** `preprocess_respiration_wholetrace` → `R` (chosen `rsp` channel,
  smoothed, Hilbert phase, a deflection metric `testSig`); `detect_sniffs_from_TTLs(R,P,outDat)`
  finds a sniff onset near each TTL (windows from `respThresh/cuedBackBuff/adjWin`);
  `refine_onsets_with_phase` snaps onsets to a respiration‑phase zero‑crossing.
- **breathing:** `process_respiration_breathing` derives the per‑breath **`bmObj`** matrix
  over the whole recording; `alignTargetBreathingTraceSimplify` aligns each paced block to
  its target trace and appends a `targTrace` channel; `processECG` band‑passes the `ECG`
  channels 5–40 Hz and z‑scores them (`buildECGz`), detects beats (`detectBeats` driven by
  `beatSpec`, min separation `fs/20`), physiologically filters the inter‑beat intervals, and
  appends an interpolated **`RRint`** (RR‑interval / HRV) channel; `flagBadBreaths` adds
  per‑breath QC.

**Behavior table:** `build_behavior_table_<task>` joins raw behavior to detected onsets
(§6.4). **Save** (§2) + `writeParams` + `writePreProcX`.

**`guess` rows** force interactive verification before any save: `paramCheck` (rsp + macro
choices), `paramCheckECG` (breathing ECG beats), and a deliberate **onset‑gate `error`** so
the user inspects figures, then promotes the row to `curated` and re‑runs.

---

## 6. THE PREPROCESSED DATA STRUCTURE  (what's in the final `.mat`)

The final variable (`outDat`, or `chanDat` for breathing) is a **struct**. The signal lives
in `outDat.data` (rows = channels, addressed by `outDat.labels`); everything else is
metadata/behavior. **Sampling rate is `fs = 500 Hz`** for all final signals.

### 6.1 Common fields (every task)

| Field | Type | Meaning |
|---|---|---|
| `data` | `double [nChan × nSamp]` | all channels; row *i* ↔ `labels{i}`. (Raw load is 2‑D ch×time; epoching happens downstream, not here.) |
| `labels` | `cell` of char/string | channel labels parallel to `data` rows. **Index channels by label string, never by fixed position** (except the EEG block, see 6.3). |
| `fs` | `500` | sample rate (Hz). Also `origFS`, `downsampled=1`. |
| `sessID` | char | e.g. `250623_Dupi_NMH_KS_2` (= raw session folder name) |
| `task` | char | `breathingTask` / `cueTask` / `threshTask` / `O15` (legacy breathing files: `breathing`) |
| `type` | char | `Dupi` / `OBE` / `EEG` (study/cohort) |
| `figs` | char | path to this session's figure folder |
| `rspIDX` | int | which `rsp`‑labelled channel is the respiration trace |
| `rspFlip` | ±1 | polarity to multiply the respiration trace by (so inhale ↑) |
| `TTL` | table or vector | task‑specific event markers in **fs=500 samples** (see 6.5) |
| `moreThan1` | 0/1 | "more than one sniff per trial" flag; **also the done‑sentinel** (1 = breathing/O15, 0 = cue/thresh) |

**O15‑only common extras:** `CSClist` (NCS channel labels), `OGdataDir`, `loadFile`
(the `*LoadData*.m` in the session dir), `preProcScript`.

### 6.2 Channel taxonomy in `outDat.data` (find rows by label)

| Channel(s) | Label match | Present when | Added by |
|---|---|---|---|
| 32‑ch scalp EEG montage | exact 10‑20 names, **rows 1–32** | `hasEEG` | raw load; validated in `preprocess_eeg` |
| macro (raw depth/strip) | `contains('macro')` | recorded | raw load |
| respiration | `contains('rsp')` | always | raw load (select with `rspIDX`, flip with `rspFlip`) |
| ECG (≈3 ch) | `contains('ECG')` | breathing | raw load (→ `buildECGz`) |
| photodiode / events | `contains('event')` | O15 (others vary) | raw load (drives `detect_ttls_O15`) |
| `blinkIndicator` | `=="blinkIndicator"` | `hasEEG` **and** >10 good EEG ch | `preprocess_eeg` |
| `badTS` | `=="badTS"` | `hasEEG` | `preprocess_eeg` (bad‑time‑sample mask) |
| `interpChan` | `=="interpChan"` | `hasEEG` | `preprocess_eeg` (interpolation mask) |
| `macBP1`…`macBP5` | `contains('macBP')` | macros run | `preprocess_macros` (bipolar pairs; count = nMacro−1) |
| `spikeCleanVec` | `=="spikeCleanVec"` | macros run | `preprocess_macros` (spike mixing vector; all‑ones if `spikeClean` off) |
| `targTrace` | `=="targTrace"` | breathing | `alignTargetBreathingTraceSimplify` (flattened paced target) |
| `RRint` | `contains('RRint')` | breathing | `processECG` (interpolated RR‑interval / HRV series, seconds) |

The 32 EEG montage order (validated, rows 1–32): `Fp1 Fz F3 F7 FT9 FC5 FC1 C3 T7 TP9 CP5
CP1 Pz P3 P7 O1 Oz O2 P4 P8 TP10 CP6 CP2 Cz C4 T8 FT10 FC6 FC2 F4 F8 Fp2`.

> The non‑EEG raw channel set (intracranial contacts, etc.) **varies by session/montage** —
> always locate channels by `cellfun(@(x) contains(x,'…'), outDat.labels)`

### 6.3 EEG‑derived fields (only if `hasEEG`)

| Field | Type | Meaning |
|---|---|---|
| `eegLocs` | table `[32×8]` | `labels, X, Y, Z, theta, phi, X_flat, Y_flat` (3‑D + flattened 2‑D coords) |
| `dataLap` | `[32×nSamp]` or `[]` | surface Laplacian (Perrin spline). `[]` if ≤10 good channels. |
| `dataLapFromInterp` | 0/1 | Laplacian was computed on a bad‑channel‑interpolated copy → carries no independent info at `badChans`; mask them for per‑channel Laplacian analysis. |
| `badChans` | cell of labels | channels flagged noisy (the broadband `data(1:32)` keeps their real, blink‑cleaned values) |
| `EEGInterpolation`, `EEGCleaning` | 1 | flags |
| `blinkRemoval` | 0/1 | 1 ⇒ blinks removed and `blinkIndicator` appended; 0 ⇒ ≤10 good channels, skipped |

### 6.4 `behDat` — the behavior table (task‑specific)

**`behDat` is a MATLAB `table`** in current finals. Access by **named column**, not numeric
index (older analysis code that does `behDat(:,13)` predates the table layout — see §8).

**cue / thresh / O15 — per‑sniff table.** Shared first 6 columns from `behDatFromSniffs`:

| Column | Meaning |
|---|---|
| `sniffOnset` | coarse sniff onset (sample @500 Hz) |
| `n` | trial number |
| `wiTriali` | sniff index within the trial |
| `TTLoffSet` | offset of the onset from its TTL (samples) |
| `sniffType` | integer sniff class |
| `sniffLabel` | readable label (cue/thresh: `cued`; O15: `start`/`free`/`confirm`) |

Then the **task‑specific** columns (raw behavior broadcast onto each sniff row):

- **cue:** `cue` (odor‑cue id 1–10), `odor` (presented‑odor id 1–10), `response` (1/2),
  `respString` (`"Yes"`/`"No"`/`"SKIP"`), `type` (signal‑detection outcome `hit`/`miss`/`cr`/`fa`)
- **thresh:** `odor` (1–3), `pleasantness`, `intensity` (rating‑slider values), `type` (`air`/`low`/`med`)
- **O15:** `target` (string odor name), `response` (string free identification), `expScore` (`0`/`0.5`/`1`)

Finally `refine_onsets_with_phase` appends the **last two** columns: **`adjust`** (signed
sample offset) and **`finalOnset`** (phase‑refined onset sample — *use this for epoching*).

> Verified against real Dupi files: cue `behDat` = `[40 sniffs × 13]`, thresh `[45 × 12]`,
> O15 `[81 × 11]`. O15 has `wiTriali` 1–8 and `sniffLabel` ∈ {`start`,`free`,`confirm`};
> cue/thresh are one sniff per trial (`wiTriali==1`, `sniffLabel=="cued"`, `moreThan1==0`).

**breathing — per‑breath table** (from `build_behavior_table_breathingTask`, derived from
`bmObj`; one row per detected breath):

| Column | Meaning |
|---|---|
| `sniffOnset`, `finalOnset` | inhale‑onset sample (equal for breathing), `finalOnset` is more accurate and should be used for alignment |
| `condition` | block/condition id (from `bmObj` col 12, tagged by `TaskBreaks`) |
| `Yonset`, `inhaleMax`, `Yend`, `exhaleMin` | respiration amplitudes (a.u.) at onset / inhale peak / breath end / exhale trough |
| `inMaxTim`, `endTim`, `exMinTim` | corresponding sample indices |
| `length` | breath duration (s) · `amp` amplitude · `index` breath index |
| `task`, `noseMouth`, `shadowFile`, `warp` | per‑block stimulus metadata (e.g. `task='audio'`, `noseMouth='nose'`) |
| `<question>_<category>` | dynamic: one column per emotion question × rating category, e.g. `calm_affective`, `tense_affective`, `emoAware_mindfulness`, `tht_come_go_mindfulness`; value = that block's rating (broadcast to every breath in the block) |
| `goodBreath` | 1/0 quality flag (from `flagBadBreaths`) |
| `maxRR`, `minRR`, `RR_max_min` | within‑breath HRV (s) |

> Empirically (a 32‑ch‑EEG Dupi breathing session) `behDat` is `[389 breaths × 33 vars]`;
> `baseEmotion` is a 1‑row table of the same `<question>_<category>` columns (the `order==0`
> baseline) plus `task`/`noseMouth`/`shadowFile`/`warp`. The emotion column set varies by
> the questions a session asked.

Breathing also stores **`baseEmotion`** (a 1‑row table: the baseline `order==0` emotion
ratings) and **`bmObj`** and **`heartBeats`** (see 6.6).

### 6.5 `TTL` conventions (samples @500 Hz)

- **O15** — `outDat.TTL` is a **table `[15 × 20]`**: `trialStart, buttonPress, confirmSniff,
  free1 … free17` (NaN‑padded). Built by `detect_ttls_O15` from the z‑scored `event`
  channel: pulse widths split `trialMarks` (<`trialSplitSamp`) from `sniffMarks`; expects
  **30** trialMarks (15 trials × start+button); `ttlRemoveIdx` drops aberrant ones.
- **thresh** — `outDat.TTL` is a **table `[45 × 3]`**: `start (=sniff−1000)`, `trial`,
  `sniff` (rebuilt in `_main` before sniff detection).
- **cue** — `outDat.TTL` is a **table `[nTrial × 3]`**: `trialStart, response, sniff`
  (one sniff per trial; e.g. `[40 × 3]`).
- **breathing** — `outDat.TTL` is a **vector** of block‑boundary samples
  (`round(intermediate.TTL/4)`, or a 5‑min fallback `0:600000:end`; e.g. `[0 1.5e5 … 7.5e5]`).

### 6.6 Breathing‑only structures

- **`bmObj`** `[nBreaths × 14]` per‑breath marker matrix:
  `1`=onset Y, `2`=onset time (s), `3`=inhale‑peak Y, `4`=peak time (s), `5`=end Y,
  `6`=end time (s), `7`=length (s), `8`=amplitude, `9`=peak idx, `10`=exhale‑peak Y,
  `11`=exhale‑peak time (s), `12`=**condition/block**, `14`=index. (col 13 unused.)
- **`heartBeats`** — ECG beat sample indices (from `detectBeats`/`beatSpec`).
- **`baseEmotion`** — baseline emotion ratings (1‑row table).

---

## 7. Using the preprocessed files downstream

**Get the respiration trace** (the canonical idiom, used everywhere):
```matlab
isRsp  = cellfun(@(x) contains(x,'rsp'), outDat.labels);
rsp    = outDat.data(isRsp,:);
rsp    = rsp(outDat.rspIDX,:) .* outDat.rspFlip;     % chosen channel, correct polarity
```
**Get analysis channels** (EEG + macro bipolar) and HRV/target (breathing):
```matlab
eeg32  = outDat.data(1:32,:);                         % only if hasEEG
isMac  = cellfun(@(x) contains(x,'macBP'), outDat.labels);
isRR   = cellfun(@(x) contains(x,'RRint'), outDat.labels);    % breathing
```


---

## 8. Known quirks & gotchas (read before analysing)

- **Index `behDat` by name, not number.** Current `behDat` is a **table**; some legacy
  analysis does `behDat(:,13)==0` assuming a numeric matrix
  with col 13 = bad‑breath flag. In the current table col 13 is `index`; the quality flag is
  **`behDat.goodBreath`**. Don't trust positional indices.
- **`task` value drift:** new files = `breathingTask`; some older finals = `breathing`.
  Handle both.
- **Top‑level var name varies:** `outDat` (cue/thresh/O15) vs `chanDat`/`out` (breathing).
  Load via `fieldnames` (§2).
- **EEG channels are exactly rows 1–32** and only when `hasEEG`; everything else is
  label‑addressed and montage‑dependent.
- **`spikeRemoval` is set to 1 in both `preprocess_macros` branches** — to know whether
  spikes were *actually* removed, inspect `spikeCleanVec` (all‑ones ⇒ none removed).
- **`spikeClean` sheet caveat (breathing):** for `KS_1, KS_2, AS(250908)` and the EEG/wave
  breathing rows the sheet records `false` but the legacy runtime used `true`; impact is nil
  where `hasMacros=false`. See `REFACTOR_NOTES.md`.
- **Filename case / naming quirks (left intentionally):** `_breathingPreProc` vs
  `_breathingPreproc`; the `250904_OBE_NWU_TI` (raw) vs sheet `…_TI_1`; `250811_Dupi_NMH_TPB_1`
  is remapped to `…_TB_1` in some downstream code and skips target‑trace alignment.
- **Memory (16 GB):** a single raw/continuous session is multi‑GB; process one session per
  MATLAB process if you hit OOM (the `_dev\run_*` harnesses isolate one session and clear big
  vars each iteration). The final `.mat` files are `-v7.3` (HDF5) — use `h5info`/`matfile` to
  inspect/partial‑load without pulling the whole `data` matrix into RAM.

---

## 9. Editing or extending the pipeline

### Hard constraints (do not violate)

- **`datPrei` index order is load‑bearing.** `*_makeOutDat.m` branches on
  `datPrei==1/2/3` to set `outDat.type`; order **must** stay `1=Dupi, 2=OBEControl,
  3=EEGbreathing`. `applyParams` returns `datPre` with those three first, extras appended.
- **"Unchanged after assemble."** In each `_main.m`, everything from the first
  `downsample_data(...)` onward is the shared pipeline — keep it byte‑for‑byte across tasks.
  O15 keeps `outDat.TTL = TTL;` right after assemble (so `assemble_outDat_all` **returns**
  `TTL`); cue/thresh/O15 use `raw` after assemble (so it **returns** `raw`).
- **Matching is case‑insensitive, whitespace‑trimmed.** Per‑session root comes from the
  sheet `datPre`, falling back to the Type→root map.
- **Prefer strict over flexible:** if the data isn't the expected shape, let it error — a
  loud failure surfaces a real data problem.


### Validating pipeline changes

Legacy oracles are kept in‑repo until retired: `getSessionParams_*.m`,
`assemble_outDat_*.m` (and the per‑task `assembleRaw_*`/`preproc\assemble_outDat_*`), plus
the `_dev\` test/harness scripts (`test_detectBeats`, `test_beatSpec_map`,
`test_integration`, `run_*`). Re‑run the relevant parity check after any change and confirm
numerically identical output where a reference exists.

---

## 10. File map (quick index)

| Kind | Files |
|---|---|
| **Config / loaders** | `labPaths.m`, `applyParams.m`, `writeParams.m`, `writePreProcX.m`, `assemble_outDat_all.m`, `preprocessAll.m` |
| **Shared pipeline** | `assembleOutDat.m`, `loadIntermediateRaw.m`, `downsample_data.m`, `preprocess_eeg.m`, `preprocess_macros.m`, `preprocess_respiration_wholetrace.m`, `detect_sniffs_from_TTLs.m`, `refine_onsets_with_phase.m`, `behDatFromSniffs.m`, `paramCheck.m`, `plot_sniff_epochs.m`, `parSave.m`, `interpolate_perrinX.m`, `laplacian_perrinX.m`, `removeNoiseChansVolt.m`, `blinkRemoveWrapper.m`, `ica_flag_spikes_targeted.m`, `detect_spikes.m` |
| **Task‑specific** | `<task>_makeOutDat.m` (`preproc\threshPreProc_makeOutDat.m`), `assembleRaw_<task>.m`, `build_behavior_table_<task>.m`; O15: `assembleRaw_O15.m`, `detect_ttls_O15.m`, `assembleOutDat_O15extras.m`; breathing: `assembleRaw_breathingTask.m`, `process_respiration_breathing.m`, `alignTargetBreathingTraceSimplify.m`, `processECG.m`/`buildECGz.m`/`paramCheckECG.m`, `detectBeats.m`, `flagBadBreaths.m`, `plotBreathLengths.m`; cue: `assembleRaw_cueTask.m`, `outMat_to_table.m`; thresh: `assembleRaw_threshTask.m` |
| **Deliverable scripts** | `breathingTaskPreProc_main.m`, `cueTaskPreProc_main.m`, `threshPreProc_main.m`, `O15PreProc_main.m` + the 3 `_makeOutDat` |
| **Docs** | `preprocessingReadme.md` (reference), `tutorialPreprocessing.md` (how‑to), `REFACTOR_NOTES.md` / `taskList.md` / `SimplifyStandardize*.md` (history) |
| **Dev / oracles** | `_dev\` (tests, `run_*` harnesses, `subfunction_catalogue.json`), `getSessionParams_*.m`, `assemble_outDat_*.m` (legacy, kept as oracles) |

Adding a participant is a **sheet edit** (a row + its parameter columns), no code change.


---

## Remote compute (lab desktop)

When a task needs more cores/RAM than the local machine or involves long-running batch work, offload to a remote Windows lab desktop over SSH. (A Quest HPC cluster + Globus transfers are also available, set up separately — see the end.)

### Prerequisite — VPN
All remote access requires Northwestern **GlobalProtect VPN** connected (`vpn-connect.northwestern.edu`; NetID + Duo). The local machine is off-domain, so without the VPN every remote connection times out. The user connects the VPN (Duo is interactive). If remote commands suddenly start timing out, assume the tunnel dropped and ask the user to reconnect — jobs already running on the remote keep going.

### Running commands — `ssh labdesktop`
- `ssh labdesktop "<command>"` runs non-interactively and returns real stdout/stderr + exit code. Key-based auth (configured in `~/.ssh/config`); no password prompt.
- Hardware: 8-core/8-thread Intel i7-10700, ~96 GB RAM, **no usable GPU**. Good for CPU/RAM-bound parallel work; not GPU compute.
- **Default remote shell is cmd.exe.** Plain cmd tools (`net use`, `robocopy`, etc.) run directly. For non-trivial **PowerShell**, base64-encode the script:
  ```powershell
  $script = @'
  $ProgressPreference = "SilentlyContinue"   # avoids CLIXML stderr noise
  # ... your PowerShell here ...
  '@
  $enc = [Convert]::ToBase64String([Text.Encoding]::Unicode.GetBytes($script))
  ssh labdesktop "powershell -NoProfile -NonInteractive -EncodedCommand $enc"
  ```

### Storage — where to read/write
- **`E:\` is the workspace.** ~880 GB free, read/write from SSH. Scratch, intermediates, large outputs go here.
- **`C:\` has only ~4 GB free — never use it for job I/O.**
- **Lab server (`R:` → `\\fsmresfiles.fsm.northwestern.edu\fsmresfiles`):** not reachable by default from an SSH session (key-based logons carry no network credentials). **`cmdkey`/Windows Credential Manager does NOT work over key-based SSH** — *tested:* a domain-share credential won't store in a passwordless logon and the no-credential mount fails, so don't retry it. The verified pattern is **explicit credentials supplied inline, in a single `ssh` invocation** (the mapping only lives in its own session):
  - **Work directly on R:** (existing `R:\...` code runs unchanged):
    ```
    ssh labdesktop "net use R: /delete /y 2>nul & net use R: \\fsmresfiles.fsm.northwestern.edu\fsmresfiles /user:fsm\dtf8829 <NetID pw> & <job that uses R:\...> & net use R: /delete /y"
    ```
  - **Heavy/repeated I/O:** robocopy a working set to local `E:` first, crunch there, write back — same single-invocation rule.
  - **Credential** = `fsm\dtf8829` + **NetID/SSO password** (sensitive, rotates). Never type it in chat or write it here. Stash it **once** as a DPAPI blob on the client (`Read-Host -AsSecureString | ConvertFrom-SecureString | Set-Content $env:USERPROFILE\.fsmcreds\netid.sec`); at runtime decrypt it **locally** and interpolate via base64 `-EncodedCommand` so the literal never appears in chat or a tracked file.
- `E:` is BitLocker-encrypted but normally already unlocked. If a check shows it **Locked**, ask the user for the E: BitLocker password (a local secret — never write it here) and unlock with `Unlock-BitLocker`.
- Clean temporary files from `E:` at the end of all jobs. Transfer updated preproc files back to `R:` and overwrite the old versions in place. 
- Keep data files and overall mirrored file structure on `E:` for future jobs. That is, keep preproc data files on `E:`


### Monitoring long jobs
Redirect output to a file on `E:` and poll it, or check process state with `ssh labdesktop "..."`. The SSH session is a separate logon from any RDP session, so mapped drives and some profile state differ.

### Code: 
- On the home local machine, code is in C:\Users\Adam\Documents\GitHub\
- On the labdesktop, code is in G:\My Drive\GitHub, but migration is currently ongoing to E:\GitHub\. For any repo in both G: and E:, prefer the copy in E:, edit in E:, utilize git pointing at E:. Only utilize G: code if no copy of the repo exists on E:. 
- Utilize git to move code between machines

### Related (separate channels) Quest is not ready for Claude use yet. This functionality is coming soon.
- **Quest HPC** — `ssh dtf8829@login.quest.northwestern.edu`, Slurm batch jobs; same VPN prerequisite, key auth set up the same way.
- **Globus** — bulk `R:` ↔ HPC transfers via the Globus service (independent of the login nodes).

> **Do not commit secrets.** This file is tracked in git. Never paste the NetID password, BitLocker password, or any credential here — supply them at runtime and keep them in a password manager or untracked location.