# Preprocessing — how it works

Reference for the respiration/EEG preprocessing pipelines in this repo. For a
step-by-step "I've never used this" walkthrough, see **tutorialPreprocessing.md**.

---

## 1. The four tasks

The same family of scripts preprocesses four experiment types. Each has its own
raw layout, photodiode/TTL structure, and behavioral file, but they share almost
all downstream signal processing.

| Task           | sheet `Task` value(s)            | what's special |
|----------------|----------------------------------|----------------|
| `breathingTask`| `breathingTasks`, `waveBreathing`| paced-breathing blocks; **ECG/HRV** and **target-trace alignment** (unique to this task) |
| `cueTask`      | `odorCueTask`                    | odor cue/sniff/response TTLs, hit/miss/cr/fa behavior |
| `threshTask`   | `Threshold`                      | PEA intensity/pleasantness threshold; 45 single-sniff trials |
| `O15`          | `O15`                            | loads genuinely raw data directly; photodiode TTLs parsed by `detect_ttls_O15` |

---

## 2. Data flow

```
RAW (Neuralynx .mat + behavioral file)
   │   <task>_makeOutDat.m      (breathing / cue / thresh only)
   │     • parse photodiode -> TTLs, load behavior, stitch blocks
   ▼
<root>/<id>/preProc/<id>_<task>PreProc.mat        ("intermediate")
   │   <task>PreProc_main.m
   │     • applyParams -> P            (per-session params from the sheet)
   │     • assembleRaw_<task>(S) + assembleOutDat(raw,S,P)  (each main calls its own loader)
   │     • shared pipeline + task-specific behavior/onset detection
   ▼
<root>/<id>/preProc/<id>_<task>preproc.mat        ("final")
```

**O15 has no `_makeOutDat`** — `O15PreProc_main.m` loads `raw_O15.mat` directly
and `assemble_outDat_all` runs `detect_ttls_O15` inside its O15 branch.

> Windows note: on a case-insensitive filesystem the intermediate
> (`..._PreProc.mat`) and the final (`..._preproc.mat`) are the **same file** for
> breathing/cue/thresh, so the main overwrites the intermediate in place. This is
> intentional and known.

---

## 3. Single source of truth — `dataTracking.xlsx`

Session lists and every per-session parameter live in
`R:\Neurology\Zelano_Lab\Lab_Common\Admin\dataTracking.xlsx`. **No session lists
or parameters are hard-coded in the scripts.** Everything is read through one
function, `applyParams.m`:

```matlab
cfg = applyParams(task, 'makeOutDat'|'main')   % Mode A: session list for a loop
P   = applyParams(task, sessID)                % Mode B: one session's parameters
```

A row is used only if its `Task` maps to one of the four pipelines **and** `Raw
Data Extracted` is non-blank and not `INCOMPLETE`.

### Parameter columns (read by `applyParams`, written by `writeParams`)

| Column | Tasks | Meaning |
|---|---|---|
| `datPre` | all | per-session data root (else the Type→root default) |
| `rspIDX` / `rspFlip` | all | which respiration channel, and ±1 polarity |
| `hasEEG` | all | run `preprocess_eeg` (default true; O15 false) |
| `spikeClean` | all | ICA spike-clean macros (default true; O15 false) |
| `spikeThresh` / `spikeWin` | all | spike detector params |
| `macroRemove` | all | macro channels to drop (`""`=none) |
| `hasMacros` | breathing | run `preprocess_macros` |
| `respThresh` / `cuedBackBuff` / `adjWin` | cue/thresh/O15 | sniff-detection windows |
| `beatSpec` | breathing | ECG beat-detection spec for `detectBeats` |
| `ttlRemoveIdx` / `ttlNote` | O15 | aberrant TTL indices to drop |
| `isNewStd` | breathing/cue/thresh | selects the "new standard" ingestion branch |
| `paramSource` | target rows | `curated` (trusted) or `guess` (review before trusting) |
| `Data Preprocessed` | all | set to `X` by `writePreProcX` when a final is saved |

---

## 4. Running on different machines — `labPaths.m`

Every machine-specific path comes from `labPaths.m`. It auto-detects the machine
(by Windows `USERNAME`, with `COMPUTERNAME` as a tiebreaker) and returns one
struct of base paths plus everything derived from them (repo, eeglab, lab-common
data roots, the Admin sheet, behavioral-results dirs, figure dir, target traces).

To run on a **new machine**: add one `case` to the switch in `labPaths.m` (or drop
in a git-ignored `labPaths_local.m`). Nothing else needs editing — unknown
machines error with a copy-pasteable template. The scripts self-bootstrap the repo
onto the path, so they run regardless of the current folder.

---

## 5. Parameter verification — `curated` vs `guess`

New rows are often pre-filled by carrying a previous same-study session forward;
those are marked `paramSource = guess`. The pipeline forces a human to verify a
guess before any result is saved:

1. **Input check** — `paramCheck(outDat, P)` plots the respiration channels (accept
   the `rspIDX`/`rspFlip`, or pick another) and the macro channels (choose any to
   remove and whether to spike-clean).
2. **ECG check (breathing only)** — `paramCheckECG(outDat, P)` shows a short ECG
   segment with the beats the current `beatSpec` detects; accept, or type a new
   spec and re-check. Uses the same `buildECGz` as `processECG`, so what you accept
   is what runs.
3. **Onset gate** — after detection, the main deliberately **errors** on a guess
   session ("inspect the figures…") so nothing is saved until you've reviewed the
   onset/beat figures in the session's figure folder.
4. **Promote + persist** — once happy, set the row's `paramSource` to `curated`
   (or call `writeParams(P, S.id)`), then re-run. The curated run skips the checks,
   `writeParams` persists the verified params, the final is saved, and
   `writePreProcX` marks `Data Preprocessed = X`.

Curated sessions run start-to-finish with no prompts.

---

## 6. What the preprocessing actually does

### Shared pipeline (identical across all four tasks — do not edit per task)

1. **`downsample_data`** — resample to `fs_target` (500 Hz), high-pass at 0.03 Hz,
   low-pass just under Nyquist, and notch at 60/120/180 Hz (line noise + harmonics).
2. **`preprocess_eeg`** (if `hasEEG`) — validate the 32 EEG electrode labels against
   `eegLocs_standard_coords.csv`, attach standard 3-D/flat coordinates, detect and
   interpolate noisy channels, remove blinks on the good channels (when >10 survive),
   and compute a surface Laplacian. EEG-specific; appends EEG-derived channels.
3. **`preprocess_macros`** — find the `macro` channels, bipolar re-reference adjacent
   pairs into `macBP1..5`; if `spikeClean`, split off the high-frequency component
   and remove spike artifacts via targeted ICA; append a `spikeCleanVec` mask.

### Respiration / onsets

- **cue / thresh / O15**: `preprocess_respiration_wholetrace` picks the chosen `rsp`
  channel, smooths it, and computes a Hilbert respiration phase; then
  `detect_sniffs_from_TTLs` finds the sniff onset near each TTL and
  `refine_onsets_with_phase` snaps onsets to the respiration phase.
- **breathing**: `process_respiration_breathing` derives per-breath metrics (the
  `bmObj` matrix) over the whole recording; `alignTargetBreathingTraceSimplify`
  aligns each block to its paced target trace; `processECG` band-passes/z-scores the
  ECG (`buildECGz`), detects beats (`detectBeats` driven by `beatSpec`), and writes an
  interpolated RR-interval/HRV channel; `flagBadBreaths` flags low-quality breaths.

### Behavior table

- **`build_behavior_table_<task>`** joins the raw behavior with the detected onsets.
  The shared first six columns (sniff onset, trial, within-trial index, TTL offset,
  sniff type + label) come from `behDatFromSniffs`; each task adds its own columns
  (O15 target/response/expScore; cue cue/odor/response/type; thresh
  odor/pleasantness/intensity/type). breathing builds a per-breath table instead.

### Save

The final `outDat` is saved to `<root>/<id>/preProc/<id>_<task>preproc.mat`,
`writeParams` persists any verified params, and `writePreProcX` marks the row.

---

## 7. Running a task

Open MATLAB with this repo reachable and run the deliverable script for the task
(curated sessions run unattended; guess sessions prompt as in §5):

```matlab
% intermediate from raw (breathing / cue / thresh only):
breathingTask_makeOutDat            % or cueTask_makeOutDat, or preproc/threshPreProc_makeOutDat
% final:
breathingTaskPreProc_main           % or cueTaskPreProc_main / threshPreProc_main / O15PreProc_main
```

Each `_main` loops over all of its task's sessions and skips ones already finished.
Memory note (16 GB machines): a single raw session can be multi-GB — process one
session per MATLAB process if you hit "Out of memory" (the `_dev/run_*` harnesses
isolate one session and clear big vars each iteration).

---

## 8. File map

| Kind | Files |
|---|---|
| **Config / loaders (shared)** | `labPaths.m`, `applyParams.m`, `writeParams.m`, `writePreProcX.m`, `assemble_outDat_all.m` (backward-compat one-call wrapper for the `_dev` harnesses — the mains don't use it) |
| **Shared pipeline (do not edit per task)** | `assembleOutDat.m`, `resolveFigDir.m`, `loadIntermediateRaw.m`, `downsample_data.m`, `preprocess_eeg.m`, `preprocess_macros.m`, `preprocess_respiration_wholetrace.m`, `detect_sniffs_from_TTLs.m`, `refine_onsets_with_phase.m`, `plot_sniff_epochs.m`, `paramCheck.m`, `behDatFromSniffs.m`, `parSave.m` |
| **Task-specific (rewrite for a new task)** | `assembleRaw_<task>.m`, `<task>_makeOutDat.m`, `build_behavior_table_<task>.m`; O15: `assembleRaw_O15.m`, `detect_ttls_O15.m`, `assembleOutDat_O15extras.m`; breathing: `process_respiration_breathing.m`, `alignTargetBreathingTraceSimplify.m`, `processECG.m`/`buildECGz.m`/`paramCheckECG.m`, `flagBadBreaths.m`, `plotBreathLengths.m` |
| **Deliverable scripts** | `breathingTaskPreProc_main.m`, `cueTaskPreProc_main.m`, `threshPreProc_main.m`, `O15PreProc_main.m` and the three `_makeOutDat` scripts |

Adding a participant is a sheet edit (a row + its parameter columns), no code
changes. Adding a whole new task means a new `<task>_makeOutDat.m`, an
`assembleRaw_<task>.m` (load its raw into a `raw` struct — reuse
`loadIntermediateRaw` if it loads a `<task>PreProc.mat`) that the new `_main.m`
calls directly, a `build_behavior_table_<task>.m`, and a `_main.m` modeled on an
existing one — see the tutorial. The shared assembler (`assembleOutDat`) and
downstream pipeline are untouched.
