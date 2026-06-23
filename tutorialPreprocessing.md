# Tutorial — preprocessing from scratch

For someone who has never used this codebase. It walks from "a new participant
was recorded" to "preprocessed data on disk." For the architecture/reference, see
`preprocessingReadme.md`.

The golden rule: **`dataTracking.xlsx` is the single source of truth.** You add a
participant by editing a spreadsheet row, not by editing code.

Contents
1. [Set up your machine](#1-set-up-your-machine-once)
2. [Add a participant to dataTracking.xlsx](#2-add-a-participant-to-datatrackingxlsx)
3. [Extract the raw data (placeholder)](#3-extract-the-raw-data-placeholder)
4. [Run preprocessing](#4-run-preprocessing)
5. [Use preprocessAll](#5-use-preprocessall)
6. [Write a pipeline for a brand-new task](#6-write-a-pipeline-for-a-brand-new-task)

---

## 1. Set up your machine (once)

Every machine-specific path lives in `labPaths.m`. If your machine isn't known, the
first call errors with a template. Add one `case`:

```matlab
case 'yourusername'        % from getenv('USERNAME')
    L.codePre   = 'C:\path\to\GitHub\';    % holds ZelanoLabScripts, slowBreathing, closed-loop-respiration
    L.eeglab    = 'C:\path\to\eeglab2026.0.0';
    L.labCommon = 'R:\Neurology\Zelano_Lab\Lab_Common\';   % wherever the lab share is mapped
    L.gdrive    = 'G:\My Drive\';           % '' if Google Drive isn't here (breathing target traces)
```

Prefer not to commit your paths? Create a git-ignored `labPaths_local.m` returning
just those four fields. Verify with:

```matlab
L = labPaths()          % should print sensible paths; the data roots should exist
```

---

## 2. Add a participant to dataTracking.xlsx

Open `R:\Neurology\Zelano_Lab\Lab_Common\Admin\dataTracking.xlsx`, `Sheet1`
(header row = row 2, data from row 3). Add **one row per task** the participant did.

### Identity + tracking columns
- `Subject ID` — the session folder name, e.g. `260420_Dupi_NMH_XY_1`.
- `Type` — study (`Dupi`, `OBE…`, `EEG…`); sets the default data root.
- `Task` — one of `breathingTasks`/`waveBreathing`, `odorCueTask`, `Threshold`, `O15`.
- `Data On Server`, `Raw Data Extracted` — leave `Raw Data Extracted` blank/`INCOMPLETE`
  until the raw `.mat` exists; the pipelines ignore rows that aren't extracted.
- `Data Preprocessed` — leave blank; `writePreProcX` sets it to `X` when a final is saved.

### Parameter columns and recommended starting guesses
Fill these for the row's task. When you're unsure, use the guess and set
`paramSource = guess` so the pipeline forces you to verify (see §4).

| Column | Tasks | Recommended guess |
|---|---|---|
| `datPre` | all | the row's data root (blank = use the Type default) |
| `rspIDX` | all | `1` (which respiration channel) — verified interactively |
| `rspFlip` | all | `1` (or `-1` if the trace is upside down) — verified interactively |
| `hasEEG` | all | `TRUE` (O15: `FALSE` unless EEG was recorded) |
| `spikeClean` | all | `TRUE` (O15: `FALSE`) |
| `spikeThresh` / `spikeWin` | all | `20` / `11` |
| `macroRemove` | all | `""` (none) — verified interactively |
| `hasMacros` | breathing | `TRUE` if macro electrodes were recorded |
| `respThresh` / `cuedBackBuff` / `adjWin` | cue/thresh/O15 | `500` / `150` / `500` |
| `beatSpec` | breathing | carry the most recent same-study value, else `1,0,gt,3` |
| `ttlRemoveIdx` / `ttlNote` | O15 | `""` unless a session has aberrant TTLs |
| `isNewStd` | breathing/cue/thresh | `TRUE` for recent recordings (new ingestion branch) |
| `paramSource` | all | **`guess`** for a new row → forces interactive verification |

The fastest way to fill a new row well: copy the most recent **curated** row of the
**same study and task** and adjust. That's exactly what a `guess` represents — values
carried forward, to be reviewed before you trust them.

---

## 3. Extract the raw data (placeholder)

> **This step is project-specific and is intentionally left blank here.** Between
> the Neuralynx/recording output and this pipeline there is a per-study "load data"
> step that produces the raw `.mat` (`curDat` with `rawData`, `outLabs`,
> `ncslabels`) in the session's `raw\` folder. Document your lab's loadData
> procedure in this section. Once the raw `.mat` exists, mark `Raw Data Extracted`
> in the sheet and continue to §4.

---

## 4. Run preprocessing

Open MATLAB with the repo reachable. Each task has (for breathing/cue/thresh) a
`makeOutDat` (raw → intermediate) and a `main` (intermediate → final). O15 has only
a main. The scripts loop over all of the task's sessions and skip finished ones.

```matlab
% breathing
breathingTask_makeOutDat
breathingTaskPreProc_main

% cue
cueTask_makeOutDat
cueTaskPreProc_main

% thresh
run('preproc/threshPreProc_makeOutDat.m')
threshPreProc_main

% O15 (no makeOutDat)
O15PreProc_main
```

### What happens for a `guess` session (interactive verification)
A curated session runs unattended. A `guess` session pauses for you to verify, in
this order:

1. **Respiration channel** (`paramCheck`) — a figure overlays all respiration
   channels with your chosen `rspIDX`/`rspFlip` highlighted. Enter `1` to accept, or
   `0` to step through channels and pick/flip the right one.
2. **Macros / spike cleaning** (`paramCheck`) — a figure shows the macro channels.
   Enter channels to remove as `[# #]` (or `[]`), and `0/1` for spike cleaning.
3. **ECG beat detection — breathing only** (`paramCheckECG`) — a figure shows a
   short ECG segment with the beats your `beatSpec` detects (magenta dashes). Enter
   `1` to accept, or `0` to type a new `beatSpec` (e.g. `1,0,gt,3 & 2,0,gt,4`) and
   re-check until the beats land on the R-peaks.
4. **Onset gate** — after onset/sniff detection the script **deliberately errors**
   ("inspect the figures…"). Open the session's figure folder
   (`…\Adam\Dupi_processing\<id>\`) and confirm the sniff/breath/ECG figures look right.
5. **Commit** — when satisfied, set that row's `paramSource` to `curated` in the
   sheet (or run `writeParams(P, S.id)` to persist what you just verified), then
   re-run. The curated run skips the prompts, saves the final, and marks
   `Data Preprocessed = X`.

### Per-task notes
- **breathing** — richest task; uses `rspIDX/rspFlip`, macro params, and `beatSpec`
  (ECG). Verifies respiration, macros, and ECG beats. Also does target-trace
  alignment and per-breath QC.
- **cue** — uses `rspIDX/rspFlip`, macro params, and `respThresh/cuedBackBuff/adjWin`
  for sniff detection. Verifies respiration + macros.
- **thresh** — like cue, plus it rebuilds a 45-trial sniff-TTL table before onset
  detection.
- **O15** — `hasEEG`/`spikeClean` default FALSE; uses `respThresh/cuedBackBuff/adjWin`
  and, when a session has aberrant photodiode pulses, `ttlRemoveIdx`. TTLs are parsed
  by `detect_ttls_O15` (a TTLs.jpg is written to the figure folder to check).

Memory (16 GB): a raw session can be multi-GB. If you hit "Out of memory", run one
session per MATLAB process, or use the `_dev/run_*` harnesses (they isolate one
session and clear big variables each iteration).

---

## 5. Use preprocessAll

`preprocessAll.m` finds every raw-extracted session not yet marked preprocessed and
can run the right pipelines for them.

```matlab
preprocessAll                                  % REPORT ONLY: lists pending per task,
                                               % flagging [GUESS] rows
setenv('PREPROCESS_RUN','1'); preprocessAll    % actually run the pending tasks
```

Curate any `[GUESS]` rows interactively first (§4) — they would otherwise halt a
batch run. The report-only mode is safe to run anytime to see the backlog.

---

## 6. Write a pipeline for a brand-new task

The shared pipeline (downsample → EEG → macros → respiration/onsets → behavior
table → save) is reused unchanged. You only write the task-specific pieces. Use an
existing task as a template (cue is the simplest end-to-end).

1. **Sheet** — teach `applyParams.m` your task: add it to `taskKey`, `canonTask`,
   and `taskCallerKey`, and add any new parameter columns + their defaults in the
   Mode-B switch. Add the parameter columns to `dataTracking.xlsx`.
2. **Ingestion** — write `<task>_makeOutDat.m`: parse the photodiode/TTLs and
   behavior for your task and save `<id>_<task>PreProc.mat` with at least
   `.data .labels .fs .behDat .TTL`. The whole file is task-specific; copy an
   existing makeOutDat's header (labPaths + the `applyParams('<task>','makeOutDat')`
   cfg block) verbatim.
3. **Assembly** — add a `case '<task>'` to `assemble_outDat_all.m` that loads your
   intermediate (or raw) into `raw`/`outDat`. Everything after the switch is shared.
4. **Behavior table** — write `build_behavior_table_<task>.m`. If your task is
   sniff-based, start from `behDatFromSniffs(sniffs, sniffTypes)` for the shared
   first six columns and add only your task's broadcast loop.
5. **Main** — copy a `_main.m`, change the task name in the `applyParams` calls and
   the load/save filenames, and swap the one task-specific block (onset detection +
   `build_behavior_table_<task>`). Keep the guess flow (`paramCheck`, the onset
   `error` gate, `writeParams`, `writePreProcX`) — it's shared.
6. **Mark it** — update the main's TASK-SHARED/TASK-SPECIFIC legend so the next
   person can see at a glance what to edit.

Conventions to keep:
- Read every session list/parameter through `applyParams`; never hard-code lists.
- Read every path through `labPaths`; never hard-code a machine path.
- Prefer simple, strict code: if the data isn't the expected shape, **let it error**
  — a loud failure surfaces a real data problem.
- Keep task-specific code in clearly-marked sections and in `<task>_*` files; leave
  the shared pipeline functions untouched.
