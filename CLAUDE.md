# CLAUDE.md — Respiration/EEG preprocessing refactor

## What this project is

MATLAB preprocessing pipelines for intracranial/EEG + respiration data across **four
tasks**: `breathingTask`, `cueTask`, `threshTask`, `O15`. Sessions are tracked in
`dataTracking.xlsx`. Today, session lists and per-session parameters are hard-coded and
duplicated across many `.m` files; adding a participant means editing parallel arrays by
hand in up to seven scripts plus a large per-session `switch`. The refactor makes
`dataTracking.xlsx` the single source of truth and routes everything through a small set
of loader functions.

## Data flow (current)

```
RAW (Neuralynx .mat + behavioral CSV)
   │   *_makeOutDat.m         (breathing/cue/thresh only; heavy per-session ingestion)
   ▼
<root>/<id>/preProc/<id>_<task>PreProc.mat   (the "outDat")
   │   *_main.m  → getSessionParams_<task>(S)  → loads that .mat back in as `raw`, sets P
   │            → (O15 only) detect_ttls_O15(raw,P)
   │            → assemble_outDat_<task>(raw,S,P)  → fresh outDat
   ▼   ...shared pipeline (downsample_data → preprocess_eeg → preprocess_macros → ...)
final <id>_<task>preproc.mat
```

O15 has **no** `_makeOutDat`: ingestion lives inside `O15PreProc_main.m` via
`detect_ttls_O15` and `getSessionParams_O15` loads genuinely raw `raw_O15.mat`.

## Refactor target (what to build)

1. **`applyParams.m`** — single source of truth. Reads `dataTracking.xlsx`. Two modes:
   session list (for the loops) and P-struct (per session). Replaces the hard-coded
   arrays and the entire `getSessionParams_*` parameter `switch`.
2. **`detectBeats.m`** — one generic ECG beat detector parameterised by a spec string;
   replaces the 25 bespoke `getBeats_*` local functions in
   `getSessionParams_breathingTask.m`.
3. **`assemble_outDat_all.m`** — `[outDat, raw, TTL] = assemble_outDat_all(S, P)`.
   Combines the raw-loading half of `getSessionParams_*` (+ `detect_ttls_O15` for O15)
   with the four `assemble_outDat_*` assemblers, branching on `P.task`. Future direction:
   this becomes the only loader (P comes in from `applyParams`, raw is loaded here).
4. **Rewrite the 7 scripts** to call `applyParams` at the top and
   `assemble_outDat_all` in place of `getSessionParams_* + detect_ttls + assemble_*`:
   `breathingTask_makeOutDat.m`, `cueTask_makeOutDat.m`, `threshPreProc_makeOutDat.m`,
   `breathingTaskPreProc_main.m`, `cueTaskPreProc_main.m`, `O15PreProc_main.m`,
   `threshPreProc_main.m`.
5. **`dataTracking.xlsx`** — already regenerated with the parameter columns (see schema
   below). It is the input the new functions read.

## Key files

| File | Role | Touch? |
|---|---|---|
| `dataTracking.xlsx` | session tracker + (new) parameter columns; the source of truth | regenerated; keep updating at the Excel level |
| `applyParams.m` | NEW loader: session lists + P struct | **create** |
| `detectBeats.m` | NEW generic ECG detector | **create** |
| `assemble_outDat_all.m` | NEW combined raw-load + assemble | **create** |
| `*_makeOutDat.m` (×3) | raw ingestion → `_PreProc.mat` | **edit head only** |
| `*_main.m` (×4) | shared pipeline | **edit head + loader calls only** |
| `getSessionParams_*.m` (×4) | OLD param/raw loaders | keep as **test oracle**; eventually delete |
| `assemble_outDat_*.m` (×3/4) | OLD assemblers | keep as oracle; eventually delete |
| `detect_ttls_O15.m` | O15 photodiode/TTL parser | unchanged; called inside `assemble_outDat_all` |
| `downsample_data / preprocess_eeg / preprocess_macros / process_respiration_breathing / alignTargetBreathingTraceSimplify / build_behavior_table_* / processECG / flagBadBreaths / preprocess_respiration_wholetrace / refine_onsets_with_phase / detect_sniffs_from_TTLs / plot_*` | shared pipeline helpers | unchanged |

## Hard constraints / conventions (do not violate)

- **`datPrei` index order is load-bearing.** `*_makeOutDat.m` and `threshPreProc_main.m`
  branch on `datPrei(...)==1/2/3` to set `outDat.type`. The order **must** stay:
  `1 = Dupi`, `2 = OBEControl`, `3 = EEGbreathing`. `applyParams` returns `datPre`
  with those three first (always present), extras appended after.
- **Roots:**
  - Dupi → `R:\Neurology\Zelano_Lab\Lab_Common\Dupi\`
  - OBEControl → `R:\Neurology\Zelano_Lab\Lab_Common\OBEControl\`
  - EEGbreathing → `R:\Neurology\Zelano_Lab\Lab_Common\AllStudyData\EEGbreathing\`
  - Per-session root now lives in the `datPre` column of the sheet (fall back to the
    Type→root map only if blank). This lets new roots be added at the Excel level.
- **"Unchanged after assemble."** In each `_main.m`, everything from the first
  `outDat = downsample_data(...)` onward is byte-for-byte unchanged. In O15 that means the
  line `outDat.TTL = TTL;` (immediately after the old assemble) stays — so
  `assemble_outDat_all` must **return** `TTL`.
- **`raw` is used after assemble** in cue/thresh/O15 (`build_behavior_table_*(sniffs,
  raw.beh)`), so `assemble_outDat_all` must **return** `raw`.
- **`_makeOutDat.m` bodies are unchanged below the array block.** The bespoke per-session
  ingestion (`if/elseif/switch` on `sessionIDs{sessi}`, the `newList`/`newSet` "standard"
  branch, `datPrei==1/2` type-setting, the save) all stay. `applyParams` only replaces the
  top arrays with drop-in variables of the **same names** (`sessionIDs`, `datPre`,
  `datPrei`, `newList`/`newSet`, `rspIDX`, `rspFlip`).
- **Matching is case-insensitive and whitespace-trimmed.** Sheet IDs are canonical at
  runtime (e.g. `…_KS_2` vs the old switch's `…KS_2` in caps are the same session).
- **Windows / case-insensitive FS.** Files use CRLF. Some legacy filename quirks exist
  (`_breathingPreproc.mat` loaded vs `_breathingPreProc.mat` written; `redoAGAIN` prefix;
  the `250904_OBE_NWU_TI` vs sheet `250904_OBE_NWU_TI_1` mismatch). These are **known and
  intentionally left as-is** for now — do not "fix" them as part of this work.

## New spreadsheet schema

Parameter columns appended to `Sheet1` (header row = row 2; data from row 3). Populated
only for rows whose `Task` maps to one of the four pipelines **and** that have raw data
extracted; `datPre` is filled for every row with a known study.

| Column | Type | Tasks | Notes |
|---|---|---|---|
| `datPre` | text path | all | per-session root |
| `rspIDX` | int | all | respiration channel |
| `rspFlip` | ±1 | all | |
| `hasEEG` | true/false | all | default true (O15 false) |
| `spikeClean` | true/false | all | default true (O15 false) |
| `spikeThresh` | int | all | default 20 |
| `spikeWin` | int | all | default 11 |
| `macroRemove` | "a,b,c" | all | "" = none → `[]` |
| `hasMacros` | true/false | breathing | |
| `respThresh` | num | cue/thresh/O15 | default 500 |
| `cuedBackBuff` | num | cue/thresh/O15 | default 150 |
| `adjWin` | num | cue/thresh/O15 | default 500 |
| `beatSpec` | text | breathing | detectBeats spec; default `1,0,gt,3` |
| `ttlRemoveIdx` | "a,b" | O15 | aberrant TTL indices to drop |
| `ttlNote` | text | O15 | |
| `isNewStd` | true/false | breathing/cue/thresh | drives `newList`/`newSet` membership |
| `paramSource` | `curated`/`guess` | target rows | `guess` = carried forward, review before trusting |

"Guess" rows were filled by carrying forward the most recent **same-study** curated
session's values (including `beatSpec`). `getSessionParams_*.m` remain the authority for
how curated values were derived — use them as the test oracle.

## Function contracts (summary; full detail in taskList.md)

```matlab
% Mode A — session list for a loop
cfg = applyParams(task, stage [, xlsxPath])   % stage = 'makeOutDat' | 'main'
%   cfg.sessionIDs (n×1 cell)  .root (n×1 cell)  .datPre (1×k cell, fixed order)
%   cfg.datPrei (1×n)  .isNewStd (1×n logical)  .newIDs (cell)  .rspIDX/.rspFlip (1×n)
%   cfg.paramSource (1×n cell)

% Mode B — parameters for one session
P = applyParams(task, sessID [, xlsxPath])

% task ∈ {'breathingTask','cueTask','threshTask','O15'}
% default xlsxPath = R:\Neurology\Zelano_Lab\Lab_Common\Admin\dataTracking.xlsx
```

```matlab
heartBeats = detectBeats(ECGz, beatSep, spec)
%   spec = 'chan,lag,op,thr & chan,lag,op,thr & ...'  (op = gt|lt)
%   must reproduce the legacy getBeats_* exactly (see taskList.md for the algorithm)
```

```matlab
[outDat, raw, TTL] = assemble_outDat_all(S, P)
%   S has .id .root .fig (and .figPath for O15); branches on P.task
%   TTL is the O15 TTL table (from detect_ttls_O15); [] for the other three tasks
```

## How to validate (always run after changes)

1. **detectBeats vs legacy.** For each curated breathing session, compare
   `detectBeats(ECGz, beatSep, spec)` against the matching legacy
   `getBeats_<name>(ECGz, beatSep)` on the same `ECGz` (real or synthetic) — outputs must
   be identical index-for-index.
2. **P parity.** For sessions present in the old switches, build P via `applyParams` and
   compare every field to `getSessionParams_<task>(S)` (ignoring `raw`). Expect a match.
3. **assemble parity.** Compare `assemble_outDat_all` output fields to the legacy
   `assemble_outDat_*` for one session per task; for O15 compare `TTL` to `detect_ttls_O15`.
4. **End-to-end.** Run one already-processed session per task through the rewritten
   make+main and confirm no errors and equivalent output where a reference exists.

Keep the legacy `getSessionParams_*.m` and `assemble_outDat_*.m` in the repo until all
validations pass; only then retire them.
