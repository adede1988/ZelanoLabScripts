## Quick orientation

This repository is a collection of MATLAB scripts for preprocessing and analyzing intracranial EEG and respiration data (Dupi / O15 pipelines). Primary entry points and patterns an agent should know:

- Run / orchestration: `O15PreProc_main.m` — builds `outDat` per session, calls `getSessionParams_O15`, `detect_ttls_O15`, `assemble_outDat_O15`, `preprocess_eeg`, `preprocess_respiration_wholetrace`, `detect_sniffs_from_TTLs`, `build_behavior_table_O15`, `refine_onsets_with_phase`, and saves `[sessID '_O15preproc.mat']` under the session `preProc` folder.
- Single-channel pipeline: `singleChanPipeline_general.m` — loads per-channel `.mat` (expects `chanDat`), computes TF (`getChanTrialTF_highPrec`), FOOOF (`fooof_basic`), HFB, ISPC/PPC, and saves into `CHANDAT_processed/` and `finished/` folders.
- EEG preprocessing: `preprocess_eeg.m` — expects `standardEEGlocs` from `myEEGcoords_thetaPhi.csv` and validates channel labels; performs noise channel detection (`removeNoiseChansVolt`), blink removal (`blinkRemoveWrapper`), interpolation (`interpolate_perrinX`), perrinX Laplacian (`laplacian_perrinX`), and appends `blinkIndicator` and `badTS` channels.

## Key repository conventions and invariants

- Data structs: two primary saved structs are `outDat` (session-level) and `chanDat` (channel-level). Scripts expect specific fields; do not rename fields without updating all callers.
- Channel count & labels: many routines assume 32 EEG channels (see `preprocess_eeg.m`) and strict label order — `preprocess_eeg` will error if `standardEEGlocs.Label` doesn't match `outDat.labels`.
- Hard-coded paths: `O15PreProc_main.m` contains `codePre` and `datPre` arrays and drive-letter network shares (e.g., `R:\Neurology\Zelano_Lab\...`). Agents should prefer editing those variables rather than changing file-level path concatenations.
- Save locations & formats: large outputs use `save(..., '-v7.3')`. Channel pipeline outputs saved under `CHANDAT_processed/` and `finished/` subfolders in the source folder of the channel file.
- Sampling/downsampling: TF computations are downsampled to ~50 Hz for storage (see `singleChanPipeline_general.m`), and some operations assume `fs` fields exist. Respect existing downsampling to avoid huge MAT files.
- Parallelism: `getChanTrialTF_highPrec(..., 'UseParfor', true)` is used — code expects MATLAB Parallel Toolbox availability. Start a `parpool` when running heavy TF loops.

## External dependencies

- EEGLAB (added via addpath in `O15PreProc_main.m`; commit shows `eeglab2025.0.0` path). Ensure a compatible EEGLAB is on MATLAB path.
- FOOOF (MATLAB port used via `fooof_basic.m`) and Signal Processing Toolbox functions like `pwelch`, `hilbert`, `smoothdata`.
- Parallel Computing Toolbox for `parfor`/`parpool`.

## Common developer workflows (how I run things)

- Start MATLAB, add repo and tools to path (example):

  addpath(genpath('G:/My Drive/GitHub/ZelanoLabScripts'));
  addpath(genpath('C:/Users/<you>/Documents/eeglab2025.0.0'));

- To preprocess a subject end-to-end: edit `O15PreProc_main.m` `datPre` / `sessionIDs` or call `O15PreProc_main` after adapting `codePre`/`datPre` to your environment. Then run `O15PreProc_main` in MATLAB.
- To process channel TF and HFB: ensure channel `.mat` files are present, then call `singleChanPipeline_general(chanFiles, idx, subFiles, type)` from a wrapper script or adapt `splitToSingleChan_allTasks.m` outputs.

## Patterns and examples to follow when editing or adding code

- Preserve struct shapes: callers rely on `outDat.data`, `outDat.labels`, `outDat.fs`, `chanDat.trialDat`, `chanDat.fs`, `chanDat.tim`, `chanDat.task`. When adding fields, use names like `XYZ_flag` or `XYZ_meta` to avoid colliding with existing names.
- Use existing helper functions: e.g., `assemble_outDat_O15`, `getSessionParams_O15`, `detect_ttls_O15`, `preprocess_macros`, `preprocess_respiration_wholetrace`. Prefer reuse over copy/paste.
- When adding heavy computations, mirror existing behavior: downsample outputs (the pipeline expects reduced-size TF variables), and save intermediate results into the `CHANDAT_processed/` or session `preProc` folders.
- Respect boolean flags: many scripts check `if ~isfield(..., 'TF')` or `if ~isfield(..., 'HFB')` to skip completed stages. New tools should follow this pattern.

## Common gotchas an agent should check for

- Channel label mismatch: `preprocess_eeg` will throw an error if labels don't match `myEEGcoords_thetaPhi.csv` — update label generation or the CSV instead of bypassing the check.
- Large MAT files: if saving fails or memory is exhausted, use `-v7.3` and consider partial saves. Many scripts already use that for big outputs.
- Hard-coded indices: e.g., `rspIDX == 999` indicates one style of respiration channel handling. Search for `rspIDX` usage when changing respiration code.
- Drive-letter and network paths: tests or CI need mapped shares or path overrides; prefer local mirrors for unit testing.

## Where to look for examples

- `O15PreProc_main.m` — session orchestration and path conventions.
- `preprocess_eeg.m` — EEG cleaning, interpolation, laplacian pipeline and label checks.
- `singleChanPipeline_general.m` — TF/FOOOF/HFB pipeline, downsampling and saving conventions.
- `myEEGcoords_thetaPhi.csv` and `standardEEGlocs.csv` — canonical electrode location files.
- `noiseFigs/` — example noise files and visual outputs used by preprocessing diagnostics.

## If you need to change or add tests

- There is no test harness in the repo; suggested quick checks: small-run scripts that load one sample `raw` file, run `getSessionParams_O15` -> `assemble_outDat_O15` -> `preprocess_eeg` and verify `outDat` fields and that files are written to `preProc` with the expected names.

---
If anything above is unclear or you'd like additional examples (e.g., a small runnable smoke-test script or instructions for adding CI), tell me which area to expand and I will update the file.
