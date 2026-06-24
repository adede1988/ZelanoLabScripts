# SimplifyStandardize.md 

Work the tasks in order. Validate after each. Keep legacy
versions of scripts as oracles available until the end

---

## Task 0 — Orient

- [ ] Read `CLAUDE.md`. Much of the work outlined in CLAUDE.md has already been accomplished, so gain awareness of what has already been done
- [ ] Read the 3 `*_makeOutDat.m`, and the 4 `*_main.m` scripts and understand what data structures they make. Map their dependencies within this repository
- [ ] Find `dataTracking.xlsx` on the R: drive and understand how its parameters feed into the *_main and *_makeOutDat scripts  — Admin master now carries all 29 cols (param cols synced); it is the live source of truth
- [ ] Catalogue what aspects of data processing is the same and what aspects are different for each of the four task types that are handled by this family of processing scripts
- [ ] Catalogue redundancies and inefficiences
- [ ] Note: cue task and threshold task data have been freshly preprocessed. Values found in preprocessed data for these tasks should be replicated precisely after any edits. 

---

## Task 1 Splitting apart task-specific subfunction elements
- [x] Read the 3 `*_makeOutDat.m`, and the 4 `*_main.m` scripts, follow all subfunctions, catalogue all occurences where there is an if statement or switch statement that implements different behavior for data from different tasks
  — The only shared *processing* subfunction with a cross-task switch was `assemble_outDat_all` (the per-task raw-load `switch` + two `if O15` blocks). The shared pipeline fns (downsample/eeg/macros/wholetrace/detect_sniffs/refine_onsets) carry NO task branch. `applyParams`'s switch is the config dispatcher (the prior refactor deliberately unified the per-task `getSessionParams_*` into it). The `behDat.task` checks in `alignTargetBreathingTraceSimplify`/`flagBadBreaths`/`build_behavior_table_breathingTask` are within-task *condition* labels (cyclicSigh/focus), not cross-task.
- [x] identify all elements of subfunctions that are not task specific and shared across all tasks
  — In `assemble_outDat_all`: the figure-dir resolution, the common `raw` fields, the common `outDat` assembly, and the TTL passthrough are shared; the per-task raw load and the O15 extras are task-specific.
- [x] for any subfunction that includes both shared and task-specific elements, make new subfunctions for task-specific elements and remove all task-specific elements from non-task-specific subfunctions
  — Split into: task-specific `assembleRaw_breathingTask` / `assembleRaw_cueTask` / `assembleRaw_threshTask` / `assembleRaw_O15` and `assembleOutDat_O15extras`; shared `assembleOutDat`. `assemble_outDat_all` is now a thin DISPATCHER (the single routing point) — the mains call it unchanged.
- [x] verify that all functions are either completely task specific or completely shared.
  — Every assemble *work* function is now wholly task-specific (`assembleRaw_*`, `assembleOutDat_O15extras`, `detect_ttls_O15`) or wholly shared (`assembleOutDat`, `loadIntermediateRaw`, `resolveFigDir`). The only task-aware function is the documented `assemble_outDat_all` dispatcher (routing only, no task-specific processing). If you want zero task-switches even there, the mains can call the loaders directly — say the word.
- [x] perform an end to end test to verify that changes still yield numerically identical outputs
  — Froze the pre-split `assemble_outDat_all` as `_dev/assemble_outDat_all_oracle.m`; `isequaln(outDat)`, `isequaln(raw)`, `isequaln(TTL)` vs oracle = **1/1/1 for all four tasks** on real data.

---

## task 2 searching for task-shared elements

- [x] Read the 3 `*_makeOutDat.m`, and the 4 `*_main.m` scripts, follow all subfunctions, catalogue all occurences where there are task specific subfunctions being used
  — Task-specific subfunctions: `assembleRaw_<task>`, `build_behavior_table_<task>`, `detect_ttls_O15` (O15), `outMat_to_table` (cue), and the breathing-only `process_respiration_breathing`/`processECG`/`buildECGz`/`paramCheckECG`/`flagBadBreaths`/`alignTargetBreathingTraceSimplify`/`plotBreathLengths`.
- [x] compare task specific subfunctions across tasks. Assess if there is any overlap in code that can be combined and shared across tasks without need for an if or switch statement
  — The 3 intermediate loaders (`assembleRaw_breathingTask`/`cueTask`/`threshTask`) shared their load-and-extract core (differing only by filename + per-task TTL handling). `build_behavior_table_*` already share `behDatFromSniffs` (prior pass). The remaining task-specific fns are single-task with no cross-task peer (no shareable overlap without a switch).
- [x] for any overlap code, find an appropriate existing shared subfunction to package it into ... make a new shared subfunction
  — Made `loadIntermediateRaw.m` (new shared subfunction): load the `<task>PreProc.mat`, extract the stored struct (`outDat`/`chanDat`/`out`), set the common `raw` fields. The 3 loaders now call it and keep only their task-specific filename + TTL handling — no if/switch on task.
- [x] perform an end to end test to verify that changes still yield numerically identical outputs
  — Re-ran the oracle comparison after the dedup: identical (1) for all four tasks.

## task 3 -- Updating documentation
- [x] update the tutorial and comments as appropriate
  — `tutorialPreprocessing.md` §6 (new-task recipe, assembly step) and `preprocessingReadme.md` (data-flow line + §8 file map) now describe the `assembleRaw_*` / `assembleOutDat` / `loadIntermediateRaw` split; `CLAUDE.md` §10 file map updated; every new function has a doc header.

