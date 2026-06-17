# Refactor notes — dataTracking.xlsx as the single source of truth

This implements the plan in `taskList.md` / `CLAUDE.md`: route session lists and
per-session parameters through `dataTracking.xlsx` via a small set of loaders, replacing
the hard-coded arrays and the per-session `switch` statements.

## New files

| File | Role |
|---|---|
| `detectBeats.m` | Generic ECG beat detector driven by a `beatSpec` string (`'chan,lag,op,thr & ...'`). Replaces the 25 bespoke `getBeats_*` local functions. |
| `applyParams.m` | Single source of truth. Mode A (`cfg`, session list for a loop) / Mode B (`P`, per-session params). Reads `dataTracking.xlsx`. |
| `assemble_outDat_all.m` | `[outDat, raw, TTL] = assemble_outDat_all(S, P)` — combined raw-load + assemble, branching on `P.task`; calls `detect_ttls_O15` for O15. |

## Edited scripts (surgical: loader swap + path fixes)

- `breathingTaskPreProc_main.m`, `cueTaskPreProc_main.m`, `threshPreProc_main.m`, `O15PreProc_main.m`
- `breathingTask_makeOutDat.m`, `cueTask_makeOutDat.m`, `preproc/threshPreProc_makeOutDat.m`

The legacy `getSessionParams_*.m` and `assemble_outDat_*.m` are kept as test oracles.

## Validation (all green) — see `_dev/`

- **`test_detectBeats.m`** — `detectBeats` vs all 25 legacy `getBeats_*` on synthetic ECG: **75/75** comparisons identical.
- **`test_beatSpec_map.m`** — each curated breathing session's sheet `beatSpec` reproduces the exact legacy detector the switch assigned it: **82/82**.
- **`test_integration.m`** — `applyParams` P-parity + `assemble_outDat_all` field-parity + O15 `TTL` parity vs the legacy loaders/assemblers across breathing/cue/thresh/O15: **0 mismatches** (only the documented `spikeClean` caveat below).

## Machine-specific path fixes (this machine)

- `codePre` → `C:\Users\Adam\Documents\GitHub\`
- eeglab → `C:\Users\Adam\Documents\eeglab2026.0.0`
- breathmetrics addpath removed (toolbox unused / not installed).
- `ft_defaults` in `threshPreProc_makeOutDat.m` commented out (FieldTrip not installed and unused by the script).

## Fixes beyond the pure loader swap (required to run)

1. **EEGLOC file.** `preprocess_eeg.m` requires `theta_deg`/`phi_deg` columns, which only
   `eegLocs_standard_coords.csv` has. The breathing main already used that file; the
   cue/thresh/O15 mains still referenced the stale `myEEGcoords_thetaPhi.csv` (only
   `Theta`/`Phi`), so every EEG session failed with *"Unrecognized table variable 'phi_deg'"*.
   The cue/thresh/O15 mains now read `eegLocs_standard_coords.csv` (matching breathing).

## Default spreadsheet path

`applyParams` resolves its default sheet as: the lab Admin master
(`R:\...\Admin\dataTracking.xlsx`) **if it already has the parameter columns**, else the
repo-local param-enriched `dataTracking.xlsx` next to `applyParams.m`. Today the Admin
master has only the original 27 columns, so the repo copy is used automatically.

**TODO (needs your action):** add the parameter columns (AB:AR, i.e. `datPre`…`paramSource`)
into the Admin master `Sheet1` so the production scripts read from there. A ready, safe
sync (row-aligned, backs up first, writes only the param columns) is in
`_dev/syncAdmin.m` — I could not run it because writing to the shared Admin master is
guarded. Run it yourself, or just keep using the repo copy.

## Known caveat — breathing `spikeClean`

For breathing rows where the legacy `getSessionParams_breathingTask` *commented out* its
`spikeClean` override, the runtime value defaulted to `true`, but the sheet records
`false` (it captured the commented intent). Affected: `KS_1, KS_2, AS(250908)` and the
EEG/wave breathing sessions. **Impact is nil** where `hasMacros=false` (EEG/wave — spike
cleaning never runs) and only matters for `KS_1/KS_2/AS` on a breathing re-run (and the
existing breathing finals were made with `true`). I did **not** edit the sheet — decide
whether you want those cells set to `true` (exact legacy repro) or kept `false` (new intent).

## Batch preprocessing runs & a memory gotcha (16 GB machine)

Running a whole task's sessions inside **one** MATLAB process accumulates the per-session
`outDat`/`raw` (each can be multi-GB at raw sampling) across loop iterations. On this
16 GB machine that exhausts RAM after ~4 heavy sessions → *"Out of memory"* for the rest,
and the OS then swap-thrashes (process startup crawls until MATLAB is killed). Two fixes,
both applied to the `_dev/run_*` harnesses:
- `clear outDat raw R sniffs TTL P ...` at the end of **every** loop iteration (so peak
  memory is one session, not the running sum).
- `maxNumCompThreads(3)` leaves a core for the OS so the box stays responsive.

If OOM still bites on a fresh run, process **one session per MATLAB invocation** (a shell
loop calling `run_O15([idx], true)` per index) so each session gets a clean process.

The double-file cue sessions (`230611_OBE_NMH_AZ`, `241017_OBE_NMH_AS`, `250310_OBE_NMH_FS`)
concatenate two full raw recordings (`comboDat`) and are the heaviest single allocations —
run them when nothing else is competing for RAM.

## Demonstrated end-to-end from raw (new pipeline)

Confirmed working on real data via the new loader: O15 (`250623_Dupi_NMH_KS_1`,
`250818_Dupi_NMH_JH_1`, `250818_Dupi_NMH_JH_2`, `250811_Dupi_NMH_TPB_1`, incl. EEG) and
cue (`250818_Dupi_NMH_JH_1`, incl. EEG). Expected non-regression failures seen in both the
new and legacy paths: missing `raw_O15.mat` for some sessions, and `detect_ttls_O15`
"wrong trial count" for a few O15 sessions with TTL anomalies (these need a curated
`ttlRemoveIdx`, same as before). Corrupt sheet rows (`Subject ID = 0`) fail fast and are
skipped.

## Machine state at end of session (IMPORTANT)

The first O15 batch OOM-ed partway (the un-cleared-loop bug above) and left Windows in a
heavily swap-committed state: afterward even a trivial `matlab -batch "disp"` took >9 min
to start, and a 12-min idle did not clear it. This is an OS/pagefile condition, not a code
problem. **Reboot the machine** (or leave it idle long enough for Windows to reclaim the
pagefile) before running the batches. The memory fixes above are already in place, so the
runs should then complete cleanly.

## Reproducing / finishing the runs (after the machine is healthy)

Deliverable scripts run directly (they self-resolve the sheet). For batch throughput with
per-session error isolation use the `_dev` harnesses (each mirrors its deliverable exactly,
wraps each session in try/catch, and now clears big vars per iteration):

```matlab
% one MATLAB session per task (serial — don't run two heavy batches at once on 16 GB):
addpath('C:\Users\Adam\Documents\GitHub\ZelanoLabScripts\_dev');
genMakeRunners      % regenerate run_cue_make/run_thresh_make WITH the per-iteration clear
                    % (this regen was interrupted by the machine hang; rerun it once)
run_O15([], true)                       % O15 from raw (overwrite)
run_cue_make;  run_cue_main([])         % cue: regenerate intermediates from raw, then main
run_thresh_make; run_thresh_main([])    % thresh: same
```

> The `run_O15` / `run_cue_main` / `run_thresh_main` functions already clear big vars per
> iteration. The two `*_make` *scripts* only get that clear after you re-run `genMakeRunners`
> above — without it, an all-sessions in-process make can OOM on 16 GB (use the per-session
> shell loop instead, which is immune since each session is its own process).

If memory is still tight, process one session per *process* instead (cleanest isolation),
e.g. a shell loop: `for i in 1..N: matlab -batch "addpath(_dev); run_O15([i],true)"`, and for
cue/thresh `matlab -batch "setenv('RUNSUBSET','i'); run_cue_make; run_cue_main(i)"`.

Done so far (new pipeline, from raw): O15 KS_1, JH_1, JH_2, TPB_1; cue JH_1.
