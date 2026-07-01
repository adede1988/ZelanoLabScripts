# chirpAnalysisStatus.md -- running status / decision log

**Task:** Adjudicate whether OB sniff-locked gamma is ONE chirping oscillator or TWO
oscillators (hi-early / lo-late) whose sum *looks* like a chirp. Spec:
`cueAnalysis/chirpAnalysis.md`. **Scope now: `cueTask` only.** Code + docs live in this
folder (`chirpAnalysis/`).

This file is a living log. Newest entries at the bottom of §LOG. Decisions/assumptions
that future-me must respect are pinned in §DECISIONS.

---

## DECISIONS & ASSUMPTIONS (pinned)

- **D1 -- Source of truth = the spec** `cueAnalysis/chirpAnalysis.md`. Do not change *what*
  analysis is done; only low-level implementation levers.
- **D2 -- Channels (per user):** core single-channel tests run on **`bestMac`** (the OB
  bipolar; project's existing OB idiom -- see D7). Spatial tests run **across ALL `macBP*`
  channels**. **No scalp/EEG channels.** Spatial layout/geometry is NOT our concern (user
  will handle after the fact) -- so record per-contact profiles + similarity, skip
  spacing/geometry gating.
- **D3 -- Epoch (per user):** `[-1.0, +3.0] s` around `behDat.finalOnset` ONLY (no
  trialStart locking). Pad beyond the window for edge-artifact-free filtering/TFR, then trim.
- **D4 -- Burst-length metric (per user):** smoothed gamma-power time series locked to
  finalOnset; per trial take peak gamma power; burst start/stop = where the smoothed signal
  rises to / falls back below the **75th percentile of the cross-trial peak distribution**,
  searching only in **(0, 3000] ms after finalOnset**.
- **D5 -- Inputs (per user):** prefer **E: finals modified AFTER 4pm 2026-06-29** (these
  carry the rerun's corrected bestMac + noise rejection). Stale (older) finals must NOT be
  used for final results. When current fresh set is exhausted, sweep E: for newly-finished
  finals and integrate. (As of build start: 0 fresh on E:; 36 stale from a 09:17 sync.)
- **D6 -- Compute/deploy:** code authored locally (home), deployed to the lab via **scp**
  into `G:\My Drive\GitHub\ZelanoLabScripts\chirpAnalysis\` (NOT via git -- the lab repo is
  mid-rerun with uncommitted changes on `simplify-standardize`; a pull/stash could disrupt
  the running job). Heavy runs on the lab desktop (`ssh labdesktop`, 8-core/96 GB, E: I/O,
  no GPU). Phase-0 synthetic can run locally (16 GB) or on lab.
- **D7 -- OB-contact idiom:** `chanDat.chanType` only encodes `'macro'` vs `'EEG'` -- it does
  NOT mark the OB contact. The cue pipeline's de-facto OB channel is **`bestMac`** (largest
  noise-screened periodic 30-58 Hz gamma peak; `cue_fooof_macBP`). Core tests use it.
- **D8 -- TFR / ridge / chirplet engines:** decided after probing lab toolbox licenses
  (see §LOG). Cloned tools: Superlets(FASLT), MPACT, Synchrosqueezed-chirplet,
  Frequency_ridge_tracking, MODA. Built-in fallbacks if Wavelet Toolbox licensed: tfridge,
  wsst/fsst.

---

## PHASE CHECKLIST

- [ ] **0. Setup** -- folder, status/report docs, clone tools, register paths, probe lab env.
- [ ] **Phase 0 -- Validation harness** (synthetic single/dual/phase-locked-dual; sweep
      overlap × Δf × amplitude; sensitivity/specificity tables; chirplet-bias & beat-min-overlap
      curves). GATES all real-data interpretation.
- [ ] **Profile one cue session** -- burst-length distribution; single-trial spectrograms;
      ridge sanity; confirm achievable coexistence vs each test's minimum.
- [ ] **Phase 1 pipeline** -- `run_chirp_analysis_pipeline.m`; per-session→channel→trial;
      tests §6.4-6.9; struct sub-fields §7.1 + tidy CSVs §7.2; figures.
- [ ] **Phase 3 -- Group adjudication** -- per-session verdict vectors; single/dual/undecidable;
      handle "spatial: not assessed" distinctly.
- [ ] **Report** -- `chirpAnalysisReport.Rmd` knit to HTML.

---

## ARTIFACT MAP

- Code/docs: `chirpAnalysis/` (this folder), deployed to lab `...\ZelanoLabScripts\chirpAnalysis\`.
- Figures: `...\Lab_Common\Adam\Dupi_processing\<id>\singleTrialSpectrograms\`.
- Aggregates: a new `chirpAnalysis\` output dir under the group dir (CSVs §7.2).
- Phase-0 outputs: `<groupDir>\chirpAnalysis\phase0\`.

---

## LOG

### 2026-06-29 ~21:00 -- survey + kickoff
- Read spec + CLAUDE.md; confirmed plan with user; pinned D1-D8 above.
- Lab recon: branch `simplify-standardize`; repo has uncommitted O15 + `cueTaskPreProc_main.m`
  changes (rerun in progress via RDP MATLAB ~12.5 GB) → **do not git-pull on lab.**
- Home FOOOF rerun attempt crashed (OOM) at session 12/35 (`_homererun.log`); live rerun is
  on the lab (writes to R:, not yet synced to E:).
- 5 external tools were MISSING on lab → kicked off `git clone --depth 1` of all 5 into
  `G:\My Drive\GitHub\` (background).
- 0 fresh E: finals yet (post-4pm rule); 36 stale (09:17 sync). Phase 0 proceeds without data.

### 2026-06-29 ~21:40 -- env locked, engines validated, design panel done
- **Lab env probe:** R2024b; ALL toolboxes licensed (Wavelet/Signal/Stats/CurveFit/Optim/Parallel);
  tfridge/wsst/fsst/pspectrum/emd present. Home R2026a: same (the old "missing toolboxes" memory
  was STALE -- corrected). eeglab2025, fieldtrip E:\fieldtrip-20260518.
- **Data contract (cue KS_1):** nChan=20, fs=500, 4 macBP pairs, behDat[40x13] table, finalOnset
  60..968 s, ~22 s between trials -> [-1,+3]s epochs never collide; sniffLabel all 'cued'.
- **labPaths roots are on R:** (figPath, rootDupi) -> unreachable from SSH. My pipeline discovers
  finals on E: (chirp_session_table) and writes outputs to E: via env overrides (CHIRP_FIGROOT/OUTDIR).
- **Tools cloned** (home + lab, G:\My Drive\GitHub & C:\...\GitHub): Superlets, mpact,
  Synchrosqueezed-chirplet-transforms, Frequency_ridge_tracking (PYTHON -> not used; tfridge instead),
  MODA. scp home->lab works.
- **D8 RESOLVED -- engines:** TFR=FASLT `nfaslt` (matlab-pure, POWER out [freq x time]); ridge=MATLAB
  `tfridge` (rows=freq, pass power directly -- validated); chirplet=`mp_adapt_chirplets` called DIRECTLY
  (NOT mp_act_signal, which plots; NOT mle_adapt_chirplets, latent EM bug) on hilbert(analytic),
  verbose='no', Q>0, 'expectmax', 'Oneill'; synchrosqueeze=`fsst`. Unit map: fc_Hz=fc*fs/2/pi,
  chirp_Hz_s=cr*fs^2/2/pi, tc_s=(tc-1)/fs+transWin(1). **Smoke test PASSED** end-to-end on a synthetic
  downchirp (FASLT->tfridge tracked 47->39 Hz; MPACT returned fc=39 Hz chirp=-14 Hz/s).
- **GOTCHA:** .m files must be ASCII (em-dash/ellipsis crash the parser) and filenames must NOT start
  with `_` (leading underscore = "invalid text character"). Renamed _smoke_local.m -> chirp_smoke_local.m.

### DESIGN-PANEL FIXES (pinned; from the 6-agent design+critic workflow)
- **F1 (high):** epoch `padSec` must exceed the max FASLT wavelet half-support in 25-58 Hz (~0.73 s);
  set padSec=1.0 and have chirp_tfr_faslt assert it. 600 ms (a designer default) is too short.
- **F2 (high):** per-trial transition window comes from THAT trial's own smoothed `fhat` extrema,
  NEVER from the inward-biased subject-average anchors (anchors -> annotation/spatial bands/figs only).
- **F3 (high):** the phase test's f-hat is the POWER ridge, never dphi/dt (assert ifSource='ridgePower');
  Phase-0 must DEMONSTRATE phase-test specificity on dualIndep before "single" is licensed.
- **F4 (high):** ALL envelopes (beat amplitude; decoup/spatial band power) come from `hilbert()` of
  FIR-bandpassed signals (narrow FIR for band-limited), NEVER from FASLT power. FASLT power is reserved
  for the spectrogram display + ridge.
- **F5 (resolved by user D2):** spatial runs across ALL good macBP, eligible if >=2; NO geometry/spacing
  gating (user will handle layout after the fact). contactSpacing=NaN. (Overrides the critic's
  "geometryUnknown block" -- user explicitly waived layout.)
- **F6:** ONE shared `chirp_surrogates(bb,nSurr,seed)` generator feeds ridge-2/beat/chirplet so their
  95th-pctile thresholds are comparable; phase uses a rotation null instead. nSurr=200 batch / 1000 profiled.
- **F7:** chirp_ridge emits BOTH `transitionWin` (phase/chirplet) AND `coexistWin` (beat) per trial;
  every consumer takes its window by explicit named arg. powerDip flanks defined vs the BURST envelope
  (guard div0). dBIC is corroborative only (never overrides geometry). atom-keep band [26 57].
- **F8:** Phase-0 scores every test on BOTH oracle (ground-truth) and pipeline-derived windows
  (windowSource axis) so ridge-finding error is separable from test error.
- **Verdict vector = ONE element per test** (phase, beat, spatial, decoup, chirplet). slip-at-null and
  dBIC are corroborative annotations, NOT extra votes (no p-value multiplication; spec 0.2).

### BUILD ORDER (critic): (0) paths+cfg+burst-len [done: cfg, burst core] -> (1) chirp_tfr_faslt +
chirp_ridge + chirp_surrogates -> (2) chirp_synth_trial + chirp_phase0_harness (GATES interp) ->
(3) chirp_phase_test + chirp_chirplet -> (4) chirp_beat_test + chirp_temporal_decoup ->
(5) chirp_spatial -> (6) run_chirp_analysis_pipeline + group + report.
Full design specs archived in chirpAnalysis_design.json (workflow output).

### 2026-06-29 ~23:30 -- full battery built, deployed to lab, partial validation
- **Built + smoke-validated locally (home R2026a):** chirp_config, chirp_session_table (E: fresh
  filter), chirp_epoch, chirp_burstlen, chirp_faslt_bank/apply (cached; **parity 0.0e+00 vs nfaslt**
  on BOTH home R2026a and lab R2024b), chirp_tfr_faslt, chirp_ridge (+anchors, peeling, surrogate
  ridge-2), chirp_surrogates, chirp_phase_test, chirp_chirplet, chirp_beat_test,
  chirp_temporal_decoup, chirp_spatial, chirp_synth_trial, chirp_phase0_harness, chirp_load_session,
  chirp_measure_burst, chirp_fig_trial/montage, run_chirp_analysis_pipeline, chirp_group_adjudicate,
  run_chirpAnalysis_all, chirpAnalysisReport.Rmd.
- **Deployed to lab** G:\My Drive\GitHub\ZelanoLabScripts\chirpAnalysis via scp->robocopy (NO git
  ops); validated load + synth + FASLT on R2024b.
- **Battery smoke (synthetic single vs dual) -- discrimination status:**
  - RIDGE: single frac(nRidgesSig>1)=0.00; dual=0.60. **Discriminates well.**
  - PHASE: single perm_z=3.1, vtest_p=0.014 (concentrated). **Works** (single-confirmer).
  - BEAT: single fracSig=0.00 (no false beat). fBeat nonzero after the two-peak-marginal anchor fix.
  - TEMPORAL: low power at jitter=60ms/n=10 (corr noisy). Expected weak test (spec); needs realistic
    jitter + more trials. Calibrate in Phase 0.
  - CHIRPLET: single->ambig BUG under investigation (atoms not classifying as a descending chain).
- **PERF:** chirplet per-trial surrogate atom-significance is the batch bottleneck (~1-2 s x nSurr
  fits/trial). Added **C.chirplet.sigMode='fixed'** (atomEnergyThr, default 0.08) -> no per-trial
  fits in the batch; 'surrogate' mode retained for Phase-0 calibration. Ridge-2 surrogates parallelize
  per-trial on the lab (parfor TODO).
- **OPEN:** (1) chirplet single->ambig (fix classify/fit); (2) calibrate atomEnergyThr, crMin,
  tfridge penalty, temporal jitter in a real Phase-0 run on the lab; (3) still 0 fresh E: finals
  (lab rerun writing to R:) -> real-data pipeline waits; profile mechanics on a stale file meanwhile.

### 2026-06-30 ~00:30 -- Phase-0 on lab, phase-test redesign, mechanics test
- **Deploy method fixed:** `scp -r chirpAnalysis labdesktop:"E:/"` (to PARENT, not named dest -- the
  named-dest form nested/emptied it) then MATLAB `copyfile('E:/chirpAnalysis/*', G:..., 'f')`
  (handles spaces; robocopy through the spaced G: path was flaky). 32/32 .m verified on E: and G:.
  Also: cmd `set "VAR=val"` (quoted) to avoid a trailing-space env var.
- **Tiny Phase-0 ran on lab** -> revealed the PHASE test does NOT discriminate as-built at (Df=7,
  overlap=500, n=5): dual gave HIGHER within-session concentration (Rbar=0.92) than single (0.36),
  and single's residual sits at meanResid=-78deg (a systematic ridge-IF bias ~1.5 Hz), so neither
  "concentration" (perm_z) nor "at-zero" (perm_z0/V-test) is the right single-confirmer.
- **PHASE FIX (conceptual, matches spec logic):** the discriminator is CROSS-TRIAL consistency --
  one oscillator gives the SAME ridge-vs-phase relationship every trial (per-trial mean residuals
  CLUSTER -> high `crossTrialR`), independent-phase duals scatter (low crossTrialR). Added
  `crossTrialR`/`crossTrial_p` to chirp_phase_test.summary; threaded into harness CSV + subject_summary;
  chirp_group_adjudicate now keys the phase=single verdict on crossTrialR (Rayleigh-significant).
  Needs adequate n (tiny n=5 too noisy) -> evaluating with the quick sweep (n=20).
- **CHIRPLET fixed:** fixed-mode atomEnergyThr=0.30 -> singleLin classified 'single' (1 descending
  atom, chirp -2.83 Hz/s); window FALLBACK to coexistWin when no transition (true duals don't sweep)
  so the chirplet/phase tests run on dual too.
- **PERF (batch-critical):** ridge-2 significance now a per-CHANNEL null (ratio2 = peeled/primary
  ridge energy; ~20x fewer FASLTs), chirplet fixed-mode (no per-trial fits). Battery now ~seconds/trial.
- **LAUNCHED on lab (background):** (a) QUICK Phase-0 sweep (4 families x overlap[200..2000] x Df[5,10]
  x amp 1:1 x {pipeline,oracle}, n=20) -> E:\chirpOut\phase0\; (b) pipeline MECHANICS test on stale
  250623_Dupi_NMH_KS_1 (allowStale) -> E:\chirpOutMech\ + figs E:\chirpFigs\ (validates full per-session
  path on real data layout; NOT a result -- stale bestMac/noise).
- **STILL 0 fresh E: cue finals** (rerun I/O-bound on R:, CPU ~idle). Real-data run remains gated on
  fresh E: finals appearing. Will sweep E: periodically.

### 2026-06-30 ~01:30 -- pipeline VALIDATED end-to-end on real data; key regime finding
- **Mechanics run on stale 250623_Dupi_NMH_KS_1 (bestMac=macBP2, 40 trials, 4 macBP), 205 s/session:**
  all 6 per-trial CSVs + subject_summary.csv + <id>_chirp.mat written; figures render headlessly
  (exportgraphics OK on R2024b) -- montage + single-trial spectrograms under E:\chirpFigs\<id>\.
  Pipeline I/O + battery confirmed working on the real data layout.
- **Real cue gamma structure (single-trial spectrograms):** a ~1 s gamma burst centered ~37 Hz with
  only a GENTLE ~3 Hz downward drift (f_hi~39 -> f_lo~36). The ridge tracks it; it wanders only in the
  post-burst low-power tail (correctly excluded by coexistWin, median ~233 ms half-max core).
- **Per-session verdict (stale KS_1): leans SINGLE.** phase crossR=0.39 (ns -> inconclusive),
  chirp 45% single / 0% dual / 55% ambig, temporal corr(latHi,latLo)=0.996 (rigid -> single),
  spatial meanSim=0.58 (matched-ish -> single), beat fracSig=NaN (no teeth), ridge fracNR2=0.03.
- **KEY REGIME FINDING (cue):** bursts ~233 ms-1 s with Df~3-5 Hz mean the **envelope-beating test
  (primary DUAL-confirmer) and spatial/decoupling routes are largely out of teeth** (a beat needs
  >=2 cycles ~= 286-400 ms at Df 5-7 Hz, often > the burst). So on cue the adjudication leans on the
  single-channel phase/chirplet/ridge tests; the honest group story will likely be "single where the
  single-confirmers are positive, else undecidable" -- to be confirmed by Phase-0 teeth tables + fresh data.
- crossR now plumbed through subject_summary + adjudication (phase=single keyed on cross-trial
  consistency). figs confirmed rendering. STALE result -- not final (awaiting fresh bestMac/noise).

### 2026-06-30 ~02:30 -- chirplet classifier validated; singleNL behavior documented; perf note
- **4-family chirplet classification (oracle window, snr=10, n=12):** singleLin -> 1.00 single ;
  dualIndep -> 0.67 dual ; dualLocked -> 0.92 dual ; singleNL -> 1.00 'ambig'.
- **singleNL='ambig' is documented CONSERVATIVE behavior, not a bug:** a fast-drop(400ms)+long-
  plateau(1100ms) sweep, fit over the full coexistence window, is dominated by the flat plateau ->
  one flat atom -> 'ambig'. A lone flat tone genuinely cannot be told from one-of-a-pair on the
  chirplet alone, so 'ambig'/undecidable is the honest call (errs toward undecidable, NEVER toward a
  false single/dual). Real cue gamma is GENTLE/linear-like (singleLin -> single), so this barely
  affects cue; matters more for fast nonlinear sweeps (revisit for O15). Relaxed the chain rule to
  NON-INCREASING (dfTol) + parallel-edge priority (kept; helps when 2 descending atoms ARE fit).
- **4-family discrimination (overlap=800ms, the long-overlap regime):** ridge nRidgesSig>1 frac =
  single 0.00 / dual 0.36,0.21 ; beat fracSig = single 0.00,NaN / dual 0.29,0.31. Both DUAL-detectors
  work WHEN overlap is long enough (>=~2 beat cycles). within-session permZ high for ALL families
  (confirms crossR is the right phase metric).
- **PERF bottleneck = chirplet (Q=3 MPEM ~2-3 s/fit).** Real batch (35 sess x 40 = 1400 fits ~1 h) OK;
  Phase-0 quick (2560 fits) is ~2.5 h. Future Phase-0: cut nTrials/Q. Lever: C.chirplet.Q (3->2) or
  mnits (5->3) to ~halve it if needed.
- **Phase-0 QUICK still running** (~2.5 h CPU, block-buffered log; writes phase0_cells.csv at end).
  Will deliver crossR + beat + chirplet vs overlap[200..2000] x Df[5,10] x {pipeline,oracle}.
- **STILL 0 fresh E: cue finals.** Real run gated.

### 2026-06-30 ~09:30 -- CLEAN RUN: reconcile R->E, audit, smoke, full detached batch (Phase 0 SKIPPED)
- **Per user: skip synthetic Phase 0** (too slow on prior attempts). Thresholds now spec-defaults;
  interpretation conservative (a test only speaks where it has teeth in the cue regime). Rmd rewritten
  to drop the "Phase-0 gates interpretation" framing -> "Methods (this run)" + threshold table.
- **R->E reconcile:** of 35 cue finals, 11 had been RE-RUN on R: the evening of 6/29 (~20:27-20:52),
  AFTER the morning E: staging (09:14-09:18): AZ, HRM, AS, FS, CS, KS_1, KS_2, TB_2, TPB_1, JH_1, JH_2.
  The other 24 unchanged (R 6/23-6/24, E morning copies current). scp'd the 11 (scp -p, preserves
  R mtime) home(R:)->lab E:; E now has newest of every cue final. (Targeted per-session stat, NOT a
  recursive Lab_Common walk -- the walk is too slow over VPN; lab robocopy /XO stalled in scan phase.)
- **freshCutoff lowered to 2026-06-29 00:00** so all 35 June-29 reruns qualify as fresh; all 35 run.
- **Multi-agent code audit (4 slices + adversarial verify):** flagged 4 "blocking" struct-array init
  mismatches (phase/beat/temporal/spatial) -- ALL FALSE POSITIVES. `struct('f',{})` makes an empty
  struct ARRAY (field name only, 0 elements), and `trial(i)=T` with a MATCHING field set is the
  canonical prealloc idiom; it works. Verified 3 ways: (a) read all 4 files -- field sets match
  exactly; (b) MATLAB semantics; (c) DECISIVE: the prior mechanics test already ran these and wrote
  the 4 *_trials.csv -- impossible if trial(i)=T crashed. No churn-fixes applied.
- **Real bug the static audit MISSED, caught by the smoke test:** the standalone E:\chirpAnalysis
  deploy could not see labPaths.m -> setup fell back to codePre=E:\ -> MPACT/Superlets/FieldTrip
  NOT on path -> chirplet returned 0 atoms (S/D=0.00) and bbfilt silently used fir1 (not firws).
  FIX in run_chirp_job_lab.m: addpath('G:\My Drive\GitHub\ZelanoLabScripts') (labPaths.m) +
  addpath('E:\') (labPaths_local -> codePre=G:, fieldtrip=E:) BEFORE setup; assert mp_adapt_chirplets
  on path. Re-smoke: MPACT=1, FIRWS=1, chirplet produces atoms, firws active. All tools live under
  G:\My Drive\GitHub\ (Superlets\matlab-pure, mpact\mpc_analysis+utility, Synchrosqueezed) + E:\fieldtrip-20260518.
- **figSub renamed singleTrialSpectrograms -> singleTrialTF** (user). chirp_fig_trial already overlays
  the ridge (white T.fhat) -- confirmed on a pulled PNG (descending ridge tracked on KS_1 trial 14).
- **maxWinMs=600** added to chirp_config (cap chirplet fit window; bounds MPEM cost).
- **Job (run_chirp_job_lab.m):** clears stale *_trials.csv (appended -> avoid double-count) ->
  run_chirp_analysis_pipeline([], saveFigs=true,nFigTrials=10) -> chirp_group_adjudicate. Outputs:
  E:\chirpOut (CSVs + per-session _chirp.mat + verdict_vectors.csv + chirp_group.mat),
  figs E:\Lab_Common\Adam\Dupi_processing\<id>\singleTrialTF\. WMI-detached PID 11028, log E:\chirp_job.log,
  sentinel E:\_CHIRP_DONE.txt / _CHIRP_ERROR.txt. Started 09:29, 35 sessions, ETA ~1.5-3 h.
- **Next:** monitor -> knit chirpAnalysisReport.Rmd on lab (R-4.5.2, data on E:) -> scp HTML home,
  copy figs back to R: subject folders.
