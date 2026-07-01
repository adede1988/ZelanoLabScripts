# chirpAnalysis.md — OB sniff-locked gamma: one chirping oscillator, or two?

**Analysis-stage spec for Claude Code.** Sits alongside `CLAUDE.md` (the pipeline +
data-structure reference). This document defines a new downstream analysis that *consumes*
the preprocessed finals / per-channel `chanDat` and adjudicates a single scientific question.
It does not touch the preprocessing pipeline.

---

## 0. The question and why it is hard

In OB sniff-locked single-trial electrophysiology, the sniff drives a gamma power increase
that exhibits a **downchirp**: instantaneous frequency starts high and descends. The question:

> Is this **one oscillator modulating its frequency** down over the trial, or **two
> oscillators** — a higher-frequency one earlier and a lower-frequency one later — whose
> aggregate *looks* like a chirp but is actually two separate generators?

### 0.1 The identifiability problem (read before building anything)

A single descending sweep and two fixed-frequency components that overlap in time and
phase-couple can produce nearly identical short signal segments. That degeneracy is *why* the
aggregate reads as a chirp. So the bare "one vs two generators" question has a region where no
fit can separate the hypotheses. We therefore answer the sharper, decidable question:

> **During the frequency transition, is there a single spectral component whose center
> frequency moves, or two simultaneously-present components at different (slowly-varying)
> frequencies?**

### 0.2 The tests are asymmetric — exploit that

No single test is symmetric; each confirms one hypothesis cleanly and is inconclusive on the
other. The design combines them so their *positive* directions cover both answers. Do **not**
multiply their p-values as if independent — they are the same data viewed from different
angles. Treat the result as triangulation across a per-session verdict vector.

| Test | Positive result | Null result |
|---|---|---|
| Phase continuity | concentrated residual → **strong for single** | dispersed → inconclusive (dual, or single w/ variable timing/noise) |
| Envelope beating | beat at \|f_hi−f_lo\| → **strong for dual** | absent → inconclusive (single also has slow AM) |
| Spatial profile (intracranial) | distinct hi/lo topographies → **strong for dual** | matched → inconclusive (one generator, *or* probe can't resolve) |
| Temporal decoupling | hi/lo latencies decorrelate across trials → **leans dual** | rigid coupling → leans single |
| Chirplet trajectory | one connected descending trajectory → **single**; two parallel overlapping ~flat atoms → **dual** | (geometry-based; see §6.7) |

**Defensible single** = concentrated phase + matched spatial profiles + rigid latency coupling
+ single connected trajectory + no beat. **Defensible dual** = beat present (long overlap) +
distinct spatial profiles + decoupled latencies + two parallel trajectories. Where the
evidence is mixed or the regime is undecidable (see §0.3), **report undecidable** rather than
defaulting to whichever model a flexibility-biased test happens to favor.

### 0.3 The regime: small Δf, long bursts

- **Frequency separation is small:** the hi/lo gap is likely only **~5–10 Hz**.
- **Bursts are long:** hundreds of ms, **up to ~3 s** at the longest.

Frequency resolution scales with the observation window (Rayleigh Δf_min ≈ 1/T):

| Window T | Resolvable Δf |
|---|---|
| 150 ms | ~7 Hz (marginal at 5 Hz) |
| 200 ms | ~5 Hz |
| 500 ms | ~2 Hz |
| 1 s | ~1 Hz |
| 3 s | ~0.3 Hz |

So a 5–10 Hz gap that would be unresolvable in a short burst is **comfortably resolvable** once
the window is hundreds of ms to seconds. **The binding constraint is the coexistence (overlap)
duration of the two putative components, not the total burst length** — but in a long burst the
non-overlap portions are themselves long and resolvable, so the worst case (close + comparable
amplitude + *brief* overlap) is now rare: you'd still see two distinct ~horizontal ridge
segments rather than a smooth diagonal. Phase 0 (§5) quantifies the minimum coexistence each
test needs against your real burst lengths — that number bounds the project.

---

## 1. Scope & current targets

- **Now:** `O15` and `cueTask`. **Eventually:** all four tasks (`breathingTask`, `cueTask`,
  `threshTask`, `O15`).
- Code must be **task-agnostic** in its signal processing; branch only where the data layout
  forces it (sniff selection, epoch guarding). Mirror how the existing pipeline keeps the
  shared signal path identical across tasks.

---

## 2. Integration contract (data, paths, idioms)

### 2.1 Input

- **Preferred input: the per-channel `chanDat`** produced by `splitToSingleChan_allTasks.m`
  under `…\QuestMirror\CHANDAT\`. It already carries `.data .behDat .rsp .fs .task .chi
  .chanType`, is per-channel (loop channels directly), and is what
  `singleChanPipeline_general.m` / `run_gamma_figures_pipeline.m` consume.
- **Fallback: the session finals** `<root>\<id>\preProc\<id>_<task>preproc.mat`. Load robustly:
  ```matlab
  s = load(finalPath); fn = fieldnames(s); outDat = s.(fn{1});   % outDat → chanDat → out
  ```
- **`fs = 500 Hz` for all final signals** (2 ms/sample). Gamma 25–58 Hz is well within Nyquist.

### 2.2 Channel selection

- The OB signal lives among the **intracranial bipolar pairs `macBP*`** (depth/strip macros,
  bipolar re-referenced; count = nMacro−1). Select the OB contact(s) via `chanDat.chanType` /
  the project's existing OB-contact identification. **Never hard-code a row index** — montages
  vary by session.
- Keep the **full `macBP*` set** available in parallel for the spatial analysis (§6.8).
- **Scalp / Laplacian channels are excluded from the spatial discriminator** — the OB does not
  project meaningfully to the 32-ch montage, so the surface Laplacian topography carries no
  usable OB spatial information. (Scalp data remain valid for other analyses, just not as a
  spatial test of OB generators.)

### 2.3 Epoching

- Build sniff-locked epochs around **`behDat.finalOnset`** (phase-refined onset, samples @500).
  Use `finalOnset`, not `sniffOnset`.
- **Window:** extend post-onset to accommodate up to **~3 s** (e.g. `[-1.0, +3.0] s`,
  configurable). Include a pre-sniff segment for the beat baseline (§6.6).
- **Sniff selection (task-specific config, not signal-path branching):**
  - `cueTask`: one sniff/trial (`sniffLabel=="cued"`, `wiTriali==1`).
  - `O15`: multiple sniffs/trial (`sniffLabel ∈ {start,free,confirm}`, `wiTriali` 1–8). Make
    the selected sniff(s) a config field (default to the task-relevant sniff).
- **O15 long-window guard:** truncate each epoch at the **next** `finalOnset` so a long window
  doesn't run into the following sniff; set a `burstTruncated` flag where the window is cut
  short, so a truncated coexistence window cannot masquerade as a genuinely brief one.

### 2.4 Trial QC — defer to the existing pipeline

Noise/sharp-transient trial rejection is **already done upstream**. Consume the cleaned trials
and respect `spikeCleanVec` (all-ones ⇒ none removed), `badTS`, `blinkIndicator`.

> **Do not add a shape-based trial filter.** Selecting "clean-looking chirps" would bias the
> whole analysis toward the single-oscillator conclusion. Select on response *existence*, never
> on response *shape*. The only permitted addition is a minimal "is there any gamma response"
> existence check so phase is not computed on pure noise — and even that is optional; a
> response-free trial simply adds noise to the phase residual, which is acceptable.

### 2.5 Figures

Save under the session figure folder (`chanDat`/`outDat.figs`, i.e.
`…\Lab_Common\Adam\Dupi_processing\<id>\`) in a new subfolder **`singleTrialSpectrograms\`**.
Filename: `sub-<sessID>_ch-<chi>_trial-<NNN>.png`.

### 2.6 Driver & code layout

- New driver **`run_chirp_analysis_pipeline.m`**, structured like
  `run_gamma_figures_pipeline.m`. Session list from `applyParams(task,'main')`. Loop
  `task ∈ {'O15','cueTask'}` now; all four later.
- New analysis functions in **`helperFuncs\`**. Outputs appended to each session's struct
  under new sub-fields (§7); aggregates to CSV under a new `chirpAnalysis\` output dir.

### 2.7 Compute

- Finals are **`-v7.3` (HDF5)** — partial-load with `matfile`/`h5info`; never pull the whole
  `data` matrix into RAM (a session is multi-GB; local box is 16 GB).
- The chirplet MP across thousands of trials × ~35 cue + 38 O15 sessions is the bottleneck.
  **Run it on the lab desktop over SSH** (`ssh labdesktop`; 8-core / ~96 GB; no GPU). `parfor`
  over trials within a channel. Write outputs to `E:\` (never `C:\`). Code synced via git to
  `G:\My Drive\GitHub`. Profile **one** session locally before launching the batch.

---

## 3. External tools (clone + register paths)

Clone into the GitHub dir on each machine (`C:\Users\Adam\Documents\GitHub\` home /
`G:\My Drive\GitHub` lab), sync via git, and register paths in `labPaths.m` (a new block beside
the eeglab root) or a `setup_chirpAnalysis_paths.m` so both machines resolve identically. All
MATLAB except where noted.

| Tool | URL | Role |
|---|---|---|
| **Superlets** | `https://github.com/TransylvanianInstituteOfNeuroscience/Superlets` | single-trial TFRs; use the **FASLT** (fractional adaptive) variant for a smooth sweep |
| **MPACT** | `https://github.com/jiecui/mpact` | adaptive chirplet transform via matching pursuit + EM refinement (MPEM); noise-robust, built for biosignals — the chirplet engine, do **not** hand-roll MP |
| **Synchrosqueezed chirplet** | `https://github.com/ziyuchen7/Synchrosqueezed-chirplet-transforms` | Chen & Wu 2023; disentangles modes with curved/crossing instantaneous frequencies — the right tool for the converging-IF read over long bursts |
| **Frequency_ridge_tracking** | `https://github.com/DavidBondesson/Frequency_ridge_tracking` | lightweight dynamic-programming ridge tracker on the superlet TFR (recommended default ridge extractor) |
| **MODA** (+ **PyMODA**) | `https://github.com/luphysics/MODA` · `https://github.com/luphysics/PyMODA` | principled ridge fallback — pull Iatsenko `ecurve`/path-optimization from `scripts/`, ignore the GUI in `allguis/` |
| **ssqueezepy** *(Python)* | `https://github.com/OverLordGoldDragon/ssqueezepy` | synchrosqueezing fallback if not using the MATLAB route |

> If the MATLAB Wavelet Toolbox is licensed, its built-in `tfridge` + `wsst`/`fsst` are a
> zero-clone option for ridge extraction and synchrosqueezing.

---

## 4. The analysis battery & evidence hierarchy

Two tiers, by data dependency:

- **Single-channel core — runs on every session with a usable OB signal:** phase continuity
  (§6.5), envelope beating (§6.6), chirplet trajectory (§6.7), temporal decoupling (§6.9), plus
  the visual triage (§6.3) and ridge extraction (§6.4). **This tier must carry the
  participants who lack multiple good contacts.**
- **Spatial — runs only where the OB probe has adequate coverage (§6.8):** the strongest
  discriminator at small Δf *when available*, but available only for a subset of participants.

The group-level story therefore splits into "what the single-channel battery concludes across
**all** participants" and "what the spatial test adds in the **subset** with adequate
coverage." Where the two agree in the overlap subset, the single-channel inference is
cross-validated against the spatial one — which then licenses extending the single-channel
conclusion to the no-spatial participants with some confidence.

---

## 5. Phase 0 — Validation harness (build and run FIRST; gates all interpretation)

Synthesize trials at `fs=500` matching empirical SNR and **real burst lengths** (hundreds of ms
to 3 s). Sweep, with **overlap (coexistence) duration as the primary axis**:

- **overlap** ∈ {brief … full-burst}, e.g. {50, 100, 200, 500, 1000, 2000 ms}
- **Δf** ∈ {3, 5, 7, 10} Hz
- **amplitude ratio** ∈ {1:1, 2:1, 4:1}

Families: (i) **single** linear *and* mildly nonlinear downchirp (fast drop → low-frequency
plateau, since a 3 s sweep is unlikely to be perfectly linear); (ii) **dual** fixed/slow tones
with logistic on/off envelopes and **independent** phases; (iii) a small **phase-locked dual**
set (the degenerate case). Push all through the full Phase-1 pipeline.

**Required deliverables (these define what you may interpret on real data):**

1. Sensitivity / specificity of each test as a function of overlap and Δf.
2. **Beat:** the minimum overlap at which ≥2 beat cycles are detected (expect recovery once
   overlap ≫ 1/Δf).
3. **Chirplet:** confirm the single-bias (Model A reflexively winning on true-dual trials) is
   **gone** once overlap comfortably exceeds 1/Δf; and confirm that a true *single* nonlinear
   sweep is captured as **one connected trajectory** even when it needs 2–3 linked atoms (so the
   discriminator is geometry, not atom count — §6.7).
4. **Phase:** confirm segmentation (§6.5) controls cumulative IF-bias error and that
   beating produces phase slips localized to beat nulls.
5. The minimum coexistence each single-channel test needs, compared against the real burst-length
   distribution (confirm on a few real O15/cue trials).

Spatial eligibility is a real-data property, not something to simulate — keep Phase-0 effort on
the single-channel tests, since for a meaningful fraction of participants that tier **is** the
whole analysis. Reuse the existing simulation / power-analysis machinery and `lme4` workflow.

---

## 6. Phase 1 — Per-session → per-channel → per-trial pipeline

Loop sessions (from `applyParams`) → OB channel(s) → trials. Keep the `macBP*` set in scope for §6.8.

### 6.0 Analysis-window definitions (used by several tests)

- **Transition window** = the frequency-changing (or two-ridge-coexisting) region read from the
  ridge (§6.4). **Use this for the chirplet and phase tests** — a long burst may end in a
  stationary low-frequency tail, and forcing a model across that tail just adds noise and atoms.
- **Coexistence window** = the full span where two components are (or would be) simultaneously
  present. **Use this for the beat test** to bank as many beat cycles as possible.
- **Baseline window** = pre-sniff segment (beat control).
- **Matched-power window** = a later within-burst segment of similar total gamma power but less
  frequency modulation (beat control).

### 6.1 Preprocess (per epoch)

Demean, detrend. Keep: broadband-gamma (zero-phase FIR **25–58 Hz**, `ft_preproc_bandpassfilter`
with `'firws'`) and its `hilbert()` (analytic signal → phase + envelope). **No in-band notch**
(58 Hz ceiling already clears line noise).

### 6.2 Time-frequency (FASLT)

FASLT, **20–70 Hz** (pad beyond 25–58 for context), base cycles ~3, adaptive order ~3–30,
frequency step ≤1 Hz. Cache the superlet kernel set across trials within a channel. Optionally
also compute the **synchrosqueezed chirplet TFR** for the converging-IF read (§6.7).

### 6.3 Single-trial spectrograms (elevated — often directly diagnostic now)

Over 1–3 s with Δf=5–10 Hz, the picture itself adjudicates much of the question: a single chirp
is a slowly descending **diagonal** ridge; two tones are two near-**horizontal** ridges (or a
high segment → low segment with a **power dip** at a sequential handover). Save one PNG per
trial (§2.5): FASLT TFR with the ridge, the f_hi/f_lo anchors, and the transition window marked;
once chirplet is done, overlay the dominant atom trajectory. Plus one **per-channel montage**:
≈12 example trials + the trial-averaged **power** TFR (power averaging preserves induced gamma;
a time-domain average would destroy it). Look at many — single trials are noisy.

### 6.4 Ridge extraction

Run the DP ridge tracker (Frequency_ridge_tracking; fallback MODA `ecurve`) on the FASLT TFR
within **25–58 Hz**. Output primary **f̂(t)** and **ridge-power(t)**. Iterative **peeling**: remove
the primary ridge's TF support, re-extract a second ridge, keep it only if its energy exceeds the
**95th percentile of ridges from per-trial phase-randomized surrogates** (now viable given long
bursts). Record: f̂(t), ridge-power(t), `nRidgesSig`, ridge overlap, and the **ridge-power dip
statistic** (min along-ridge power in the transition ÷ flanking maxima — the sequential-handover
signature). Derive the **transition window** here.

**Anchors:** trace the ridge of the **subject-average power TFR** once; read **f_hi** (early,
~20th-percentile transition time) and **f_lo** (late, ~80th), clamped to **[27, 56]**. Rough by
construction (onset jitter smears the average inward); they feed only the phase-test crossing
times, the spatial bands, and figure annotations. The chirplet never uses them.

### 6.5 Phase continuity (PRIMARY single-channel test) — **segmented**

Logic: for one coherent oscillator the broadband Hilbert phase accumulates consistently with the
ridge-IF integral; two independent generators inject independent phase offsets that disperse the
residual across trials.

> **Independence is mandatory:** f̂ must be the **power-based ridge** frequency, never dφ/dt of
> the broadband phase — using dφ/dt makes the residual identically zero and the test vacuous.

> **Why segment (long-burst correction):** φ_pred = 2π∫f̂ dt integrates ridge-IF error over the
> window. A 0.5 Hz IF bias over 3 s = 1.5 cycles of accumulated phase error → wraps and destroys
> the test (even 0.1 Hz over 3 s ≈ 100°). **Do not** run one t_hi→t_lo increment across a
> multi-second transition.

Procedure: tile the transition window into short consecutive sub-windows (~100–200 ms; ~4–8
gamma cycles, short enough that IF-bias accumulation ≪ 1 cycle). Per sub-window: `phiObsInc`
(unwrapped broadband Hilbert phase increment) vs `phiPredInc = 2π∫f̂ dt` over that sub-window;
`resid = wrap(obs − pred)`. A single oscillator is concentrated in **every** sub-window; two
beating oscillators show phase **slips localized to the beat nulls** — so segmentation both
controls cumulative error and *localizes* any discontinuity in time. Store per-sub-window
residuals and a per-trial summary (mean resultant length across sub-windows).

### 6.6 Envelope beating (confirms dual; strong now over long overlaps) — **corrected**

At Δf=5–10 Hz the beat period is 100–200 ms, so a 1 s overlap gives 5–10 cycles, a 2 s overlap
10–20.

> **Drop the 10 Hz envelope high-pass entirely** — it was a short-window patch and would
> *destroy* a 5–10 Hz beat. Over a long window the envelope spectrum has ~1/T resolution and the
> onset/offset shape sits at ≪ Δf, so a plain envelope spectrum isolates the Δf bin cleanly. If
> you filter the envelope at all, the cutoff must be **below the minimum expected Δf (≤3 Hz)**,
> never 10.

Procedure: `env = abs(hilbert(broadband))` over the **coexistence window**; subtract a fitted
smooth envelope (low-order polynomial, or the Model-A reconstruction envelope) to remove the
onset/offset shape; spectrum of the residual. `modSNR` = power at `fBeat = |f_hi−f_lo|` ÷ mean of
flanking bins; `modIndex` = depth at fBeat; require **≥2 beat cycles** in the window. Also use the
more specific signature available at long overlap: comparable-amplitude components drive the
envelope to **deep, periodic nulls** at period 1/fBeat — detect those directly. Repeat in the
**baseline** and **matched-power** windows and store the contrasts.

### 6.7 Chirplet matching pursuit (MPACT) — judge **geometry, not count**

> **Long-burst readout change:** with long overlap the two-atom model is well-conditioned and
> the BIC single-bias is gone (Phase 0 confirms). **But** a 3 s sweep is rarely perfectly linear,
> so a single oscillator may need 2–3 linked linear atoms to trace a curved sweep — meaning
> "number of significant atoms" is **no longer** the discriminator. Judge **trajectory geometry**.

Run MPACT (MPEM coarse→refine). Atom: Gaussian chirplet with center time t₀, width α, chirp rate
β, center frequency ω₀ (maps to t_c, σ, c, f_c). Coarse dictionary: t_c on a ~10 ms grid over the
**transition window**, f_c 25–58 Hz (~3 Hz), c biased negative (downchirp), σ at burst-relevant
scales. Atom significance via per-trial phase-randomized surrogates (95th percentile).

Classify each trial by how the significant atoms arrange:

- **Single:** atoms **chain into one continuous descending trajectory** (centers connect end-to-
  end in time and frequency; chirp rates same sign).
- **Dual:** **two parallel, temporally-overlapping, near-constant-frequency atoms** (|c| small,
  distinct f_c, coexisting in time).

Corroborate with the optional structural comparison where Phase 0 licenses it: Model A (single
atom, free c) vs Model B (two atoms, |c| ≤ c_tol, distinct f_c), BIC on **Neff = number of gamma
cycles in the window** (not raw samples — oscillatory residuals are autocorrelated; raw-N BIC is
wildly overconfident). Record atoms, trajectory-connectivity metric, c1, errA/errB, dBIC.

**Lean on the synchrosqueezed chirplet transform** for the trajectory read over long bursts
rather than forcing linear atoms across 3 s — it is built to track curved/crossing IF.

### 6.8 Spatial profile (intracranial only; **gated** on contact availability)

> **Scalp excluded** (OB doesn't project to the montage). This runs **only across the OB-probe
> `macBP*` contacts**.

> **Eligibility gate (per session, first-class output):** count good OB-probe contacts that pass
> the existing QC, sit on the bulb, and are adequately spaced. Run the test only with **≥2 well-
> separated good contacts** (ideally **≥3** for a profile worth correlating). Record
> `spatial.eligible`, `spatial.nGoodContacts`, and the **contact geometry/spacing**.

> **Bipolar caveat:** `macBP` referencing is itself a spatial high-pass, so a hi and lo generator
> separated by less than the bipolar spacing can be attenuated or mixed within a pair. A **null**
> spatial result on coarse/few contacts is therefore **weaker** evidence than a null on a dense,
> well-covered probe — which is exactly why geometry is recorded alongside the verdict.

Procedure: extract **early-high-band power** (around f_hi, early window) and **late-low-band
power** (around f_lo, late window) on each good `macBP*` contact; form per-trial spatial vectors;
compute the hi-vs-lo **profile similarity** (correlation / cosine across contacts). Two generators
⇒ low similarity. A single oscillator predicts **identical** spatial profiles for its early and
late phases. Store per-trial profiles + similarity; per-session mean.

### 6.9 Cross-trial temporal decoupling (single-channel; pairs with phase)

Under one chirping oscillator, the latency of peak power at f_hi and at f_lo are **rigidly
coupled** across trials (the whole burst shifts together → correlation ≈ 1). Two independently-
timed generators **decouple** them (correlation < 1; the hi→lo gap varies more than one
trajectory allows). Per trial: band-limited envelope at f_hi and f_lo (from the TFR); record peak
**latency** of each and their **difference**. Group-level: across-trial Pearson `corr(latHi,
latLo)` and `var(latDiff)`, compared against the single-oscillator prediction calibrated from
Phase-0 single-chirp sims.

---

## 7. Phase 2 — Outputs

### 7.1 Per-session struct sub-fields (append to the session's `.mat`)

```
sub.ridge.trial(i)          : { fhat[], ridgePower[], nRidgesSig, ridgeOverlap,
                                powerDipStat, transitionWin, burstTruncated, included }
sub.phaseTest.params        : { f_hi, f_lo, subWinLen, filterSpec, ifSource='ridgePower' }
sub.phaseTest.trial(i)      : { subWinResid[], meanRbar, slipLatencies[], included }
sub.phaseTest.summary       : { nUsed, Rbar, rayleigh_z, rayleigh_p, vtest_p_vs0,
                                perm_z, perm_p }
sub.beating.params          : { coexistWin, envFitOrder, hpCutoff(≤3 or none), surrogateSpec }
sub.beating.trial(i)        : { fBeat, modSNR, modIndex, nBeatCycles, deepNullPeriodicity,
                                modSNR_baseline, modSNR_matchedPower, included }
sub.beating.summary         : { nUsed, medianLogSNR, fracSig, sepBeatSlope }
sub.chirplet.params         : { dictRanges, c_tol, Neff='cycles', stoppingRule, ssChirpUsed }
sub.chirplet.trial(i)       : { atoms[tc,fc,c,sigma,amp,phase,energyFrac], nAtomsSig,
                                trajConnectivity, classification('single'|'dual'|'ambig'),
                                c1, errA, errB, dBIC, f1, f2, tSep, overlap, included }
sub.chirplet.summary        : { medianC1, fracDownchirp, fracSingleTraj, fracDualTraj }
sub.spatial.params          : { fHiBand, fLoBand, earlyWin, lateWin }
sub.spatial.trial(i)        : { profileHi[over macBP], profileLo[], profileSim, included }
sub.spatial.summary         : { eligible(bool), nGoodContacts, contactSpacing, meanProfileSim,
                                nUsed }
sub.temporalDecoup.trial(i) : { latHi, latLo, latDiff, included }
sub.temporalDecoup.summary  : { corrHiLo, varLatDiff, nUsed }
```

### 7.2 Aggregate CSVs (tidy; one row per trial; feed `lmer`/`ggplot` directly)

`ridge_trials.csv`, `phaseTest_trials.csv`, `beating_trials.csv`, `chirplet_trials.csv`,
`spatial_trials.csv`, `temporalDecoup_trials.csv` — each keyed on
`session, channel, trial, task, sniffLabel`. Plus **`subject_summary.csv`** (one row per
session×channel with all `.summary` fields, **including `spatial.eligible` and
`spatial.nGoodContacts`**). Everything needed for group stats lives in both the struct and the
CSVs; the group script never recomputes single-trial quantities.

---

## 8. Phase 3 — Group-level adjudication

Template throughout: **per-trial stat → per-session summary → random-effects across sessions**
(session = sampling unit). With unequal n, a mixed model `metric ~ 1 + (1|session)` testing the
intercept is cleaner than two-stage; `lme4`.

- **Phase (confirms single):** per-session resultant length of {sub-window residuals} with a
  per-session permutation null (rotate residuals to uniform → `perm_z`, absorbs n-dependence);
  group one-sample t of `perm_z > 0` + a **V-test against 0** (you predict concentration *at zero
  residual*). Non-significant ⇏ dual — inconclusive.
- **Beat (confirms dual; long-overlap sessions only):** per-session Wilcoxon of modSNR vs the
  neighbor null → `medianLogSNR`; group one-sample t > 0 + `fracSig`. Gate on the Phase-0 minimum
  overlap. Confirmatory mixed model `envPeakFreq ~ Δf + (1|session)`, positive slope ⇒ the
  modulation really is beating. Plus baseline / matched-power contrasts.
- **Spatial (confirms dual; subset only):** **only over eligible sessions.** Test that mean hi-vs-
  lo profile similarity is below a within-component noise ceiling (estimate the ceiling by
  splitting trials and comparing hi-to-hi across halves). Weight by contact count/spacing.
- **Temporal decoupling (leans dual):** group test that `corr(latHi,latLo)` is below the single-
  oscillator prediction and `var(latDiff)` above it.
- **Chirplet (geometry):** group `fracSingleTraj` vs `fracDualTraj`; the chirp-real check
  (`medianC1 < 0`); structural dBIC only where Phase 0 licenses it.

### 8.1 Convergence logic — handle missing spatial correctly

> A per-session verdict vector must distinguish **"spatial: not assessed"** (ineligible session)
> from **"spatial: profiles matched"** (eligible, single-like). **Never average absent-spatial
> together with negative-spatial** — they are different states. The group analysis splits into
> (a) the single-channel battery across **all** sessions and (b) the spatial test over the
> **eligible subset**; agreement in the overlap subset cross-validates (a) against (b) and
> licenses extending (a)'s conclusion to ineligible sessions.

Read the convergent pattern from the verdict vectors (not one omnibus p): defensible single vs
defensible dual per §0.2; where evidence is mixed or Phase 0 marks the regime undecidable,
**report undecidable**.

---

## 9. Definition of done

1. Phase-0 harness exists, sweeps overlap × Δf × amplitude, and emits the power/specificity
   tables + the chirplet-bias and beat-minimum-overlap curves; real burst lengths confirmed
   against it.
2. `run_chirp_analysis_pipeline.m` runs end-to-end on `O15` and `cueTask` from
   `applyParams(task,'main')`, task-agnostic in the signal path.
3. Per-trial single-trial spectrograms + per-channel montages saved under
   `…\<id>\singleTrialSpectrograms\`.
4. All struct sub-fields (§7.1) populated; all CSVs (§7.2) written, with `spatial.eligible` /
   `nGoodContacts` present.
5. Group adjudication (§8) produces per-session verdict vectors with a distinct
   "spatial: not assessed" state and a documented single/dual/undecidable call.
6. Heavy chirplet batch validated on the lab desktop with `parfor`, outputs on `E:\`.

---

## Appendix — parameter defaults (tune in Phase 0)

| Parameter | Default |
|---|---|
| `fs` | 500 Hz (fixed by pipeline) |
| Broadband filter | zero-phase FIR, 25–58 Hz (`firws`) |
| FASLT range / order | 20–70 Hz, base cycles ~3, order ~3–30, freq step ≤1 Hz |
| Epoch window | `[-1.0, +3.0] s` around `finalOnset` (O15: truncate at next `finalOnset`) |
| f_hi / f_lo anchors | from subject-average TFR ridge; clamp `[27, 56]` Hz |
| Phase sub-window | ~150 ms, consecutive (≪ 1-cycle IF-bias accumulation) |
| Beat envelope filter | none, or high-pass ≤3 Hz (must be < min Δf); **never 10 Hz** |
| Beat minimum cycles | ≥2 within the coexistence window |
| Chirplet dictionary | t_c ~10 ms grid (transition window); f_c 25–58 Hz @~3 Hz; c biased negative; σ at burst scales |
| Chirplet stationary tol `c_tol` | small \|c\| (e.g. ≤10–20 Hz/s) for Model-B atoms |
| Chirplet BIC `Neff` | number of gamma cycles in the window (not raw samples) |
| Surrogates | per-trial FFT phase-randomization, ~200–1000, 95th-pctile threshold |
| Spatial bands | f_hi ± ~2–3 Hz (early window), f_lo ± ~2–3 Hz (late window) |
| Spatial eligibility | ≥2 well-separated good `macBP*` contacts (≥3 preferred); record spacing |
| Group test | per-session summary → one-sample t / `lmer` intercept across sessions; permutation null for circular (phase) |

---

### Items for the user to confirm before/with the build
- **OB contact identification idiom** — what `chanDat.chanType` actually encodes for the OB-probe `macBP*` contacts (so selection matches the existing splits).
- **Real burst-length distribution** on a few O15/cue trials (sets the achievable coexistence and thus which tests have teeth).
- **Exact epoch window** and the O15 sniff(s) of interest (`start`/`free`/`confirm`).
