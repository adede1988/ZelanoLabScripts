# cueTask_analysis — analysis plan 


---

## 0. Scope, inputs, and shared conventions

### 0.1 Data source (what each file gives us)
- **Files:** every cue final `<root>\<id>\preProc\<id>_cueTaskPreproc.mat`, top-level var
  **`outDat`** (load via `fieldnames`; cue is always `outDat`). Session list from
  `applyParams('cueTask','main')` (≈35 sessions), filtered to those that actually exist on
  disk **and** have ≥1 `macBP` channel.
- **Signal:** `outDat.data` `[nChan × nSamp]` @ **`outDat.fs = 500 Hz`**; channels addressed
  by `outDat.labels`.
- **macBP channels:** `find(cellfun(@(x) contains(x,'macBP'), outDat.labels))` → `macBP1..5`
  (bipolar macro pairs; **not all files have all 5**; count = nMacro−1).
- **Events (all in 500 Hz samples):**
  - `outDat.TTL.trialStart` — trial-onset sample per trial (cue `TTL` table is `[nTrial×3]`
    `{trialStart, response, sniff}`).
  - `outDat.behDat.finalOnset` — phase-refined sniff onset per trial (cue = one sniff/trial,
    `wiTriali==1`); also `outDat.behDat.n` (trial number) and `type` (`hit/miss/cr/fa`).
- **Identity:** `outDat.sessID` (e.g. `250623_Dupi_NMH_KS_2`), `outDat.type`
  (`Dupi`/`OBE`/`EEG`). Parse `subID`/`sessNum` from the sessID like
  `splitToSingleChan_allTasks.m`: `bits=strsplit(sessID,'_')` → `subID=bits{4}`,
  `sessNum=str2double(bits{5})` (default 1 if absent).
- **Figure folders:** per-subject = recompute portably as
  `fullfile(labPaths().figPath, sessID, 'cueTask')` (≡ stored `outDat.figs`); create if
  absent. Group = `R:\Neurology\Zelano_Lab\Lab_Common\Adam\Dupi_processing\groupStatFigs`.

### 0.2 Toolboxes & coding rules (per request)
- Use **EEGLAB** (`C:\Users\Adam\Documents\eeglab2026.0.0`) and **FieldTrip**
  (`C:\Users\Adam\Documents\fieldtrip-20260518`) functions; write only glue/formatting code.
- **Pipe to FieldTrip *through* EEGLAB** (`eeglab2fieldtrip`) so formatting stays consistent.
- Prefer EEGLAB unless a FieldTrip function is specifically needed.
- **Resolved tool choices (from clarifying Q&A):**
  - **FOOOF → FieldTrip** `ft_freqanalysis` (`cfg.output='fooof_peaks'` / `'fooof_aperiodic'`;
    uses FieldTrip's bundled brainstorm `process_fooof`, no Python).
  - **Spectrograms → EEGLAB `newtimef`** (log-spaced freqs, dB baseline); group means by
    averaging the per-subject dB ERSP matrices.
- **Path hygiene (gotcha):** EEGLAB and FieldTrip share many function names. Plan to
  `ft_defaults` for the FOOOF step and run EEGLAB (`eeglab nogui`) for the spectrogram step,
  managing `addpath` order so neither shadows the other (and `eeglab2fieldtrip` lives in the
  EEGLAB `dipfit` plugin). Document the working path order in the implementation.

### 0.3 Time/sample conventions
- `1 sample = 2 ms` (fs=500). Helper: `ms→samp = round(ms/1000*fs)`.
- Spectrogram **display** window: **−1000 … +3000 ms**. **Epoch** (extraction) window is
  padded to **−1750 … +3750 ms** to avoid wavelet edge effects at 2 Hz (≈3 cycles = 1500 ms);
  display is cropped to −1000…+3000 ms after TF estimation.
- **Baseline** (both lockings): **−700 … −200 ms relative to `trialStart`**.
- Drop any event whose padded epoch would run off either end of the recording.

### 0.4 Groups (updated to FOUR groups)
Assign each cue session by `type` + `sessNum`:
| Group | Rule |
|---|---|
| `DupiS1` | `type=='Dupi'` & `sessNum==1` |
| `DupiS2` | `type=='Dupi'` & `sessNum==2` |
| `DupiS3` | `type=='Dupi'` & `sessNum==3` |
| `Control` | **all non-Dupi** cue sessions (OBEControl + EEG studies pooled) |
Dupi sessions with `sessNum>3` (if any) are reported in the data table but excluded from the
4-group means unless you say otherwise.

### 0.5 Trial inclusion (assumption to confirm)
Use **all trials** for spectrograms (not split by `hit/miss/cr/fa`). Flagged in §6.

---

## Tasks: 

*Shared pipeline (order matters — analyses reuse intermediates)
Run per subject, use the pre-determined best mac channel.
---
Load & preprocess. Per trial: demean, detrend. Keep three copies: raw (for TF), broadband-gamma (zero-phase FIR 25–58, ft_preproc_bandpassfilter with 'firws'), and the analytic signal of the broadband-gamma copy (Hilbert). Your 58 Hz ceiling already dodges 60 Hz line — good; don't notch inside band.
- Per-trial time-frequency (also produces the spectrogram figures). Primary: superlets (Get code from Moca et al. 2021 at https://github.com/TransylvanianInstituteOfNeuroscience/Superlets.git) for short-burst resolution; fallback: FieldTrip ft_freqanalysis mtmconvol, freq-dependent window ~5–7 cycles, 20–70 Hz, 5 ms steps; optional reassigned spectrogram for the figure overlay. From this extract the ridge f̂(t) = power-weighted peak frequency per time bin over the burst.
- Chirplet MP (below). Its dominant atoms give independent estimates of component frequencies, which feed the phase anchors and the beat frequency — this keeps all three analyses mutually consistent rather than each inventing its own f_hi/f_lo.
- Phase continuity test.
- Beating test.
- QC / inclusion. Define if each trial is included using the noise flagging that has already been implemented. All analyses run on included trials only; store the flag so group n is honest.
Save struct + figures, append trial rows to aggregate output .csv file

Anchors, per subject: f_hi = ridge frequency at ~20th percentile of time-in-burst, f_lo at ~80th percentile, clamped to [27, 56] so filter skirts stay in band. f_beat, per trial: |f1−f2| from the two dominant chirplet atoms (model B below), or ridge early-minus-late as fallback.
Output schema
sub.phaseTest

.params: f_hi, f_lo, window, filter spec, IF source ('ridge')
.trial(i): t_hi, t_lo, dur, phiObsInc, phiPredInc, resid (wrapped), included
.summary: nUsed, Rbar, meanAngle, rayleigh_z, rayleigh_p, vtest_p (against 0), perm_z, perm_p

sub.beating

.params: env method, overlap window, surrogate spec
.trial(i): fBeat, envPeakFreq, modSNR, modIndex, p_vs_neighbor, included
.summary: nUsed, medianLogSNR, fracSig, sepBeatCorr (see stats), test stats

sub.chirplet

.params: dictionary ranges, stopping rule, refinement, A/B constraints, Neff definition
.trial(i): atoms (array of tc,fc,c,sigma,amp,phase,energyFrac), nAtomsSig, c1, absC1, signC1, errA, errB, kA, kB, dBIC (B−A; >0 favors single), f1,f2,tSep,overlap, included
.summary: medianC1, fracDownchirp, meanDBIC, fracFavorSingle, test stats

CSVs (tidy, one row per trial — feeds lmer/ggplot directly): phaseTest_trials.csv, beating_trials.csv, chirplet_trials.csv, each keyed on subject, channel, trial, plus subject_summary.csv (one row per subject×channel with the .summary fields). Everything needed to redo group stats lives in the struct and the CSVs, so the group script never recomputes single-trial quantities.
Figures: …/<subjectFigPath>/singleTrialSpectrograms/sub-XX_ch-YY_trial-NNN.png — TF with ridge overlaid and, once chirplet is done, the dominant atom's instantaneous-frequency line overlaid. Plus one per-subject montage: a grid of ~12 example trials + the trial-averaged power TF (the power average is fine and shows the mean chirp; just remember the time-domain average would kill induced gamma).

Task list, by analysis
A. Single-trial spectrograms
Mostly covered above. The deliverable is visual triage: in the transition window, does any trial show two ridges (one fading high, one growing low) versus a single moving ridge? Annotate each figure with the f_hi/f_lo anchors and the burst window so they're comparable across trials. Look at many; single trials are noisy.
B. Phase continuity (the duration-confound–corrected version)
Per included trial:

From the ridge, find t_hi (first descending crossing of f_hi) and t_lo (crossing of f_lo); dur = t_lo − t_hi.
phiObsInc = unwrapped broadband-gamma Hilbert phase at t_lo minus at t_hi.
phiPredInc = 2π · ∫_{t_hi}^{t_lo} f̂(t) dt, integrating the ridge IF (trapezoid).
resid = wrap(phiObsInc − phiPredInc).

Critical pitfall to hard-code against: f̂ must come from the spectrogram ridge (power-based), not from the derivative of the Hilbert phase. If you set f̂ = dφ/dt the residual is identically zero and the test is vacuous. The whole point is that ridge frequency (power) and analytic phase are independent estimators that coincide for one coherent oscillator and diverge (because of the independent φ₁, φ₂) for two.
Why this answers your critique: variable chirp duration changes both phiObsInc and phiPredInc together via the measured trajectory, so it cancels in the residual. What does not cancel is a true phase discontinuity through the transition. Interpretation stays asymmetric: concentrated resid → single; dispersed → inconclusive (now also because IF estimation is noisy, not only because of timing). The one residual confound is two oscillators that are phase-locked to each other — their fixed φ₁−φ₂ also yields concentration — which is the "tightly coupled = partly semantic" case; flag it, don't pretend the test beats it.
C. Envelope beating
Per included trial, within the overlap window only (the transition region where both putative components would coexist; exclude onset/offset transients):

Envelope = abs(hilbert(broadband_gamma)), lightly smoothed; high-pass the envelope (~10 Hz) to kill the slow onset/offset shape that otherwise masquerades as low-frequency modulation.
Spectrum of the envelope over the window. modSNR = power at fBeat ÷ mean power in flanking bins (the neighbor-bin null is simpler and more robust than surrogates here; AAFT as a secondary check). modIndex = modulation depth at fBeat. p_vs_neighbor from comparing the fBeat bin to the flanking distribution.

Constraints to respect: you need at least one beat period inside the overlap window (≈50 ms at 20 Hz separation), and unequal component amplitudes shrink modulation depth — so a null here is genuinely uninformative, as you said. A nice confirmatory check at group level: the detected envelope-modulation frequency should track the measured |f_hi−f_lo| across trials/subjects (sepBeatCorr); real beating scales with separation, artifact doesn't.
D. Chirplet matching pursuit (Mann & Haykin)
Atom (Mann & Haykin 1991/1995): Gaussian-windowed linear-FM,

g(t) ∝ exp(−(t−t_c)²/2σ²) · exp(i[2π f_c(t−t_c) + π c(t−t_c)² + φ]), instantaneous frequency f_c + c(t−t_c), five real structural params (t_c, f_c, c, σ) plus complex amplitude.

Coarse dictionary over constrained ranges: t_c spanning the burst (~10 ms grid), f_c ∈ [25,58] (~3 Hz), c ∈ [−400, +100] Hz/s biased negative for downchirp (~50 Hz/s), σ at 2–3 burst-relevant scales. Keep it coarse — it's only for initialization.
MP loop: find atom maximizing |⟨residual, g⟩|, refine its continuous params with local nonlinear optimization (Newton/EM à la Bultan 1999 — lsqnonlin with the linear amplitude projected out via VARPRO so you only search t_c, f_c, c, σ), subtract the projection, repeat.
Stopping / significance: stop on residual-energy fraction; declare an atom "significant" (nAtomsSig) only if it captures more energy than the 95th percentile of MP run on phase-randomized surrogates of that trial. This gives an honest atom count rather than a fixed number.

Then the clean two-model comparison inside the chirplet framework (this is what collapses your original step 2–4):

Model A (single): best 1 atom, free chirp rate c → 6 real params; errA.
Model B (dual): best 2 atoms constrained to |c| ≤ c_tol (≈stationary) at distinct f_c → 10 params; errB. Record f1, f2, tSep, overlap.
Compare with BIC on effective df: Neff ≈ number of gamma cycles in the window, not raw samples — oscillatory residuals are autocorrelated and raw-N BIC will be wildly overconfident. dBIC = BIC_B − BIC_A with the right k (6 vs 10). dBIC > 0 favors single. The rigorous alternative is cross-trial CV of the structural params (c for A; f1,f2 for B) — fit on a training split, evaluate held-out reconstruction — which shrinks the flexibility gap to near zero; offer both, default to BIC-with-Neff for tractability.

Readouts per trial: c1 (dominant-atom chirp rate; sign and magnitude), nAtomsSig, dBIC. A strongly negative c1 in one dominant atom is the single signature; two near-zero-c atoms at separated f_c with offset t_c is the dual signature.

Group-level statistical adjudication
The recurring template — and it fits your lmer workflow — is per-trial statistic → per-subject summary → random-effects across subjects, with subject as the sampling unit so it generalizes to the population. Either two-stage (summary-statistic RFX) or single-stage mixed model metric ~ 1 + (1|subject) testing the intercept; for thousands of trials with unequal n per subject the mixed model is cleaner and you already use lme4.
Phase continuity (circular, asymmetric). Per subject: mean resultant length Rbar and angle of {resid}, plus a per-subject permutation null (rotate/shuffle resid to uniform, recompute Rbar, → perm_z/perm_p) — this absorbs the n-dependence of Rbar. Group: one-sample t-test of perm_z > 0 across subjects (random-effects), plus proportion of subjects individually significant. Use a V-test against 0 (you predict concentration at zero residual, not just non-uniformity) for sensitivity. The gold-standard alternative is a Bayesian von Mises mixed model (subject random effects on mean and concentration κ) — mention it, but permutation→RFX-t is the practical primary. Crucially: a non-significant group result cannot be read as "dual." It's inconclusive; equivalence testing doesn't help because dispersion is the genuinely ambiguous regime.
Beating (asymmetric, mirror image). Per subject: distribution of modSNR at fBeat; per-subject Wilcoxon of observed vs neighbor-null, summarized as medianLogSNR. Group: one-sample t of medianLogSNR > 0 across subjects, plus fracSig. Confirmatory: sepBeatCorr — a mixed model envPeakFreq ~ freqSeparation + (1|subject) with a positive slope is real corroboration that the modulation is beating. Significant → dual; null → inconclusive, symmetric to the phase logic.
Chirplet (symmetric — does the real adjudicating). Two separate group questions:

Is the chirp real? Per-subject medianC1; group one-sample t of medianC1 < 0 (downchirp), plus fracDownchirp. This validates the phenomenon before interpreting model preference.
Single vs dual? Per-subject meanDBIC; group one-sample t of meanDBIC vs 0 (positive → single wins on average), plus fracFavorSingle. Because this one is two-directional, it's the test that can actually reject single in favor of dual, unlike the other two.
Optional heterogeneity probe: if you suspect some trials are single and some dual, test bimodality of per-trial c1 or dBIC (Hartigan dip test, or a 2-component mixture). A bimodal dBIC across trials within subjects would itself be an interesting finding (mixed generators across sniffs).

Combining. Build a per-subject (or per-channel) verdict vector and read the convergent pattern from the top table, rather than a single omnibus p. The strongest claims come when chirplet's symmetric verdict agrees with whichever asymmetric test is positive: chirp-real + concentrated phase + no beat = single you can defend; chirp-real + beat present + dBIC<0 = dual you can defend. If you have multiple contacts, add the spatial check (do model-B f1/f2 localize to different contacts?) as a separate, possibly decisive line.
