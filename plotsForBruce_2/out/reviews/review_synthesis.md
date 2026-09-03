Three independent expert reviews were obtained (personas modeled on **Christina Zelano**, **Bradley Voytek**, and **Michael X Cohen**; each preceded by a deep-research pass on that scientist's methods and priors — full reviews in `out/reviews/review_*.md`). The reviews **converged** on four load-bearing issues, all of which were then tested directly. The results are sobering and are reported honestly below.

### What the reviewers agreed on

1. **The ridge often tracks 1/f, not an oscillation.** The FOOOF-selected "best gamma" channel had a *detectable periodic peak in only 39% of recordings* (reported in bold in the macBP section). A forward–backward ridge tracker always returns a ridge, so ridge-frequency/burst measures on the other 61% describe the aperiodic background. (Voytek, Cohen, Zelano.)
2. **The lowest-variance "winners" were bounded-tracker artifacts.** Ridge-frequency is confined to 25–58 Hz, so its mean sits near band-center (~42 Hz) with mechanically tiny CV; the original `goodness = |sep| / (1+CV)` *rewarded* that uninformative stability. (Cohen, Voytek.)
3. **The responder/non-responder axis rests on n=2.** There are exactly two non-responders (JA, PC), both of whom started near control-level composite and could only regress downward. The original ranking was driven by this `sep_resp_nonresp ≈ 2.0` on two participants. (All three.)
4. **Confounds were carried but not removed:** airflow covariates, aperiodic exponent, and the pre-inhale baseline z were all computed but not used to *test* the effects. (Zelano airflow; Voytek 1/f; Cohen z-window.)

### What was done in response (new analyses)

- **Bootstrap CIs on every control↔Dupi separation (session-level).** Result: **0 of 460** measure×task cells (across raw, peak-gated, and airflow-adjusted analyses) have a 95% CI excluding zero. The apparent large effects were underpowered (controls n=5–19) and/or artifact-driven. (`gamma_goodness_v2.csv`, `gamma_headline_controlDupi.csv`.)
- **Peak-gated re-ranking** (restrict to breaths with a FOOOF gamma peak): the separations do not strengthen; ridge-frequency raw≈gated, confirming Voytek's point that gating barely changes the ridge. (`gamma_headline_peakgated.csv`.)
- **Airflow regression** (residualize each metric on inhale volume/duration/peak-flow, re-rank): still 0 CI-significant. (`analysis=="airflow_adjusted"` in `gamma_goodness_v2.csv`.)
- **Aperiodic exponent per group + 1/f-controlled separation.** Dupi have a *steeper* aperiodic exponent than controls (d≈0.28–0.79 by task) — a genuine E/I-type difference (Gao/Voytek). The absolute-power separations (Dupi < control) **partially survive** residualizing on the exponent (e.g., threshTask w5_rpowZ d = −1.30 → −1.00; O15 gammaBumpDb −0.94 → −1.07), so the power gap is *not purely* aperiodic. (`gamma_aperiodic_separation.csv`, `gamma_1f_controlled_separation.csv`, `ap_apExp.png`.)
- **Whole-window vs pre-inhale baseline z concordance.** peakZ correlates only r≈0.11–0.30 with its baseline-z twin (the z-window matters, per Cohen), but *burst presence* agrees 99% and time-above agrees moderately (r≈0.33–0.42). (`gamma_zwindow_concordance.csv`.)
- **Broadband-transient simulation + superlet FWHM.** Passing a purely broadband, aperiodic inhale transient (no oscillation) through the identical pipeline yields a ridge tracker that wanders mid-band (~43 Hz) and an onset-locked broadband response — confirming ridge/onset measures need cautious interpretation. Superlet effective resolution: temporal FWHM ≈ 64–80 ms, spectral FWHM ≈ 4 Hz across 25–58 Hz. (`sim_broadband_transient.png`, `superlet_FWHM.csv`.)

### Not addressable with the current data (acknowledged limitations)

Electrode localization / olfactory-ROI specificity (no MNI here); OBE as the control cohort (different implant rationale); evoked-vs-induced separation (no ERP subtraction); surrogate nulls for the coupling MI (the continuous gamma-amplitude series was not retained per session); more than two non-responders. These are the priorities for a follow-up.
