# Review — Christina Zelano (in character)

**Overall:** Careful first pass with the right instincts (inhale-onset locking, per-breath single-trial z, FOOOF, session-level aggregation, validated superlet). But the central biomarker claim is not yet supported because the confounds I care about are carried as *columns* rather than *removed*. The descending inhale-locked gamma in control/responder spectrograms is real and encouraging, but nothing rules out airflow, electrode location, or two subjects.

**Major concerns (ranked):**
1. **Airflow/sniff confound named but never tested.** Inhale volume/duration/peak-flow are "carried as covariates" but never used. The ridge descends ~50→30 Hz over ~1500 ms — exactly the time course of decaying inspiratory flow. `w3_rfreqGated__mean` (#1, control CV 0.049) may read peak flow, not neural gamma. *Did the sniff change or the brain?*
2. **Responder/non-responder rests on n=2 + regression to mean.** Only JA, PC are non-responders; both started at near-control S1 composite (0.456, 0.448) with nowhere to go but down; responders were selected from the floor. d≈2 on n=2 is a coin.
3. **Wrong controls, no localization.** Controls = OBE (different disease/site); no MNI, no olfactory-ROI, no inside>outside specificity. Control 14–20 dB flattened gamma vs Dupi 2–5 dB is as consistent with placement/SNR as with olfaction.
4. **Gamma metric computed where there is no gamma** (39% peak rate); `w3_rfreqGated`≈41.9 Hz sits at band center → low CV may mean "pinned to center, uninformative."
5. **Evoked vs induced not separated; coupling has no null.** Ridge highest at onset = where the sniff ERP leaks in. MI/resultant/preferred-phase computed but no surrogate + no non-sinusoidality correction.
6. **Selection without inference** — ~40 metrics × 5 tasks ranked, no permutation/FDR/CV/CIs (winner's curse). Recovery (goal #3) barely present (rho 0.13–0.52).

**Actionable:** (1) regress airflow out per breath, re-rank; (2) airflow-matched spectrograms; (3) drop/gate the n=2 contrast, per-subject dots + LOO, rank on control-vs-Dupi + recovery; (4) MNI + olfactory-ROI + inside>outside; restrict to the 59 peak-present recordings; (5) evoked/induced split (subtract ERP or PLV); (6) 200+ phase-shuffled surrogates for MI + non-sinusoidality control; (7) permutation+FDR, bootstrap CIs, nested CV; (8) within-subject test-retest reliability (ICC) + per-patient sensitivity to change.

**Strong:** methodological hygiene (inhale-onset epoching, single-trial z, FOOOF, session aggregation, both z-windows, superlet 4e-14, goodBreath QC); `inhExhRatio` very low control CV; reproducible control+responder inhale-locked gamma across tasks.

**Verdict:** Trust no single measure as a biomarker yet. If forced: watch **O15 `w3_rfreq*` family** (low control CV ~0.045–0.05) and **`coup_inhExhRatio`** — but only after airflow regression and reframed as control-vs-Dupi, not the two-subject split. Convince me: ridge separation survives regressing out inhale volume/flow; contacts anatomically in the olfactory pathway with inside>outside; within-subject test-retest reliability + per-patient recovery tracking.
