# olfactoryHRV — respiration → HRV coupling vs olfactory recovery (Dupi cohort)

Does the cardiac response to breathing change as olfactory function recovers on dupilumab?
This folder holds the whole analysis: extraction from the preprocessed `.mat` finals, per-session
coupling measures, change scores against olfactory outcome, and the two generated reports.

**Headline.** Over the matched session 1 → 2 interval (~1 month, n = 7 patients), change in the
cardiac response to breath *depth* tracks change in olfactory function (r = 0.85, ρ = 0.75). The
same holds for RSA measured at a fixed breath size (ρ = 0.89) and for bias-calibrated odor intensity
sensitivity (ρ = 0.71). General heart-rate variability measures do **not** track recovery — the
effect is specific to respiratory–cardiac coupling.

---

## Run it

```powershell
powershell -ExecutionPolicy Bypass -File run_all.ps1
```

Edit the two interpreter paths at the top of `run_all.ps1` if MATLAB or Python live elsewhere.
Python needs `h5py` and `numpy`; MATLAB needs the Signal Processing and Statistics toolboxes.
Pass `-keep` to reuse cached intermediates instead of re-extracting.

All paths come from **`ohrv_config.m`** — nothing else is machine-specific. Set `OHRV_WORK` to move
the scratch directory, or edit `P.dataRoot` to point at a different set of preprocessed sessions.

> **Always re-extract after a preprocessing run.** The `.mat` finals get rewritten in place. Reading
> them while a batch is still writing silently mixes two data versions — this happened once and
> invalidated a full set of results. `run_all.ps1` clears the cache by default for that reason.

---

## Pipeline

| Step | Script | Produces |
|---|---|---|
| 1 | `rsa_extract.m` | `work/<sess>_slim.mat` — per-breath `behDat` + the `RRint` channel |
| 2 | `get_beats.py` | `work/<sess>_beats.npz` — `heartBeats` sample indices, read via HDF5 without loading the data matrix |
| 3 | `rsa_olf.m` | `work/olfactory_scores.csv` — cue d′, threshold, O15 per session |
| 4 | `rsa_analyze.m` | `work/session_gains.csv` — per-session slopes with both surrogate nulls |
| 5 | `rsa_panel.m` | `work/panel_pairs.csv` — vagal metric panel over spontaneous blocks |
| 6 | `rsa_report.m` | `work/report_changes.csv` — change scores, reliability gate, inference |
| 7 | `thresh_analysis.m` | `work/thresh_changes.csv` — bias-calibrated intensity sensitivity |
| 8 | `rsa_figs.m`, `fig_grant*.m`, `fig_duration.m`, `fig_raw.m` | `figures/*.png` |
| 9 | `mkreport.py`, `mkgrant.py` | `reports/*.html` |

`qc/` holds diagnostics that are not part of the main run: `diag_excl.m` (why sessions drop out),
`diag_thresh.m` (rating-scale check), `probe_ecg.py` (which finals carry ECG), `rsa_secheck.m`
(standard errors three ways), `rsa_stab.m` (rank stability), `rsa_robust.m` (RMSSD vs MASD).

---

## The measures, and why

**Duration slope — do not use as a result.** Regressing within-breath heart-period range on breath
*duration* has a built-in positive slope: a longer breath is a longer window, more beats fall in it,
and the range of more samples is larger. Measured across 25 sessions the surrogate null runs **+0.44
to +1.03** — the size of the raw effect. It also fails the reliability gate (κ ≈ 1). It is computed
and reported for completeness only.

**Depth slope — the well-posed version.** Breath depth does not set the window length, so the
artifact does not apply. The measured null is **−0.011 to +0.022**, i.e. zero. Interpretation is
unambiguous: how much cardiac modulation a deeper breath recruits.

**`adjLogRSA`** — fitted log RSA at a fixed reference breath. A *level* rather than a slope, anchored
by a precisely estimated intercept. Arguably the cleanest single measure.

**Two surrogate nulls.** The circular shift keeps the RR series intact and breaks only its alignment
to the breaths — it retains real RSA, so it is conservative. The i.i.d. permuted-NN null rebuilds the
series from the same interval distribution with all temporal structure destroyed, giving the pure
order-statistics floor. Their difference separates "mechanical" from "oscillation present but
misaligned".

**Reliability gate (κ).** Before any outcome is touched, each metric is screened on
κ = mean SE² / var(Δ). κ > 0.5 means the between-subject spread is no larger than the measurement
noise, and the metric is dropped without ever being correlated with olfaction.

---

## Data handling that matters

- **Ectopic filter.** Three stages: physiological range, deviation from an 11-wide local median
  (wide enough that genuine RSA averages out of the reference), then a successive-difference rule at
  5 robust SDs. That third stage is essential — an ectopic beat gives a short interval followed by a
  compensatory long one, each within tolerance individually while their *difference* is impossible.
  RMSSD squares differences, so a handful of these dominate a session. One session had
  RMSSD/MASD = 5.99 and manufactured a spurious perfect correlation. `rmssdRatio` is carried as a QC
  readout; above ~3 means outlier-driven.
- **`noseMouth` is often empty**, not `"mouth"` — no session in this cohort carries an explicit
  `"mouth"`. Since respiration is recorded with a nasal cannula, a detected breath with a valid RSA
  value was necessarily nasal, so empty means unrecorded metadata. Filtering on `== "nose"` discards
  31–46% of breaths in 11 sessions. Filter on `~= "mouth"` instead.
- **Threshold intensity is a screen pixel coordinate** on an ~800 px track (floor 335, ceiling
  ~1135), not a rating value. Subtracting the blank-air trial removes both the pixel offset and any
  response bias; `thresh_analysis.m` expresses the result as a proportion of the track. Of the two
  contrasts, `high − air` carries the dynamic range; `med − air` sits near floor for most patients.
- **Vagal metrics use spontaneous blocks only** (`audio`, `focus`, `naturalFocus`), with block types
  matched within each subject's session pair and durations matched. Slow paced breathing inflates RSA
  and HF power through baroreflex resonance, and block composition differs across sessions.
- **Known bad file:** `250811_Dupi_NMH_TB_2_PEA_threshold_preproc.mat` fails to open with an HDF5
  `inflate()` error on repeated attempts — genuine corruption, needs re-preprocessing.

---

## Reports

- `reports/grant_report.html` — preliminary-data write-up, descriptive, minimal inference.
  Has a `CITATIONS TO ADD` marker where the literature grounding belongs.
- `reports/rsa_report.html` — full analysis with the reliability gate, multiplicity correction and
  every caveat. The max-statistic permutation across the metric family gives p = 0.159; the pooled
  analysis across all 13 intervals gives p = 0.729. Neither is significant, and the grant report is
  written accordingly as preliminary data rather than a result.

`work/` is git-ignored — it holds regenerable intermediates and multi-hundred-MB extracts.
