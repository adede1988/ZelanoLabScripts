# Dupi Olfactory Gamma — Review Packet

## Design

Cohorts: Dupi intervention (S1/2/3) vs OBE controls; EEG excluded (no nasal electrode). Best macBP channel per recording = max peak flattened power (FOOOF) in 25-58 Hz. Superlet TFR (fractional adaptive, c1=3, order [3 30]), validated to 4e-14 vs reference faslt. Ridge via Frequency_ridge_tracking forward-backward tracker (ported to MATLAB) on per-breath within-frequency z (myChanZscore, no baseline). Per-breath FOOOF-lite (aperiodic exponent/offset + gamma peak presence). Session-level respiration-gamma coupling (Tort MI, preferred phase, resultant length, inhale/exhale ratio). Group stats aggregate to session level first.

## Inventory

On-disk finals:



|task          |cohort |  n|
|:-------------|:------|--:|
|O15           |Dupi   | 28|
|O15           |OBE    | 11|
|breathingTask |Dupi   | 29|
|breathingTask |OBE    |  7|
|cueTask       |Dupi   | 30|
|cueTask       |OBE    | 11|
|threshTask    |Dupi   | 30|
|threshTask    |OBE    |  5|

Missing (sheet->disk): 260504_Dupi_NMH_JA_2


## Behavioral & responder split

Per-session composite means by X-position:




|xpos    |  n| mean_comp|
|:-------|--:|---------:|
|Control | 19|     0.644|
|S1      | 12|     0.091|
|S2      | 12|     0.310|
|S3      |  6|     0.262|


Responder table:



|participant | nSess| firstSess| lastSess| composite_S1| composite_final|      delta|class                         |
|:-----------|-----:|---------:|--------:|------------:|---------------:|----------:|:-----------------------------|
|JA          |     2|         1|        2|    0.4559237|       0.3412058| -0.1147179|non-responder                 |
|PC          |     2|         1|        2|    0.4475183|       0.2803149| -0.1672034|non-responder                 |
|GH          |     3|         1|        3|    0.0967874|       0.1515851|  0.0547977|non-responder                 |
|DB          |     2|         1|        2|   -0.0489492|       0.1360609|  0.1850101|non-responder                 |
|TB          |     3|         1|        3|   -0.1037268|       0.0563964|  0.1601232|non-responder                 |
|JL          |     2|         1|        2|   -0.0749739|      -0.0094785|  0.0654954|non-responder                 |
|KS          |     3|         1|        3|    0.1965322|      -0.0730248| -0.2695570|non-responder                 |
|JN          |     2|         2|        3|    0.2407440|      -0.1043481| -0.3450921|non-responder                 |
|JH          |     3|         1|        3|    0.0370648|       0.8497849|  0.8127202|responder                     |
|AB          |     3|         1|        3|   -0.0816880|       0.6936784|  0.7753664|responder                     |
|BS          |     2|         1|        2|    0.0630853|       0.6190203|  0.5559350|responder                     |
|PD          |     2|         1|        2|    0.2263226|       0.4213628|  0.1950402|responder                     |
|DL          |     1|         1|        1|   -0.1160737|      -0.1160737|  0.0000000|unclassified (single session) |


## macBP gamma selection

Selected channel had FOOOF-detectable gamma peak in 59/150 recordings (39%). Mean non-selected channels with a detectable peak: 0.59.

Selected-channel spikeFrac summary (noise proxy): median 0.0001, >0.02 in 4 recordings.


## Group spectrograms

5 task-rows (cueTask, threshTask, O15, audiobook, focusedBreathing) x 5 columns (control, S2/3 responder, S1 responder, S2/3 non-responder, S1 non-responder). See figs/spectrograms_5x5.png. Group map = mean over breaths of per-breath single-trial within-frequency z (no baseline).


## Gamma measures — goodness ranking (KEY RESULT)

Goodness = |max separation effect| / (1 + |control CV|). Top measures:




|task             |metric                 | control_CV| sep_control_dupi| sep_resp_nonresp| recovery_rho_resp| goodness|
|:----------------|:----------------------|----------:|----------------:|----------------:|-----------------:|--------:|
|focusedBreathing |p1_rpowZ__mean         |      0.017|            -1.60|             0.25|             -0.36|     1.58|
|focusedBreathing |w5_rfreqRaw__mean      |      0.044|            -1.52|             0.59|             -0.25|     1.46|
|threshTask       |w4_rfreqGated__mean    |      0.053|             1.52|             0.39|             -0.44|     1.44|
|cueTask          |w4_rfreqPW__mean       |      0.071|             1.52|            -0.15|             -0.54|     1.42|
|threshTask       |w4_rfreqPW__mean       |      0.056|             1.46|             0.30|             -0.49|     1.38|
|O15              |maxBurstMs__mean       |      0.098|            -1.51|             0.01|              0.36|     1.38|
|cueTask          |w4_rfreqGated__mean    |      0.075|             1.48|            -0.07|             -0.51|     1.38|
|threshTask       |freqSpan__mean         |      0.068|            -1.45|            -0.26|             -0.76|     1.36|
|O15              |timeAboveMs_base__mean |      0.056|            -1.41|             0.60|              0.03|     1.33|
|threshTask       |peakFreq__mean         |      0.053|             1.40|             0.17|             -0.58|     1.33|
|O15              |dutyCycle__mean        |      0.123|            -1.46|             0.34|              0.39|     1.30|
|O15              |timeAboveMs__mean      |      0.123|            -1.46|             0.34|              0.39|     1.30|
|O15              |w4_rfreqPW__mean       |      0.033|             1.32|             0.00|             -0.42|     1.28|
|cueTask          |w4_rfreqRaw__mean      |      0.070|             1.35|            -0.25|             -0.54|     1.27|
|focusedBreathing |w3_rpowDb__mean        |      0.434|            -1.81|            -0.82|              0.60|     1.26|
|threshTask       |w4_rfreqRaw__mean      |      0.056|             1.32|             0.30|             -0.49|     1.25|
|focusedBreathing |w2_rpowDb__mean        |      0.437|            -1.79|            -0.92|              0.38|     1.24|
|threshTask       |w3_rfreqRaw__mean      |      0.030|             1.28|             0.30|             -0.54|     1.24|
|threshTask       |w3_rfreqPW__mean       |      0.038|             1.29|             0.44|             -0.49|     1.24|
|focusedBreathing |maxBurstMs__mean       |      0.061|            -1.30|             0.13|              0.45|     1.23|
|threshTask       |w3_rfreqGated__mean    |      0.043|             1.28|             0.49|             -0.39|     1.23|
|O15              |w4_rfreqRaw__mean      |      0.033|             1.27|             0.02|             -0.49|     1.23|
|threshTask       |w5_aucZ__mean          |      0.067|            -1.30|             0.02|              0.54|     1.22|
|threshTask       |w5_rpowZ__mean         |      0.067|            -1.30|             0.02|              0.54|     1.22|
|focusedBreathing |w1_rpowDb__mean        |      0.425|            -1.71|            -0.85|              0.38|     1.20|
|threshTask       |p3_rpowZ__mean         |      0.091|            -1.30|            -0.13|              0.47|     1.19|
|focusedBreathing |p1_rpowDb__mean        |      0.435|            -1.69|            -0.91|              0.38|     1.18|
|focusedBreathing |w5_rfreqPW__mean       |      0.044|            -1.22|             0.59|             -0.25|     1.17|
|focusedBreathing |p4_rpowDb__mean        |      0.422|            -1.66|            -0.85|              0.66|     1.17|
|threshTask       |apExp__mean            |      0.332|             0.79|             1.53|             -0.27|     1.15|


## Coupling summary

Session-level respiration-gamma coupling per recording in out/gamma/coupling/*.csv (coup_MI, coup_prefPhaseRad, coup_resultantLen, coup_rayleighP, coup_inhExhRatio).

