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



|participant | composite_S1| composite_S2|usedAll3 |      delta|class                                     |
|:-----------|------------:|------------:|:--------|----------:|:-----------------------------------------|
|BS          |    0.0630853|    0.6190203|TRUE     |  0.5559350|responder                                 |
|TB          |   -0.1037268|    0.4231614|FALSE    |  0.5268882|responder                                 |
|JH          |    0.0370648|    0.3706548|TRUE     |  0.3335900|responder                                 |
|AB          |   -0.0816880|    0.2375768|TRUE     |  0.3192648|responder                                 |
|GH          |    0.0967874|    0.3088261|TRUE     |  0.2120387|responder                                 |
|PD          |    0.2263226|    0.4213628|FALSE    |  0.1950402|responder                                 |
|DB          |   -0.0489492|    0.1360609|TRUE     |  0.1850101|responder                                 |
|KS          |    0.1965322|    0.3478550|TRUE     |  0.1513228|responder                                 |
|JL          |   -0.0749739|   -0.0094785|TRUE     |  0.0654954|responder                                 |
|JA          |    0.4559237|    0.3412058|FALSE    | -0.1147179|non-responder                             |
|PC          |    0.4475183|    0.2803149|TRUE     | -0.1672034|non-responder                             |
|DL          |   -0.1160737|           NA|FALSE    |         NA|unclassified (missing S1 or S2 composite) |
|JN          |           NA|    0.2407440|FALSE    |         NA|unclassified (missing S1 or S2 composite) |


## macBP gamma selection

Selected channel had FOOOF-detectable gamma peak in 59/150 recordings (39%). Mean non-selected channels with a detectable peak: 0.59.

Selected-channel spikeFrac summary (noise proxy): median 0.0001, >0.02 in 4 recordings.


## Group spectrograms

5 task-rows (cueTask, threshTask, O15, audiobook, focusedBreathing) x 5 columns (control, S2/3 responder, S1 responder, S2/3 non-responder, S1 non-responder). See figs/spectrograms_5x5.png. Group map = mean over breaths of per-breath single-trial within-frequency z (no baseline).


## Gamma measures — goodness ranking (KEY RESULT)

Goodness = |max separation effect| / (1 + |control CV|). Top measures:




|task             |metric                 | control_CV| sep_control_dupi| sep_resp_nonresp| recovery_rho_resp| goodness|
|:----------------|:----------------------|----------:|----------------:|----------------:|-----------------:|--------:|
|O15              |w3_rfreqGated__mean    |      0.049|             0.94|             2.10|              0.24|     2.00|
|O15              |w3_rfreqRaw__mean      |      0.045|             1.09|             2.07|              0.22|     1.98|
|audiobook        |w5_aucZ__mean          |      0.048|            -0.41|            -2.03|              0.27|     1.93|
|audiobook        |w5_rpowZ__mean         |      0.048|            -0.41|            -2.03|              0.27|     1.93|
|O15              |w3_rfreqPW__mean       |      0.048|             1.10|             2.00|              0.18|     1.90|
|audiobook        |timeAboveMs_base__mean |      0.049|            -0.60|            -1.99|              0.52|     1.89|
|audiobook        |w3_rfreqRaw__var       |      0.253|             0.71|             2.22|             -0.33|     1.77|
|audiobook        |w3_rfreqPW__var        |      0.260|             0.76|             2.11|             -0.35|     1.68|
|audiobook        |nBursts__mean          |      0.158|            -0.02|            -1.94|             -0.10|     1.67|
|threshTask       |w5_aucZ__mean          |      0.067|            -1.30|            -1.74|              0.15|     1.63|
|threshTask       |w5_rpowZ__mean         |      0.067|            -1.30|            -1.74|              0.15|     1.63|
|focusedBreathing |p1_rpowZ__mean         |      0.017|            -1.63|             0.16|              0.15|     1.61|
|O15              |p2_rfreqRaw__mean      |      0.042|             0.26|             1.65|              0.14|     1.58|
|focusedBreathing |p2_aucZ__mean          |      0.201|             0.59|            -1.86|              0.32|     1.55|
|O15              |p2_rfreqPW__mean       |      0.047|             0.28|             1.62|              0.13|     1.55|
|audiobook        |p2_aucZ__mean          |      0.172|            -0.07|            -1.81|              0.04|     1.54|
|O15              |p2_rfreqGated__mean    |      0.051|             0.27|             1.61|             -0.01|     1.53|
|cueTask          |peakFreq__var          |      0.181|            -0.44|             1.80|             -0.05|     1.52|
|audiobook        |w3_rfreqGated__var     |      0.304|             0.78|             1.97|             -0.25|     1.51|
|audiobook        |dutyCycle__mean        |      0.141|            -0.70|            -1.72|              0.28|     1.50|
|audiobook        |timeAboveMs__mean      |      0.141|            -0.70|            -1.72|              0.28|     1.50|
|O15              |peakFreq__var          |      0.190|            -0.36|             1.79|             -0.09|     1.50|
|audiobook        |w4_rfreqGated__var     |      0.336|             0.21|             1.96|             -0.37|     1.47|
|focusedBreathing |nBursts__mean          |      0.024|             0.23|            -1.50|              0.21|     1.46|
|cueTask          |w2_rfreqGated__mean    |      0.042|            -0.36|             1.51|              0.10|     1.45|
|threshTask       |w4_rfreqGated__mean    |      0.053|             1.52|             0.46|             -0.12|     1.44|
|focusedBreathing |w5_rfreqPW__mean       |      0.044|            -1.19|             1.50|              0.08|     1.43|
|audiobook        |w3_aucZ__mean          |      0.095|            -0.99|            -1.56|              0.16|     1.43|
|audiobook        |w3_rpowZ__mean         |      0.095|            -0.99|            -1.56|              0.16|     1.43|
|focusedBreathing |w5_rfreqRaw__mean      |      0.044|            -1.49|             1.25|              0.08|     1.42|


## Coupling summary

Session-level respiration-gamma coupling per recording in out/gamma/coupling/*.csv (coup_MI, coup_prefPhaseRad, coup_resultantLen, coup_rayleighP, coup_inhExhRatio).

