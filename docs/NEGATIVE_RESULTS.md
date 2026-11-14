# Negative Results, Failure Modes & Observational Confounding Audit

1. **Where Empirical-Bayes Shrinkage Under-Corrected (`TCGA-GBM`)**: In Glioblastoma (`N=596`, median follow-up `12.0` months, 5-year cancer mortality `78%`, held-out `O/E = 1.38`), shrinking toward the pan-cancer prior (`gamma = 1.00`) produced a 5-year calibration error of `3.4` percentage points compared to `2.9` percentage points for unregularized MLE.
2. **Small-Sample Spline Instability (`TCGA-KICH`, `N=113`)**: With fewer than 15 primary cancer deaths (`EPV < 10`), `TCGA-KICH` has a wide held-out 5-year `O/E` 95% confidence interval (`0.41 to 1.48`) and cannot support standalone spline fitting without pooling.
3. **Observational Confounding by Indication (`aSHR = 0.75`)**: Adjuvant therapy (`aSHR = 0.75`, `aRMTLd = -4.15` months) is derived from observational multivariable regression adjusted for age, T/N stage, and grade, not a randomized trial.
