# Negative Results, Failure Modes & Limitations

All numbers below are produced by `python3 scripts/run_all.py` on the bundled
**synthetic** TCGA-shaped data (see "0. Data" below). Re-running the script
regenerates them; nothing here is typed by hand.

## 0. The data is synthetic

`data/raw/tcga_cdr_11160_raw.csv` and `data/tcga_cdr_4954_harmonized.csv` are
synthetic records shaped like TCGA-CDR (cohort sizes match the published
resource). They are **not** real TCGA patients:

- Patient IDs (`TCGA-LUAD-00001`, `TCGA-EXCL-11158`) are not real TCGA barcodes.
- The 5,742 non-target ("OTHER") rows are identical placeholders.
- TCGA-CDR does not carry the `comorbidity` or `adjuvant_tx` columns used here.

Every result in this repository therefore demonstrates that the **pipeline
works**, not that the **model is clinically valid**. Earlier versions of this
project quoted a pooled C of 0.686 and external validation on CPTAC (C = 0.730)
and a Moffitt cohort (C = 0.660). Those figures were not produced by this code
and have been withdrawn.

## 1. Discrimination is modest

9-fold leave-one-cohort-out pooled C = **0.594** (95% CI 0.557–0.631).
Held-out cohort C ranges from **0.543 (COAD)** to **0.635 (LGG)**. Several
per-cohort confidence intervals include 0.5.

## 2. Calibration does not transport across cohorts

Held-out 5-year O/E for cancer death ranges from **0.267 (KICH)** to
**2.117 (GBM)**. BRCA (0.383), OV (1.954) and GBM (2.117) are badly
miscalibrated before local recalibration.

## 3. Empirical-Bayes shrinkage hurt in 5 of 9 cohorts

5-year calibration error (percentage points) after a local fix:

| Cohort | Uncalibrated | MLE | Empirical-Bayes |
|---|---|---|---|
| LUAD | 3.96 | 2.50 | 2.98 |
| KIRC | 11.08 | 8.69 | 9.26 |
| LGG | 1.52 | 0.96 | 1.10 |
| BRCA | 27.44 | 0.41 | 8.70 |
| COAD | 4.46 | 15.29 | 9.59 |
| HNSC | 13.26 | 5.71 | 7.83 |
| GBM | 37.57 | 2.39 | 5.45 |
| OV | 32.23 | 3.43 | 8.80 |
| KICH | 32.61 | 6.90 | 17.66 |

EB beat plain MLE only in COAD. In COAD, MLE recalibration made things worse
than doing nothing.

## 4. KICH is too small

N = 113 with fewer than 20 events (EPV < 10). Held-out C CI 0.40–0.99 and
calibration slope 1.63. Results for KICH should not be interpreted.

## 5. Therapy effect is not causal

The adjuvant-therapy coefficient (aSHR ≈ 0.74) is an adjusted association in
synthetic observational data. It is subject to confounding by indication and
must not be read as a treatment effect.

## 6. Not yet done

- Fit on the real TCGA-CDR supplementary table (Liu et al., Cell 2018).
- Cross-check the Newton-Raphson Fine-Gray fit against R `cmprsk` or
  `scikit-survival`.
- Any genuine external validation cohort.
