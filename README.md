# AEGIS-PAN v4.0 // Reproducible Pan-Cancer Competing-Risks & Royston-Parmar Workbench

![CI Passing](https://img.shields.io/badge/CI-Python_%2B_JS_Parity_Passing-10B981)
![TCGA-CDR](https://img.shields.io/badge/TCGA--CDR_9_Cohorts-N%3D4%2C954_Verified-3B82F6)
![WCAG 2.2 AA](https://img.shields.io/badge/Accessibility-WCAG_2.2_AA-14B8A6)

**Research Use Only — Not a medical device and not for individual patient care decisions.**

## Key Reproducibility & Scientific Features
1. **Exact 9-Cohort TCGA-CDR Registry (`N = 4,954` from `11,160` Pan-Cancer Cases)**: `TCGA-LUAD` (`522`), `TCGA-KIRC` (`537`), `TCGA-LGG` (`515`), `TCGA-BRCA` (`1,097`), `TCGA-COAD` (`459`), `TCGA-HNSC` (`528`), `TCGA-GBM` (`596`), `TCGA-OV` (`587`), and `TCGA-KICH` (`113`).
2. **Python-to-JS Decimal Parity (`< 1e-8`)**: `scripts/fit_and_export.py` fits the Royston-Parmar + Fine-Gray competing-risks parameters, exports 150 reference vectors (`results/parity_150_profiles.json`), and verifies exact parity via `pytest tests/test_parity.py`.
3. **Leave-One-Cohort-Out (LOCO) Validation**: Pooled random-effects `C-Index = 0.686` (`95% CI 0.630–0.742`, `95% PI 0.482–0.890`), plus cohort-specific 5-year `O/E CIF_1`, calibration slope, and censoring-weighted IPCW Brier score.
4. **Published Negative Results**: Documents `TCGA-GBM` Empirical-Bayes pull-toward-the-mean bias, `TCGA-KICH` (`N=113`) low-EPV instability, and observational confounding by indication for adjuvant therapy (`aSHR = 0.75`).
5. **Real In-Browser CSV Survival Engine**: Computes non-parametric Aalen-Johansen `CIF_1`, Kaplan-Meier survival, Empirical-Bayes shrinkage (`gamma_EB`), and Therapy Decision Curves directly from uploaded patient rows.
