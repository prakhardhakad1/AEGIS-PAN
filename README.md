# AEGIS-PAN // Pan-Cancer Competing-Risks & Royston-Parmar Survival Workbench (v3.0)

**Multi-Organ NCI GDC TCGA-CDR Survival Transportability, Fine-Gray Competing-Risks (`sHR`), Adjuvant Therapy Counterfactuals & Empirical-Bayes Calibration**

## 100% Distinct Dataset, Model & Scope
- **Dataset**: Curated **NCI GDC Pan-Cancer Clinical Data Resource (TCGA-CDR, Liu et al., *Cell* / PMC6066282)** across **8 Solid-Tumor Organ Systems (`N = 4,841` primary cases from the `11,160`-case TCGA Pan-Cancer resource)**: `TCGA-LUAD` (`N=522`), `TCGA-KIRC` (`N=537`), `TCGA-LGG` (`N=515`), `TCGA-BRCA` (`N=1,097`), `TCGA-COAD` (`N=459`), `TCGA-HNSC` (`N=528`), `TCGA-GBM` (`N=596`), and `TCGA-OV` (`N=587`), plus external validation benchmarks from **CPTAC (`C = 0.730`)** and **Moffitt Cancer Center (`n = 108, C = 0.660, p < 0.05`, PMC12347193)**.
- **Mathematical Model**: **Royston-Parmar Flexible Parametric Survival Model (FPM, PMC7953595)** coupled with **Fine-Gray Subdistribution Hazard (`sHR`) Competing-Risks Decomposition** (`CIF_1` Primary Cancer Death vs. `CIF_2` Competing Non-Cancer Death) and **Restricted Mean Time Lost (`RMTL`, PMC8889377)**.
- **Treatment Counterfactuals**: Models **Adjuvant Radiotherapy / Systemic Therapy (`aSHR = 0.75`, 95% CI `0.66–0.85`; `aRMTLd = -4.15` months saved, CLEP.S521309)**.
- **State-Space Verification**: In-browser **43,200-Profile Competing-Risks Partition-of-Unity & Monotonicity Auditor**.
- **User-Friendly Dual Theme UI**: Explicit `☀️ Light Theme` and `🌙 Dark Theme` segmented switcher + `✨ Friendly Clinical Guide` vs. `📐 Biostatistical Formulas` toggle.
