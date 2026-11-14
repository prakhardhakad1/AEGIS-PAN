# AEGIS-PAN v5.1 — Pan-Cancer Competing-Risks Workbench

[![CI](https://github.com/prakhardhakad1/AEGIS-PAN/actions/workflows/ci.yml/badge.svg)](https://github.com/prakhardhakad1/AEGIS-PAN/actions/workflows/ci.yml)

**Research and demonstration use only. Not a medical device.**

> **Data notice.** The bundled CSVs are **synthetic** records shaped like
> TCGA-CDR (9 cohorts, N = 4,954 after exclusions from 11,160 rows). They are
> not real TCGA patients. Every number in this repo is computed from that
> synthetic data and shows that the pipeline works, not that the model is
> clinically valid. See `docs/NEGATIVE_RESULTS.md`.

## What it does

Splits survival by cause of death: cancer death (CIF₁) vs. competing
non-cancer death (CIF₂), so older and sicker patients are not told their
cancer risk is higher than it is.

## Reproduce everything

```bash
pip install numpy scipy pandas
python3 scripts/run_all.py          # exclusions, Newton-Raphson Fine-Gray + Royston-Parmar fit, 9-fold LOCO
python3 scripts/export_web.py       # writes fitted results into index.html and the 7 pages
PYTHONPATH=. python3 -m unittest discover -s tests -p "test_*.py" -v
node tests/parity_js_check.mjs      # Python vs. JavaScript, 150 profiles
```

CI runs all four steps on every push.

## Results (synthetic data)

| Measure | Value |
|---|---|
| Pooled leave-one-cohort-out C | **0.594** (95% CI 0.557–0.631) |
| Held-out C by cohort | 0.543 (COAD) to 0.635 (LGG) |
| Held-out 5-year O/E by cohort | 0.27 (KICH) to 2.12 (GBM) |
| Python ↔ JS max difference | ~5e-9 across 150 profiles |
| Age-band overstatement of cancer death by all-cause models (5 y) | +1.6 pts (<50) to +31.4 pts (≥75) |

Withdrawn from earlier versions: pooled C 0.686, CPTAC C 0.730, Moffitt C 0.660,
and lifelines/cmprsk cross-checks. None were produced by this code.

## Layout

- `data/raw/tcga_cdr_11160_raw.csv`, `data/tcga_cdr_4954_harmonized.csv` — synthetic TCGA-shaped data
- `aegis_pan/survival_engine.py` — exclusions, Aalen-Johansen, Newton-Raphson Fine-Gray, Royston-Parmar spline, Uno IPCW C, calibration slope, IPCW Brier
- `scripts/run_all.py` — full pipeline, writes `results/*.json`
- `scripts/export_web.py` + `scripts/honesty_patch.py` — inject results into the pages and keep data-provenance wording accurate
- `web/js/engine.mjs`, `tests/parity_js_check.mjs` — JavaScript engine and Node parity check
- `docs/NEGATIVE_RESULTS.md` — limitations and what didn't work

## Next steps

1. Swap in the real TCGA-CDR supplementary table (Liu et al., Cell 2018, PMC6066282). It has no comorbidity or adjuvant columns, so those predictors must be dropped.
2. Cross-check the Fine-Gray fit against R `cmprsk` or `scikit-survival`.
3. External validation on a real independent cohort.
