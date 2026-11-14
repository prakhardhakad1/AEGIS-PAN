# What changed in v5.1 (mentor fixes)

1. **Export works in CI.** `scripts/export_web.py` pointed at a file that isn't in the repo
   (`artifacts/file_generation/ttl=63d/output/oncostrata_prism.html`) and crashed. It now uses
   `index.html` as its template and is idempotent (re-running gives byte-identical output).
2. **CI is real.** The static "CI passing" badge is replaced with the live GitHub Actions badge.
   CI also fails if committed pages or results are stale.
3. **Synthetic data is labelled as synthetic** in the README, every page (top banner), the exclusion
   flow, and `docs/NEGATIVE_RESULTS.md`.
4. **Withdrawn claims the code doesn't back:** pooled C 0.686, CPTAC 0.730, Moffitt 0.660, the
   "Moffitt-like external validation" button (now "Synthetic Demo Cohort"), lifelines/cmprsk
   cross-checks, the r = 0.9994 concordance line, and multimodal C values now marked as published
   SeNMo figures not run here.
5. **Every page shows the computed result:** pooled LOCO C = 0.594 (95% CI 0.557–0.631).
6. **NEGATIVE_RESULTS.md rewritten** from the pipeline's actual output (EB shrinkage hurt in 5 of 9
   cohorts, O/E 0.27–2.12, KICH too small, therapy effect not causal).

Still to do before claiming clinical results: fit on the real TCGA-CDR table, cross-check against
`cmprsk`/`scikit-survival`, and validate on a real external cohort. Also update the Devpost write-up,
which still quotes 0.730 / 0.660.
