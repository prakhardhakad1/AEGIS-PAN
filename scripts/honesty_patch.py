"""
Honesty patch applied to every generated HTML page.

The bundled CSV files are SYNTHETIC records shaped like TCGA-CDR. They are not
downloaded TCGA data. This module makes every page say so, removes claims the
code does not back (lifelines/cmprsk cross-checks, "Moffitt" external
validation), and writes the LOCO numbers the pipeline actually computed.
"""

BANNER = (
    '<div role="note" style="background:#fef3c7;color:#78350f;border-bottom:2px solid #f59e0b;'
    'padding:10px 16px;font:600 13px/1.4 system-ui,sans-serif;text-align:center">'
    'Data notice: this demo runs on SYNTHETIC patient records shaped like TCGA-CDR '
    '(9 cohorts, N = 4,954). They are not real TCGA patients. All model fits and '
    'leave-one-cohort-out results shown here are computed from that synthetic data '
    '(pooled C = {c:.3f}, 95% CI {lo:.3f}&ndash;{hi:.3f}). Research and demonstration use only.'
    '</div>'
)


def apply(html: str, loco: dict) -> str:
    c, (lo, hi), (plo, phi) = loco["pooled_c"], loco["ci"], loco["pi"]

    reps = [
        # Pooled LOCO summary that was still hardcoded
        ("Pooled C-index 0.686 (95% PI 0.482–0.890)",
         f"Pooled C-index {c:.3f} (95% PI {plo:.3f}–{phi:.3f})"),
        # Synthetic demo cohort is not external validation
        ("'Moffitt-Like External Validation Cohort'",
         "'Synthetic Demo Cohort (not external validation)'"),
        # Cross-checks that the code does not perform
        ("Python/lifelines Cross-Check &amp; 150-Profile Decimal Parity",
         "Python &amp; Node.js 150-Profile Decimal Parity"),
        ("Fixed Reference Vectors + lifelines/cmprsk Cross-Check",
         "Fixed Reference Vectors + Node.js Parity Check"),
        (", cross-checked against <code class=\"font-mono\">lifelines</code> Aalen-Johansen / CoxPH and "
         "<code class=\"font-mono\">cmprsk</code> subdistribution hazards,",
         ","),
        ("fits Royston-Parmar & Fine-Gray on TCGA-CDR and cross-checks against `lifelines`/`scipy`",
         "fits Royston-Parmar & Fine-Gray on the synthetic TCGA-shaped data with a hand-written Newton-Raphson solver"),
        # Data provenance wording
        ("9 Real TCGA-CDR Cohorts", "9 Synthetic TCGA-Shaped Cohorts"),
        ("Grounded in NCI GDC TCGA-CDR (", "Cohort sizes modelled on NCI GDC TCGA-CDR, data synthetic ("),
        ("N = 4,954 TCGA-CDR Analysed Cohort", "N = 4,954 Synthetic TCGA-Shaped Cohort"),
        ("Exact accounting from total TCGA-CDR resource to the 9 solid-tumor cohorts",
         "Exclusion steps applied to the synthetic 11,160-row file"),
        ("Total TCGA Pan-Cancer Clinical Data Resource (TCGA-CDR, 33 tumor types)",
         "Synthetic TCGA-CDR-shaped file (11,160 rows)"),
        ("of 4,954 TCGA-CDR)", "of 4,954 synthetic)"),
        ("1.00× (TCGA-CDR Ref)", "1.00× (Pooled Ref)"),

        ("Fine-Gray sHR vs R cmprsk / lifelines concordance r = 0.9994.",
         "No external cmprsk or lifelines comparison has been run yet."),
        (" Cross-checked against lifelines/cmprsk.`;", "`;"),
        ("'Moffitt-Like External Validation CSV (n=240)'",
         "'Synthetic Demo CSV (n=240, not external validation)'"),
        ("Pooled LOCO C = 0.686 (95% PI 0.482 to 0.890)",
         f"Pooled LOCO C = {c:.3f} (95% PI {plo:.3f} to {phi:.3f})"),
        ("${x(0.686)}", f"${{x({c:.3f})}}"),
        ('<td class="py-2 px-2 text-right text-teal-500">0.686 (Pooled)</td>',
         f'<td class="py-2 px-2 text-right text-teal-500">{c:.3f} (Pooled)</td>'),
        ("{ label: '1. Clinical Staging + Age', c: 0.686,",
         f"{{ label: '1. Clinical Staging + Age (computed here)', c: {c:.3f},"),
        ("'2. + Transcriptomics (3-Modal)'", "'2. + Transcriptomics (published SeNMo, not run here)'"),
        ("'3. + Proteomics + Mutations (5-Modal)'", "'3. + Proteomics + Mutations (published SeNMo, not run here)'"),
        ("'4. Full 6-Modal Pan-Cancer Ensemble'", "'4. 6-Modal Ensemble (published SeNMo, not run here)'"),
        ("(Pooled C = 0.686, 95% PI 0.482–0.890)",
         f"(Pooled C = {c:.3f}, 95% PI {plo:.3f}–{phi:.3f}, synthetic data)"),
    ]
    for old, new in reps:
        html = html.replace(old, new)

    # Any remaining stale numbers from earlier versions
    html = html.replace("N = 4,854", "N = 4,954").replace("displayed as 4,954 instead of 4,954", "displayed as 4,854 instead of 4,954")

    banner = BANNER.format(c=c, lo=lo, hi=hi)
    if 'role="note"' not in html.split("<body", 1)[1][:2000]:
        i = html.index("<body")
        j = html.index(">", i) + 1
        html = html[:j] + "\n" + banner + "\n" + html[j:]
    return html
