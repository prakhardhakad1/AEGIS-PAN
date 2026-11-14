"""
Injects the exact fitted JSON artifacts (`results/model_artifact.json`,
`results/loco_validation.json`, `results/age_band_bias.json`,
`results/consort_counts.json`, and `results/parity_150_profiles.json`)
directly into the self-contained web application (`index.html`)
and generates the 7 multi-page HTML entry points so zero coefficients or LOCO metrics
are ever hand-typed.
"""
import json, math, os, shutil, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import honesty_patch

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RESULTS_DIR = os.path.join(REPO_ROOT, "results")

model_art = json.load(open(os.path.join(RESULTS_DIR, "model_artifact.json")))
loco_val = json.load(open(os.path.join(RESULTS_DIR, "loco_validation.json")))
age_bias = json.load(open(os.path.join(RESULTS_DIR, "age_band_bias.json")))
consort = json.load(open(os.path.join(RESULTS_DIR, "consort_counts.json")))
parity_150 = json.load(open(os.path.join(RESULTS_DIR, "parity_150_profiles.json")))

# Build the FPM parameter table directly from the Newton-Raphson fitted coefficients & standard errors
b1 = model_art["cause1_beta"]
se1 = model_art["cause1_se"]
rp1 = model_art["cause1_rp_gamma"]
b2 = model_art["cause2_beta"]
se2 = model_art["cause2_se"]
rp2 = model_art["cause2_rp_gamma"]

def fmt_hr_ci(beta, se):
    hr = math.exp(beta)
    lo = math.exp(beta - 1.96 * se)
    hi = math.exp(beta + 1.96 * se)
    return round(hr, 3), f"{lo:.2f} to {hi:.2f}", f"{beta:+.4f} (SE {se:.4f})"

fpm_rows = [
    {"term": f"Royston-Parmar Spline γ₀ (Intercept, knots={model_art['cause1_rp_knots']})", "shr": round(math.exp(rp1[0]), 4), "ci": f"γ₀ = {rp1[0]:.4f}", "chr": round(math.exp(rp2[0]), 4), "rmtl": f"γ₀,₂ = {rp2[0]:.4f}"},
    {"term": "Royston-Parmar Spline γ₁ · v₁(log t) (Log-Time Linear Basis)", "shr": round(math.exp(rp1[1]), 3), "ci": f"γ₁ = {rp1[1]:+.4f}", "chr": round(math.exp(rp2[1]), 3), "rmtl": f"γ₁,₂ = {rp2[1]:+.4f}"},
    {"term": "Royston-Parmar Spline γ₂ · v₂(log t) (Restricted Cubic Basis)", "shr": round(math.exp(rp1[2]), 3), "ci": f"γ₂ = {rp1[2]:+.4f}", "chr": round(math.exp(rp2[2]), 3), "rmtl": f"γ₂,₂ = {rp2[2]:+.4f}"},
]

term_labels = [
    ("Age per +10 yrs above Organ Mean", 0, True),
    ("Tumor Extent T2+ (step above T1)", 1, False),
    ("Tumor Extent T3+ (step above T2)", 2, False),
    ("Tumor Extent T4 (step above T3)", 3, False),
    ("Regional Nodes N1+ (step above N0)", 4, False),
    ("Extensive Nodes N2+ (step above N1)", 5, False),
    ("Histologic Grade (per step above G2)", 6, True),
    ("Adjuvant Therapy (Observational · aSHR)", 7, False),
]
for label, idx, has_c2 in term_labels:
    hr, ci_s, b_s = fmt_hr_ci(b1[idx], se1[idx])
    if idx == 0:
        hr2, _, _ = fmt_hr_ci(b2[0], se2[0])
    elif idx == 6:
        hr2, _, _ = fmt_hr_ci(b2[1], se2[1])
    else:
        hr2 = 1.000
    fpm_rows.append({
        "term": label,
        "shr": hr,
        "ci": ci_s,
        "chr": hr2,
        "rmtl": b_s,
    })

print("Constructed fitted parameter rows from Newton-Raphson output:", len(fpm_rows))

# Now patch `index.html` so it uses
# the exact fitted MODEL_ARTIFACT, LOCO_VALIDATION, AGE_BAND_BIAS, CONSORT_FLOW_STEPS,
# FPM_PARAMETERS, and all 150 PYTHON_REFERENCE_VECTORS, and evaluates the exact
# Royston-Parmar spline basis `roystonParmarBasisJS(logT, knots)`.

# The repo's own index.html is the template; every marker below survives a re-run,
# so the export is idempotent and works in CI.
html_path = os.path.join(REPO_ROOT, "index.html")
html = open(html_path, "r", encoding="utf-8").read()

# 1. Replace PAN_CANCER_COHORTS with the exact `loco_val["cohorts"]` computed by `scripts/run_all.py`
cohorts_js = "const PAN_CANCER_COHORTS = " + json.dumps(loco_val["cohorts"], indent=2) + ";\n"
cohorts_js += "    const LOCO_SUMMARY = " + json.dumps({k: v for k, v in loco_val.items() if k != "cohorts"}, indent=2) + ";\n"
cohorts_js += "    const MODEL_ARTIFACT = " + json.dumps(model_art, indent=2) + ";\n"
cohorts_js += "    const AGE_BAND_BIAS = " + json.dumps(age_bias, indent=2) + ";\n"

start_marker = "const PAN_CANCER_COHORTS = ["
end_marker = "const HORIZON_MONTHS = [0, 6, 12, 24, 36, 48, 60, 84, 96, 120];"
idx_s = html.index(start_marker)
idx_e = html.index(end_marker)

consort_formatted = [
    {"step": s["step"], "n": s["n"], "drop": None if i == 0 else (0 if i == len(consort) - 1 else s["excluded"])}
    for i, s in enumerate(consort)
]

replacement_data_block = (
    cohorts_js
    + "\n    const CONSORT_FLOW_STEPS = " + json.dumps(consort_formatted, indent=2) + ";\n"
    + "\n    const PYTHON_REFERENCE_VECTORS = " + json.dumps(parity_150, indent=2) + ";\n"
    + "\n    const FPM_PARAMETERS = " + json.dumps(fpm_rows, indent=2) + ";\n    "
)

html = html[:idx_s] + replacement_data_block + html[idx_e:]

# 2. Replace `computeCompetingRisksProfile` with the exact Royston-Parmar natural cubic spline + Fine-Gray evaluator
old_fn_start = "function computeCompetingRisksProfile(pt, customKappa) {"
old_fn_end = "// =========================================================================\n    // 4. 2-MINUTE GUIDED JUDGE DEMO FLOW"
# On a re-run the helpers from a previous export sit just above the function; replace them too.
idx_fs = html.index("function roystonParmarBasisJS(") if "function roystonParmarBasisJS(" in html else html.index(old_fn_start)
idx_fe = html.index(old_fn_end)

new_fn_block = """function roystonParmarBasisJS(logT, knots) {
      const kMin = knots[0], kInt = knots[1], kMax = knots[2];
      const lam = (kMax - kInt) / (kMax - kMin);
      const posInt = Math.max(0.0, logT - kInt);
      const posMin = Math.max(0.0, logT - kMin);
      const posMax = Math.max(0.0, logT - kMax);
      const v1 = logT;
      const v2 = (Math.pow(posInt, 3) - lam * Math.pow(posMin, 3) - (1.0 - lam) * Math.pow(posMax, 3)) / Math.pow(kMax - kMin, 2);
      return [v1, v2];
    }

    function evalRoystonParmarH0JS(m, rpGamma, rpKnots) {
      if (m <= 0) return 0.0;
      const basis = roystonParmarBasisJS(Math.log(m), rpKnots);
      const logH0 = rpGamma[0] + rpGamma[1] * basis[0] + rpGamma[2] * basis[1];
      return Math.exp(Math.max(-20.0, Math.min(10.0, logH0)));
    }

    function computeCompetingRisksProfile(pt, customKappa) {
      const kVal = typeof customKappa === 'number' ? customKappa : state.kappa;
      const cohort = getCohortMeta(pt.organ);
      const art = MODEL_ARTIFACT;
      const ageMean = art.organ_age_means[pt.organ];
      const ob = art.organ_baselines[pt.organ];
      const b1 = art.cause1_beta;
      const b2 = art.cause2_beta;

      const x1 = [
        (pt.age - ageMean) / 10.0,
        pt.T >= 2 ? 1.0 : 0.0,
        pt.T >= 3 ? 1.0 : 0.0,
        pt.T >= 4 ? 1.0 : 0.0,
        pt.N >= 1 ? 1.0 : 0.0,
        pt.N >= 2 ? 1.0 : 0.0,
        pt.grade - 2.0,
        pt.adjuvant ? 1.0 : 0.0
      ];
      const x2 = [
        (pt.age - 60.0) / 10.0,
        pt.comorb - 1.0
      ];

      const logShrAge = b1[0] * x1[0];
      const logShrT = b1[1] * x1[1] + b1[2] * x1[2] + b1[3] * x1[3];
      const logShrN = b1[4] * x1[4] + b1[5] * x1[5];
      const logShrGrade = b1[6] * x1[6];
      const logShrTherapy = b1[7] * x1[7];

      let eta1 = Math.log(kVal);
      let eta1NoTx = Math.log(kVal);
      for (let i = 0; i < b1.length; i++) {
        eta1 += b1[i] * x1[i];
        if (i < b1.length - 1) eta1NoTx += b1[i] * x1[i];
      }
      const eta2 = b2[0] * x2[0] + b2[1] * x2[1];

      const shr1 = Math.exp(eta1);
      const hr2 = Math.exp(eta2);

      const curve = art.horizons.map(m => {
        if (m === 0) return { m: 0, os: 1.0, cif1: 0.0, cif2: 0.0, cif1NoTx: 0.0, allCauseMort: 0.0 };
        const h1Rp = evalRoystonParmarH0JS(m, art.cause1_rp_gamma, art.cause1_rp_knots) * ob.h1_scale;
        const h2Rp = evalRoystonParmarH0JS(m, art.cause2_rp_gamma, art.cause2_rp_knots) * ob.h2_scale;

        const h1 = h1Rp * shr1;
        const h1NoTx = h1Rp * Math.exp(eta1NoTx);
        const h2 = h2Rp * hr2;
        const hTotal = h1 + h2;

        const os = Math.exp(-hTotal);
        const totalDeath = 1.0 - os;
        const cif1 = totalDeath * (h1 / Math.max(1e-12, hTotal));
        const cif2 = totalDeath * (h2 / Math.max(1e-12, hTotal));

        const osNoTx = Math.exp(-(h1NoTx + h2));
        const cif1NoTx = (1.0 - osNoTx) * (h1NoTx / Math.max(1e-12, h1NoTx + h2));

        return {
          m,
          os,
          cif1,
          cif2,
          cif1NoTx,
          allCauseMort: totalDeath
        };
      });

      let rmst10y = 0, rmtlCancer10y = 0, rmtlComp10y = 0;
      for (let i = 0; i < curve.length - 1; i++) {
        const dt = curve[i + 1].m - curve[i].m;
        rmst10y += 0.5 * (curve[i].os + curve[i + 1].os) * dt;
        rmtlCancer10y += 0.5 * (curve[i].cif1 + curve[i + 1].cif1) * dt;
        rmtlComp10y += 0.5 * (curve[i].cif2 + curve[i + 1].cif2) * dt;
      }

      return {
        pt,
        cohort,
        shr1,
        hr2,
        factors: [
          { label: `Age (${pt.age}y vs ${ageMean}y)`, val: logShrAge },
          { label: `Tumor Extent (T${pt.T})`, val: logShrT },
          { label: `Lymph Nodes (N${pt.N})`, val: logShrN },
          { label: `Biology Grade (G${pt.grade})`, val: logShrGrade },
          { label: `Adjuvant Tx (${pt.adjuvant ? 'Yes' : 'No'})`, val: logShrTherapy },
          { label: `Hospital Prior (γ=${kVal.toFixed(2)})`, val: Math.log(kVal) }
        ],
        curve,
        rmst10y,
        rmtlCancer10y,
        rmtlComp10y
      };
    }

    """

html = html[:idx_fs] + new_fn_block + html[idx_fe:]

# 3. Also update setCompLabComorb to use the computed `AGE_BAND_BIAS` from `results/age_band_bias.json`
old_comp_start = "function setCompLabComorb(horizonKey) {"
old_comp_end = "function setMcParam(key, val) {"
idx_cs = html.index(old_comp_start)
idx_ce = html.index(old_comp_end)

new_comp_block = """function setCompLabComorb(horizonKey) {
      state.compLabLevel = horizonKey;
      document.querySelectorAll('.clab-btn').forEach(b => b.classList.toggle('active', b.dataset.clab === horizonKey));

      const ageBands = (AGE_BAND_BIAS[horizonKey] || AGE_BAND_BIAS['5y']).map(b => ({
        band: b.band,
        n: b.n,
        cif1: b.cif1,
        cif2: b.cif2,
        allCause: b.allCause,
        bias: b.overstatement
      }));

      document.getElementById('comp-bias-tbody').innerHTML = ageBands.map(b => `
        <tr>
          <td class="py-2 font-bold">${b.band}</td>
          <td class="py-2 text-right">${b.n.toLocaleString()}</td>
          <td class="py-2 text-right text-teal-600 dark:text-teal-400 font-bold">${b.cif1.toFixed(2)}%</td>
          <td class="py-2 text-right text-amber-500 font-semibold">${b.cif2.toFixed(2)}%</td>
          <td class="py-2 text-right text-rose-500 font-bold">${b.allCause.toFixed(2)}%</td>
          <td class="py-2 text-right font-extrabold text-rose-500">+${b.bias.toFixed(2)} pts Overstated</td>
        </tr>
      `).join('');

      const w = 540, h = 215, padL = 175, padR = 55;
      const scale = (w - padL - padR) / 90.0;
      document.getElementById('comp-bias-svg').innerHTML = `
        <svg viewBox="0 0 ${w} ${h}" class="w-full h-full">
          ${ageBands.map((b, i) => {
            const y = 12 + i * 48;
            return `
              <text x="${padL - 8}" y="${y + 15}" text-anchor="end" class="fill-slate-500 font-mono text-[10px] font-bold">${b.band.slice(0, 22)}</text>
              <rect x="${padL}" y="${y}" width="${b.cif1 * scale}" height="12" rx="2" fill="#14B8A6" />
              <rect x="${padL}" y="${y + 15}" width="${b.allCause * scale}" height="12" rx="2" fill="#F43F5E" />
              <text x="${padL + b.allCause * scale + 6}" y="${y + 24}" class="fill-rose-500 font-mono text-[10px] font-bold">+${b.bias.toFixed(2)} pts overstatement</text>
            `;
          }).join('')}
        </svg>
      `;
    }

    """

html = html[:idx_cs] + new_comp_block + html[idx_ce:]

# 4. Update LOCO pooled summary display in HTML to match `loco_val`
html = html.replace("LOCO Pooled C: 0.686 (95% CI 0.630–0.742)", f"LOCO Pooled C: {loco_val['pooled_c']:.3f} (95% CI {loco_val['ci'][0]:.3f}–{loco_val['ci'][1]:.3f})")
html = html.replace("95% Prediction Interval: 0.482 to 0.890", f"95% Prediction Interval: {loco_val['pi'][0]:.3f} to {loco_val['pi'][1]:.3f}")
html = html.replace("Pooled RE (95% PI: 0.482–0.890)", f"Pooled RE (95% PI: {loco_val['pi'][0]:.3f}–{loco_val['pi'][1]:.3f})")
html = html.replace("<td class=\"py-1.5 text-right\">0.686</td>\n            <td class=\"py-1.5 text-right\">0.630–0.742</td>", f"<td class=\"py-1.5 text-right\">{loco_val['pooled_c']:.3f}</td>\n            <td class=\"py-1.5 text-right\">{loco_val['ci'][0]:.3f}–{loco_val['ci'][1]:.3f}</td>")

html = honesty_patch.apply(html, loco_val)
open(html_path, "w", encoding="utf-8").write(html)

for p in ["workbench", "transport", "censoring", "calibration", "gdc-live", "methodology", "judge-audit"]:
    page_html = html.replace(
        "currentPage: 'workbench',",
        f"currentPage: '{p}',"
    ).replace(
        "syncSidebarControls();\n      updateFromUI();",
        f"syncSidebarControls();\n      updateFromUI();\n      switchPage('{p}');"
    )
    with open(os.path.join(REPO_ROOT, f"{p}.html"), "w", encoding="utf-8") as out:
        out.write(page_html)

print("Injected fitted artifacts into HTML and synced 7 multi-page files.")
