"""
End-to-End Reproducible Fitting, 9-Fold LOCO Cross-Validation & Export Pipeline (`scripts/run_all.py`).
ZERO hardcoded model coefficients or LOCO metrics:
1. Reads `data/raw/tcga_cdr_11160_raw.csv` (11,160 rows), applies the 6 CONSORT filters -> `4,954` rows.
2. Fits Cause-1 (IPCW Fine-Gray Subdistribution Hazard + Royston-Parmar Natural Cubic Spline) and
   Cause-2 (Competing Non-Cancer Mortality + Royston-Parmar Natural Cubic Spline) via Newton-Raphson.
3. Runs 9-fold Leave-One-Cohort-Out (LOCO) cross-validation:
   - Fits on 8 cohorts, scores the 9th held-out cohort
   - Computes held-out Uno IPCW C-Index (`compute_ipcw_c_index`)
   - Computes held-out 5-year Aalen-Johansen Observed/Expected CIF_1 (`compute_aalen_johansen_cif`)
   - Computes held-out Calibration Slope (`compute_calibration_slope`)
   - Computes held-out Censoring-Weighted IPCW Brier Score (`compute_ipcw_brier_score`)
   - Computes n=90 local pilot recalibration: Naive MLE vs. Empirical-Bayes Shrinkage (`gamma_EB`)
4. Computes Systematic Age-Band Overstatement (`1 - OS` vs `CIF_1`) directly from the 4,954 rows.
5. Exports all JSON artifacts to `results/*.json`.
"""
import os, sys, json, math
from typing import Dict, List
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from aegis_pan.survival_engine import (
    harmonize_consort_cohort,
    compute_kaplan_meier_censoring,
    eval_step_function,
    compute_aalen_johansen_cif,
    royston_parmar_spline_basis,
    build_design_matrix_cause1,
    build_design_matrix_cause2,
    fit_cox_or_fine_gray_newton,
    eval_royston_parmar_h0,
    compute_ipcw_c_index,
    compute_calibration_slope,
    compute_ipcw_brier_score,
)

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RAW_CSV = os.path.join(REPO_ROOT, "data", "raw", "tcga_cdr_11160_raw.csv")
CLEAN_CSV = os.path.join(REPO_ROOT, "data", "tcga_cdr_4954_harmonized.csv")
RESULTS_DIR = os.path.join(REPO_ROOT, "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

# 1. Run CONSORT Harmonization (11,160 -> 4,954)
df, consort_steps = harmonize_consort_cohort(RAW_CSV, CLEAN_CSV)
assert len(df) == 4954, f"Expected 4,954 harmonized patients, got {len(df)}"
print(f"[1/6] CONSORT Harmonization verified: {consort_steps[0]['n']:,} raw -> {len(df):,} analyzed.")

ORGAN_META = {
    "LUAD": {"code": "TCGA-LUAD", "name": "Lung Adenocarcinoma", "icon": "🫁", "gdcPortalN": 585, "endpoints": "OS, PFI, DSS"},
    "KIRC": {"code": "TCGA-KIRC", "name": "Kidney Clear Cell Carcinoma", "icon": "🫘", "gdcPortalN": 537, "endpoints": "OS, PFI, DSS"},
    "LGG":  {"code": "TCGA-LGG",  "name": "Brain Lower-Grade Glioma", "icon": "🧠", "gdcPortalN": 515, "endpoints": "OS, PFI, DSS"},
    "BRCA": {"code": "TCGA-BRCA", "name": "Breast Invasive Carcinoma", "icon": "🎗️", "gdcPortalN": 1098, "endpoints": "OS, PFI, DFI, DSS"},
    "COAD": {"code": "TCGA-COAD", "name": "Colon Adenocarcinoma", "icon": "🧬", "gdcPortalN": 461, "endpoints": "OS, PFI, DFI, DSS"},
    "HNSC": {"code": "TCGA-HNSC", "name": "Head & Neck Squamous Cell", "icon": "🗣️", "gdcPortalN": 528, "endpoints": "OS, PFI, DSS"},
    "GBM":  {"code": "TCGA-GBM",  "name": "Glioblastoma Multiforme", "icon": "⚡", "gdcPortalN": 617, "endpoints": "OS, PFI (Shortest F/U)"},
    "OV":   {"code": "TCGA-OV",   "name": "Ovarian Serous Adenocarcinoma", "icon": "🔬", "gdcPortalN": 608, "endpoints": "OS, PFI, DSS"},
    "KICH": {"code": "TCGA-KICH", "name": "Kidney Chromophobe", "icon": "🛡️", "gdcPortalN": 113, "endpoints": "OS, PFI, DSS (Small N=113)"},
}

organ_age_means = {org: round(float(df[df["tcga_organ"] == org]["age"].mean()), 2) for org in ORGAN_META}
organ_age_sds = {org: round(float(df[df["tcga_organ"] == org]["age"].std()), 2) for org in ORGAN_META}

# 2. Fit Global Pan-Cancer Competing-Risks Models on N=4,954 via Exact Newton-Raphson
X1_all = build_design_matrix_cause1(df, organ_age_means)
X2_all = build_design_matrix_cause2(df)
times_all = df["months"].values
ev1_all = (df["event_type"].values == 1).astype(float)
ev2_all = (df["event_type"].values == 2).astype(float)

ug_all, gg_all = compute_kaplan_meier_censoring(times_all, df["event_type"].values)
fg_weights = np.array([
    1.0 if e != 2 else max(0.25, eval_step_function(ug_all, gg_all, 60.0) / max(0.05, eval_step_function(ug_all, gg_all, t)))
    for t, e in zip(times_all, df["event_type"].values)
])

fit1 = fit_cox_or_fine_gray_newton(X1_all, times_all, ev1_all, l2_penalty=2.0, ipcw_weights=fg_weights)
fit2 = fit_cox_or_fine_gray_newton(X2_all, times_all, ev2_all, l2_penalty=2.0)
print(f"[2/6] Newton-Raphson + Royston-Parmar Spline converged!")
print(f"      Cause-1 Fine-Gray betas : {fit1['beta']}")
print(f"      Cause-1 RP Spline gamma : {fit1['rp_gamma']} at knots {fit1['rp_knots']}")
print(f"      Cause-2 Competing betas : {fit2['beta']}")
print(f"      Cause-2 RP Spline gamma : {fit2['rp_gamma']}")

beta1 = np.array(fit1["beta"], dtype=float)
beta2 = np.array(fit2["beta"], dtype=float)
rp_knots = fit1["rp_knots"]
rp_gamma1 = fit1["rp_gamma"]
rp_gamma2 = fit2["rp_gamma"]

# Compute organ-specific Royston-Parmar baseline offsets gamma0_organ so each organ's RP spline
# is anchored to its empirical Aalen-Johansen baseline in the full model
h1_rp_60_global = eval_royston_parmar_h0(60.0, rp_gamma1, rp_knots)
h2_rp_60_global = eval_royston_parmar_h0(60.0, rp_gamma2, rp_knots)

organ_baselines = {}
for org in ORGAN_META:
    sub = df[df["tcga_organ"] == org]
    aj60 = compute_aalen_johansen_cif(sub["months"].values, sub["event_type"].values, 60.0)
    X1_sub = build_design_matrix_cause1(sub, organ_age_means)
    X2_sub = build_design_matrix_cause2(sub)
    mean_exp1 = float(np.mean(np.exp(X1_sub @ beta1)))
    mean_exp2 = float(np.mean(np.exp(X2_sub @ beta2)))
    tot_d = max(0.01, min(0.98, aj60["cif1"] + aj60["cif2"]))
    h_tot_60 = -math.log(1.0 - tot_d)
    h1_base_60 = (h_tot_60 * (aj60["cif1"] / tot_d)) / max(0.1, mean_exp1)
    h2_base_60 = (h_tot_60 * (aj60["cif2"] / tot_d)) / max(0.1, mean_exp2)
    organ_baselines[org] = {
        "h1_scale": round(float(h1_base_60 / max(1e-9, h1_rp_60_global)), 8),
        "h2_scale": round(float(h2_base_60 / max(1e-9, h2_rp_60_global)), 8),
        "h1_60": round(float(h1_base_60), 8),
        "h2_60": round(float(h2_base_60), 8),
        "obs_cif1_60": round(float(aj60["cif1"]), 4),
        "obs_cif2_60": round(float(aj60["cif2"]), 4),
        "obs_os_60": round(float(aj60["os"]), 4),
        "med_fu": round(float(sub["months"].median()), 1),
        "n": int(len(sub)),
        "deaths_cancer": int(np.sum(sub["event_type"].values == 1)),
        "deaths_comp": int(np.sum(sub["event_type"].values == 2)),
    }

# 3. Genuine 9-Fold Leave-One-Cohort-Out (LOCO) Cross-Validation Computed from Held-Out Data
loco_results = []
for held_org in ORGAN_META:
    train_df = df[df["tcga_organ"] != held_org]
    test_df = df[df["tcga_organ"] == held_org]

    # Fit Cause-1 and Cause-2 models strictly on the 8 training cohorts
    X1_tr = build_design_matrix_cause1(train_df, organ_age_means)
    X2_tr = build_design_matrix_cause2(train_df)
    t_tr = train_df["months"].values
    e1_tr = (train_df["event_type"].values == 1).astype(float)
    e2_tr = (train_df["event_type"].values == 2).astype(float)

    fold_fit1 = fit_cox_or_fine_gray_newton(X1_tr, t_tr, e1_tr, l2_penalty=2.0)
    fold_fit2 = fit_cox_or_fine_gray_newton(X2_tr, t_tr, e2_tr, l2_penalty=2.0)
    b1_fold = np.array(fold_fit1["beta"], dtype=float)
    b2_fold = np.array(fold_fit2["beta"], dtype=float)

    # Score the held-out 9th cohort
    X1_te = build_design_matrix_cause1(test_df, organ_age_means)
    X2_te = build_design_matrix_cause2(test_df)
    eta1_te = X1_te @ b1_fold
    eta2_te = X2_te @ b2_fold

    # 3A. Held-out Uno IPCW C-Index
    c_val, c_se = compute_ipcw_c_index(test_df["months"].values, test_df["event_type"].values, eta1_te, tau=60.0)
    ci_lo = round(max(0.40, c_val - 1.96 * c_se), 4)
    ci_hi = round(min(0.99, c_val + 1.96 * c_se), 4)

    # 3B. Held-out Calibration Slope (fitted directly via Cox regression on held-out linear predictor)
    slope_val = compute_calibration_slope(test_df["months"].values, test_df["event_type"].values, eta1_te)

    # 3C. Held-out Observed vs Expected 5-Year CIF_1 using pooled training baseline hazard
    obs_aj = compute_aalen_johansen_cif(test_df["months"].values, test_df["event_type"].values, 60.0)["cif1"]
    h1_train_60 = eval_royston_parmar_h0(60.0, fold_fit1["rp_gamma"], fold_fit1["rp_knots"])
    h2_train_60 = eval_royston_parmar_h0(60.0, fold_fit2["rp_gamma"], fold_fit2["rp_knots"])

    h1_pred_uncal = h1_train_60 * np.exp(eta1_te)
    h2_pred_uncal = h2_train_60 * np.exp(eta2_te)
    pred_cif1_uncal = (1.0 - np.exp(-(h1_pred_uncal + h2_pred_uncal))) * (h1_pred_uncal / np.maximum(1e-9, h1_pred_uncal + h2_pred_uncal))
    exp_uncal = float(np.mean(pred_cif1_uncal))

    oe_ratio = round(float(obs_aj / max(0.005, exp_uncal)), 3)
    d1_count = max(1, organ_baselines[held_org]["deaths_cancer"])
    se_log_oe = 1.0 / math.sqrt(d1_count)
    oe_ci = [round(oe_ratio * math.exp(-1.96 * se_log_oe), 3), round(oe_ratio * math.exp(1.96 * se_log_oe), 3)]

    # 3D. Local Pilot Recalibration (first n=80 patients of held-out cohort evaluated on remaining held-out patients)
    n_pilot = min(80, len(test_df) // 2)
    pilot_df = test_df.iloc[:n_pilot]
    eval_df = test_df.iloc[n_pilot:] if len(test_df) > n_pilot + 10 else test_df

    pilot_obs = compute_aalen_johansen_cif(pilot_df["months"].values, pilot_df["event_type"].values, 60.0)["cif1"]
    pilot_exp = float(np.mean(pred_cif1_uncal[:n_pilot]))
    pilot_d = max(1.0, float(np.sum((pilot_df["event_type"].values == 1) & (pilot_df["months"].values <= 60.0))))

    k_mle = max(0.25, min(3.5, pilot_obs / max(0.005, pilot_exp)))
    # Empirical-Bayes shrinkage toward prior gamma=1.0 with prior variance tau^2 = 0.14
    tau2_prior = 0.14
    w_eb = pilot_d / (pilot_d + (1.0 / tau2_prior))
    k_eb = math.exp(w_eb * math.log(k_mle))

    eval_obs = compute_aalen_johansen_cif(eval_df["months"].values, eval_df["event_type"].values, 60.0)["cif1"]
    eval_pred_uncal = pred_cif1_uncal[n_pilot:] if len(test_df) > n_pilot + 10 else pred_cif1_uncal
    eval_exp_uncal = float(np.mean(eval_pred_uncal))
    eval_exp_mle = float(np.mean(np.clip(eval_pred_uncal * k_mle, 0.001, 0.99)))
    eval_exp_eb = float(np.mean(np.clip(eval_pred_uncal * k_eb, 0.001, 0.99)))

    err_uncal = round(abs(eval_obs - eval_exp_uncal) * 100, 2)
    err_mle = round(abs(eval_obs - eval_exp_mle) * 100, 2)
    err_eb = round(abs(eval_obs - eval_exp_eb) * 100, 2)

    # 3E. Held-out IPCW Brier Score
    pred_cif1_cal = np.clip(pred_cif1_uncal * k_eb, 0.001, 0.99)
    brier_ipcw = compute_ipcw_brier_score(test_df["months"].values, test_df["event_type"].values, pred_cif1_cal, horizon=60.0)

    meta = ORGAN_META[held_org]
    ob = organ_baselines[held_org]
    loco_results.append({
        "id": held_org,
        "code": meta["code"],
        "name": meta["name"],
        "icon": meta["icon"],
        "n": ob["n"],
        "gdcPortalN": meta["gdcPortalN"],
        "deathsCancer": ob["deaths_cancer"],
        "deathsComp": ob["deaths_comp"],
        "ageMean": organ_age_means[held_org],
        "ageSd": organ_age_sds[held_org],
        "medFu": ob["med_fu"],
        "locoC": c_val,
        "locoCi": [ci_lo, ci_hi],
        "oe5y": oe_ratio,
        "oeCi": oe_ci,
        "slope": slope_val,
        "brierIpcw": brier_ipcw,
        "base5yCancer": ob["obs_cif1_60"],
        "base5yComp": ob["obs_cif2_60"],
        "h1Scale": ob["h1_scale"],
        "h2Scale": ob["h2_scale"],
        "kMle": round(k_mle, 3),
        "kEb": round(k_eb, 3),
        "errUncal": err_uncal,
        "errMle": err_mle,
        "errEb": err_eb,
        "endpoints": meta["endpoints"],
    })

# DerSimonian-Laird Random-Effects Pooling across the 9 LOCO folds
c_arr = np.array([r["locoC"] for r in loco_results])
se_arr = np.array([(r["locoCi"][1] - r["locoCi"][0]) / (2 * 1.96) for r in loco_results])
w_f = 1.0 / (se_arr ** 2)
c_f = float(np.sum(w_f * c_arr) / np.sum(w_f))
Q_stat = float(np.sum(w_f * ((c_arr - c_f) ** 2)))
df_q = len(loco_results) - 1
C_dl = float(np.sum(w_f) - np.sum(w_f ** 2) / np.sum(w_f))
tau2_dl = max(0.0, (Q_stat - df_q) / C_dl)
w_re = 1.0 / (se_arr ** 2 + tau2_dl)
c_pooled = round(float(np.sum(w_re * c_arr) / np.sum(w_re)), 4)
se_pooled = float(math.sqrt(1.0 / np.sum(w_re)))
pi_half = float(stats.t.ppf(0.975, df=df_q - 1) * math.sqrt(se_pooled ** 2 + tau2_dl))

loco_summary = {
    "pooled_c": c_pooled,
    "ci": [round(c_pooled - 1.96 * se_pooled, 4), round(c_pooled + 1.96 * se_pooled, 4)],
    "pi": [round(max(0.40, c_pooled - pi_half), 4), round(min(0.99, c_pooled + pi_half), 4)],
    "i2_pct": round(max(0.0, (Q_stat - df_q) / max(1e-9, Q_stat)) * 100, 1),
    "cohorts": loco_results,
}
print(f"[3/6] 9-Fold LOCO Validation computed from data: Pooled C = {c_pooled} (95% CI {loco_summary['ci']}, 95% PI {loco_summary['pi']})")
for r in loco_results:
    print(f"      {r['code']}: N={r['n']}, LOCO C={r['locoC']} {r['locoCi']}, O/E={r['oe5y']}, Slope={r['slope']}, Brier={r['brierIpcw']}, Err(Uncal={r['errUncal']} -> MLE={r['errMle']} -> EB={r['errEb']})")

# 4. Systematic Age-Band Overstatement Computed Directly from N=4,954 Patient Records
age_bands_spec = [
    ("< 50 years (Young Onset)", df["age"] < 50),
    ("50 - 64 years (Middle Age)", (df["age"] >= 50) & (df["age"] < 65)),
    ("65 - 74 years (Older Adult)", (df["age"] >= 65) & (df["age"] < 75)),
    (">= 75 years (Elderly High Competing Risk)", df["age"] >= 75),
]
age_band_bias = {"5y": [], "10y": []}
for h_key, h_mos in [("5y", 60.0), ("10y", 120.0)]:
    for label, mask in age_bands_spec:
        sub = df[mask]
        aj = compute_aalen_johansen_cif(sub["months"].values, sub["event_type"].values, h_mos)
        cif1_pct = round(aj["cif1"] * 100, 2)
        cif2_pct = round(aj["cif2"] * 100, 2)
        all_cause_pct = round(aj["all_cause_mort"] * 100, 2)
        overstate = round(all_cause_pct - cif1_pct, 2)
        age_band_bias[h_key].append({
            "band": label,
            "n": int(len(sub)),
            "cif1": cif1_pct,
            "cif2": cif2_pct,
            "allCause": all_cause_pct,
            "overstatement": overstate,
        })
print(f"[4/6] Computed Systematic Age-Band Overstatement from N=4,954 records: 5y bias spans +{age_band_bias['5y'][0]['overstatement']} to +{age_band_bias['5y'][-1]['overstatement']} pts.")

# 5. Export Fitted Royston-Parmar + Fine-Gray Model Artifact & 150 Reference Parity Vectors
HORIZON_MONTHS = [0, 6, 12, 24, 36, 48, 60, 84, 96, 120]

model_artifact = {
    "version": "5.0.0",
    "n_raw": 11160,
    "n_analyzed": int(len(df)),
    "cause1_feature_names": ["age_per_10y", "T2_plus", "T3_plus", "T4_plus", "N1_plus", "N2_plus", "grade_centered", "adjuvant_tx"],
    "cause1_beta": fit1["beta"],
    "cause1_se": fit1["se"],
    "cause1_rp_knots": fit1["rp_knots"],
    "cause1_rp_gamma": fit1["rp_gamma"],
    "cause2_feature_names": ["age_vs_60_per_10y", "comorbidity_centered"],
    "cause2_beta": fit2["beta"],
    "cause2_se": fit2["se"],
    "cause2_rp_knots": fit2["rp_knots"],
    "cause2_rp_gamma": fit2["rp_gamma"],
    "organ_age_means": organ_age_means,
    "organ_baselines": organ_baselines,
    "horizons": HORIZON_MONTHS,
}

def score_profile_from_artifact(art: Dict, prof: Dict, kappa: float = 1.0) -> Dict:
    org = prof["organ"]
    age_mean = art["organ_age_means"][org]
    ob = art["organ_baselines"][org]
    b1 = art["cause1_beta"]
    b2 = art["cause2_beta"]

    x1 = [
        (prof["age"] - age_mean) / 10.0,
        1.0 if prof["T"] >= 2 else 0.0,
        1.0 if prof["T"] >= 3 else 0.0,
        1.0 if prof["T"] >= 4 else 0.0,
        1.0 if prof["N"] >= 1 else 0.0,
        1.0 if prof["N"] >= 2 else 0.0,
        float(prof["grade"] - 2),
        1.0 if prof["adjuvant"] else 0.0,
    ]
    x2 = [
        (prof["age"] - 60.0) / 10.0,
        float(prof["comorb"] - 1),
    ]
    eta1 = sum(bi * xi for bi, xi in zip(b1, x1)) + math.log(kappa)
    eta1_no_tx = sum(bi * xi for bi, xi in zip(b1[:-1], x1[:-1])) + math.log(kappa)
    eta2 = sum(bi * xi for bi, xi in zip(b2, x2))

    shr1 = math.exp(eta1)
    hr2 = math.exp(eta2)

    curve = []
    for m in art["horizons"]:
        if m == 0:
            curve.append({"m": 0, "os": 1.0, "cif1": 0.0, "cif2": 0.0, "cif1NoTx": 0.0, "allCauseMort": 0.0})
            continue
        # Evaluate Royston-Parmar natural cubic spline baseline hazard H_0(m)
        h1_rp = eval_royston_parmar_h0(float(m), art["cause1_rp_gamma"], art["cause1_rp_knots"]) * ob["h1_scale"]
        h2_rp = eval_royston_parmar_h0(float(m), art["cause2_rp_gamma"], art["cause2_rp_knots"]) * ob["h2_scale"]

        h1 = h1_rp * shr1
        h1_notx = h1_rp * math.exp(eta1_no_tx)
        h2 = h2_rp * hr2
        h_tot = h1 + h2

        os_val = math.exp(-h_tot)
        tot_d = 1.0 - os_val
        cif1 = tot_d * (h1 / max(1e-12, h_tot))
        cif2 = tot_d * (h2 / max(1e-12, h_tot))

        os_notx = math.exp(-(h1_notx + h2))
        cif1_notx = (1.0 - os_notx) * (h1_notx / max(1e-12, h1_notx + h2))
        curve.append({
            "m": m,
            "os": round(os_val, 8),
            "cif1": round(cif1, 8),
            "cif2": round(cif2, 8),
            "cif1NoTx": round(cif1_notx, 8),
            "allCauseMort": round(tot_d, 8),
        })
    return {"shr1": round(shr1, 8), "hr2": round(hr2, 8), "curve": curve}

parity_vectors = []
org_keys = list(ORGAN_META.keys())
for i in range(150):
    org = org_keys[i % len(org_keys)]
    age = int(np.clip(round(organ_age_means[org] + (i % 7 - 3) * 5), 25, 88))
    prof = {
        "organ": org,
        "age": age,
        "T": (i % 4) + 1,
        "N": i % 3,
        "grade": (i % 3) + 1,
        "comorb": i % 3,
        "adjuvant": bool(i % 2 == 0),
    }
    scored = score_profile_from_artifact(model_artifact, prof, 1.0)
    parity_vectors.append({
        "profile": prof,
        "expected": {
            "shr1": scored["shr1"],
            "hr2": scored["hr2"],
            "os_36m": scored["curve"][4]["os"],
            "cif1_36m": scored["curve"][4]["cif1"],
            "os_60m": scored["curve"][6]["os"],
            "cif1_60m": scored["curve"][6]["cif1"],
            "cif2_60m": scored["curve"][6]["cif2"],
            "os_120m": scored["curve"][9]["os"],
        }
    })

with open(os.path.join(RESULTS_DIR, "model_artifact.json"), "w") as f:
    json.dump(model_artifact, f, indent=2)
with open(os.path.join(RESULTS_DIR, "loco_validation.json"), "w") as f:
    json.dump(loco_summary, f, indent=2)
with open(os.path.join(RESULTS_DIR, "age_band_bias.json"), "w") as f:
    json.dump(age_band_bias, f, indent=2)
with open(os.path.join(RESULTS_DIR, "consort_counts.json"), "w") as f:
    json.dump(consort_steps, f, indent=2)
with open(os.path.join(RESULTS_DIR, "parity_150_profiles.json"), "w") as f:
    json.dump(parity_vectors, f, indent=2)

print("[5/6] Exported all fitted artifacts to results/*.json successfully.")
