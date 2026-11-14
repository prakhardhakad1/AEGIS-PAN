"""Automated Python & JS Parity Verification Suite (Tolerance < 1e-7)."""
import json, math
from aegis_pan.models import COHORTS, score_profile_py

def test_cohort_sum_is_4954():
    assert sum(c["n"] for c in COHORTS) == 4954

def test_partition_of_unity_and_parity_vectors():
    with open("results/parity_150_profiles.json") as f:
        vectors = json.load(f)
    assert len(vectors) == 150
    max_diff = 0.0
    for vec in vectors:
        p = vec["profile"]
        res = score_profile_py(p["organ"], p["age"], p["T"], p["N"], p["grade"], p["comorb"], p["adjuvant"], 1.0)
        c60 = res["curve"][6]
        assert abs((c60["os"] + c60["cif1"] + c60["cif2"]) - 1.0) < 1e-7
        max_diff = max(
            max_diff,
            abs(res["shr1"] - vec["expected"]["shr1"]),
            abs(c60["os"] - vec["expected"]["os_60m"]),
            abs(c60["cif1"] - vec["expected"]["cif1_60m"]),
            abs(c60["cif2"] - vec["expected"]["cif2_60m"]),
        )
    assert max_diff < 1e-7
