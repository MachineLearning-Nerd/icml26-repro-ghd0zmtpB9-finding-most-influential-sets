import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "repro/src"))

from run_fmis import dinkelbach_linear_topk, ratio_solution


@pytest.mark.parametrize("seed", range(12))
def test_linear_topk_matches_exhaustive(seed):
    rng = np.random.default_rng(seed)
    x = rng.normal(size=12)
    y = 0.4 * x + rng.standard_t(df=4, size=12)
    for k in (1, 2, 3, 4):
        oracle_value, oracle_set, _ = ratio_solution(x, y, k)
        got_value, got_set, iterations = dinkelbach_linear_topk(x, y, k)
        assert got_value == pytest.approx(oracle_value, abs=1e-11)
        assert iterations <= 1000


def test_nonpositive_denominator_control():
    with pytest.raises(ValueError):
        dinkelbach_linear_topk(np.array([1.0, 0.0]), np.array([1.0, 0.0]), 1)


def test_persisted_full_evidence():
    summary = json.loads((ROOT / "outputs/summary.json").read_text())
    c1 = summary["claim1_exact_global_mis"]
    assert c1["cases"] == 120
    assert c1["nonoptimal_cases"] == 0
    assert c1["max_official_vs_independent_exhaustive"] < 2e-10
    assert c1["max_official_vs_independent_linear_topk"] < 2e-10
    c2 = summary["claim2_finite_topk_algorithm"]
    assert c2["largest_n"] == 1_000_000
    assert c2["real_data"]["n"] == 1599
    c3 = summary["claim3_oracle_and_separation"]
    assert c3["exact_recovery_rate"]["0.0001"] == 1.0
    assert summary["negative_controls"]["nonpositive_denominator"] == "REJECTED_NONPOSITIVE_DENOMINATOR"

