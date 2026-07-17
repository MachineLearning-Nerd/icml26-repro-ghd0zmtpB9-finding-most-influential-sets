#!/usr/bin/env python3
"""Independent certificates for Finding Most Influential Sets (ghd0zmtpB9)."""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import subprocess
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
IMAGE = "icml26-findingmis-r:4.5.1"


def parse_set(text: str) -> tuple[int, ...]:
    return tuple(int(x) for x in text.split(";") if x)


def ratio_solution(x: np.ndarray, y: np.ndarray, k: int):
    """Independent exhaustive oracle for the exact finite deletion identity."""
    beta = float(np.dot(x, y) / np.dot(x, x))
    residual = y - beta * x
    w, c = x * residual, x * x
    total = float(c.sum())
    best_ratio, best_set = -np.inf, None
    second = -np.inf
    for subset in itertools.combinations(range(len(x)), k):
        denom = total - float(c[list(subset)].sum())
        if denom <= 0:
            continue
        ratio = float(w[list(subset)].sum() / denom)
        if ratio > best_ratio:
            second, best_ratio, best_set = best_ratio, ratio, subset
        elif ratio > second:
            second = ratio
    return beta - best_ratio, tuple(i + 1 for i in best_set), best_ratio - second


def dinkelbach_linear_topk(x: np.ndarray, y: np.ndarray, k: int):
    """Independent Algorithm 1 using expected-linear np.argpartition top-k."""
    beta = float(np.dot(x, y) / np.dot(x, x))
    residual = y - beta * x
    w, c = x * residual, x * x
    total = float(c.sum())
    singleton_denominators = total - c
    if np.any(singleton_denominators <= 0):
        raise ValueError("nonpositive denominator")
    eta = float(np.max(w / singleton_denominators))
    previous = None
    for iteration in range(1, 1001):
        scores = w + eta * c
        selected = np.argpartition(-scores, k - 1)[:k]
        selected.sort()
        denom = total - float(c[selected].sum())
        if denom <= 0:
            raise ValueError("nonpositive denominator")
        updated = float(w[selected].sum() / denom)
        chosen = tuple(int(i + 1) for i in selected)
        if chosen == previous or abs(updated - eta) < 1e-12:
            return beta - updated, chosen, iteration
        eta, previous = updated, chosen
    raise RuntimeError("Dinkelbach did not terminate")


def run_official() -> None:
    command = [
        "docker", "run", "--rm", "-v", f"{ROOT}:/work", "-w", "/work",
        IMAGE, "Rscript", "repro/src/run_official.R",
    ]
    subprocess.run(command, check=True)


def load_csv(name: str):
    with (ROOT / "outputs" / name).open(newline="") as handle:
        return list(csv.DictReader(handle))


def validate() -> dict:
    run_official()
    rows = load_csv("official_cases.csv")
    values = defaultdict(list)
    for row in load_csv("official_case_data.csv"):
        values[int(row["case_id"])].append((int(row["index"]), float(row["x"]), float(row["y"])))

    max_official_vs_enum = 0.0
    max_independent_error = 0.0
    max_linear_topk_error = 0.0
    set_mismatches = 0
    max_iterations = 0
    max_linear_iterations = 0
    for row in rows:
        case_id = int(row["case_id"])
        ordered = sorted(values[case_id])
        x = np.array([z[1] for z in ordered])
        y = np.array([z[2] for z in ordered])
        independent_value, independent_set, _ = ratio_solution(x, y, int(row["k"]))
        linear_value, linear_set, linear_iterations = dinkelbach_linear_topk(x, y, int(row["k"]))
        official_value = float(row["official_value"])
        enum_value = float(row["enumerated_value"])
        max_official_vs_enum = max(max_official_vs_enum, abs(official_value - enum_value))
        max_independent_error = max(max_independent_error, abs(official_value - independent_value))
        max_linear_topk_error = max(max_linear_topk_error, abs(official_value - linear_value))
        if abs(official_value - independent_value) > 1e-10:
            set_mismatches += 1
        max_iterations = max(max_iterations, int(row["iterations"]))
        max_linear_iterations = max(max_linear_iterations, linear_iterations)

    scaling = load_csv("official_scaling.csv")
    scaling_out = [{
        "n": int(r["n"]), "k": int(r["k"]), "seconds": float(r["seconds"]), "iterations": int(r["iterations"])
    } for r in scaling]
    real = load_csv("official_realdata.csv")[0]

    rng = np.random.default_rng(260605919)
    linear_scaling = []
    for n in (10_000, 100_000, 1_000_000):
        x = rng.normal(size=n)
        y = x + rng.normal(size=n)
        start = time.perf_counter()
        _, _, iterations = dinkelbach_linear_topk(x, y, 100)
        linear_scaling.append({"n": n, "seconds": time.perf_counter() - start, "iterations": iterations})
    separation = load_csv("official_separation.csv")
    recovery = defaultdict(list)
    for row in separation:
        recovery[float(row["alpha"])].append(row["exact_recovery"].upper() == "TRUE")
    recovery_rate = {str(alpha): float(np.mean(hits)) for alpha, hits in sorted(recovery.items())}
    control = (ROOT / "outputs" / "official_control.txt").read_text().strip()

    summary = {
        "paper": "ghd0zmtpB9",
        "official_commit": "12fdd4775a3df694926bfabe45421aebaa281eb1",
        "claim1_exact_global_mis": {
            "cases": len(rows),
            "max_official_vs_official_enumeration": max_official_vs_enum,
            "max_official_vs_independent_exhaustive": max_independent_error,
            "max_official_vs_independent_linear_topk": max_linear_topk_error,
            "nonoptimal_cases": set_mismatches,
            "max_dinkelbach_iterations": max_iterations,
            "max_linear_topk_iterations": max_linear_iterations,
        },
        "claim2_finite_topk_algorithm": {
            "official_scaling": scaling_out,
            "independent_expected_linear_topk_scaling": linear_scaling,
            "largest_n": max(row["n"] for row in scaling_out),
            "max_iterations": max(row["iterations"] for row in scaling_out),
            "real_data": {"n": int(real["n"]), "k": int(real["k"]), "iterations": int(real["iterations"])},
        },
        "claim3_oracle_and_separation": {
            "trials_per_alpha": len(next(iter(recovery.values()))),
            "exact_recovery_rate": recovery_rate,
            "minimum_positive_separation_gap": min(float(r["separation_gap"]) for r in separation),
        },
        "negative_controls": {"nonpositive_denominator": control},
    }
    (ROOT / "outputs" / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs")
    parser.parse_args()
    validate()
