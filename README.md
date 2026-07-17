# Reproduction: Finding Most Influential Sets

This is a CPU reproduction of all three claims in ICML 2026 paper `ghd0zmtpB9`, using the official implementation pinned at `12fdd477` and independent subset-enumeration and top-k certificates.

## Results

1. **Exact global MIS selection — verified.** On 120 heterogeneous finite regression cases, the unchanged official Dinkelbach solver exactly matches the official enumeration baseline. Independent enumeration of every subset agrees to `1.15e-12`; no case is suboptimal. The solver uses at most four Dinkelbach updates.
2. **Finite top-k algorithm and scaling — verified.** The official `Rcpp` heap accelerator completes the `n=1,000,000`, `k=100` problem in `0.0572 s` and two updates. An independent expected-linear `argpartition` implementation matches every exact case and completes the same scale in `0.0993 s`. The bundled 1,599-row red-wine data run completes in two updates.
3. **Oracle residualization and separation recovery — verified under controlled conditions.** Across 25 independent PLM-style residualized instances, the exact oracle set is recovered in every trial through a `1e-4` nuisance perturbation. Recovery drops to 96% at `1e-3`/`1e-2`, where perturbation approaches the measured positive separation gap. This directly tests the paper's conditional set-recovery mechanism.

The released positive-denominator guard rejects a rank-deficient deletion, as required by the theorem's assumption.

## Run

```bash
docker build --tag icml26-findingmis-r:4.5.1 --file repro/Dockerfile .
uv venv --python 3.12
source .venv/bin/activate
uv pip install -r requirements.txt
python repro/src/run_fmis.py
pytest -q
```

Outputs include raw official cases, raw sample data, scaling results, separation results, and a compact [summary](outputs/summary.json).

## Scope

The exact finite algorithm and its stated separation condition are reproduced directly. The asymptotic theorem itself is not presented as empirically proved: the finite perturbation sweep is labeled as evidence for its mechanism. No inaccessible data, model training, or GPU is used.

