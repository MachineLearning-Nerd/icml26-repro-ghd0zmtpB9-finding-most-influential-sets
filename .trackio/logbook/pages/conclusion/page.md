# Conclusion


---
<!-- trackio-cell
{"type": "markdown", "id": "cell_a364f02d1271", "created_at": "2026-07-17T05:05:18+00:00", "title": "Executive summary", "pinned": true, "pinned_at": "2026-07-17T05:05:19+00:00"}
-->
**Outcome: all three claims are verified within their stated scopes.** The released fractional solver agrees with complete enumeration across every finite test; an independent expected-linear top-k implementation reaches the same global optima; and the PLM-style separation sweep demonstrates exact oracle set recovery until perturbations approach the measured gap.

## Scope & cost

| | Scope | Hardware | Time | Cost | Outcome |
|---|---|---|---|---|---|
| This reproduction | 120 exhaustive finite cases, official + linear top-k through 1e6, 25 separation trials | 4-vCPU CPU, R Docker | ~9 s/run plus image build | $0 | 3/3 verified |
| Full empirical replication | all paper simulations and applications | CPU cluster optional | broader benchmark suite | not incurred | distinct empirical extension |

The asymptotic selection-consistency theorem is represented honestly: the finite sweep is evidence for its specified generated-score/separation mechanism, not a universal proof. Artifacts include raw official cases, vectors, timing CSVs, separation data, source, and 14 passing tests.
