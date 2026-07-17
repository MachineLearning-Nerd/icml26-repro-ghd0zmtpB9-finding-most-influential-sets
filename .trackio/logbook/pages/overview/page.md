# Overview


---
<!-- trackio-cell
{"type": "markdown", "id": "cell_776e0f6e62d7", "created_at": "2026-07-17T05:04:27+00:00", "title": "Outcome"}
-->
All three claims are **verified within their stated scopes**. The unchanged official R implementation is pinned at `12fdd477` and run in a reproducible R 4.5.1 container. Independent enumeration and expected-linear top-k certificates agree to `1.15e-12`.

- Claim 1: 120 exact finite instances, no suboptimal result.
- Claim 2: released Rcpp solver and independent linear top-k run through one million rows.
- Claim 3: exact oracle set recovery measured through an explicit separation sweep.
