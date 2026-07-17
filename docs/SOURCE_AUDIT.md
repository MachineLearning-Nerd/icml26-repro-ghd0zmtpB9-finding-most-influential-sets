# Source and scope audit

- Paper: *Finding Most Influential Sets*, OpenReview `ghd0zmtpB9`, arXiv:2606.05919.
- Official repository: `nk027/findingMIS`.
- Vendored source commit: `12fdd4775a3df694926bfabe45421aebaa281eb1`.
- Executed official paths: `R/10_miss-frac.R`, `R/11_miss-enum.R`, `R/50_order-partial.R`, and `R/90_helpers.R`.

The release implements the exact Dinkelbach/top-k solver, complete enumeration baseline, optional compiled `Rcpp` top-k routine, simulation code, and the red-wine dataset. The host did not have R, so `repro/Dockerfile` creates a pinned `rocker/r-ver:4.5.1` image with only the release's optional `Rcpp` dependency. The source itself is not modified.

Independent verification enumerates every subset on 120 finite instances in Python and separately implements Algorithm 1 using `numpy.argpartition`, an expected-linear top-k selection. The official accelerated implementation uses a size-k heap, which is `O(n log k)`; that is a valid implementation choice noted by the paper. The independent algorithm tests the stated expected-`O(n)` selection version directly.

The first-order PLM claim is validated as a finite controlled instance: true residualized inputs define an oracle MIS and estimated residualized inputs receive a known perturbation. Exact set recovery is measured across a separation sweep. This is direct evidence for the stated separation mechanism, not a claim that a finite simulation proves the asymptotic theorem universally.

