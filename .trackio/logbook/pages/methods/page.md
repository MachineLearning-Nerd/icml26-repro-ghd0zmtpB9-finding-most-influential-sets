# Methods


---
<!-- trackio-cell
{"type": "markdown", "id": "cell_f1fbd81a7223", "created_at": "2026-07-17T05:05:15+00:00", "title": "Method and independent certificates"}
-->
The released code is vendored unchanged and executed in a Docker image based on `rocker/r-ver:4.5.1` with only `Rcpp`, the release’s optional compiled dependency. The official Rcpp route uses a size-k heap.

Independent checks are deliberately different:

1. complete Python subset enumeration verifies every released finite case;
2. a separate Dinkelbach implementation uses expected-linear `numpy.argpartition` top-k;
3. an oracle/estimated residualized PLM sweep measures exact set recovery against a known separation gap.

All fixed seeds, raw data, commands, and CSV outputs are retained.
