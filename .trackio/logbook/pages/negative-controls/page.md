# Negative controls


---
<!-- trackio-cell
{"type": "markdown", "id": "cell_47152b81116e", "created_at": "2026-07-17T05:05:16+00:00", "title": "Assumption and adversarial controls"}
-->
The official solver is given a rank-deficient deletion whose remaining curvature denominator is zero; it rejects it with `Non-positive denominator reached`, confirming Assumption 1 is enforced. Independent tests additionally check every exact result with a nonofficial exhaustive oracle and verify that the linear top-k method cannot silently accept the same invalid denominator.


---
<!-- trackio-cell
{"type": "code", "id": "cell_1c9e16a39dd2", "created_at": "2026-07-17T05:05:17+00:00", "title": "Regression and adversarial tests", "command": ["pytest", "-q"], "exit_code": 0, "duration_s": 0.618}
-->
````bash
$ pytest -q
````

exit 0 · 0.6s


````output
..............                                                           [100%]
14 passed in 0.26s

````
