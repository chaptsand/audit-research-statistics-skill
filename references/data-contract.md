# Data Contract

## Contents

1. [Canonical long-format scores](#1-canonical-long-format-scores)
2. [Alternative raw formats](#2-alternative-raw-formats)
3. [Source manifest](#3-source-manifest)
4. [Validation requirements](#4-validation-requirements)
5. [Analysis-plan and helper inputs](#5-analysis-plan-and-helper-inputs)
6. [Result contract](#6-result-contract)

## 1. Canonical long-format scores

Use UTF-8 CSV with a header, stable identifiers, and unrounded numerical values. Leave `fold` empty for a fixed-test repeated-run design.

| Field | Definition |
| --- | --- |
| `dataset` | Dataset or network identifier. |
| `task` | Prediction task or subgroup identifier. |
| `protocol` | Training and information-flow protocol. |
| `method` | Method identifier, including feature mode when relevant. |
| `metric` | Exact metric definition, such as `AUROC`, `PR_AUC_trapezoid`, or `AP`. |
| `run` | Stable repetition identifier. |
| `fold` | Stable cross-validation fold identifier, when applicable. |
| `score` | Original evaluation score at full available precision. |
| `split_id` | Verifiable split identifier, preferably linked to evaluation-sample IDs. |
| `evaluation_set` | `test` or `validation`. |
| `seed` | Training seed; leave empty and explain when unknown. |

Declare a unique observation key. A typical key is

$$
(\mathrm{dataset},\mathrm{task},\mathrm{protocol},\mathrm{method},
\mathrm{metric},\mathrm{run},\mathrm{fold},\mathrm{evaluation\_set}).
$$

Add identifiers for nested resampling, repeated model instances, or other experimental levels. Preserve task IDs for cross-task pairing. Arbitrary file order is not a matching key.

## 2. Alternative raw formats

Accept TXT, CSV, NPY, and complete raw logs when their mappings are documented.

| Design | Numerical structure | Required mapping |
| --- | --- | --- |
| Repeated cross-validation | $R\times K$ matrix | Rows are runs; columns are folds. |
| Fixed-test repeated training | Length-$R$ vector | Each score maps to a run and evaluation set. |
| Multiple tasks | One $R_t\times K_t$ matrix per task | Task, run, and fold identities. |

Document every reshape, transpose, or aggregation. Infer neither orientation nor pairing from performance trends. Keep original files read-only and store normalized derivatives separately.

Published Mean ± SD, images, and rounded manuscript tables do not identify the paired score distribution. They cannot recover signed-rank tests, paired effect sizes, or HL intervals.

## 3. Source manifest

Record:

- Source paths, SHA-256 checksums, data versions, original precision, and metric scales.
- Parser version and any normalization or orientation transformations.
- Model configuration, feature mode, training protocol, seeds, and code version.
- Split-file hashes, evaluation-sample identities, and run/fold mapping.
- Missing, failed, undefined, or overwritten experiments.
- Exact Python, NumPy, SciPy, and other relevant package versions.

A checksum identifies a file version; it does not independently establish that two methods evaluated the same samples. A manually assigned `split_id` also requires supporting evidence.

When deriving scores from predictions, retain sample IDs, true labels, prediction scores, the positive-class convention, and the exact metric calculation and aggregation rules.

## 4. Validation requirements

Check finite values, declared ranges, duplicate keys, missing pairs, metric definitions, and evaluation-set identities. If the metric has a declared range $[a,b]$, verify every score lies within it. Do not assume every metric is bounded by $[0,1]$. Average precision and trapezoidal PR area are distinct quantities.

Document failed runs and undefined metrics, including single-class test sets. Do not silently drop observations, impute scores, or duplicate folds. If a complete-pair analysis is authorized, record the original and retained counts, exclusion reasons, and implications for the analysis and correction family.

Validate pairing from actual split evidence. Equal seeds do not prove common splits; different seeds do not invalidate independently confirmed common splits.

Distinguish

$$
n_{\mathrm{pairs}}=n
\quad\text{from}\quad
n_{\mathrm{nonzero}}
=\sum_{i=1}^{n}\mathbf{1}(\widetilde d_i\neq0),
$$

where $\widetilde d_i$ is the difference after any explicitly declared rank-only precision handling. Retain the original $d_i$, including zeros, for mean difference, HL, and bootstrap calculations.

With missing folds or unequal fold counts, the mean of all available folds need not equal the equally weighted mean of run means.

## 5. Analysis-plan and helper inputs

Complete `assets/analysis-plan.template.json` with comparison IDs, left and right methods, metric definitions, evaluation sets, units, alternatives, descriptive conventions, rank settings, complete families, intervals, and authorized outputs.

Prepare helper input only after verifying these requirements:

| `unit` | Input shape | `pair_ids` |
| --- | --- | --- |
| `run_mean` | One-dimensional paired run means | Run IDs. |
| `fixed_run` | One-dimensional paired repeated-run scores | Run IDs. |
| `task_mean` | One-dimensional paired task means | Task IDs. |
| `fold_supplement` | Paired two-dimensional run-by-fold matrices | Run-block IDs. |

Supply a meaningful `pairing_evidence` reference. The helper checks input structure and numerical values; it cannot authenticate the referenced evidence or validate the study design.

The bundled helper implements `zero_method=wilcox`, `p_method=exact_signflip` or `asymptotic`, and two-sided 95% percentile-bootstrap HL intervals. For a different declared protocol, use and validate an appropriate implementation rather than relabeling these outputs.

## 6. Result contract

Retain machine-readable fields for:

| Category | Required information |
| --- | --- |
| Identity and provenance | `comparison_id`, task, metric, left/right methods, difference direction, sources, hashes. |
| Descriptive statistics | Unit, both means and SDs, `ddof`, and observation counts. |
| Inferential design | Paired unit, `n_pairs`, `n_nonzero`, pairing evidence, and alternative. |
| Test implementation | Method, zero rule, rank precision, tolerance, continuity correction, $W_+$, $W_-$, statistic, and $p$. |
| Multiple testing | `family_id`, membership, $m$, primary/supplementary role, and BH $q$. |
| Effects and intervals | $r_{\mathrm{rb}}$, mean difference, HL difference, CI limits, target, method, level, resampling unit, resample count, seed, and multiplicity status. |
| Verification | Environment, integrity status, numerical status, protocol status, and evidence. |

Label different levels explicitly. Across-task SD describes dispersion among task means; it is not run-level SD. Preserve full precision in analysis files and format manuscript tables as separate derivatives.
