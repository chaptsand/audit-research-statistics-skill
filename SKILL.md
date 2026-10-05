---
name: audit-research-statistics
description: "Specify, reproduce, and audit statistical comparisons of research models. Use for repeated cross-validation, repeated training on fixed test sets, cross-task macro comparisons, ablation studies, keyword hit-rate analyses, reviewer-requested statistics, and statistical CSV/Excel quality assurance. Covers raw-data contracts, pairing evidence, Wilcoxon signed-rank tests, Benjamini-Hochberg adjustment, matched-pairs rank-biserial correlation, paired Hodges-Lehmann estimates, and confidence intervals. Require an explicit analysis plan; do not hardcode project-specific datasets, sample sizes, or correction families."
---

# Research Statistics Audit

## 1. Establish scope and authority

Default to `audit-only`: permit local numerical verification and create separate QA outputs without overwriting source data, existing results, or scripts. If the user prohibits recalculation or file creation, inspect existing evidence and report in the conversation. Apply corrections or regenerate results only within the explicitly authorized scope. Obtain separate authorization for model changes, retraining, or changes to evaluation sets.

Read [Data contract](references/data-contract.md) and [Statistical protocol](references/statistical-protocol.md) before analysis. For existing deliverables, also read [Audit checklist](references/audit-checklist.md).

Limit this workflow to scientifically justified paired model-score comparisons and independent-group binary hit-rate analyses. Establish another analysis plan for survival outcomes, causal questions, complex clustering, or incompatible study designs.

## 2. Specify an analysis plan

Complete `assets/analysis-plan.template.json` before calculating significance. Declare:

- Source versions, evaluation sets, metric definitions, and comparison direction.
- Pairing keys, inferential units, estimands, and test alternatives.
- Descriptive units, SD denominator, zero handling, ties, and numerical precision.
- Complete multiple-testing families, primary and supplementary roles, and CI methods.
- Authorized outputs and the treatment of missing or failed experiments.

Use the supplied scientific protocol when available. Record deviations and their rationale. Distinguish a prospective plan from a retrospective standardization; do not describe a retrospective plan as preregistered. Keep datasets, run counts, task counts, and family sizes configurable.

## 3. Validate provenance and pairing

Use raw score arrays or metrics reconstructed from complete sample-level predictions. Distinguish test scores from validation scores. Record source hashes, parsing rules, original precision, software versions, and the mapping of runs and folds.

Verify pairing through evaluation-sample identities and actual split mappings. Equal shapes, matching file order, or equal seeds alone are insufficient. Confirmed common splits may justify pairing when training seeds differ; disclose the additional training randomness. Mark unsupported pairing as `UNVERIFIED`.

Never reconstruct paired observations from published means and SDs or generate synthetic observations that match a summary.

## 4. Separate descriptive and inferential units

For repeated $K$-fold cross-validation, aggregate each method $M$ within a run:

$$
\bar{x}^{(M)}_r=\frac{1}{K}\sum_{k=1}^{K}x^{(M)}_{rk},
\qquad r=1,\ldots,R.
$$

Use the $R$ paired run means for the primary comparison. Performance tables may display all $RK$ fold scores if their footnotes identify the descriptive unit and the inferential sample size $n=R$.

For fixed-test repeated training, pair the actual runs. For cross-task macro comparisons, aggregate each task and pair task means. For a single-task analysis, pair run means within that task.

Treat direct fold-level tests as explicitly labeled supplementary analyses. Run aggregation reduces within-run pseudoreplication but does not create independent datasets. State the remaining dependence from shared data and the scope of inference.

## 5. Calculate under an explicit protocol

Fix the direction throughout:

$$
d_i=L_i-R_i.
$$

Use a two-sided alternative for exploratory or ablation comparisons unless a directional hypothesis was specified before inspecting results. Record zero handling, average ranks for ties, precision adjustments, the exact or asymptotic method, continuity correction, and software versions.

Compute matched-pairs rank-biserial correlation from signed rank sums. Compute the paired Hodges-Lehmann estimate from all Walsh averages. Test sidedness does not change these point estimates. Mean difference, Hodges-Lehmann difference, and rank-biserial correlation can have different signs; verify their definitions rather than forcing agreement.

For the default HL interval, bootstrap paired inferential units, recompute HL in each resample, and use a two-sided percentile interval with at least $10{,}000$ resamples and a fixed seed. For fold-level supplements, resample complete run blocks. Label intervals as approximate, marginal, and unadjusted for multiplicity; do not call them exact Wilcoxon-inversion intervals.

Apply BH to the complete predefined family. Record its membership and size. Do not redefine families based on significance, export filenames, or successful rows.

Use `scripts/paired_stats.py` for the supported numerical operations. Run `python scripts/test_paired_stats.py` when validating its implementation. The helper calculates statistics; it does not authenticate provenance, pairing, or scientific assumptions.

## 6. Audit and repair minimally

Assign `PASS`, `FAIL`, `UNVERIFIED`, or `LEGACY` with evidence. Report file integrity, numerical reproducibility, and protocol compliance separately. Preserve historical records and identify the current primary analysis. Relabeling an old table cannot resolve a missing compliant analysis.

Prefer an existing correct supplementary table or correction column when appropriate. After authorization, repair only affected outputs and their generators. Preserve raw data and verify that future reproduction retains the corrected protocol.

## 7. Report without overstating evidence

Include both methods' descriptive Mean ± SD, inferential unit and $n$, Wilcoxon $p$, primary BH $q$, $r_{\mathrm{rb}}$, HL difference, and its CI. Retain full precision in machine-readable results and round only for presentation.

Interpret a significant two-sided result using its stated effect direction. A one-sided result supports only the prespecified direction. For $q\geq\alpha$, report that no statistically significant difference was detected. Such a result does not establish equivalence, noninferiority, or absence of information leakage.

## Resources

| Resource | Purpose |
| --- | --- |
| [Data contract](references/data-contract.md) | Input schemas, provenance, pairing, missingness, and output fields. |
| [Statistical protocol](references/statistical-protocol.md) | Assumptions, equations, test implementation, intervals, BH, and hit rates. |
| [Audit checklist](references/audit-checklist.md) | Evidence requirements, QA status, counting rules, and minimal corrections. |
| [Analysis-plan template](assets/analysis-plan.template.json) | Configurable protocol specification. |
| [Long-format metrics template](assets/metrics-long.template.csv) | Canonical score-table columns. |
| [Paired-input template](assets/pairs-input.template.json) | Explicit inputs to the numerical helper. |
| [Numerical helper](scripts/paired_stats.py) | Paired comparisons and complete-family BH; writes to stdout. |
| [Numerical tests](scripts/test_paired_stats.py) | Synthetic cases and independent numerical checks. |
