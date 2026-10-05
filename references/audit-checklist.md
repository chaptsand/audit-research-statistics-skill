# Audit Checklist and Minimal Corrections

## 1. Evidence required for each comparison

| Check | Required evidence |
| --- | --- |
| Provenance | Original source, version/hash, parsing rules, original precision, and evaluation-set identity. |
| Pairing | Stable IDs, actual split or matching evidence, run/fold mapping, and missing-pair handling. |
| Units | Separate descriptive and inferential units, counts, aggregation rules, and SD denominator. |
| Test | Declared alternative, zero and tie rules, precision handling, numerical method, continuity correction, and environment. |
| Multiplicity | Complete family membership, scientific rationale, size, and primary/supplementary role. |
| Effects | Declared Left minus Right direction, correct rank-biserial formula, and paired Walsh HL calculation. |
| Intervals | Correct target, paired or block resampling unit, method, level, seed, resample count, and multiplicity label. |
| Presentation | Accurate labels, full-precision source values, consistent primary tables, and qualified interpretation. |
| Reproduction | Executed numerical checks, logs, and the exact coverage of verification scripts. |

Do not classify a sign difference between mean difference, HL, and rank-biserial correlation as an error without checking their definitions and calculations.

## 2. Status definitions

| Status | Definition |
| --- | --- |
| `PASS` | All applicable checks have adequate evidence and satisfy the declared protocol. |
| `FAIL` | A demonstrated data, numerical, or current-protocol error exists; identify the failed check. |
| `UNVERIFIED` | Required provenance, pairing, environment, or execution evidence is unavailable. |
| `LEGACY` | An explicitly isolated historical result with its original unit, sidedness, family, and an identified current primary analysis. |
| `N/A` | A particular check is inapplicable, with a reason; it is not a comparison-level success. |

Possession of raw arrays alone does not establish valid pairing. Mark an unexecuted CI verification as unverified even if its point estimate matches. Historical labeling cannot substitute for a missing compliant primary analysis.

## 3. QA output and counting

Record per comparison:

- Comparison ID, file/sheet, task, metric, methods, sources, and hashes.
- Descriptive unit, `ddof`, inferential unit, counts, and pairing evidence.
- Alternative, numerical test method, precision settings, and correction family.
- Expected and observed quantities at full precision.
- Separate integrity, numerical, protocol, and interpretation statuses.
- Evidence, identified issue, and the smallest justified correction.

Declare numerical tolerances and their rationale. For stochastic intervals, distinguish reproducibility under the same seed and algorithm from expected Monte Carlo differences across implementations.

Count unique statistical comparisons separately from files and rendered rows. The same comparison in CSV, Excel, and JSON is one unique comparison. State the denominator for every percentage. Do not include `LEGACY` entries in a current-protocol `PASS` rate.

Inspect what a verification command actually checks. File existence, successful generation, numerical agreement, and methodological compliance are different verification levels. Do not report checks that were not executed.

## 4. Audit-only boundary

In `audit-only` mode, create separate QA artifacts and preserve existing data, results, and generators. If the user also prohibits calculations or file creation, restrict the audit accordingly.

Suggest promotion of an already correct supplementary table or correction column when that resolves the issue. Report missing compliant analyses as outstanding work. Do not remove outstanding work by changing status labels.

## 5. Authorized minimal corrections

After authorization:

1. Preserve the raw inputs and historical outputs.
2. Modify only the affected results and their generation logic.
3. Validate revised outputs in an independent destination.
4. Compare source hashes before and after the correction.
5. Verify that reproduction retains the corrected primary and supplementary roles.

For an unchanged statistical definition and unchanged numerical inputs, verify that

$$
(p,q,r_{\mathrm{rb}},\widehat{\theta}_{\mathrm{HL}},
\mathrm{CI}_{\mathrm{low}},\mathrm{CI}_{\mathrm{high}})
_{\mathrm{before}}
=
(p,q,r_{\mathrm{rb}},\widehat{\theta}_{\mathrm{HL}},
\mathrm{CI}_{\mathrm{low}},\mathrm{CI}_{\mathrm{high}})
_{\mathrm{after}},
$$

subject to the declared numerical tolerance and identical stochastic settings. A change to family membership can legitimately alter $q$; a change to the inferential unit can alter the entire analysis. Do not describe either as a purely cosmetic repair.
