# Statistical Protocol

## Contents

1. [Scope and assumptions](#1-scope-and-assumptions)
2. [Descriptive statistics and inferential units](#2-descriptive-statistics-and-inferential-units)
3. [Wilcoxon signed-rank test](#3-wilcoxon-signed-rank-test)
4. [Effect estimates](#4-effect-estimates)
5. [Confidence intervals](#5-confidence-intervals)
6. [Benjamini-Hochberg adjustment](#6-benjamini-hochberg-adjustment)
7. [Binary keyword hit rates](#7-binary-keyword-hit-rates)
8. [Reporting and primary references](#8-reporting-and-primary-references)

## 1. Scope and assumptions

Define paired observations $(L_i,R_i)$ and their differences:

$$
d_i=L_i-R_i,\qquad i=1,\ldots,n.
$$

For a signed-rank location interpretation, use the model

$$
D_i=\theta+\varepsilon_i,
\qquad
\varepsilon_i\overset{\mathrm{iid}}{\sim}F_0,
\qquad
F_0\text{ is symmetric about }0.
$$

Under this model, the signed-rank test addresses the location parameter $\theta$. It does not directly test equality of arithmetic means or, under unrestricted distributions, equality of medians.

Establish the observational unit, scientific pairing, symmetry assumption, and scope of inference. Repeated cross-validation shares data and overlapping training sets. Run aggregation reduces within-run pseudoreplication; it does not remove all dependence. Fixed-test repeated training principally reflects algorithmic variation conditional on that evaluation set.

If the design contains additional clustering, strong asymmetry, or a different target population, document the limitation and establish an appropriate alternative protocol before changing methods.

## 2. Descriptive statistics and inferential units

For $N$ displayed observations, define

$$
\bar{x}=\frac{1}{N}\sum_{i=1}^{N}x_i,
\qquad
s_{\delta}=
\sqrt{\frac{1}{N-\delta}\sum_{i=1}^{N}(x_i-\bar{x})^2},
\qquad
\delta\in\{0,1\}.
$$

Record `ddof` $=\delta$. Permit `ddof=0` for compatibility with existing model summaries, or a declared `ddof=1` convention. The latter gives the usual unbiased sample variance under IID sampling; its square root is not generally an unbiased SD estimate. Do not label SD as standard error.

For a complete equally weighted $R\times K$ score matrix:

$$
\bar{x}_r=\frac{1}{K}\sum_{k=1}^{K}x_{rk},
\qquad
\bar{x}_{\mathrm{fold}}
=\frac{1}{RK}\sum_{r=1}^{R}\sum_{k=1}^{K}x_{rk}
=\frac{1}{R}\sum_{r=1}^{R}\bar{x}_r.
$$

The overall means coincide, but the SDs describe different variation. With `ddof=0`,

$$
s_{\mathrm{fold},0}^{2}
=
\frac{1}{R}\sum_{r=1}^{R}
\left[
\frac{1}{K}\sum_{k=1}^{K}(x_{rk}-\bar{x}_r)^2
\right]
+s_{\mathrm{run},0}^{2}.
$$

Use paired differences between run means for the primary repeated-CV test:

$$
d_r=\bar{x}^{(L)}_r-\bar{x}^{(R)}_r,
\qquad n=R.
$$

| Design | Primary paired unit | Descriptive unit |
| --- | --- | --- |
| Repeated cross-validation | Run mean | All folds or run means, explicitly labeled. |
| Fixed-test repeated training | Actual run | Runs. |
| Cross-task macro comparison | Task mean | Task means. |
| Single-task repeated-CV comparison | Run mean within that task | Task-specific folds or run means. |
| Direct fold comparison | Supplementary fold-level analysis | Folds, with dependence disclosed. |

For cross-task macro analysis, average within each task before comparison and give each task the declared weight. Do not pool heterogeneous task/run/fold levels into an unlabeled sample. Missing folds or unequal weights invalidate the simple equality above.

## 3. Wilcoxon signed-rank test

### 3.1 Alternative hypothesis

Under the symmetric location model, specify:

$$
\begin{array}{lll}
\text{Two-sided:} & H_0:\theta=0 & H_1:\theta\neq0,\\
\text{Greater:} & H_0:\theta\leq0 & H_1:\theta>0,\\
\text{Less:} & H_0:\theta\geq0 & H_1:\theta<0.
\end{array}
$$

Use the boundary $\theta=0$ to calculate directional null probabilities. Choose sidedness before examining results. Default exploratory and ablation comparisons to two-sided tests when no directional plan exists. Do not obtain a desired one-sided result by arbitrarily halving a two-sided $p$-value.

### 3.2 Zeros, precision, and ties

Use `zero_method=wilcox` by default: exclude zero differences from rank calculation. Keep all original differences for mean difference, HL, and bootstrap estimation.

If numerical precision handling is justified, define its transformation before analysis:

$$
\widetilde d_i=\mathcal{T}(d_i),
\qquad
\mathcal{I}=\{i:\widetilde d_i\neq0\},
\qquad
n_0=\lvert\mathcal{I}\rvert.
$$

Record any rank rounding and zero tolerance. Default to no rounding and zero tolerance $0$. Use source precision or floating-point evidence to restore theoretical ties; never lower precision to obtain a favorable result.

Rank the absolute nonzero differences with average ranks for ties:

$$
a_i=\operatorname{rank}_{\mathrm{avg}}(\lvert\widetilde d_i\rvert),
\qquad i\in\mathcal{I}.
$$

Then

$$
W_+=\sum_{i\in\mathcal{I}}a_i\mathbf{1}(\widetilde d_i>0),
\qquad
W_-=\sum_{i\in\mathcal{I}}a_i\mathbf{1}(\widetilde d_i<0).
$$

Report $n_{\mathrm{pairs}}=n$, $n_{\mathrm{nonzero}}=n_0$, both rank sums, and the test statistic. The two-sided statistic is $\min(W_+,W_-)$; the directional statistic is $W_+$.

### 3.3 Conditional exact calculation

For small primary samples, condition on the observed absolute ranks and enumerate independent null signs:

$$
S_i\overset{\mathrm{iid}}{\sim}\operatorname{Bernoulli}(1/2),
\qquad
W_+^{*}=\sum_{i\in\mathcal{I}}a_iS_i,
\qquad
c=\frac{1}{2}\sum_{i\in\mathcal{I}}a_i.
$$

Calculate inclusive tail probabilities:

$$
\begin{aligned}
p_{\mathrm{greater}}&=\Pr(W_+^{*}\geq W_+),\\
p_{\mathrm{less}}&=\Pr(W_+^{*}\leq W_+),\\
p_{\text{two-sided}}&=
\Pr\!\left(\lvert W_+^{*}-c\rvert\geq\lvert W_+-c\rvert\right).
\end{aligned}
$$

Enumerate all $2^{n_0}$ sign assignments, including repeated rank-sum outcomes. This handles tied absolute ranks under the conditional sign-symmetry null. The bundled helper limits full enumeration to $n_0\leq20$.

For larger samples, explicitly select and document an asymptotic method or a separately validated Monte Carlo sign-flip implementation with its seed and finite-resampling rule. Avoid version-dependent `auto` selection. A software option named `exact` does not by itself guarantee exact handling of ties and zeros.

### 3.4 Asymptotic and degenerate cases

For fold-level supplements, use an explicitly declared asymptotic test, `zero_method=wilcox`, and no continuity correction by default. Report the dependence limitation. Record the package version and actual option name; SciPy releases may use `approx` or `asymptotic`.

If no rankable differences remain, the helper convention is $p=1$ and $r_{\mathrm{rb}}=0$. Label this a degenerate case. If every raw difference is zero, raw HL and its bootstrap interval are also zero. If rank precision handling suppresses nonzero raw differences, calculate HL and its interval from the raw differences without forcing them to zero, and flag the discrepancy.

## 4. Effect estimates

### 4.1 Matched-pairs rank-biserial correlation

Define

$$
r_{\mathrm{rb}}=\frac{W_+-W_-}{W_++W_-},
\qquad -1\leq r_{\mathrm{rb}}\leq1.
$$

Use the same zero and rank rules as the test. Positive values indicate positive signed-rank dominance for Left minus Right. Negative values indicate the reverse. This measure is distinct from Cohen's $d$; a value of $1$ describes rank dominance among retained differences, not perfect prediction.

### 4.2 Paired Hodges-Lehmann difference

Use every original paired difference, including zeros. Form the Walsh averages

$$
w_{ij}=\frac{d_i+d_j}{2},
\qquad 1\leq i\leq j\leq n,
$$

and calculate

$$
\widehat{\theta}_{\mathrm{HL}}
=
\operatorname{median}
\{w_{ij}:1\leq i\leq j\leq n\}.
$$

For an even number of Walsh averages, use the midpoint of the two central ordered values.

The population target is the pseudomedian

$$
\theta_{\mathrm{pseudo}}
=
\operatorname{median}\!\left(\frac{D+D'}{2}\right),
$$

where $D$ and $D'$ are independent draws from the paired-difference distribution. Under symmetry it equals the distribution's location and median.

Distinguish this estimate from the arithmetic mean difference

$$
\bar d=\frac{1}{n}\sum_{i=1}^{n}d_i,
$$

the sample median of $d_i$, and the difference between the two marginal sample medians. These quantities need not coincide. Mean difference, HL, and $r_{\mathrm{rb}}$ can also have different signs. Changing test sidedness leaves the point estimates unchanged when all other settings remain fixed.

## 5. Confidence intervals

### 5.1 Primary paired-unit bootstrap

For each bootstrap replication $b=1,\ldots,B$, sample paired-unit indices with replacement:

$$
I_1^{(b)},\ldots,I_n^{(b)}
\overset{\mathrm{iid}}{\sim}\operatorname{Uniform}\{1,\ldots,n\},
\qquad
d_j^{*(b)}=d_{I_j^{(b)}}.
$$

Recompute the HL estimate $\widehat{\theta}_{\mathrm{HL}}^{*(b)}$. For a declared interval level $1-\alpha_{\mathrm{CI}}$, calculate

$$
\mathrm{CI}_{1-\alpha_{\mathrm{CI}}}
=
\left[
Q_{\alpha_{\mathrm{CI}}/2}
\!\left(\widehat{\theta}_{\mathrm{HL}}^{*}\right),
Q_{1-\alpha_{\mathrm{CI}}/2}
\!\left(\widehat{\theta}_{\mathrm{HL}}^{*}\right)
\right].
$$

Default to $B\geq10{,}000$, $\alpha_{\mathrm{CI}}=0.05$, a fixed seed, and linear empirical-quantile interpolation. Record the RNG, resample count, seed, quantile method, and software versions. Resample the paired differences or whole paired records; do not resample Left and Right independently.

The bundled helper implements the two-sided 95% interval. Another interval level or estimator requires a validated implementation consistent with the declared plan.

### 5.2 Fold-level supplementary bootstrap

Retain complete run blocks. For a paired difference matrix $d_{rk}$, sample $R$ run indices:

$$
J_1^{(b)},\ldots,J_R^{(b)}
\overset{\mathrm{iid}}{\sim}\operatorname{Uniform}\{1,\ldots,R\},
\qquad
d_{rk}^{*(b)}=d_{J_r^{(b)},k}.
$$

Recompute the fold-level HL estimate over the resampled matrix. If run-block identities are unavailable, mark the interval `UNVERIFIED` rather than inventing blocks. Block resampling preserves within-run structure; it does not eliminate all dependence from a shared dataset.

### 5.3 Interval interpretation

Label these intervals paired-unit or run-block percentile-bootstrap intervals. They are approximate marginal intervals, unadjusted for multiplicity. They are not exact intervals obtained by inverting the reported Wilcoxon test.

A clearly labeled two-sided interval may accompany a prespecified one-sided test. Do not force agreement between interval inclusion of zero and either raw $p$ or adjusted $q$: they arise from different procedures. Small-sample coverage remains approximate.

## 6. Benjamini-Hochberg adjustment

Define the complete scientific family $\mathcal{F}$ before examining significance, with $m=\lvert\mathcal{F}\rvert$. Sort its raw $p$-values:

$$
p_{(1)}\leq p_{(2)}\leq\cdots\leq p_{(m)}.
$$

Compute

$$
q_{(i)}
=
\min\!\left\{
1,\,
\min_{k=i,\ldots,m}\frac{m}{k}p_{(k)}
\right\},
\qquad i=1,\ldots,m,
$$

and restore the original comparison order. Implement the reverse cumulative minimum; simple multiplication and clipping alone are insufficient.

Use $q$ as shorthand for the BH-adjusted $p$-value. Do not conflate it with a separately estimated Storey $q$-value. Its calculation depends on the complete family's raw $p$-values and membership, without directly using mean difference, HL, or effect size.

Record `family_id`, all member IDs, $m$, scientific rationale, and primary or supplementary role. Different declared families answer different multiplicity questions. Do not select a favorable family after inspecting results or define it merely from the rows exported to a file.

Report incomplete or failed family members before adjustment. Do not silently reduce $m$ and claim that the original family was completed. Verify the original family before reusing historical adjusted values.

BH provides FDR control under its applicable independence or positive-dependence conditions. Do not claim an unconditional guarantee for arbitrary dependence.

## 7. Binary keyword hit rates

### 7.1 Audit unit and matching rule

Define one binary observation per sample or entity:

$$
H_g=\mathbf{1}
\{\text{at least one match in the predefined text scope for entity }g\}.
$$

If the combined scope includes several descriptions,

$$
H_g^{\mathrm{combined}}=\max_{f\in\mathcal{S}}H_{gf},
$$

where $\mathcal{S}$ is the predefined set of text fields. Multiple descriptions of the same entity are not independent observations.

Declare keyword lists, text scope, case sensitivity, word boundaries, inflections, hyphens, exact or expanded rules, and treatment of negation. Report the group denominators and hit proportions.

### 7.2 Independent-group test and odds ratio

For two disjoint independent groups, form

$$
\begin{array}{c|cc}
 & \mathrm{Hit} & \mathrm{No\ hit}\\ \hline
\mathrm{Left} & a & b\\
\mathrm{Right} & c & d
\end{array}
$$

and use a declared two-sided Fisher exact test on the original integer counts. Overlapping or matched samples require another design.

Define the odds ratio

$$
\widehat{\mathrm{OR}}=\frac{ad}{bc}.
$$

Values above $1$ indicate higher hit odds in the Left group. If any cell is zero and the plan specifies a Haldane-Anscombe correction, add $1/2$ to every cell for the point estimate:

$$
\widehat{\mathrm{OR}}_{\mathrm{corrected}}
=
\frac{(a+\tfrac12)(d+\tfrac12)}
{(b+\tfrac12)(c+\tfrac12)}.
$$

Use the unmodified integer table for Fisher $p$. This small-cell OR correction is distinct from BH adjustment. Declare any OR interval separately, including its estimator and approximation.

Construct BH families from the actual prespecified keyword and text-scope comparisons. Do not assume a fixed number of tests.

## 8. Reporting and primary references

Report a significant two-sided comparison using the declared difference direction and effect estimate. If mean, HL, and rank directions disagree, state that explicitly. A significant directional test supports only its prespecified direction.

For $q\geq\alpha$, state that no statistically significant difference was detected. Equivalence and noninferiority require prespecified margins and appropriate tests. A nonsignificant masking experiment alone cannot establish absence of semantic leakage.

Use the following primary documentation to verify implementation details. Record the software version actually used because live documentation can describe a different release.

- [SciPy: Wilcoxon signed-rank test](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.wilcoxon.html)
- [R: Wilcoxon tests and pseudomedian estimation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/wilcox.test.html)
- [SciPy: false discovery control](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.false_discovery_control.html)
- [SciPy: bootstrap intervals](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.bootstrap.html)
- [SciPy: Fisher exact test](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html)
