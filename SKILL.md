---
name: audit-research-statistics
description: "Standardize and audit research model-comparison statistics and raw metric data. Use for repeated cross-validation, fixed-test repeated runs, cross-task macro comparisons, ablation studies, keyword hit-rate audits, reviewer-requested statistics, and statistical CSV/Excel QA. Covers data contracts, pairing evidence, Wilcoxon signed-rank tests, BH FDR families, rank-biserial effect sizes, Hodges-Lehmann differences and confidence intervals. Use project-supplied analysis plans rather than fixed datasets or hardcoded test counts."
---

# 科研统计检验与审计

## 确定范围

默认采用 audit-only：允许局部复算核对，只生成新的 QA 文件，不覆盖原始数据、现有结果或脚本。用户明确禁止复算或写文件时，只检查已有证据并在对话中报告，不运行计算helper。用户明确授权计算、修复或重建时再处理相应文件。修改模型、重新训练和选择新评价集需要另行授权。

先完整阅读 [数据契约](references/data-contract.md)、[统计协议](references/statistical-protocol.md)。做现有结果检查时另读 [审计检查表](references/audit-checklist.md)。本 skill 针对模型指标的配对比较和独立基因/样本组的词汇命中率；不将其强行用于生存分析、因果推断或任意研究设计。

## 固定计划和来源

使用 `assets/analysis-plan.template.json`，在结果检验前声明：数据版本、评价集、左右方法、配对键、统计单位、alternative、描述统计口径、完整 BH families、CI 方法及输出范围。使用预先确定的研究/论文协议，记录与推荐方法的区别。不要固定为某个项目的样本数、方法数或 family 大小。

仅读取原始指标数组或逐样本预测重建的准确指标。区分测试与验证数据，记录解析步骤、SHA256、源代码版本、指标定义和原始精度。仅有 Mean ± SD 时不能反推出配对分布。

按实际 ID 和划分证据配对。相同数组形状、文件顺序或随机种子不能单独证明配对。核查 Run/Fold 映射和评价样本集；同一划分在种子不同的情况下可构成配对依据，但记录额外训练随机性。缺失证据时标记 UNVERIFIED。

## 区分描述和推断单位

重复 K 折：先对每个 Run 内 K 折平均，再比较 R 个配对 Run 均值。展示可使用全部 R×K 折 Mean ± SD，但必须注明推断单位和 n=R。完整等权矩阵两种总体均值相同，SD不同。

固定测试重复运行：根据真实重复实验设计按 Run 配对，n=R。跨任务总体：每个任务聚合为一个任务均值，按任务配对，n=T。单任务分析：在该任务内用 R 个 Run 均值。不得用排序得分或虚构 ID 恢复配对。

全部折值直接检验仅作为明确约定的补充敏感性分析，不能视为独立重复。Run 聚合减少折内伪重复，但重复 CV 仍共享样本；条件于同一数据集的波动不能代表新数据集不确定性。发现更复杂聚类、强不对称差值或假设不适用时报告限制，获确认后再改方法。

## 计算和检验

固定 Left 减 Right 的差值方向。探索/消融默认双侧；单侧须有事前方向依据。不得根据结果选择侧数。按统计协议记录零差值、ties、浮点精度、精确/渐近方法、连续性校正及环境，避免版本相关的 auto 默认。

用正负秩和计算 r_rb；用配对差值全部 i≤j Walsh averages 的中位数计算 HL。单侧与双侧不改变点估计。Mean difference、HL 和 r_rb 可具有不同符号，检查方向定义与公式，不强求同号。

主要 CI 对配对单位 bootstrap 重算 HL，至少10000次、固定种子，使用双侧95% percentile区间。折级补充 CI 按 Run 区组保留全部折重抽。说明该 CI 未多重校正，且不是当前 Wilcoxon 的反演区间。单侧 p 可搭配明确标记的双侧 CI。

BH 只使用预先定义 family 的完整原始 p 列表，记录 family_id、成员和 m。主要与补充校正范围都需声明，不能按显著性或导出文件随意重组。缺失 family 成员时不能直接对已成功的行校正并声称原 family 已完成。

复用 `scripts/paired_stats.py` 的计算函数和显式配置接口；运行 `python scripts/test_paired_stats.py` 验证实现。该脚本接受显式组织的数值对和单位，只负责计算，不能认证来源、配对或科学假设。

## 审计、解释和交付

逐项标记 PASS、FAIL、UNVERIFIED、LEGACY，分开统计文件完整性、数值复现及协议合规。保留真正历史记录并隔离当前主分析；不能靠改名把缺少主分析的问题清零。优先复用正确副表或调整主列，获得修复授权后再局部重算与改生成逻辑。

论文简表保留两组 Mean ± SD、配对单位及n、p、主要BH q、r_rb、HL和95% CI。完整精度表可保留Mean difference及Conclusion。脚注明确描述与推断的层级。不要用格式化后的数值检验。

双侧显著结果结合 HL/秩方向解释；单侧只支持预设方向。q≥alpha写未检测到显著差异。等效、非劣效和彻底排除泄露需要其他预设设计，不能由不显著检验推出。

## 资源

- `references/data-contract.md`：格式、ID、矩阵、缺失及来源记录。
- `references/statistical-protocol.md`：公式、假设、CI、BH、命中率及官方参考。
- `references/audit-checklist.md`：逐项QA、计数、历史口径和最小修复。
- `assets/analysis-plan.template.json`：可填写的通用分析计划。
- `assets/metrics-long.template.csv`：标准长表字段。
- `assets/pairs-input.template.json`：计算helper输入字段。
- `scripts/paired_stats.py`：单比较与完整BH family计算，默认stdout不覆盖文件。
- `scripts/test_paired_stats.py`：合成数据与独立数值核验。
