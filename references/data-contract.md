# 数据契约

## 首选 long CSV

使用 UTF-8、逗号分隔、完整精度。按实际设计填写字段，固定测试时 fold 为空，不制造虚拟折：

| 字段 | 定义 |
| --- | --- |
| dataset | 数据集或网络ID |
| task | 实际预测任务或子组ID |
| protocol | 训练/信息流协议 |
| method | 方法及特征模式ID |
| metric | 精确指标定义，如 AUROC、PR_AUC_trapezoid、AP |
| run | 稳定重复实验ID |
| fold | 稳定CV折ID，无CV则为空 |
| score | 原始测试指标值 |
| split_id | 划分ID，优先使用评价样本集合校验和 |
| evaluation_set | test 或 validation |
| seed | 随机种子；未知填写空并在manifest解释 |

声明唯一键，通常为 dataset/task/protocol/method/metric/run/fold/evaluation_set。如果一个Run还有重复采样、模型副本等层级，增加相应ID。跨任务分析保留明确task_id，不能用任意文件顺序作任务配对。

## 文件与来源

记录源路径、SHA256、版本、数据尺度、原始数值精度、解析代码、模型参数、特征模式、种子、划分文件/哈希及Run/Fold映射。SHA256识别版本，不能单独证明两个方法评价集相同。CSV标准化副本应能无损追溯原始输入。

兼容TXT/CSV/NPY数组及完整原始日志。重复CV格式为R行×K列，行Run、列Fold；固定测试为R个分数加ID；跨任务为每个任务一份R×K数组加task ID。展平/转置必须有来源映射依据，不能凭性能趋势猜测。

仅有汇总Mean±SD、图片或论文格式化表不能恢复Wilcoxon、配对效应量或HL/CI。不得造出满足汇总数字的人工原始数据。

## 质量和配对

检查有限值、指标声明的范围/尺度、重复键、配对缺失、指标计算定义、真实测试集合及重跑覆盖。AUROC/AUPRC可声明[0,1]；其他指标按实际范围。AP和梯形PR AUC不能改名后直接视为相同。

缺失Run、失败折和无阳性测试集均需记录。未经批准不能删行、填均值或复制折；批准完整配对分析后记录原n、保留n和删除原因，说明family处理。缺折时，不保证折均值等于等权Run均值。

种子相同不证明划分一致；种子不同不自动取消同一划分的配对。需要实际split证据和对应映射。只提供同形矩阵则配对可信性标记UNVERIFIED。

原始文件只读。新解析副本、QA、汇总和论文排版表分开保存。记录环境的Python、NumPy、SciPy及所有实际包确切版本。

## 分析计划

使用 `assets/analysis-plan.template.json`，逐项声明comparison_id、Left/Right、metric、评价集、推断单位、alternative、描述单位与ddof、完整family成员、CI、零值与精度规则。计划不能从结果显著性倒推。

计算helper的输入只能在上述验证后创建。`unit=run_mean`与`fixed_run`输入一维配对值；`task_mean`输入任务均值；`fold_supplement`输入有真实Run区组映射的R×K数组。helper验证数值但不证明配对设计。

## 输出字段

完整结果至少有comparison_id、Left/Right、差值定义、source/hash、descriptive_unit、mean/sd/ddof、paired_unit、n_pairs/n_nonzero、alternative、test_method、zero_method、rank_precision、W_plus/W_minus、statistic、p、family_id/m、q、r_rb、mean_diff、hl_diff、ci_low/high、ci_method、bootstrap_n/seed、status。

描述统计和推断可用不同单位，必须标注。跨任务SD按任务均值计算，不能把混合层级的列全部标成Run SD。保留full_precision，最后单独格式化论文表。
