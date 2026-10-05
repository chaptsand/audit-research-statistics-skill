# 统计协议

## 适用范围与假设

针对有科学配对依据的连续模型指标比较。先确认实验单位与目标推断总体。Wilcoxon signed-rank检验配对差分布在零附近的对称性/位置差，依赖独立配对差等假设；不直接检验算术均值相等，也不在任意分布下等于中位数检验。

同一Run的K折训练集重叠，直接把R×K折当独立重复会夸大样本量。Run聚合是本工作流的主方法，但重复CV仍共享样本。固定数据的重复训练主要描述算法随机性，不能视为独立外部数据集。复杂聚类、不对称差值或总体泛化要求需说明限制，确认后另选方法。

## 描述统计

为兼容常见模型输出可默认ddof=0，但分析前允许明确选择ddof=1，全程标明。ddof=0分母N；ddof=1分母N−1且是样本方差的无偏修正，不能直接称其平方根为无偏SD。

完整等权CV展示全部R×K折mean/sd、检验R个Run均值；固定测试展示R次；任务宏比较展示T个任务均值。两组SD可以不同层级，只要准确标注，不混淆SD与SE。

## Wilcoxon方向和数值

d_i=Left_i−Right_i。默认探索/消融two-sided。greater/less需事前方向依据，结果后不得改变。一般不能将双侧p直接除2得到想要的单侧p。

默认zero_method=wilcox，删除真实零差值再排绝对值平均秩；p及r_rb遵循同一处理。报告n_pairs及n_nonzero。均值差、HL、bootstrap保留全部配对差。全零按约定输出p=1、r_rb=0、HL=0、CI=[0,0]，注明退化数据。

显式声明浮点误差策略。先计算差值，必要时按保存精度恢复理论ties，不能为了结果有利降精度。zero_tolerance及rank_decimals必须有科学或数值依据；默认0和null不修改差值。

小样本采用条件于观测绝对平均秩的sign-flip全枚举：穷举2^n符号。greater用P(W_plus≥observed)，less用P(W_plus≤observed)，two-sided用离中心距离双尾。helper限制n_nonzero≤20；更大主样本需明确渐近或带种子的permutation方案。

折级补充用渐近Wilcoxon，默认zero_method=wilcox、correction=False，并注明依赖局限。SciPy方法名可能为approx/asymptotic，记录版本和实际方法，避免auto漂移。单侧统计量W_plus；双侧min(W_plus,W_minus)。

## 效应量与CI

r_rb=(W_plus−W_minus)/(W_plus+W_minus)。范围[−1,1]，为配对rank-biserial correlation，不能混作Cohen's d。1表示配对差同向秩优势，不表示模型完美。

HL=median{(d_i+d_j)/2 : i≤j}。估计差分布pseudomedian，在对称条件下解释为位置差。一般不等于均值差、差中位数或两组中位数之差。单/双侧不改变HL或r_rb点估计。

均值差、HL、r_rb可不同号。这不是公式错误的充分证据；核验各自计算和差值方向。近零时保留足够小数位。

默认对配对单位的完整d进行至少10000次有放回重抽，每次重算HL，取2.5%和97.5% quantile。固定seed，记录RNG、quantile定义和版本。称paired bootstrap percentile 95% CI，不能称exact。

折级补充CI按真实Run区组抽取全部K折，再算总体折差HL。缺少区组映射则无法计算该区组CI。区组bootstrap不能消除共享数据的一切依赖。

双侧95%CI可搭配单侧p但需说明。此CI未多重校正，且不反演当前Wilcoxon检验；不能强制跨零与p/q阈值严格一致。小样本覆盖率近似。

## BH families

事前按科学问题声明完整family。给m个p排序，q_(i)=min(1,min_{k≥i}(m*p_(k)/k))，再恢复原ID。反向单调累积最小不可省略。q与mean_diff、HL、CI无直接计算关系。

记录family_id、全部comparison_ids、m及primary/supplementary角色。分指标、分协议、合并family可以有不同有效解释，但不能在看到结果后选择有利q，不能以导出文件名临时定义family。

失败/缺失成员先报告，不直接缩小m。BH有依赖条件，在未知强依赖场景不要宣称无条件FDR保证。复用历史q必须核实完整family。

双侧q<alpha结合HL与秩方向解释。单侧显著只支持预设方向，不能据单侧结果宣布反向显著。q≥alpha写未检测到显著差异。等效/非劣效需要预设界值及适当检验。词汇遮蔽不显著不能证明无语义泄露。

## 词汇命中率

每个基因/样本做hit/no_hit。预设关键词、文本范围、大小写、词界、复数/连字符、exact/expanded及否定语境。任一文本命中可定义combined=1，不将同一基因多个描述当独立样本。

互不重叠独立两组用2×2表[[左hit,左no_hit],[右hit,右no_hit]]的双侧Fisher exact。重叠或配对样本需其他设计。

OR=ad/(bc)，>1左组hit odds高。若任一格零，按预设Haldane-Anscombe四格各加0.5计算corrected OR，注明规则；Fisher p始终原整数表。该corrected与BH无关。OR CI单独声明方法。

按实际关键词×文本范围列明完整BH family，不预设固定数量。

## 官方参考

- SciPy Wilcoxon：https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.wilcoxon.html
- R Wilcoxon/HL：https://stat.ethz.ch/R-manual/R-devel/library/stats/html/wilcox.test.html
- statsmodels BH：https://www.statsmodels.org/stable/generated/statsmodels.stats.multitest.multipletests.html
- SciPy Fisher：https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html
- SciPy bootstrap：https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.bootstrap.html
