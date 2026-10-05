# 审计与最小修复

## 每项比较

检查source/hash与来源、真实测试指标、完整精度、解析映射；检查唯一ID、split配对证据及缺失；分开核查描述单位/ddof与推断单位/n；验证alternative、ties/zeros、精确/渐近方法；确认BH family成员完整；验证r_rb、Walsh HL和CI的目标/抽样单位；检查主表与副表口径和结论。均值/秩/HL不同号不自动判错。

状态定义：

- PASS：全部适用检查有证据且通过。不适用项为N/A。
- FAIL：确定存在数据、计算或当前主分析的协议冲突。具体列出失败检查。
- UNVERIFIED：缺乏来源、配对、环境或运行证据。只有原始数组不自动PASS；未复算CI则标记CI未验证。
- LEGACY：明确隔离的历史论文结果，列明实际单位、方向、family与当前主分析入口。不能用LEGACY掩盖缺少符合协议主分析的事实。

## QA输出

逐项记录comparison_id、file/sheet、task/metric、source/hash、descriptive_unit/ddof、paired_unit/n、alternative、p_method、family_id/m、数值子状态、方向/语言子状态、总状态、evidence、问题和最小修复。

文件级与唯一comparison级分别汇总。CSV/Excel/JSON的同一比较不重复计入unique comparisons，另计rendered rows时明确命名。报告PASS/FAIL/UNVERIFIED/LEGACY各自分母，不把LEGACY算当前标准PASS。

分开报告文件完整性、数值复现和协议合规。先检查verify脚本实际覆盖范围，不把文件存在或生成成功称为所有统计审核通过。没有执行不能声称已验证。

## 修复边界

audit-only只生成新QA文件。若已存在正确副表，优先建议复用及调整主列/Sheet顺序。用户批准后再局部修改结果和生成脚本，避免下次运行恢复旧口径。保留原数据与旧结果版本。

修复后在独立输出目录验证，检查raw哈希不变、变更范围受限，未改计算定义时p/效应量/HL/CI保持不变。缺少正确主分析时报告仍需工作，不能通过状态重分类清零。
