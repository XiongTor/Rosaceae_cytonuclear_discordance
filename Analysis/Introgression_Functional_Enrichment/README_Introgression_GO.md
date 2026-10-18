# 六组 GO 富集韦恩图说明

- 分析对象：6 个 SRR 开头的组合文件夹。
- 阈值：Top 3%（f_dM_abs 的 97% 分位数）和 Top 5%（95% 分位数）。
- 富集方法：沿用各文件夹现有脚本的 clusterProfiler::enrichGO 参数。
- 韦恩集合：每个组合在 BP、CC、MF 中按 p.adjust、pvalue 排序后最多前 10 个 GO 条目。
- 韦恩图交集数字：只属于该精确集合组合的 GO 条目数。
- 完整复核：GO_Venn_exact_intersections.csv 与 GO_Venn_term_membership.csv。
- 跨阈值稳定条目：GO_stable_enriched_functions.csv。
