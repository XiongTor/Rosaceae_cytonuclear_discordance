library(dplyr)

# 合并SNP和QUIBL的结果

# ========== 设定要处理的trio ==================
trio = "SRR14099059_SRR15691176_SRR17214826"

# ========== 读取要处理的来自SNP分析方法的所有直系同源基因的Fdm信息文件 ======
snp <- read.csv(paste0(trio,"/SNP_",trio,"_allgene_Fdm.csv"))
snp <- snp %>%
  group_by(gene_id) %>%
  mutate(window_span = windowEnd - windowStart) %>%
  slice_max(window_span, n = 1, with_ties = FALSE) %>%
  mutate(
    D        = mean(snp$D[snp$gene_id == cur_group()$gene_id], na.rm = TRUE),
    f_d      = mean(snp$f_d[snp$gene_id == cur_group()$gene_id], na.rm = TRUE),
    f_dM     = mean(snp$f_dM[snp$gene_id == cur_group()$gene_id], na.rm = TRUE),
    d_f      = mean(snp$d_f[snp$gene_id == cur_group()$gene_id], na.rm = TRUE),
    f_dM_abs = mean(snp$f_dM_abs[snp$gene_id == cur_group()$gene_id], na.rm = TRUE)
  ) %>%
  select(-window_span) %>%
  ungroup() %>%
  arrange(desc(f_dM_abs)) %>%
  mutate(SNP_number = row_number())
# ========== 读取QUIBL结果 ============================
quibl <- read.csv(paste0(trio,"/final_",trio,"_ih_result.csv"))
quibl$gene <- gsub("_noout.*","",quibl$gene)

# ========= 结果合并 ================
final <- merge(snp,quibl,by.x="gene_id",by.y="gene",all.x = T)

final_plot <- final %>%
  filter(!is.na(intro_prob)) %>%
  arrange(desc(intro_prob)) %>%
  mutate(QUIBL_number = row_number())

write.csv(final_plot,paste0(trio,"/final_snp_quibl_",trio,".csv"))

# ========= 画图，x轴按照quibl的intro_prob的大小来排序 =======
p <- ggplot(final_plot,aes(x = QUIBL_number, y = f_dM_abs))+
  geom_col()+
  theme_bw()

ggsave(paste0(trio,"/Fdm_snp_quibl_",trio,".pdf"), p, height = 10, width = 20)
# 
# library(ggplot2)
# library(dplyr)
# 
# final_plot <- final %>%
#   filter(!is.na(intro_prob)) %>%
#   arrange(desc(intro_prob)) %>%
#   mutate(group = cut(intro_prob,
#                      breaks = quantile(intro_prob, probs = seq(0, 1, 0.5)),
#                      include.lowest = TRUE,
#                      labels = paste0("G", 2:1)),
#          QUIBL_number = row_number())  # G1=最高组
# 
# ggplot(final_plot, aes(x = group, y = f_dM_abs)) +
#   geom_boxplot(fill = "#4C9BE8", color = "black", outlier.size = 0.8, outlier.alpha = 0.5) +
#   labs(x = "intro_prob group (low → high)", y = expression(f[dM])) +
#   theme_classic(base_size = 12) +
#   theme(
#     axis.text.x  = element_text(size = 11, color = "black"),
#     axis.text.y  = element_text(size = 11, color = "black"),
#     axis.title.x = element_text(size = 14, color = "black", margin = margin(t = 8)),
#     axis.title.y = element_text(size = 14, color = "black", margin = margin(r = 8)),
#     axis.line    = element_line(color = "black", linewidth = 0.5),
#     axis.ticks   = element_line(color = "black")
#   )
