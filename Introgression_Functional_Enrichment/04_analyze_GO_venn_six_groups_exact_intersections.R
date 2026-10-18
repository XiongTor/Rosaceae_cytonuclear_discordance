#!/usr/bin/env Rscript

# 本脚本重新执行六个组合在 Top 3% 与 Top 5% 阈值下的 GO 富集，
# 并以每个组合、每个本体最终排名前 10 的 GO 条目构建六集合韦恩图。

# 定位脚本所在目录，确保从任意工作目录执行时都能找到输入文件。
command_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", command_args, value = TRUE)
if (length(file_arg) > 0L) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/")
  root_dir <- dirname(script_path)
} else {
  root_dir <- normalizePath(getwd(), winslash = "/")
}
setwd(root_dir)

# 优先使用项目本地的六集合韦恩图库，并沿用现有 GO 流程的 R 包库。
.libPaths(c(
  file.path(root_dir, ".r_lib"),
  "E:/R/R/R-4.4.1/library",
  .libPaths()
))

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(clusterProfiler)
  library(dplyr)
  library(grid)
  library(gridExtra)
  library(purrr)
  library(ragg)
  library(readr)
  library(stringr)
  library(svglite)
  library(tibble)
  library(tidyr)
  library(venn)
})

# 六组定义与各文件夹中既有 GO 脚本保持一致。
group_map <- tribble(
  ~Group, ~Folder, ~DataFile, ~OrgDbPackage, ~MappingFile,
  "G1", "SRR10377315_SRR15691171_SRR24154117", "SNP_SRR10377315_SRR15691171_SRR24154117_allgene_Fdm.csv", "org.DR.eg.db", "mapping_sp.txt",
  "G2", "SRR10377315_SRR15691173_SRR24154117", "SNP_SRR10377315_SRR15691173_SRR24154117_allgene_Fdm.csv", "org.Cf.eg.db", "map_sp.txt",
  "G3", "SRR14099043_SRR18054568_Torminalis_glaberrima", "SNP_SRR14099043_SRR18054568_Torminalis_glaberrima_allgene_Fdm.csv", "org.Aa.eg.db", "mapping_sp.txt",
  "G4", "SRR14099058_SRR14099059_SRR15691176", "SNP_SRR14099058_SRR14099059_SRR15691176_allgene_Fdm.csv", "org.Cd.eg.db", "mapping_sp.txt",
  "G5", "SRR14099059_SRR15691171_Micromeles_alnifolia_trimmo", "SNP_SRR14099059_SRR15691171_Micromeles_alnifolia_trimmo_allgene_Fdm.csv", "org.Dt.eg.db", "mapping_sp.txt",
  "G6", "SRR14099059_SRR15691176_SRR17214826", "SNP_SRR14099059_SRR15691176_SRR17214826_allgene_Fdm.csv", "org.Cd.eg.db", "mapping_sp.txt"
)

# 两个阈值均纳入，以免在 Top 3% 与 Top 5% 之间武断取舍。
threshold_map <- tribble(
  ~ThresholdTag, ~ThresholdLabel, ~QuantileValue,
  "top0.03", "Top 3%", 0.97,
  "top0.05", "Top 5%", 0.95
)

# 图中使用短组名，详细映射另存为表格，避免六集合图标签拥挤。
group_colors <- c(
  G1 = "#4477AA",
  G2 = "#EE6677",
  G3 = "#228833",
  G4 = "#CCBB44",
  G5 = "#66CCEE",
  G6 = "#AA3377"
)

ontology_labels <- c(
  BP = "Biological Process",
  CC = "Cellular Component",
  MF = "Molecular Function"
)

# 输出集中放在一个新目录，避免覆盖六个子文件夹中的既有结果。
output_dir <- file.path(root_dir, "GO_Venn_six_groups")
dir.create(output_dir, showWarnings = FALSE, recursive = FALSE)

# 读取并校验每组的物种前缀、输入数据与 OrgDb 包。
group_map <- group_map %>%
  mutate(
    DataPath = file.path(root_dir, Folder, DataFile),
    MappingPath = file.path(root_dir, Folder, MappingFile),
    SpeciesPrefix = map_chr(MappingPath, function(path) {
      if (!file.exists(path)) {
        stop("缺少物种映射文件：", path)
      }
      prefix <- trimws(readLines(path, n = 1L, encoding = "UTF-8"))
      if (length(prefix) == 0L || is.na(prefix) || prefix == "") {
        stop("物种映射文件为空：", path)
      }
      prefix
    })
  )

missing_inputs <- group_map$DataPath[!file.exists(group_map$DataPath)]
if (length(missing_inputs) > 0L) {
  stop("缺少输入 Fdm 文件：\n", paste(missing_inputs, collapse = "\n"))
}

missing_packages <- group_map$OrgDbPackage[
  !map_lgl(group_map$OrgDbPackage, requireNamespace, quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("缺少 OrgDb 包：", paste(unique(missing_packages), collapse = ", "))
}

# 对一个组合和一个阈值执行与现有绘图脚本一致的 GO 富集流程。
run_one_enrichment <- function(config_row, threshold_row) {
  group_code <- config_row$Group[[1]]
  folder_name <- config_row$Folder[[1]]
  data_path <- config_row$DataPath[[1]]
  orgdb_package <- config_row$OrgDbPackage[[1]]
  species_prefix <- config_row$SpeciesPrefix[[1]]
  threshold_tag <- threshold_row$ThresholdTag[[1]]
  threshold_label <- threshold_row$ThresholdLabel[[1]]
  quantile_value <- threshold_row$QuantileValue[[1]]

  message("正在分析 ", group_code, " / ", threshold_label, "：", folder_name)

  gene_df <- read.csv(
    data_path,
    header = TRUE,
    stringsAsFactors = FALSE,
    fileEncoding = "GBK",
    check.names = FALSE
  )

  required_columns <- c("f_dM_abs", "gene_id", "gene_location")
  absent_columns <- setdiff(required_columns, names(gene_df))
  if (length(absent_columns) > 0L) {
    stop("输入文件缺少列：", paste(absent_columns, collapse = ", "), "；文件：", data_path)
  }

  fdm_cutoff <- as.numeric(quantile(gene_df$f_dM_abs, quantile_value, na.rm = TRUE))
  top_genes <- gene_df[gene_df$f_dM_abs >= fdm_cutoff, , drop = FALSE]

  expanded <- top_genes %>%
    separate_rows(gene_id, gene_location, sep = ",")

  # 与原脚本保持一致：去除基因编号中的“_1”，再添加物种前缀。
  gene_ids <- unique(gsub("_1", "", expanded$gene_id, fixed = TRUE))
  input_gene_ids <- unique(paste0(species_prefix, "-", gene_ids))

  orgdb <- getExportedValue(orgdb_package, orgdb_package)
  all_gids <- AnnotationDbi::keys(orgdb, keytype = "GID")
  matched_gene_ids <- input_gene_ids[input_gene_ids %in% all_gids]

  if (length(matched_gene_ids) == 0L) {
    stop("没有基因匹配到 OrgDb：", group_code, " / ", threshold_label)
  }

  ego <- enrichGO(
    gene = matched_gene_ids,
    universe = all_gids,
    OrgDb = orgdb,
    keyType = "GID",
    qvalueCutoff = 0.5,
    pvalueCutoff = 0.5,
    ont = "ALL",
    minGSSize = 1
  )

  result_df <- as.data.frame(ego)
  if (nrow(result_df) == 0L) {
    stop("没有得到 GO 富集条目：", group_code, " / ", threshold_label)
  }

  result_df <- result_df %>%
    filter(!is.na(ONTOLOGY), ONTOLOGY %in% names(ontology_labels)) %>%
    mutate(
      Description = as.character(Description),
      BgCount = as.numeric(sub(".*/", "", BgRatio)),
      GeneCount = as.numeric(sub(".*/", "", GeneRatio)),
      RichFactor = Count / BgCount,
      Group = group_code,
      Folder = folder_name,
      SpeciesPrefix = species_prefix,
      ThresholdTag = threshold_tag,
      ThresholdLabel = threshold_label,
      QuantileValue = quantile_value,
      FdmCutoff = fdm_cutoff,
      InputGeneCount = length(input_gene_ids),
      MatchedGeneCount = length(matched_gene_ids)
    ) %>%
    relocate(
      Group, Folder, SpeciesPrefix, ThresholdTag, ThresholdLabel,
      QuantileValue, FdmCutoff, InputGeneCount, MatchedGeneCount
    )

  ontology_counts <- result_df %>%
    count(ONTOLOGY, name = "EnrichedTermCount") %>%
    complete(ONTOLOGY = names(ontology_labels), fill = list(EnrichedTermCount = 0L)) %>%
    pivot_wider(names_from = ONTOLOGY, values_from = EnrichedTermCount, names_prefix = "Enriched_")

  diagnostics <- tibble(
    Group = group_code,
    Folder = folder_name,
    SpeciesPrefix = species_prefix,
    OrgDbPackage = orgdb_package,
    ThresholdTag = threshold_tag,
    ThresholdLabel = threshold_label,
    QuantileValue = quantile_value,
    FdmCutoff = fdm_cutoff,
    TotalInputRows = nrow(gene_df),
    SelectedRows = nrow(top_genes),
    InputGeneCount = length(input_gene_ids),
    MatchedGeneCount = length(matched_gene_ids),
    MatchRate = length(matched_gene_ids) / length(input_gene_ids)
  ) %>%
    bind_cols(ontology_counts)

  list(results = result_df, diagnostics = diagnostics)
}

# 顺序执行 12 次富集并保留每次的诊断信息。
enrichment_runs <- list()
run_index <- 1L
for (group_index in seq_len(nrow(group_map))) {
  for (threshold_index in seq_len(nrow(threshold_map))) {
    enrichment_runs[[run_index]] <- run_one_enrichment(
      group_map[group_index, , drop = FALSE],
      threshold_map[threshold_index, , drop = FALSE]
    )
    run_index <- run_index + 1L
  }
}

all_results <- bind_rows(map(enrichment_runs, "results"))
diagnostics <- bind_rows(map(enrichment_runs, "diagnostics"))

# “最终富集条目”采用现有图的口径：每组、每阈值、每个本体按校正后 P 值取前 10。
top10_results <- all_results %>%
  group_by(ThresholdTag, ThresholdLabel, Group, ONTOLOGY) %>%
  arrange(p.adjust, pvalue, ID, .by_group = TRUE) %>%
  slice_head(n = 10L) %>%
  mutate(RankWithinOntology = row_number()) %>%
  ungroup()

# 建立每个 GO 条目在六组中的 0/1 成员矩阵。
membership_flags <- top10_results %>%
  distinct(ThresholdTag, ThresholdLabel, ONTOLOGY, ID, Group) %>%
  mutate(Present = 1L) %>%
  pivot_wider(names_from = Group, values_from = Present, values_fill = 0L)

for (group_code in group_map$Group) {
  if (!group_code %in% names(membership_flags)) {
    membership_flags[[group_code]] <- 0L
  }
}

term_information <- top10_results %>%
  group_by(ThresholdTag, ThresholdLabel, ONTOLOGY, ID) %>%
  arrange(p.adjust, pvalue, RankWithinOntology, .by_group = TRUE) %>%
  summarise(
    Description = first(Description),
    BestPAdjust = min(p.adjust, na.rm = TRUE),
    BestPValue = min(pvalue, na.rm = TRUE),
    BestRank = min(RankWithinOntology, na.rm = TRUE),
    .groups = "drop"
  )

term_membership <- membership_flags %>%
  select(ThresholdTag, ThresholdLabel, ONTOLOGY, ID, all_of(group_map$Group)) %>%
  left_join(
    term_information,
    by = c("ThresholdTag", "ThresholdLabel", "ONTOLOGY", "ID")
  ) %>%
  rowwise() %>%
  mutate(
    NumberOfGroups = sum(c_across(all_of(group_map$Group))),
    Membership = paste(group_map$Group[c_across(all_of(group_map$Group)) == 1L], collapse = "&")
  ) %>%
  ungroup() %>%
  arrange(ThresholdTag, ONTOLOGY, desc(NumberOfGroups), BestPAdjust, ID)

# 精确交集表可补足六集合韦恩图在高阶交集区域较拥挤的问题。
intersection_summary <- term_membership %>%
  group_by(ThresholdTag, ThresholdLabel, ONTOLOGY, NumberOfGroups, Membership) %>%
  summarise(
    IntersectionSize = n(),
    GO_IDs = paste(ID, collapse = ";"),
    Functions = paste(Description, collapse = " | "),
    .groups = "drop"
  ) %>%
  arrange(ThresholdTag, ONTOLOGY, desc(NumberOfGroups), desc(IntersectionSize), Membership)

# 跨阈值稳定性表用于判断主要功能是否同时出现在 Top 3% 与 Top 5%。
stable_functions <- top10_results %>%
  distinct(ThresholdTag, ThresholdLabel, Group, ONTOLOGY, ID, Description, p.adjust) %>%
  group_by(ONTOLOGY, ID, Description) %>%
  summarise(
    NumberOfGroups = n_distinct(Group),
    Groups = paste(sort(unique(Group)), collapse = ";"),
    NumberOfThresholds = n_distinct(ThresholdTag),
    Thresholds = paste(sort(unique(ThresholdLabel)), collapse = ";"),
    TotalAppearances = n(),
    BestPAdjust = min(p.adjust, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(ONTOLOGY, desc(NumberOfThresholds), desc(NumberOfGroups), desc(TotalAppearances), BestPAdjust)

# 导出可追溯的数据表。
write_csv(group_map %>% select(Group, Folder, SpeciesPrefix, OrgDbPackage, DataFile), file.path(output_dir, "GO_Venn_group_mapping.csv"))
write_csv(diagnostics, file.path(output_dir, "GO_enrichment_diagnostics.csv"))
write_csv(all_results, file.path(output_dir, "GO_enrichment_all_terms.csv"))
write_csv(top10_results, file.path(output_dir, "GO_enrichment_top10_terms.csv"))
write_csv(term_membership, file.path(output_dir, "GO_Venn_term_membership.csv"))
write_csv(intersection_summary, file.path(output_dir, "GO_Venn_exact_intersections.csv"))
write_csv(stable_functions, file.path(output_dir, "GO_stable_enriched_functions.csv"))

# 根据一个阈值和一个本体构建六个 GO ID 集合。
build_venn_sets <- function(threshold_tag, ontology_code) {
  set_names(group_map$Group) %>%
    map(function(group_code) {
      top10_results %>%
        filter(
          ThresholdTag == threshold_tag,
          ONTOLOGY == ontology_code,
          Group == group_code
        ) %>%
        pull(ID) %>%
        unique()
    })
}

# 绘制一个六集合韦恩图；交集标签直接显示 GO 条目数量。
draw_one_venn <- function(threshold_tag, ontology_code, compact = FALSE) {
  venn_sets <- build_venn_sets(threshold_tag, ontology_code)
  threshold_label <- threshold_map$ThresholdLabel[match(threshold_tag, threshold_map$ThresholdTag)]
  # 在组名后直接标注集合大小；例如 G3 (0) 可明确表示该本体没有富集条目。
  set_labels <- paste0(names(venn_sets), " (", lengths(venn_sets), ")")

  venn::venn(
    x = venn_sets,
    snames = set_labels,
    ilabels = "counts",
    zcolor = unname(group_colors[names(venn_sets)]),
    opacity = if (compact) 0.18 else 0.22,
    plotsize = if (compact) 14 else 18,
    ilcs = if (compact) 0.30 else 0.43,
    sncs = if (compact) 0.52 else 0.68,
    borders = TRUE,
    box = FALSE,
    par = FALSE
  )

  title(
    main = paste0(threshold_label, " · ", ontology_code, " — ", ontology_labels[[ontology_code]]),
    line = if (compact) -1.2 else -0.7,
    cex.main = if (compact) 0.92 else 1.15,
    font.main = 2
  )

  if (!compact) {
    mtext("Each set contains up to 10 enriched GO terms; labels show exact-region counts.", side = 1, line = 0.2, cex = 0.66)
  }
}

# 分别输出六张独立矢量图和 600 dpi 位图。
for (threshold_tag in threshold_map$ThresholdTag) {
  for (ontology_code in names(ontology_labels)) {
    output_stem <- file.path(output_dir, paste0("GO_Venn_", threshold_tag, "_", ontology_code))

    grDevices::cairo_pdf(paste0(output_stem, ".pdf"), width = 7.4, height = 7.4, bg = "white")
    par(mar = c(2.1, 0.6, 2.8, 0.6), xpd = NA, family = "sans")
    draw_one_venn(threshold_tag, ontology_code, compact = FALSE)
    dev.off()

    ragg::agg_png(paste0(output_stem, ".png"), width = 7.4, height = 7.4, units = "in", res = 600, background = "white")
    par(mar = c(2.1, 0.6, 2.8, 0.6), xpd = NA, family = "sans")
    draw_one_venn(threshold_tag, ontology_code, compact = FALSE)
    dev.off()

    svglite::svglite(paste0(output_stem, ".svg"), width = 7.4, height = 7.4, bg = "white")
    par(mar = c(2.1, 0.6, 2.8, 0.6), xpd = NA, family = "sans")
    draw_one_venn(threshold_tag, ontology_code, compact = FALSE)
    dev.off()
  }
}

# 输出 2×3 总图：行对应 Top 3%/Top 5%，列对应 BP/CC/MF。
draw_composite <- function() {
  par(
    mfrow = c(2, 3),
    mar = c(1.1, 0.25, 2.2, 0.25),
    oma = c(2.4, 0.4, 3.0, 0.4),
    xpd = NA,
    family = "sans"
  )

  for (threshold_tag in threshold_map$ThresholdTag) {
    for (ontology_code in names(ontology_labels)) {
      draw_one_venn(threshold_tag, ontology_code, compact = TRUE)
    }
  }

  mtext("Six-group GO enrichment overlap", side = 3, outer = TRUE, line = 1.6, cex = 1.45, font = 2)
  mtext("Sets contain up to 10 terms per group and ontology; exact-region labels are GO-term counts.", side = 3, outer = TRUE, line = 0.3, cex = 0.82)
  mtext("Group codes G1–G6 are defined in GO_Venn_group_mapping.csv.", side = 1, outer = TRUE, line = 1.0, cex = 0.72)
}

composite_pdf <- file.path(output_dir, "GO_Venn_BP_CC_MF_top3_top5_composite.pdf")
grDevices::cairo_pdf(composite_pdf, width = 18, height = 11.8, bg = "white")
draw_composite()
dev.off()

composite_png <- file.path(output_dir, "GO_Venn_BP_CC_MF_top3_top5_composite_450dpi.png")
ragg::agg_png(composite_png, width = 18, height = 11.8, units = "in", res = 450, background = "white")
draw_composite()
dev.off()

composite_svg <- file.path(output_dir, "GO_Venn_BP_CC_MF_top3_top5_composite.svg")
svglite::svglite(composite_svg, width = 18, height = 11.8, bg = "white")
draw_composite()
dev.off()

# 另存组名映射图，便于与韦恩图中的 G1–G6 对照。
mapping_table <- group_map %>%
  transmute(
    Group,
    Folder,
    Species = SpeciesPrefix,
    OrgDb = OrgDbPackage
  )

mapping_grob <- tableGrob(
  mapping_table,
  rows = NULL,
  theme = ttheme_minimal(
    base_size = 10,
    core = list(fg_params = list(hjust = 0, x = 0.02)),
    colhead = list(fg_params = list(fontface = 2))
  )
)

mapping_pdf <- file.path(output_dir, "GO_Venn_group_mapping.pdf")
grDevices::cairo_pdf(mapping_pdf, width = 13.5, height = 4.3, bg = "white")
grid.newpage()
grid.draw(mapping_grob)
dev.off()

mapping_png <- file.path(output_dir, "GO_Venn_group_mapping_300dpi.png")
ragg::agg_png(mapping_png, width = 13.5, height = 4.3, units = "in", res = 300, background = "white")
grid.newpage()
grid.draw(mapping_grob)
dev.off()

# 写出简洁方法说明，明确韦恩图集合口径与配套表格用途。
readme_lines <- c(
  "# 六组 GO 富集韦恩图说明",
  "",
  "- 分析对象：6 个 SRR 开头的组合文件夹。",
  "- 阈值：Top 3%（f_dM_abs 的 97% 分位数）和 Top 5%（95% 分位数）。",
  "- 富集方法：沿用各文件夹现有脚本的 clusterProfiler::enrichGO 参数。",
  "- 韦恩集合：每个组合在 BP、CC、MF 中按 p.adjust、pvalue 排序后最多前 10 个 GO 条目。",
  "- 韦恩图交集数字：只属于该精确集合组合的 GO 条目数。",
  "- 完整复核：GO_Venn_exact_intersections.csv 与 GO_Venn_term_membership.csv。",
  "- 跨阈值稳定条目：GO_stable_enriched_functions.csv。"
)
writeLines(readme_lines, file.path(output_dir, "README_GO_Venn.md"), useBytes = TRUE)

message("分析完成，结果目录：", output_dir)
