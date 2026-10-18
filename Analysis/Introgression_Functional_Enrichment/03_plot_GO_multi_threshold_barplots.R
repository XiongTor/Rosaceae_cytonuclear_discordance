#!/usr/bin/env Rscript
# ==============================================================================
# Script: plot_GO_template2_multi_threshold.R
# Author: Antigravity Assistant
# Description: Clean Template 2 GO enrichment plots for Top 3% and Top 5% (title=NULL for multi-panel layout)
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggh4x)
  library(GO.db)
  library(clusterProfiler)
  library(org.DR.eg.db)
})

gene_df <- read.csv("SNP_SRR10377315_SRR15691171_SRR24154117_allgene_Fdm.csv", 
                    header = TRUE, stringsAsFactors = FALSE, fileEncoding = "GBK")
all_gids <- keys(org.DR.eg.db, keytype = "GID")

category_names <- c("BP" = "Biological Process", "CC" = "Cellular Component", "MF" = "Molecular Function")
strip_colors <- c("Biological Process" = "#4DAF4A", "Cellular Component" = "#FF7F00", "Molecular Function" = "#984EA3")

threshold_list <- list(
  "top0.05" = 0.95,
  "top0.03" = 0.97
)

for (tag in names(threshold_list)) {
  q_val <- threshold_list[[tag]]
  thresh <- quantile(gene_df$f_dM_abs, q_val, na.rm = TRUE)
  top_genes <- gene_df[gene_df$f_dM_abs >= thresh, ]
  
  expanded <- top_genes %>% separate_rows(gene_id, gene_location, sep = ",")
  gene_ids <- unique(expanded$gene_id) %>% gsub("_1", "", .) %>% unique()
  IH_gene <- unique(paste0("Dichotomanthes_tristaniicarpa-", gene_ids))
  gene_diff <- IH_gene[IH_gene %in% all_gids]
  
  ego <- enrichGO(
    gene         = gene_diff,
    universe     = all_gids,
    OrgDb        = org.DR.eg.db,
    keyType      = "GID",
    qvalueCutoff = 0.5,
    pvalueCutoff = 0.5,
    ont          = "ALL",
    minGSSize    = 1
  )
  
  res_df <- as.data.frame(ego)
  
  res_df <- res_df %>%
    mutate(
      BgCount    = as.numeric(gsub(".*/", "", BgRatio)),
      GeneCount  = as.numeric(gsub(".*/", "", GeneRatio)),
      RichFactor = Count / BgCount
    )
  
  go_data <- res_df %>%
    filter(!is.na(ONTOLOGY)) %>%
    mutate(GO_Category = factor(category_names[ONTOLOGY], 
                                levels = c("Biological Process", "Cellular Component", "Molecular Function")))
  
  top_go <- go_data %>%
    group_by(GO_Category) %>%
    arrange(p.adjust, pvalue) %>%
    slice_head(n = 10) %>%
    ungroup()
  
  top_go <- top_go %>%
    arrange(GO_Category, desc(p.adjust), desc(pvalue)) %>%
    mutate(Description = factor(Description, levels = unique(Description)))
  
  p_mid <- (max(top_go$p.adjust) + min(top_go$p.adjust)) / 2
  present_categories <- levels(droplevels(top_go$GO_Category))
  actual_strip_colors <- strip_colors[present_categories]
  
  p_template2 <- ggplot(top_go, aes(x = RichFactor, y = Description)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_colour_gradient2(
      low      = "#B2182B", 
      mid      = "#FEE08B", 
      high     = "#4575B4", 
      midpoint = p_mid,
      name     = "p.adjust"
    ) +
    guides(color = guide_colorbar(reverse = TRUE)) +
    facet_grid2(
      GO_Category ~ ., 
      scales = "free_y", 
      space  = "free_y", 
      strip  = strip_themed(background_y = elem_list_rect(fill = actual_strip_colors))
    ) +
    scale_size_continuous(range = c(3, 9), name = "Gene Count") +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(x = "Rich Factor", y = NULL, title = NULL) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(linetype = "dashed", color = "grey90"),
      axis.text.y        = element_text(color = "black", size = 10, face = "plain"),
      axis.text.x        = element_text(color = "black", size = 10),
      axis.title.x       = element_text(color = "black", size = 11, face = "bold"),
      strip.text.y       = element_text(angle = 0, face = "bold", size = 11, color = "white"),
      legend.position    = "right",
      legend.title       = element_text(face = "bold", size = 10),
      legend.key         = element_blank(),
      plot.margin        = margin(t = 5, r = 15, b = 10, l = 10, unit = "pt")
    )
  
  pdf_out <- paste0("GO_enrichment_Template2_", tag, ".pdf")
  png_out <- paste0("GO_enrichment_Template2_", tag, ".png")
  
  ggsave(pdf_out, plot = p_template2, width = 9.5, height = 8.5)
  ggsave(png_out, plot = p_template2, width = 9.5, height = 8.5, dpi = 300, bg = "white")
}
