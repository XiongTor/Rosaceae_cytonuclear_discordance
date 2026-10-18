# !/usr/bin/Rscript
# date: 2026-06-27
# author: xiongtao
# description: This script is used to perform KEGG and GO enrichment analysis on introgressed genes.
# ======================================================================================================
#### 加载包
library(AnnotationForge)
library(clusterProfiler)
library(GO.db)
library(HDO.db)
library(tidyverse)
library(ggplot2)
library(forcats)
library(viridis)
library(patchwork)
library(enrichplot)

#### 读取注释数据（无论是否构建OrgDb都需要）
annot_file <- "all_OGsgenes.emapper.annotations"
orgdb_dir <- "org.DR.eg.db"

# 读取注释文件，跳过以#开头的行
annot_data <- read.delim(annot_file, comment.char = "#", header = FALSE, stringsAsFactors = FALSE)

# 设置列名
colnames(annot_data) <- c(
  "query", "seed_ortholog", "evalue", "score", "eggNOG_OGs",
  "max_annot_lvl", "COG_category", "Description", "Preferred_name",
  "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module",
  "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy",
  "BiGG_Reaction", "PFAMs"
)

# 将"-"替换为NA
annot_data <- annot_data %>%
  dplyr::mutate(across(everything(), ~ifelse(. == "-", NA, .)))

emapper <- annot_data %>%
  dplyr::select(GID=query,Gene_Symbol=Preferred_name,
                GO=GOs,KO=KEGG_ko,Pathway =KEGG_Pathway,
                OG=eggNOG_OGs,Gene_Name =seed_ortholog)
#GID=query,这里的query就是csv文件中的表头信息

#这里共提取了gene_info,  gene2go,gene2ko,gene2pathway,gene2symbol
#你用哪些信息就可以进行相应的增减。
#如果你只是做go富集分析，其实gene_info和gene2go就足够了
#gene2symbol在这里感觉纯粹是凑数，这是参考引文2弄的。
#gene2ko，gene2pathway是用来做kegg富集分析的。

#提取GID与Gene_Name信息，参考1是将X.4作为Gene_Name，应该是自己定义的信息，这里将seed_ortholog作为Gene_Name信息，你可以根据实际再调整。
gene_info <- dplyr::select(emapper,GID,Gene_Name) %>%
  dplyr::filter(!is.na(Gene_Name)) %>%
  dplyr::distinct(GID, .keep_all = TRUE)

#提取GID与GO信息，组成goTable，建库时的goTable需要三列（GID, GO和EVIDENCE），少一列，就会出错.
gene2go <- dplyr::select(emapper,GID,GO) %>%
  separate_rows(GO, sep = ',', convert = F) %>%
  dplyr::filter(!is.na(GO)) %>%
  mutate(EVIDENCE = 'A') %>%
  dplyr::distinct(GID, GO, .keep_all = TRUE)
dim(gene2go)    #查看数据维度。


#提取GID与KO信息，这里只有2列信息
gene2ko<- dplyr::select(emapper,GID,KO) %>%
 separate_rows(KO, sep = ',', convert = F) %>%
  dplyr::filter(!is.na(KO)) %>%
  dplyr::distinct(GID, KO, .keep_all = TRUE)
dim(gene2ko)


#提取GID与Pathway信息，这里只有2列信息
gene2pathway<- dplyr::select(emapper,GID,Pathway) %>%
separate_rows(Pathway, sep = ',', convert = F) %>%
  dplyr::filter(!is.na(Pathway)) %>%
  dplyr::distinct(GID, Pathway, .keep_all = TRUE)
 dim(gene2pathway)


#提取GID与Gene_Symbol信息，Gene_Symbol是Preferred_name信息，这里只有2列信息
gene2symbol<- dplyr::select(emapper,GID,Gene_Symbol) %>%
  dplyr::filter(!is.na(Gene_Symbol)) %>%
  dplyr::distinct(GID, .keep_all = TRUE)
dim(gene2symbol)

if (!dir.exists(orgdb_dir)) {
  message("正在构建OrgDb数据库...")

  # 使用makeOrgPackage构建OrgDb - 需要给参数命名
  AnnotationForge::makeOrgPackage(gene_info=gene_info,
                                  go=gene2go,
                                  ko=gene2ko,
                                  pathway=gene2pathway,
                                  symbol=gene2symbol,
                                  maintainer='XT <xt@example.com>',
                                  author='XT',
                                  version="0.1",
                                  outputDir=".",
                                  tax_id="4097",
                                  genus="D",
                                  species="R",
                                  goTable="go")

  # 安装构建好的包
  install.packages(orgdb_dir, repos = NULL, type = "source")
  message("OrgDb数据库构建完成！")
} else {
  message("OrgDb数据库已存在，跳过构建步骤")
}


# =======================================================================================================
# 加载OrgDb数据库
library(org.DR.eg.db)

# GO富集分析

# 选取Dsuite Fdm值前5%的基因作为候选基因集
gene <- read.csv("SRR10377315_SRR15691171_SRR24154117/SNP_SRR10377315_SRR15691171_SRR24154117_allgene_Fdm.csv", header = T, 
                 stringsAsFactors = FALSE, fileEncoding = "GBK")
threshold <- quantile(gene$f_dM_abs,0.9,na.rm = TRUE)

top0.05_genes <- gene[gene$f_dM_abs >= threshold, ]
top0.05_genes_expanded <- top0.05_genes %>%
  separate_rows(gene_id, gene_location, sep = ",")

gene_id <- unique(top0.05_genes_expanded$gene_id) %>% 
	gsub("_1", "", .) %>% 
	unique()

IH_gene <- paste0("Dichotomanthes_tristaniicarpa-", gene_id)

# 去重，确保每个基因只出现一次
IH_gene <- unique(IH_gene)


# IH_gene <- read.csv("SRR10377315_SRR15691171_SRR24154117/IH_genelist.csv",header = F)
# IH_gene <- IH_gene$V1
# 确认这些基因在 OrgDb 里存在
all_gids <- keys(org.DR.eg.db, keytype = "GID")
gene_diff <- IH_gene[IH_gene %in% all_gids]

# 看看有多少基因匹配上了
cat("输入基因数：", length(IH_gene), "\n")
cat("匹配到OrgDb的基因数：", length(gene_diff), "\n")

# GO 富集分析
ego <- enrichGO(
  gene         = gene_diff,
  universe     = all_gids,          # 背景基因 = OrgDb 里所有基因
  OrgDb        = org.DR.eg.db,
  keyType      = "GID",
  qvalueCutoff = 0.2,
  pvalueCutoff = 0.2,
  ont          = "ALL",
  minGSSize    = 1         # 加上这行，避免基因数少的 term 被过滤
)

head(as.data.frame(ego))

# # 柱状图
# barplot(ego, showCategory=20)+    
#   #### 从高到低，按p值升序排列；横轴标示差异基因与通路基因的重复基因数
#   scale_y_discrete(labels=function(x) str_wrap(x, width = 100))

# # 分面的柱状图
# barplot(ego,
#         split="ONTOLOGY", 
# 				showCategory=20)+
#   facet_grid(ONTOLOGY~., scale="free")+
#   # scale_y_discrete(labels=function(x) str_wrap(x, width = 100))+
#   theme(
#     axis.text = element_text(size = 1, angle = 0, hjust = 1, vjust = 0.5),
#   )

# 分面的点图
p <- dotplot(ego, split = "ONTOLOGY",showCategory=20) +  
  facet_grid(ONTOLOGY~., scale="free")+    ##### 分面展示
  # scale_y_discrete(labels=function(x) str_wrap(x, width = 100))+
  theme(
    axis.text = element_text(size = 1, angle = 0, hjust = 1, vjust = 0.5),
  ) 
p
ggsave("gomf_0.05.pdf", p, width = 10, height = 12)
# #气泡图
# goresult <- data.frame(ego@result)
# go_select <- goresult %>% 
#   filter(ONTOLOGY == "CC" & pvalue < 0.05) %>% 
#   arrange(pvalue) %>% 
#   head(20)  %>% 
#   dplyr::mutate(Description = factor(Description, levels = rev(unique(Description))))


# p <- ggplot(go_select, aes(x = -log10(pvalue), y = Description)) +
#   geom_bar(stat = "identity", fill = "#7f4ea850", width = 0.8) +
#   geom_text(aes(x = 0.05, label = Description), hjust = 0, vjust = 0.8, size = 4.5) +
#   theme_minimal() +
#   theme(
#     axis.text.y = element_blank(),
#     axis.title.x = element_text(margin = margin(t = 10)),
#     axis.title.y = element_text(margin = margin(r = 10)),
#     plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
#     axis.line = element_line(linewidth = 0.75, color = "black"),
#     axis.ticks.x = element_line(linewidth = 0.75),
#     axis.ticks.y = element_blank(),
#     axis.text = element_text(size = 12, color = "black"),
#     axis.title = element_text(size = 12, color = "black"),
#     plot.title = element_text(size = 12, color = "black"),
#     panel.grid = element_blank()
#   ) +
#   labs(x = "-log10(pvalue)", y = "") +
#   scale_x_continuous(expand = c(0, 0)) +
#   coord_cartesian(clip = 'off')

# p

# ggsave("gomf.pdf", p, width = 4, height = 4.5)