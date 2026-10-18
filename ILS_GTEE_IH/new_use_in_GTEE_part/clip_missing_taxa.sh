#!/bin/bash

# ==================== 路径配置 ====================
PERFECT_DIR="simulated_genes"     # 实际上是 FASTA 格式的源目录
GENETREE_DIR="../../04-genetree_reroot_3outg"             # 存放真实基因树的目录
GENE_LIST="./gene_list.txt" # 基因名列表
OUTPUT_DIR="sim_genes_fasta"        # 输出的标准 FASTA 目录
TMP_DIR="clip_tmp"                  # 临时对齐目录
# ==================================================

mkdir -p $OUTPUT_DIR
mkdir -p $TMP_DIR

while read -r line;do
  realsp=$(nw_labels -I $GENETREE_DIR/${line}*.tre)
  simsp=$(grep ">" $PERFECT_DIR/${line}*.fasta|sed 's/>//g')
# comm 需要两边都排好序：-23 表示只保留仅存在于文件1（simsp）中的行
  comm -23 <(echo "$simsp" | sort) <(echo "$realsp" | sort) > taxa_to_delete.txt

  pxrms -s $PERFECT_DIR/${line}*.fasta -f taxa_to_delete.txt > $OUTPUT_DIR/${line}_rm.fasta

done < $GENE_LIST
