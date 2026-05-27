#!/bin/bash

# ==================== 路径配置 ====================
PERFECT_DIR="sim_genes_perfect"     # 实际上是 FASTA 格式的源目录
GENETREE_DIR="genetree"             # 存放真实基因树的目录
OUTPUT_DIR="sim_genes_fasta"        # 输出的标准 FASTA 目录
TMP_DIR="clip_tmp"                  # 临时对齐目录
# ==================================================

mkdir -p $OUTPUT_DIR
mkdir -p $TMP_DIR

echo "正在处理假装成 .phy 的 FASTA 文件..."

shopt -s nullglob

for phy_file in ${PERFECT_DIR}/*.phy; do
    [[ -f "$phy_file" ]] || continue
    
    gene_id=$(basename "$phy_file" .phy)
    matched_trees=(${GENETREE_DIR}/${gene_id}*.tre ${GENETREE_DIR}/${gene_id}*.tree)
    
    output_fa="${OUTPUT_DIR}/${gene_id}.fasta"
    
    if [ ${#matched_trees[@]} -gt 0 ]; then
        actual_tree="${matched_trees[0]}"
        
        # 1. 极其稳健地从树中剥离出纯物种名单
        tr -d '[:space:]' < "$actual_tree" | tr '()[],;:+' '\n' | grep -vE '^[0-9.]+$' | grep -v '^$' | sort -u > ${TMP_DIR}/taxa_${gene_id}.txt
        
        # 如果名单为空，直接整流拷贝
        if [ ! -s "${TMP_DIR}/taxa_${gene_id}.txt" ]; then
            cp "$phy_file" "$output_fa"
            continue
        fi
        
        # 2. ⚡【核心魔法】：利用 RS=">" 降维打击多行 FASTA
        awk -v taxa_file="${TMP_DIR}/taxa_${gene_id}.txt" '
            BEGIN {
                # 先把物种名单读进数组
                while ((getline < taxa_file) > 0) {
                    taxa[$1] = 1
                }
                close(taxa_file)
                # 设置记录分隔符为 >，这样每个物种块会被当成一整行处理
                RS = ">"
                ORS = ""
            }
            NR > 1 {
                # 提取当前块的第一行作为物种名
                match($0, /^[^\n]+/)
                id = substr($0, RSTART, RLENGTH)
                # 去除可能残存的尾部回车
                gsub(/\r/, "", id)
                
                # 如果这个物种存在于树的名单中，原封不动写出
                if (id in taxa) {
                    print ">" $0
                }
            }
        ' "$phy_file" > "$output_fa"
        
        rm -f ${TMP_DIR}/taxa_${gene_id}.txt
    else
        # 兜底：如果没找到树，直接拷贝过去改名即可
        cp "$phy_file" "$output_fa"
    fi
done

shopt -u nullglob
rmdir $TMP_DIR

echo "🎉 修复版转换完成！"
echo "🔍 正在抽查新生成的 FASTA 内容："
head -n 6 ${OUTPUT_DIR}/*.fasta | head -n 12