#!/bin/bash
# Author: Tao Xiong
# Date: 2026-06-08
# Description: gcf scf qcf
# ==== 主体代码开始 ====

# 全部流程的教程：https://iqtree.github.io/doc/recipes/concordance-vector

# 计算每一个基因分区文件的最佳模型
iqtree3 -s rosa_orthofinder_supermatrix.fasta -S rosa_orthofinder_partition.txt --prefix loci -T 50

# 计算scfl,注意这里的astral_species_annotated.tree是需要通过`-u 2`或者`-t 2`选项生成的，带有q1,q2,q3注释的树文件
iqtree3 -s rosa_orthofinder_supermatrix.fasta -te rosa_orthofinder_MO_treeshrink_sp_u2.tre -p loci.best_model.nex --scfl 100 --prefix scfl -T 50

# 计算gcf,注意这里的gene trees文件需要是无根的
iqtree3 -te scfl.cf.tree --gcf rosa_orthofinder_genetrees_noreroot.tre --prefix gcf -T 50

# 计算一致性分数
iqtree3 -s rosa_orthofinder_supermatrix.fasta -te rosa_orthofinder_MO_treeshrink_sp_u2.tre -blfix -p loci.best_model.nex --scfl 1 --prefix coalescent_bl -T 50

# 统计绘表
Rscript concordance_vector.R
Rscript change_labels.R

