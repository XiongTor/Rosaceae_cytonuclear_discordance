#!/usr/bin/env python3
"""
提取四棵树中对应内部节点的值，计算 GTEE、ILS、IH 对基因树不一致的贡献。

与原版 extract_node_values.py 的区别：
  - genetree_var = 100 - gCF（总不一致程度）
  - GTEE = 100 - GTEE_recovery（纯估计误差导致的不一致）
  - ILS = concatenation_bl / coalescent_bl（theta）
  - IH = reticulation_index（保持不变）

依赖: ete3 (pip install ete3)
"""

import re
import csv
import sys
from ete3 import Tree


def get_bipartition(node, all_leaves):
    """
    返回节点的 bipartition 标识：
    取该节点下叶节点集合中较小的那一半作为 frozenset
    """
    leaves_below = frozenset(node.get_leaf_names())
    leaves_above = all_leaves - leaves_below
    if len(leaves_below) <= len(leaves_above):
        return leaves_below
    else:
        return leaves_above


def extract_genetree_var(tree_path):
    """
    从 concord.cf.tree 提取 gCF 值
    返回 100 - gCF 作为 genetree_var（总不一致程度）
    """
    t = Tree(tree_path, format=0)
    all_leaves = frozenset(t.get_leaf_names())

    result = {}
    for node in t.traverse():
        if node.is_leaf() or node.is_root():
            continue
        bp = get_bipartition(node, all_leaves)
        result[bp] = 100.0 - node.support
    return result, all_leaves


def extract_gtee(tree_path):
    """
    从 GTEE_raxml_format.tre 提取 recovery %
    返回 100 - recovery 作为 GTEE（纯误差导致的不一致）
    """
    t = Tree(tree_path, format=0)
    all_leaves = frozenset(t.get_leaf_names())

    result = {}
    for node in t.traverse():
        if node.is_leaf() or node.is_root():
            continue
        bp = get_bipartition(node, all_leaves)
        result[bp] = 100.0 - node.support
    return result


def extract_ils_theta(concatenation_path, coalescent_path):
    """
    ILS theta = concatenation_bl / coalescent_bl
    """
    # 加载 concatenation 树 (sortadata)
    with open(concatenation_path) as f:
        nwk = f.read().strip()
    nwk_clean = re.sub(r'\)([\d.]+)/\d+:', r')\1:', nwk)
    t_concat = Tree(nwk_clean, format=0)

    # 加载 coalescent 树
    t_coal = Tree(coalescent_path, format=0)

    all_leaves_concat = frozenset(t_concat.get_leaf_names())
    all_leaves_coal = frozenset(t_coal.get_leaf_names())

    # 提取 concatenation 枝长
    concat_bl = {}
    for node in t_concat.traverse():
        if node.is_leaf() or node.is_root():
            continue
        bp = get_bipartition(node, all_leaves_concat)
        concat_bl[bp] = node.dist

    # 提取 coalescent 枝长
    coal_bl = {}
    for node in t_coal.traverse():
        if node.is_leaf() or node.is_root():
            continue
        bp = get_bipartition(node, all_leaves_coal)
        coal_bl[bp] = node.dist

    # 计算 theta = concatenation_bl / coalescent_bl
    theta = {}
    for bp in concat_bl:
        if bp in coal_bl and coal_bl[bp] > 0:
            theta[bp] = concat_bl[bp] / coal_bl[bp]
        else:
            theta[bp] = None
    return theta, concat_bl, coal_bl


def extract_ih(tree_path):
    """
    从 unbalanced_triples_perc_reticulation_index.tre 提取 reticulation index
    """
    t = Tree(tree_path, format=0)
    all_leaves = frozenset(t.get_leaf_names())

    result = {}
    for node in t.traverse():
        if node.is_leaf() or node.is_root():
            continue
        bp = get_bipartition(node, all_leaves)
        result[bp] = node.support
    return result


def bipartition_label(bp, max_taxa=3):
    """生成可读的 bipartition 标签"""
    taxa = sorted(bp)
    if len(taxa) <= max_taxa:
        return " + ".join(taxa)
    else:
        return f"{taxa[0]} + {taxa[1]} + ... ({len(taxa)} taxa)"


def main():
    # 文件路径
    tree_dir = "."
    gcf_path = f"{tree_dir}/concord.cf.tree"
    gtee_path = f"{tree_dir}/GTEE_raxml_format.tre"
    concatenation_path = f"{tree_dir}/rosa_orthofinder_sortadata_connect_rebranch_rt_oneoutg.tre"
    coalescent_path = f"{tree_dir}/rosa_orthofinder_MO_treeshrink_sp_rt_oneoutg.tre"
    ih_path = f"{tree_dir}/unbalanced_triples_perc_reticulation_index.tre"

    print("正在加载树文件...")

    # 1. genetree_var = 100 - gCF
    print("  提取 genetree_var (100 - gCF)...")
    genetree_var_vals, all_leaves = extract_genetree_var(gcf_path)
    print(f"    找到 {len(genetree_var_vals)} 个内部节点")

    # 2. GTEE = 100 - recovery%
    print("  提取 GTEE (100 - recovery%)...")
    gtee_vals = extract_gtee(gtee_path)
    print(f"    找到 {len(gtee_vals)} 个内部节点")

    # 3. ILS = theta (concatenation_bl / coalescent_bl)
    print("  计算 ILS theta (concatenation_bl / coalescent_bl)...")
    ils_vals, concat_bl, coal_bl = extract_ils_theta(
        concatenation_path, coalescent_path)
    print(f"    找到 {len(ils_vals)} 个内部节点")

    # 4. IH = reticulation index
    print("  提取 IH (reticulation index)...")
    ih_vals = extract_ih(ih_path)
    print(f"    找到 {len(ih_vals)} 个内部节点")

    # 5. 合并输出
    all_bps = (set(genetree_var_vals.keys()) | set(gtee_vals.keys())
               | set(ils_vals.keys()) | set(ih_vals.keys()))

    output_file = "relative_contribution_ILS_Err_Intro_count_use_theta.csv"

    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'node_id', 'bipartition',
            'genetree_var',       # 100 - gCF（总不一致，因变量）
            'GTEE',               # 100 - recovery%（误差导致的不一致）
            'ILS',                # theta = concatenation_bl / coalescent_bl
            'IH',                 # reticulation index
            'concatenation_bl',   # 串联法枝长
            'coalescent_bl'       # 溯祖法枝长
        ])

        node_id = 1
        for bp in sorted(all_bps, key=lambda x: len(x)):
            gvar = genetree_var_vals.get(bp, None)
            gtee_v = gtee_vals.get(bp, None)
            ils_v = ils_vals.get(bp, None)
            ih_v = ih_vals.get(bp, None)
            concat_v = concat_bl.get(bp, None)
            coal_v = coal_bl.get(bp, None)

            # 格式化
            gvar_out = round(gvar, 4) if gvar is not None else 'NA'
            gtee_out = round(gtee_v, 4) if gtee_v is not None else 'NA'
            ils_out = round(ils_v, 6) if ils_v is not None else 'NA'
            ih_out = round(ih_v, 8) if ih_v is not None else 'NA'
            concat_out = round(concat_v, 8) if concat_v is not None else 'NA'
            coal_out = round(coal_v, 8) if coal_v is not None else 'NA'

            label = bipartition_label(bp)
            writer.writerow([
                node_id, label,
                gvar_out, gtee_out, ils_out, ih_out,
                concat_out, coal_out
            ])
            node_id += 1

    print(f"\n完成！结果保存到: {output_file}")
    print(f"   共 {node_id - 1} 个内部节点")
    print(f"\n变量含义:")
    print(f"  genetree_var = 100 - gCF（总不一致程度，因变量）")
    print(f"  GTEE = 100 - GTEE_recovery（纯估计误差导致的不一致）")
    print(f"  ILS = theta = concatenation_bl / coalescent_bl")
    print(f"  IH = reticulation index（基因流导致的不一致）")

    # 打印前几行预览
    print(f"\n前 5 行预览:")
    print(f"{'Node':<6}{'gVar':<10}{'GTEE':<10}"
          f"{'ILS(θ)':<12}{'IH':<12}{'Bipartition'}")
    print("-" * 72)

    with open(output_file) as f:
        reader = csv.DictReader(f)
        for i, row in enumerate(reader):
            if i >= 5:
                break
            print(f"{row['node_id']:<6}"
                  f"{row['genetree_var']:<10}"
                  f"{row['GTEE']:<10}"
                  f"{row['ILS']:<12}"
                  f"{row['IH']:<12}"
                  f"{row['bipartition'][:30]}")


if __name__ == '__main__':
    main()
