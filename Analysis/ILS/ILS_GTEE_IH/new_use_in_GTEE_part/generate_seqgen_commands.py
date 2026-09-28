#!/usr/bin/env python3
"""
generate_seqgen_commands.py
============================
从 IQ-TREE 的 .iqtree 文件中提取替代模型参数，生成 Seq-Gen 模拟命令。

用法：
    python generate_seqgen_commands.py

输入文件（当前目录下）：
    gene_list.txt          - 基因 ID 列表，每行一个
    alignment_length.txt   - 制表符分隔，列1: 文件名, 列2: 序列长度

输出：
    seqgen_commands.sh     - 包含所有 seq-gen 命令的 shell 脚本
    simulated_genes/       - 模拟序列输出目录（如果指定）

依赖的外部文件：
    ../../03-genetrees/{gene_id}_noout_output.iqtree   - IQ-TREE 模型文件
    ../../03-genetrees/{gene_id}_noout_output.treefile - 树文件 (Newick, 用于参考)
"""

import os
import re
import sys
import glob
import argparse
from pathlib import Path


# ============================================================
# 1. 解析 IQ-TREE .iqtree 文件
# ============================================================

def parse_iqtree(filepath):
    """
    解析 IQ-TREE 的 .iqtree 文件，提取替代模型参数。
    
    返回 dict，包含：
        model       : 完整模型名 (如 "TVM+F+G4")
        base_model  : 基础模型名 (如 "TVM")
        has_F       : 是否使用经验碱基频率
        has_G       : 是否使用 Gamma 异质性
        has_I       : 是否使用不变位点
        gamma_cats  : Gamma 离散类别数 (如有)
        gamma_alpha : Gamma 形状参数 α (如有)
        inv_prop    : 不变位点比例 (如有)
        rates       : GTR 6个速率 [A-C, A-G, A-T, C-G, C-T, G-T] (如有)
        tstv_ratio  : 转换/颠换比 (HKY/F84 模型)
        freqs       : 碱基频率 [A, C, G, T]
    """
    with open(filepath, 'r') as f:
        content = f.read()

    params = {}

    # --- 模型名称 ---
    m = re.search(r'Model of substitution:\s*(\S+)', content)
    if not m:
        raise ValueError(f"Cannot find 'Model of substitution' in {filepath}")
    params['model'] = m.group(1)

    # 解析模型组件
    model_str = params['model']
    params['has_F'] = '+F' in model_str
    params['has_G'] = '+G' in model_str
    params['has_I'] = '+I' in model_str
    params['base_model'] = model_str.split('+')[0]

    # --- Gamma 类别数 ---
    m = re.search(r'(?:Model of rate heterogeneity|Among-site rate heterogeneity):\s*Gamma\s*(?:\+\s*Invar)?\s*(?:with\s*(\d+)\s*categor)', content)
    if not m:
        # 备用模式
        m = re.search(r'Gamma\s*(?:with\s*(\d+)\s*categor)', content)
    if m:
        params['gamma_cats'] = int(m.group(1))

    # --- Gamma 形状参数 ---
    m = re.search(r'Gamma shape alpha:\s*([\d.]+)', content)
    if m:
        params['gamma_alpha'] = float(m.group(1))

    # --- 不变位点比例 ---
    m = re.search(r'Proportion of invariable sites:\s*([\d.]+)', content)
    if m:
        params['inv_prop'] = float(m.group(1))

    # --- GTR 速率参数 (Rate parameter R) ---
    # 匹配从 "Rate parameter R:" 到下一个空行或章节标题
    m = re.search(
        r'Rate parameter R:\s*\n(.*?)(?=\n\s*\n|\n(?:State frequencies|Rate matrix|Model of rate|Base frequencies|Among-site|\Z))',
        content, re.DOTALL
    )
    if m:
        rates = {}
        for line in m.group(1).strip().split('\n'):
            line = line.strip()
            if not line:
                continue
            rm = re.match(r'\s*([ACGT])-([ACGT]):\s*([\d.]+)', line)
            if rm:
                rates[f"{rm.group(1)}-{rm.group(2)}"] = float(rm.group(3))
        # 标准顺序: A-C, A-G, A-T, C-G, C-T, G-T
        rate_order = ['A-C', 'A-G', 'A-T', 'C-G', 'C-T', 'G-T']
        if all(k in rates for k in rate_order):
            params['rates'] = [rates[k] for k in rate_order]

    # --- 转换/颠换比 (HKY, K2P, F81 等) ---
    # 可能的格式: "Transition/transversion ratio: 3.5", "TS/TV ratio: 2.0"
    m = re.search(
        r'(?:Transition/transversion ratio|TS/TV ratio|kappa)\s*:?\s*([\d.]+)',
        content
    )
    if m:
        params['tstv_ratio'] = float(m.group(1))

    # --- 碱基频率 ---
    # 尝试匹配经验频率
    freq_section = re.search(
        r'State frequencies:.*?\n(.*?)(?=\n\s*\n|\n(?:Rate matrix|Model of rate|Base frequencies|Among-site|\Z))',
        content, re.DOTALL
    )
    if freq_section:
        freqs = {}
        for line in freq_section.group(1).strip().split('\n'):
            fm = re.search(r'pi\(([ACGT])\)\s*=\s*([\d.]+)', line)
            if fm:
                freqs[fm.group(1)] = float(fm.group(2))
        if len(freqs) == 4:
            params['freqs'] = [freqs.get(b, 0.0) for b in 'ACGT']

    # 如果没有找到频率（无 +F 的模型），但 IQ-TREE 仍可能报告 "Base frequencies"
    if 'freqs' not in params:
        # 尝试查找 "Base frequencies" 或 "Equal frequencies"
        freq_section2 = re.search(
            r'(?:Base frequencies|State frequencies):.*?\n(.*?)(?=\n\s*\n|\n(?:Rate matrix|Model of rate|Among-site|\Z))',
            content, re.DOTALL
        )
        if freq_section2:
            freqs = {}
            for line in freq_section2.group(1).strip().split('\n'):
                fm = re.search(r'pi\(([ACGT])\)\s*=\s*([\d.]+)', line)
                if fm:
                    freqs[fm.group(1)] = float(fm.group(2))
            if len(freqs) == 4:
                params['freqs'] = [freqs.get(b, 0.0) for b in 'ACGT']

    # 仍然没有 → 使用等频率
    if 'freqs' not in params:
        params['freqs'] = [0.25, 0.25, 0.25, 0.25]

    return params


# ============================================================
# 2. IQ-TREE 模型 → Seq-Gen 参数映射
# ============================================================

# 属于 GTR 家族的模型（Seq-Gen 用 -mGTR + -r 6个速率）
# 注意：IQ-TREE 中同一模型可能有多个名称（如 K3P=K81, TN=TN93）
GTR_FAMILY = {
    'GTR', 'TVM', 'TVMe', 'SYM',
    'TIM', 'TIM2', 'TIM3', 'TIM3e', 'TIMe', 'TIM2e',
    'TPM2', 'TPM2u', 'TPM3', 'TPM3u',
    'TNe',
    'K81', 'K81u', 'K3P', 'K3Pu',  # K3P = K81 (Kimura 3-parameter)
    'TN',                            # TN = TN93 (Tamura-Nei)
}

# 属于 HKY 家族的模型（Seq-Gen 用 -mHKY + -t 转换/颠换比）
HKY_FAMILY = {
    'HKY', 'K2P', 'F81',
}

# 纯 JC69 模型
JC_FAMILY = {'JC', 'JC69'}


def iqtree_to_seqgen(params):
    """
    将 IQ-TREE 参数映射为 Seq-Gen 参数。
    
    返回 dict:
        seqgen_model : 'GTR' 或 'HKY'
        seqgen_args  : 命令行参数字典
    """
    base = params['base_model']
    args = {}

    # --- 确定 Seq-Gen 模型 ---
    if base in GTR_FAMILY:
        args['model'] = 'GTR'
        if 'rates' in params:
            args['rates'] = params['rates']
        else:
            print(f"  WARNING: GTR-family model '{base}' but no rate parameters found", file=sys.stderr)
    elif base in HKY_FAMILY:
        args['model'] = 'HKY'
        if 'tstv_ratio' in params:
            args['tstv'] = params['tstv_ratio']
        elif base == 'K2P':
            args['tstv'] = 2.0  # fallback
        elif base == 'F81':
            args['tstv'] = 0.5  # F81 = HKY with equal ts/tv rates
    elif base in JC_FAMILY:
        args['model'] = 'HKY'
        args['tstv'] = 0.5
        args['freqs'] = [0.25, 0.25, 0.25, 0.25]
    else:
        # 未知模型 → 默认用 GTR
        print(f"  WARNING: Unknown model '{base}', defaulting to GTR", file=sys.stderr)
        args['model'] = 'GTR'
        if 'rates' in params:
            args['rates'] = params['rates']

    # --- 碱基频率 ---
    if 'freqs' in params:
        if base not in JC_FAMILY:  # JC 已经强制等频率
            args['freqs'] = params['freqs']

    # --- Gamma 速率异质性 ---
    if 'gamma_alpha' in params and 'gamma_cats' in params:
        args['gamma_alpha'] = params['gamma_alpha']
        args['gamma_cats'] = params['gamma_cats']
    elif 'gamma_alpha' in params:
        # 有 alpha 但没有类别数，默认 4
        args['gamma_alpha'] = params['gamma_alpha']
        args['gamma_cats'] = 4

    # --- 不变位点 ---
    if 'inv_prop' in params and params['inv_prop'] > 0:
        args['inv_prop'] = params['inv_prop']

    return args


# ============================================================
# 3. 构建 Seq-Gen 命令行
# ============================================================

def build_seqgen_cmd(gene_id, args, seq_length, tree_file, output_file,
                     n_reps=1, scale=1.0, seed=None, out_format='f'):
    """
    构建 seq-gen 命令行字符串。
    
    参数:
        gene_id      : 基因 ID
        args         : iqtree_to_seqgen() 返回的参数字典
        seq_length   : 序列长度
        tree_file    : 输入树文件路径
        output_file  : 输出文件路径 (None 表示 stdout)
        n_reps       : 重复数据集数 (-n)
        scale        : 分支长度缩放因子 (-s)
        seed         : 随机种子 (-z)
        out_format   : 输出格式 'p'=PHYLIP, 'n'=NEXUS, 'f'=FASTA
    """
    parts = ['seq-gen']

    # 模型
    parts.append(f'-m{args["model"]}')

    # GTR 速率
    if 'rates' in args:
        parts.append(f'-r{",".join(str(round(r, 6)) for r in args["rates"])}')

    # HKY 转换/颠换比
    if 'tstv' in args:
        parts.append(f'-t{args["tstv"]}')

    # 碱基频率
    if 'freqs' in args:
        if all(abs(f - 0.25) < 0.0001 for f in args['freqs']):
            parts.append('-fe')
        else:
            parts.append(f'-f{",".join(str(round(f, 6)) for f in args["freqs"])}')

    # Gamma
    if 'gamma_alpha' in args:
        parts.append(f'-a{args["gamma_alpha"]}')
    if 'gamma_cats' in args:
        parts.append(f'-g{args["gamma_cats"]}')

    # 不变位点
    if 'inv_prop' in args:
        parts.append(f'-i{args["inv_prop"]}')

    # 序列长度
    parts.append(f'-l{seq_length}')

    # 重复数
    if n_reps > 1:
        parts.append(f'-n{n_reps}')

    # 缩放
    if scale != 1.0:
        parts.append(f'-s{scale}')

    # 种子
    if seed is not None:
        parts.append(f'-z{seed}')

    # 输出格式
    parts.append(f'-o{out_format}')

    # 输入重定向
    cmd = ' '.join(parts)
    cmd += f' < {tree_file}'

    # 输出重定向
    if output_file:
        cmd += f' > {output_file}'

    return cmd


# ============================================================
# 4. 主程序
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description='从 IQ-TREE .iqtree 文件生成 Seq-Gen 模拟命令'
    )
    parser.add_argument('-g', '--gene-list', default='gene_list.txt',
                        help='基因列表文件 (默认: gene_list.txt)')
    parser.add_argument('-a', '--align-length', default='alignment_length.txt',
                        help='序列长度文件 (默认: alignment_length.txt)')
    parser.add_argument('-t', '--tree-dir', default='../../03-genetrees',
                        help='基因树目录 (默认: ../../03-genetrees)')
    parser.add_argument('--species-tree', default='rosa_orthofinder_MO_treeshrink_sp_rt_oneoutg_nonodelabel_nozero_final.tre',
                        help='物种树文件，所有基因统一使用 (默认: rosa_orthofinder_MO_treeshrink_sp_rt_oneoutg_nonodelabel_nozero_final.tre)')
    parser.add_argument('-o', '--output-dir', default='simulated_genes',
                        help='输出目录 (默认: simulated_genes)')
    parser.add_argument('-n', '--nreps', type=int, default=1,
                        help='每个基因生成的重复数据集数 (默认: 1)')
    parser.add_argument('-s', '--scale', type=float, default=1.0,
                        help='分支长度缩放因子 (默认: 1.0)')
    parser.add_argument('--seed-base', type=int, default=None,
                        help='随机种子基数 (默认: 使用基因ID)')
    parser.add_argument('--format', choices=['f', 'p', 'n'], default='f',
                        help='输出格式: f=FASTA, p=PHYLIP, n=NEXUS (默认: f)')
    parser.add_argument('--no-output-dir', action='store_true',
                        help='不创建输出目录，所有输出到 stdout')
    parser.add_argument('--dry-run', action='store_true',
                        help='仅打印命令，不写入文件')
    args = parser.parse_args()

    # 读取基因列表
    with open(args.gene_list, 'r') as f:
        genes = [line.strip() for line in f if line.strip()]
    print(f"# 读取到 {len(genes)} 个基因", file=sys.stderr)

    # 读取序列长度
    align_lengths = {}
    with open(args.align_length, 'r') as f:
        header = f.readline()  # 跳过表头
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                # 文件名格式: "OG0013459_1_noout_output.fasta" 或 "OG0013471_1_noout_output.tri.fasta"
                # 提取基因 ID（"_noout_output" 之前的部分）
                gene_name = os.path.basename(parts[0])
                gene_id = gene_name.split('_noout_output')[0]
                align_lengths[gene_id] = int(float(parts[1]))
    print(f"# 读取到 {len(align_lengths)} 个基因的序列长度", file=sys.stderr)

    # 准备输出
    if not args.dry_run and not args.no_output_dir:
        os.makedirs(args.output_dir, exist_ok=True)

    # 检查物种树文件
    if not os.path.exists(args.species_tree):
        print(f"ERROR: 物种树文件不存在: {args.species_tree}", file=sys.stderr)
        sys.exit(1)
    print(f"# 参考树: {args.species_tree}", file=sys.stderr)

    commands = []
    errors = []
    model_summary = {}

    for gene_id in genes:
        # 文件路径（使用通配符匹配 .iqtree 文件）
        iqtree_candidates = glob.glob(os.path.join(args.tree_dir, f'{gene_id}*.iqtree'))
        if not iqtree_candidates:
            errors.append(f'{gene_id}: .iqtree not found')
            continue
        iqtree_file = iqtree_candidates[0]  # 取第一个匹配

        # 检查文件存在
        if not os.path.exists(iqtree_file):
            errors.append(f'{gene_id}: .iqtree not found')
            continue
        if gene_id not in align_lengths:
            errors.append(f'{gene_id}: length not in alignment_length.txt')
            continue

        try:
            # 解析参数
            params = parse_iqtree(iqtree_file)
            seqgen_args = iqtree_to_seqgen(params)

            # 统计模型
            model = params['model']
            model_summary[model] = model_summary.get(model, 0) + 1

            # 输出文件
            if args.no_output_dir:
                output_file = None
            else:
                output_file = os.path.join(args.output_dir,
                                           f'{gene_id}_simulated.fasta')

            # 种子（使用基因 ID 的哈希值确保可复现）
            if args.seed_base is not None:
                seed = args.seed_base + hash(gene_id) % 100000
            else:
                seed = hash(gene_id) % 100000

            # 构建命令
            cmd = build_seqgen_cmd(
                gene_id, seqgen_args,
                seq_length=align_lengths[gene_id],
                tree_file=args.species_tree,
                output_file=output_file,
                n_reps=args.nreps,
                scale=args.scale,
                seed=seed,
                out_format=args.format,
            )
            commands.append((gene_id, params['model'], cmd))

        except Exception as e:
            errors.append(f'{gene_id}: {e}')

    # --- 输出 ---
    if not args.dry_run:
        # 写入 shell 脚本
        with open('seqgen_commands.sh', 'w') as f:
            f.write('#!/bin/bash\n')
            f.write(f'# Seq-Gen commands for {len(commands)} genes\n')
            f.write(f'# Generated by generate_seqgen_commands.py\n')
            f.write(f'# Model: IQ-TREE .iqtree -> Seq-Gen\n')
            f.write(f'# Scale: {args.scale}, Replicates: {args.nreps}\n')
            f.write(f'#\n')
            f.write(f'# Model summary:\n')
            for model, count in sorted(model_summary.items()):
                f.write(f'#   {model}: {count} genes\n')
            f.write(f'\n')
            f.write(f'mkdir -p {args.output_dir}\n\n')
            for gene_id, model, cmd in commands:
                f.write(f'# {gene_id}  (model: {model})\n')
                f.write(f'{cmd}\n\n')

        print(f'# 命令已写入: seqgen_commands.sh', file=sys.stderr)
    else:
        # 仅打印
        for gene_id, model, cmd in commands:
            print(f'# {gene_id} ({model})')
            print(cmd)
            print()

    # 模型统计
    print(f'\n# ========== 模型统计 ==========', file=sys.stderr)
    for model, count in sorted(model_summary.items(), key=lambda x: -x[1]):
        print(f'#   {model}: {count}', file=sys.stderr)

    # 错误报告
    if errors:
        print(f'\n# ========== 错误/跳过 ({len(errors)}) ==========', file=sys.stderr)
        for err in errors:
            print(f'#   {err}', file=sys.stderr)

    print(f'\n# 成功: {len(commands)}, 失败/跳过: {len(errors)}', file=sys.stderr)
    print(f'# 运行: bash seqgen_commands.sh', file=sys.stderr)


if __name__ == '__main__':
    main()