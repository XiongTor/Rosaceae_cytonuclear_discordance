#!/usr/bin/env julia
# Author: Tao Xiong
# Date: 2026-02-12
# Description: For displaying the results of PhyloNetworks analysis
#              Single-h version: h is passed in from the command line,
#              meant to be driven by an outer bash loop, e.g.:
#                for h in $(seq 1 6); do julia ../script/plot_result_julia.R $h; done
# 注意：本脚本内容为 Julia 代码，虽然文件名按要求存为 .R，
#      建议后续改名为 plot_result_julia.jl 以避免与 R 脚本混淆。

# ==== Main Script Start ====
using PhyloNetworks
using PhyloPlots  # 注意更新到最新版 julia, >1.9 否则可能报错
using RCall

# ---- 从命令行读取 h ----
if length(ARGS) < 1
    error("用法: julia plot_result_julia.R <h>   (h 为 1~6 的整数)")
end
h = parse(Int, ARGS[1])

infile        = "net$(h)_10runs.out"        # 输入的最佳网络树文件
outtree       = "bestnet_h$(h).tre"          # 输出的 dendroscope 可读文件
imagefilename = "snaqplot_net$(h)_root.svg"  # 输出的绘图文件

println("==== Processing h = $h  (input: $infile) ====")

if !isfile(infile)
    error("找不到文件 $infile ，请检查路径或先完成该组网络分析。")
end

net = readTopology(infile)  # 将最佳网络树存为变量 net

# 尝试定根；如果失败（比如该组网络里没有指定的外群，或定根后不满足
# level-1 网络的根方向要求），则打印警告并跳过定根，继续用原始（未定根）
# 的网络画图，而不中断整个脚本。
try
    rootatnode!(net, "Outgroup")  # 如需替换外群名称，改这里
    println("  已成功定根 (Outgroup)。")
catch e
    println("  警告: h=$h 定根失败 (", e, ")，将使用未定根的网络继续画图。")
end

writeTopology(net, outtree, di=true)  # 写出 dendroscope 可读形式

# 可以直接到 dendroscope 中可视化，或者按照如下流程在 julia 中继续操作
R"svg"(imagefilename, width=8, height=6)
R"par"(mar=[0, 0, 0, 0])
plot(net,
     showgamma=true,
     tipcex=0.9,
     xlim=[0, 30])  # 横向空间扩大
R"dev.off()"

println("已生成: $outtree , $imagefilename")