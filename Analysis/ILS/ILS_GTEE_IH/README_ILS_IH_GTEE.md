#### 代码文件夹说明;
new_use_in_global: 该文件夹用于存放ILS_GTEE_IH相关的全局代码，包括所有的步骤
new_use_in_GTEE_part： 该文件夹用于存放ILS_GTEE_IH相关的GTEE部分代码

old_one_script_count_and_plot:老版本的代码，可能存在与新的代码相同的情况
old_onw_raw_script：:老版本的代码，可能存在与新的代码相同的情况

#### 新老版本代码说明
新版本代码，主要是分析流程有变化，最大的变化在于：

**GTEE的分析部分**
主要是seq-gen对应的部分，此前主要是通过超矩阵重新估算ASTRAL树的过程中，生成的.iqtree文件中包含了各个分区，即各个基因的碱基替换速率等信息。seq-gen的分析部分主要是通过这些.iqtree文件，提取各个分区的碱基替换速率等信息，进而进行seq-gen模拟。
新版本代码中，这些替换速率改成了从各个单基因树的.iqtree文件中提取，而不通过超矩阵的分区文件了
主要涉及的脚本为：generate_seqgen_commands.py 以及 clip_missing_taxa.sh

