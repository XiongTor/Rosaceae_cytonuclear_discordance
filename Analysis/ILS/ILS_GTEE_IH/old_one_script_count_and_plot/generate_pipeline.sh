#!/bin/bash

# ==================== 路径与文件配置 ====================
PARAMS_FILE="speed_distribution.txt"
PXCAT_FILE="rosa_orthofinder_partition.txt"
SPECIES_TREE="rosa_orthofinder_sptree_rebranch_rt_oneoutgroup.tre"
OUTPUT_SH="run_seq_gen.sh"
# ========================================================

echo "#!/bin/bash" > $OUTPUT_SH
echo "mkdir -p sim_genes_perfect" >> $OUTPUT_SH
echo "" >> $OUTPUT_SH

echo "正在生成带真实 Gene ID 的 seq-gen 模拟脚本..."

grep -E '^[[:space:]]*[0-9]+' "$PARAMS_FILE" | while read -r line; do
    
    part_id=$(echo "$line" | awk '{print $1}')
    model_name=$(echo "$line" | awk '{print $2}')
    speed=$(echo "$line" | awk '{print $3}')
    param_str=$(echo "$line" | awk '{for(i=4;i<=NF;i++) printf "%s", $i; print ""}')

    pxcat_line=$(sed -n "${part_id}p" "$PXCAT_FILE")
    [[ -z "$pxcat_line" ]] && continue

    # 1. 计算长度
    start_end=$(echo "$pxcat_line" | awk -F'=' '{print $2}' | tr -d '[:space:]')
    start_pos=$(echo "$start_end" | awk -F'-' '{print $1}')
    end_pos=$(echo "$start_end" | awk -F'-' '{print $2}')
    real_length=$((end_pos - start_pos + 1))
    
    # 2. ⚡【关键修正】提取真实的核心 Gene ID
    raw_fasta=$(echo "$pxcat_line" | awk -F'=' '{print $1}' | sed 's/DNA,//g' | tr -d '[:space:]' | xargs basename)
    gene_id=$(echo "$raw_fasta" | sed 's/\.tri\.fasta//g' | sed 's/\.fasta//g')

    # 3. 提取模型参数
    rates_raw=$(echo "$param_str" | grep -oE '(GTR|TPM3u|TPM2u|TIM3|TN|HKY)\{[^}]+\}' | sed -E 's/.*\{([^}]+)\}/\1/')
    freqs_raw=$(echo "$param_str" | grep -oE 'F\{[^}]+\}' | sed -E 's/.*\{([^}]+)\}/\1/')
    inv_raw=$(echo "$param_str" | grep -oE 'I\{[^}]+\}' | sed -E 's/.*\{([^}]+)\}/\1/')
    gamma_raw=$(echo "$param_str" | grep -oE 'G4\{[^}]+\}' | sed -E 's/.*\{([^}]+)\}/\1/')
    [[ -z "$gamma_raw" ]] && gamma_raw=$(echo "$param_str" | grep -oE 'R[3-5]\{[^,]+' | sed -E 's/.*\{([^,]+)/\1/')

    freqs_cmd=""; [[ ! -z "$freqs_raw" ]] && freqs_cmd="-f $(echo "$freqs_raw" | tr ',' ' ')"
    inv_cmd="";   [[ ! -z "$inv_raw" ]]   && inv_cmd="-i $inv_raw"
    gamma_cmd=""; [[ ! -z "$gamma_raw" ]] && gamma_cmd="-a $gamma_raw -g 4"

    r1=$(echo "$rates_raw" | cut -d',' -f1); r2=$(echo "$rates_raw" | cut -d',' -f2)
    r3=$(echo "$rates_raw" | cut -d',' -f3); r4=$(echo "$rates_raw" | cut -d',' -f4)
    r5=$(echo "$rates_raw" | cut -d',' -f5)

    gtr_rates="1.0 1.0 1.0 1.0 1.0 1.0"
    if [[ "$model_name" == *"GTR"* ]]; then         gtr_rates="$r1 $r2 $r3 $r4 $r5 1.0"
    elif [[ "$model_name" == *"TPM3u"* ]]; then     gtr_rates="1.0 $r1 $r2 1.0 $r1 1.0"
    elif [[ "$model_name" == *"TPM2u"* ]]; then     gtr_rates="1.0 $r1 1.0 $r2 $r1 1.0"
    elif [[ "$model_name" == *"TIM3"* ]]; then      gtr_rates="1.0 $r1 $r2 1.0 $r3 1.0"
    elif [[ "$model_name" == *"HKY"* ]]; then       gtr_rates="1.0 $r1 1.0 1.0 $r1 1.0"
    fi

    # 输出到以真实 gene_id 命名的文件中
    echo "seq-gen -mGTR -r $gtr_rates -l $real_length -s $speed $freqs_cmd $inv_cmd $gamma_cmd < $SPECIES_TREE > sim_genes_perfect/${gene_id}.phy" >> $OUTPUT_SH

done

echo "🎉 成功！请检查并运行 sh $OUTPUT_SH 生成全物种完美序列。"