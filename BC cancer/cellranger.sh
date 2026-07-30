#!/usr/bin/env bash
set -e

# ==============================================================================
# 路径与资源配置
# ==============================================================================
FQ_DIR="/home/stu_wangyanan/working/BC_cancer/01.fq.gz"
OUT_DIR="/home/stu_wangyanan/working/BC_cancer/02.cellranger"
REF_DIR="/home/stu_wangyanan/working/BC_cancer/reference_genome/refdata-gex-GRCh38-2024-A/"

THREADS=16
MEM_GB=64

mkdir -p "$OUT_DIR"

ABS_FQ_DIR=$(cd "$FQ_DIR" && pwd)
ABS_OUT_DIR=$(cd "$OUT_DIR" && pwd)

# ==============================================================================
# 步骤 1: 在 01.fq.gz 文件夹内直接原地重命名文件（符合 10x 规范）
# ==============================================================================
echo "=================================================="
echo "[1/2] Renaming files in-place inside ${FQ_DIR}..."
echo "=================================================="

for r1_file in "$ABS_FQ_DIR"/*_1.fastq.gz; do
    if [ -f "$r1_file" ]; then
        sample=$(basename "$r1_file" _1.fastq.gz)
        r2_file="${ABS_FQ_DIR}/${sample}_2.fastq.gz"

        new_r1="${ABS_FQ_DIR}/${sample}_S1_L001_R1_001.fastq.gz"
        new_r2="${ABS_FQ_DIR}/${sample}_S1_L001_R2_001.fastq.gz"

        mv "$r1_file" "$new_r1"
        mv "$r2_file" "$new_r2"
        echo "Renamed: ${sample}_1.fastq.gz -> $(basename "$new_r1")"
    fi
done

# ==============================================================================
# 步骤 2: 直接在 01.fq.gz 上批量运行 Cell Ranger Count
# ==============================================================================
echo ""
echo "=================================================="
echo "[2/2] Starting Cell Ranger batch processing..."
echo "=================================================="

cd "$ABS_OUT_DIR"

for r1_file in "$ABS_FQ_DIR"/*_S1_L001_R1_001.fastq.gz; do
    sample=$(basename "$r1_file" _S1_L001_R1_001.fastq.gz)

    echo ""
    echo "--------------------------------------------------"
    echo "Processing sample: ${sample} (Time: $(date '+%Y-%m-%d %H:%M:%S'))"
    echo "--------------------------------------------------"

    # 已成功添加 --create-bam true 参数
    cellranger count \
        --id="${sample}" \
        --transcriptome="${REF_DIR}" \
        --fastqs="${ABS_FQ_DIR}" \
        --sample="${sample}" \
        --create-bam true \
        --localcores=${THREADS} \
        --localmem=${MEM_GB}

    echo "Sample ${sample} completed!"
done

echo ""
echo "=================================================="
echo "All samples processed successfully!"
echo "Results saved in: ${ABS_OUT_DIR}"
echo "=================================================="
