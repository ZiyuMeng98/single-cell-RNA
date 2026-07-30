#!/bin/bash

RAW_DIR="/home/stu_wangyanan/working/BC_cancer/00.raw_data"
OUT_DIR="/home/stu_wangyanan/working/BC_cancer/01.fq.gz"
MAX_JOBS=4

# 自动创建输出文件夹 01.fq.gz
mkdir -p "${OUT_DIR}"

echo "========== Processing begin (Running ${MAX_JOBS} samples at the same time) =========="

for SAMPLE_DIR in "${RAW_DIR}"/*/; do
    [ -d "${SAMPLE_DIR}" ] || continue
    SAMPLE_NAME=$(basename "${SAMPLE_DIR}")

    # 后台并发执行
    (
        # 匹配找到该样本目录下的 .sra 文件
        SRA_FILE=$(ls "${SAMPLE_DIR}"/*.sra 2>/dev/null | head -n 1)

        if [ -n "${SRA_FILE}" ]; then
            SRA_NAME=$(basename "${SRA_FILE}" .sra)
            echo ">>> 开始处理样本 [${SAMPLE_NAME}] (${SRA_NAME}.sra)"

            # 1. 使用 -O 参数直接输出解压结果到 01.fq.gz 文件夹下
            fasterq-dump "${SRA_FILE}" --split-files --threads 2 -O "${OUT_DIR}"

            # 2. 仅对该样本在 01.fq.gz 中生成的 FASTQ 进行 gzip 压缩
            gzip -f "${OUT_DIR}/${SRA_NAME}"*.fastq

            echo "<<< 完成样本 [${SAMPLE_NAME}]"
        fi
    ) &

    # 控制后台同时运行的任务不超过 4 个
    while [ $(jobs -r | wc -l) -ge ${MAX_JOBS} ]; do
        sleep 2
    done
done

# 等待所有样本处理完成
wait

echo "========== Processing finished =========="
