bash srr_to_fq.sh 2>&1 | tee -a srr_to_fq.log
bash cellranger.sh 2>&1 | tee -a cellranger.log
