# 使用cellranger生成表达矩阵，生成后导入R语言Seruat库进行分析
bash srr_to_fq.sh 2>&1 | tee -a srr_to_fq.log
bash cellranger.sh 2>&1 | tee -a cellranger.log
# R语言Seruat库
library(Seurat)
library(ggplot2)

set.seed(1234)

meta <- data.frame(
  sample_id = c("SRR12603790", "SRR12603789", "SRR12603787","SRR12603786", "SRR12603785", "SRR12603784","SRR12603783", "SRR12603782", "SRR12603781","SRR12603780", "SRR12603788"),
  patient = c("P01", "P02", "P03", "P04", "P05", "P06","P07", "P08", "N01", "N02", "N03"),
  tissue = c(rep("Tumor", 8), rep("Normal", 3))
)
dir.create("03.seurat", showWarnings = FALSE)
objects <- list()
qc_summary <- data.frame()
# 读取、质控、标准化和高变基因筛选
for (i in seq_len(nrow(meta))) {
  sid <- meta$sample_id[i]
  matrix_path <- file.path("02.cellranger",sid,"outs","filtered_feature_bc_matrix")
  counts <- Read10X(data.dir = matrix_path)
  obj <- CreateSeuratObject(counts = counts,project = sid,min.cells = 3)
  obj$sample_id <- sid
  obj$patient <- meta$patient[i]
  obj$tissue <- meta$tissue[i]
  obj[["percent.mt"]] <- PercentageFeatureSet(obj,pattern = "^MT-")
  cells_before <- ncol(obj)
  obj <- subset(obj,subset = nCount_RNA >= 1000 & percent.mt <= 10 & nFeature_RNA <= 6000)
  cells_after <- ncol(obj)
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000,verbose = FALSE)
  obj <- FindVariableFeatures(obj,selection.method = "vst",nfeatures = 3000,verbose = FALSE)
  objects[[sid]] <- obj
  qc_summary <- rbind(qc_summary,data.frame(sample_id = sid,patient = meta$patient[i],tissue = meta$tissue[i],cells_before = cells_before,cells_after = cells_after,retained_percent = round(100 * cells_after / cells_before,2),median_UMI = median(obj$nCount_RNA),median_genes = median(obj$nFeature_RNA),median_percent_mt = median(obj$percent.mt))
  )
}
