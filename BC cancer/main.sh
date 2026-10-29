##############################################################################
# 通用分析
# QC → 标准化 → 整合 → PCA → 聚类 → UMAP → marker → 大类注释
##############################################################################
library(Seurat)
library(ggplot2)
library(dplyr)

set.seed(1234)
options(future.globals.maxSize = 8 * 1024^3)

dir.create("03.seurat", showWarnings = FALSE)
dir.create("r_tmp", showWarnings = FALSE)

Sys.setenv(TMPDIR = file.path(getwd(), "r_tmp"))
Sys.setenv(TMP = file.path(getwd(), "r_tmp"))
Sys.setenv(TEMP = file.path(getwd(), "r_tmp"))
##############################################################################
# 1. 样本信息
##############################################################################
meta <- data.frame(
  sample_id = c("SRR12603790", "SRR12603789", "SRR12603787", "SRR12603786", "SRR12603785", "SRR12603784", "SRR12603783", "SRR12603782", "SRR12603781", "SRR12603780", "SRR12603788"),
  patient = c("P01", "P02", "P03", "P04", "P05", "P06", "P07", "P08", "N01", "N02", "N03"),
  tissue = c(rep("Tumor", 8), rep("Normal", 3))
)
##############################################################################
# 2. 读取、QC、标准化和高变基因
##############################################################################
objects <- list()
qc_summary <- data.frame()

for (i in seq_len(nrow(meta))) {
  sid <- meta$sample_id[i]
  matrix_path <- file.path("02.cellranger", sid, "outs", "filtered_feature_bc_matrix")
  counts <- Read10X(data.dir = matrix_path)
  # 建立 Seurat 对象
  obj <- CreateSeuratObject(counts = counts, project = sid, min.cells = 3)
  # 保存样本信息
  obj$sample_id <- sid
  obj$patient <- meta$patient[i]
  obj$tissue <- meta$tissue[i]
  # 计算线粒体基因比例
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  cells_before <- ncol(obj)
  # QC 条件
  obj <- subset(obj, subset = nCount_RNA >= 1000 & percent.mt <= 10 & nFeature_RNA <= 6000)
  cells_after <- ncol(obj)
  # LogNormalize 标准化
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  # 每个样本筛选 3000 个高变基因
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  objects[[sid]] <- obj
  # 记录 QC 结果
  qc_summary <- rbind(qc_summary,
    data.frame(sample_id = sid, patient = meta$patient[i], tissue = meta$tissue[i], cells_before = cells_before, cells_after = cells_after,
      retained_percent = round(cells_after / cells_before * 100, 2),
      median_UMI = median(obj$nCount_RNA),
      median_genes = median(obj$nFeature_RNA),
      median_percent_mt = median(obj$percent.mt)
    )
  )
}
saveRDS(objects, "03.seurat/all_samples_qc_norm.rds")
write.csv(qc_summary, "03.seurat/QC_summary.csv", row.names = FALSE)
##############################################################################
# 3. 样本整合
# 论文：3000 个特征基因，使用 30 个维度
##############################################################################
features <- SelectIntegrationFeatures(object.list = objects, nfeatures = 3000)
anchors <- FindIntegrationAnchors(object.list = objects, anchor.features = features, dims = 1:30)
bladder <- IntegrateData(anchorset = anchors, dims = 1:30)
saveRDS(bladder, "03.seurat/bladder_integrated_raw.rds")
##############################################################################
# 4. PCA、邻居图、聚类和 UMAP
##############################################################################
DefaultAssay(bladder) <- "integrated"
bladder <- ScaleData(bladder, verbose = FALSE)
bladder <- RunPCA(bladder, npcs = 30,verbose = FALSE)
pdf("03.seurat/elbowplot.pdf", width = 6, height = 5)

print(ElbowPlot( bladder, ndims = 30))
dev.off()
bladder <- FindNeighbors(bladder, dims = 1:30)
bladder <- FindClusters(bladder, resolution = 0.8)
bladder <- RunUMAP(bladder, dims = 1:30)
pdf("03.seurat/umap_clusters.pdf", width = 8, height = 7)
print(DimPlot(bladder, reduction = "umap", group.by = "seurat_clusters", label = TRUE))
dev.off()
pdf("03.seurat/umap_sample.pdf", width = 8, height = 7)
print(DimPlot(bladder, reduction = "umap", group.by = "sample_id"))
dev.off()
pdf("03.seurat/umap_tissue.pdf", width = 8, height = 7)
print(DimPlot(bladder, reduction = "umap", group.by = "tissue"))
dev.off()
saveRDS(bladder, "03.seurat/bladder_clustered.rds")
##############################################################################
# 5. 查找每个 cluster 的 marker
##############################################################################
DefaultAssay(bladder) <- "RNA"
# Seurat v5 中，差异表达前合并 RNA layers
bladder <- JoinLayers(bladder, assay = "RNA")
bladder <- NormalizeData(bladder, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
markers <- FindAllMarkers(bladder, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top20_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 20
  ) %>%
  ungroup()

write.csv(
  markers,
  "03.seurat/all_cluster_markers.csv",
  row.names = FALSE
)

write.csv(
  top20_markers,
  "03.seurat/top20_markers_each_cluster.csv",
  row.names = FALSE
)
