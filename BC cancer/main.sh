# 使用cellranger生成表达矩阵，生成后导入R语言Seruat库进行分析
bash srr_to_fq.sh 2>&1 | tee -a srr_to_fq.log
bash cellranger.sh 2>&1 | tee -a cellranger.log
# R语言Seruat库
library(Seurat)
library(ggplot2)
set.seed(1234)
# 输入文件
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
# 保存质控、标准化结果
saveRDS(objects,"03.seurat/all_samples_qc_norm.rds")
write.csv(qc_summary,"03.seurat/QC_summary.csv",row.names = FALSE)

# 选择3000个用于进一步分析的基因
features <- SelectIntegrationFeatures(object.list = objects,nfeatures = 3000)
length(features) # [1] 3000
head(features) # 显示基因名
# 寻找样本间的anchors（主要用来识别和消除批次效应，或实现跨模态（如 RNA 与 ATAC）的数据对齐）
options(future.globals.maxSize = 8 * 1024^3)
anchors <- FindIntegrationAnchors(object.list = objects , anchor.features = features , dims = 1:30)
# 整合所有样本
bladder <- IntegrateData(anchorset = anchors, dims = 1:30)
# 如果出现Warning: Layer counts isn't present in the assay object; returning NULL
# 即为Seurat v5 中，Assay 采用了“图层”（Layers，如 counts、data、scale.data）机制。整合函数在运行时会尝试调取原始的 counts 图层。如果你在此之前对对象进行了数据提取、筛选，或者直接基于归一化数据创建了 Assay，导致对象中缺少原始计数图层，系统就会返回 NULL
# 这个报错不影响后续分析
# 或者Warning: Different cells in new layer data than already exists for scale.data
# 即为当前 data 图层中的细胞数量/标签与已有的 scale.data 图层不匹配。通常发生在合并（Merge）数据集或对细胞进行子集筛选（Subset）后，没有重新运行 ScaleData() 的情况
# 这个报错可能影响后续分析

# 检查整合结果
bladder
DefaultAssay(bladder)
table(bladder$sample_id)
table(bladder$tissue)
# 保存整合对象
saveRDS(bladder,"03.seurat/bladder_integrated_raw.rds")
# 用 integrated assay 做降维和聚类
DefaultAssay(bladder) <- "integrated"
bladder <- ScaleData(bladder, verbose = FALSE)
bladder <- RunPCA(bladder, npcs = 30, verbose = FALSE)
# 查看每个主成分解释的变异
ElbowPlot(bladder, ndims = 30)  # 看图在什么位置变平缓，作为PC的取值（文章中复现结果拐点在8-10，但是按照文章的选了30）
# 使用前 30 个 PC
bladder <- FindNeighbors(bladder, dims = 1:30)
bladder <- FindClusters(bladder,resolution = 0.8)
bladder <- RunUMAP(bladder,dims = 1:30)
# 查看聚类结果
DimPlot(bladder, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
# 查看样本是否被批次效应主导
DimPlot(bladder, reduction = "umap", group.by = "sample_id")
# 查看肿瘤和正常组织分布
DimPlot(bladder, reduction = "umap", group.by = "tissue")
# 保存当前对象
saveRDS(bladder, "03.seurat/bladder_clustered.rds")
# 找各聚类的标记基因。注意这里要切换回 RNA assay，因为差异表达应使用原始/标准化表达数据，不用 integrated assay
DefaultAssay(bladder) <- "RNA"
bladder <- NormalizeData(bladder, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
markers <- FindAllMarkers(bladder, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(markers)
colnames(markers)
# 查看每个 cluster 的前 20 个 marker
library(dplyr)
top20_markers <- markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 20) %>% ungroup()
top20_markers
# 只看基因名
top20_markers %>% select(cluster, gene, avg_log2FC, pct.1, p_val_adj)
# 绘制 marker 热图
DoHeatmap(bladder, features = unique(top20_markers$gene), group.by = "seurat_clusters")
# 导出
write.csv(
  top20_markers,
  "03.seurat/top20_markers_each_cluster.csv",
  row.names = FALSE
)
# 根据这些 marker 注释主要细胞类型
marker_list <- list(
  Epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "KRT20", "UPK1A", "UPK2", "UPK3A", "GATA3", "FOXA1"),
  T_cell = c("CD3D", "CD3E", "CD3G", "TRBC1", "TRBC2", "IL7R", "CD4", "CD8A", "CD8B", "NKG7", "GNLY", "CCL5"),
  B_cell = c("CD79A", "CD79B", "MS4A1", "CD37", "CD74", "HLA-DRA", "CD19", "CD22"),
  Myeloid = c("LYZ", "LST1", "FCER1G", "TYROBP", "AIF1", "CTSS", "CD68", "CD163", "C1QC", "FCN1", "S100A8", "S100A9"),
  Mast = c("TPSAB1", "TPSB2", "KIT", "MS4A2", "HDC", "CPA3"),
  Endothelial = c("PECAM1", "VWF", "EMCN", "KDR", "ESAM", "CDH5", "RAMP2", "PLVAP"),
  Fibroblast = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A1", "COL6A2", "PDGFRA", "CFD"),
  iCAF = c("PDGFRA", "CXCL12", "IL6", "CXCL14", "CXCL1", "CXCL2", "CXCL8", "CCL2", "CFD"),
  mCAF = c("RGS5", "ACTA2", "TAGLN", "MYL9", "TPM2", "PDGFRB", "MCAM", "CSPG4", "DES")
)
# 确认当前对象和聚类
DefaultAssay(bladder) <- "RNA"
table(bladder$seurat_clusters)
levels(bladder$seurat_clusters)
# 保留表达矩阵中真实存在的基因
marker_list <- lapply(
  marker_list,
  function(x) {
    intersect(x, rownames(bladder))
  }
)
sapply(marker_list, length)
# 查看各 cluster 的 marker
DotPlot(
  bladder,
  features = marker_list,
  group.by = "seurat_clusters",
  dot.scale = 6
) +
  RotatedAxis()
# 保存图
pdf(
  "03.seurat/major_celltype_markers.pdf",
  width = 18,
  height = 10
)

print(
  DotPlot(
    bladder,
    features = marker_list,
    group.by = "seurat_clusters",
    dot.scale = 6
  ) +
    RotatedAxis()
)

dev.off()
# 查看每个 cluster 的真实差异基因
DefaultAssay(bladder) <- "RNA"

markers <- FindAllMarkers(
  bladder,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
# 查看每个 cluster 前 20 个 marker
library(dplyr)

top20_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 20
  ) %>%
  ungroup()

top20_markers
# 保存文件
write.csv(
  top20_markers,
  "03.seurat/top20_markers_each_cluster.csv",
  row.names = FALSE
)
# FeaturePlot 重点核对大类细胞
FeaturePlot(
  bladder,
  features = c(
    "EPCAM",
    "PTPRC",
    "COL1A1",
    "PECAM1",
    "CD3D",
    "CD79A",
    "LYZ",
    "TPSAB1",
    "PDGFRA",
    "CXCL12",
    "RGS5"
  ),
  reduction = "umap",
  ncol = 4
)
# 主要表达	                     初步类型
# EPCAM, KRT8, KRT18, KRT19    上皮细胞
# CD3D, CD3E, TRBC1	            T 细胞
# CD79A, MS4A1, CD74	          B 细胞
# LYZ, LST1, TYROBP	            髓系细胞
# TPSAB1, KIT, CPA3	            肥大细胞
# PECAM1, VWF, EMCN	            内皮细胞
# COL1A1, COL1A2, DCN, LUM	    成纤维细胞
# COL1A1 + PDGFRA + CXCL12	    iCAF 倾向
# COL1A1 + RGS5 + ACTA2 + TAGLN	mCAF 倾向



# 建立 cluster 注释表
cluster_annotation <- data.frame(
  cluster = levels(bladder$seurat_clusters),
  celltype = NA_character_,
  stringsAsFactors = FALSE
)

cluster_annotation
# 这时不要直接把 0、1、2 等编号套用到论文。你必须根据自己的 DotPlot 和 top20_markers 判断
# 例如，如果观察到
# cluster 0：EPCAM、KRT8、KRT18、KRT19 高
# cluster 1：CD3D、CD3E、TRBC1 高
# cluster 2：LYZ、LST1、TYROBP 高
# cluster 3：COL1A1、DCN、LUM、PDGFRA 高
# 才可以这样填写
# cluster_annotation$celltype[
#   cluster_annotation$cluster == "0"
# ] <- "Epithelial"
# cluster_annotation$celltype[
#   cluster_annotation$cluster == "1"
# ] <- "T cell"
# cluster_annotation$celltype[
#   cluster_annotation$cluster == "2"
# ] <- "Myeloid"
# cluster_annotation$celltype[
#   cluster_annotation$cluster == "3"
# ] <- "Fibroblast"

# 把注释写回 Seurat 对象
cluster_to_celltype <- setNames(
  cluster_annotation$celltype,
  cluster_annotation$cluster
)

bladder$celltype <- unname(
  cluster_to_celltype[
    as.character(bladder$seurat_clusters)
  ]
)
# 检查结果
table(
  bladder$celltype,
  useNA = "ifany"
)
# 如果存在 NA，说明有 cluster 还没有完成注释，需要补充
cluster_annotation
# 绘制注释后的UMAP图
DimPlot(
  bladder,
  reduction = "umap",
  group.by = "celltype",
  label = TRUE,
  repel = TRUE
)
# 按组织来源查看细胞类型
DimPlot(
  bladder,
  reduction = "umap",
  group.by = "tissue"
)
# 保存文件
saveRDS(
  bladder,
  "03.seurat/bladder_annotated.rds"
)
