# 安装Seurat包（如果尚未安装）
if (!require("Seurat")) install.packages("Seurat")
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")

# 加载包
library(Seurat)
library(dplyr)
library(ggplot2)
# 读取表达矩阵
count_matrix <- ReadMtx(mtx = "matrix.mtx", cells = "barcodes.tsv", features = "features.tsv")

# 读取基因注释文件
features <- read.table("features.tsv", header = FALSE, sep = "\t", row.names = 1)

# 创建Seurat对象
seurat_obj <- CreateSeuratObject(counts = count_matrix, project = "Sample1")

# 安装并加载 readr 包
if (!require("readr")) install.packages("readr")
library(readr)

# 直接读取 barcodes.tsv.gz 文件
barcodes <- read_tsv("barcodes.tsv.gz", col_names = FALSE)

# 查看前几行
head(barcodes)

# 读取表达矩阵
count_matrix <- ReadMtx(mtx = "matrix.mtx", cells = "barcodes.tsv.gz", features = "features.tsv")

# 读取基因注释文件
features <- read.table("features.tsv", header = FALSE, sep = "\t", row.names = 1)

# 创建Seurat对象
seurat_obj <- CreateSeuratObject(counts = count_matrix, project = "Sample1")

# 计算线粒体基因比例
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# 绘制质量控制图
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# 过滤低质量细胞
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)

seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# 绘制高变基因图
top10 <- head(VariableFeatures(seurat_obj), 10)
plot1 <- VariableFeaturePlot(seurat_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
# 假设你的 Seurat 对象是 seurat_obj
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# 查看高变基因信息
hvf_info <- HVFInfo(seurat_obj)
head(hvf_info)

seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))

seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))

# 绘制PCA图
DimPlot(seurat_obj, reduction = "pca")

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# 绘制UMAP图
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters")

# 查找每个簇的标记基因
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# 查看每个簇的 top 标记基因
top_markers <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
print(top_markers)

# 绘制标记基因表达热图
DoHeatmap(seurat_obj, features = top_markers$gene) + NoLegend()

# 比较两个簇的差异表达基因
cluster1_vs_cluster2 <- FindMarkers(seurat_obj, ident.1 = 0, ident.2 = 1, min.pct = 0.25)

# 查看显著差异基因
head(cluster1_vs_cluster2)

# 绘制火山图
EnhancedVolcano(cluster1_vs_cluster2, lab = rownames(cluster1_vs_cluster2), 
                x = "avg_log2FC", y = "p_val_adj")

# 安装 EnhancedVolcano 包
if (!require("EnhancedVolcano")) BiocManager::install("EnhancedVolcano")

# 加载 EnhancedVolcano 包
library(EnhancedVolcano)
n
# 绘制火山图
EnhancedVolcano(cluster1_vs_cluster2, 
                lab = rownames(cluster1_vs_cluster2), 
                x = "avg_log2FC", 
                y = "p_val_adj")

# 安装并加载 clusterProfiler
if (!require("clusterProfiler")) install.packages("clusterProfiler")
library(clusterProfiler)

# 加载 clusterProfiler 包
library(clusterProfiler)

# 安装人类注释数据库
if (!require("org.Hs.eg.db")) BiocManager::install("org.Hs.eg.db")

# 设置 Bioconductor 镜像源为清华镜像
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

# 重新安装 org.Hs.eg.db
BiocManager::install("org.Hs.eg.db")

install.packages("AnnotationDbi")
install.packages("org.Hs.eg.db")
# 手动安装下载的包
install.packages("D:/org.Hs.eg.db_3.20.0.tar.gz", repos = NULL, type = "source")


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("AnnotationDbi")

packageVersion("AnnotationDbi")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("AnnotationDbi", force = TRUE)



install.packages("D:/AnnotationDbi_1.68.0.tar.gz", repos = NULL, type = "source")


BiocManager::install("org.Hs.eg.db")

force = TRUE
n

install.packages("AnnotationDbi")
# 加载人类注释数据库
library(org.Hs.eg.db)
n
# 提取显著差异基因
sig_genes <- rownames(subset(cluster1_vs_cluster2, p_val_adj < 0.05))

library(clusterProfiler)
# GO 富集分析
go_enrich <- enrichGO(gene = sig_genes, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
dotplot(go_enrich)

go_enrich <- enrichGO(gene = sig_genes, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "CC")
dotplot(go_enrich)

go_enrich <- enrichGO(gene = sig_genes, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "MF")
dotplot(go_enrich)








# 检查 sig_genes 的格式
head(sig_genes)

# 如果 sig_genes 是 Symbol，转换为 Entrez ID
if (all(grepl("^[A-Za-z]", sig_genes))) {
  library(org.Hs.eg.db)
  sig_genes <- mapIds(org.Hs.eg.db, keys = sig_genes, column = "ENTREZID", keytype = "SYMBOL")
  sig_genes <- na.omit(sig_genes)
}

# 如果 sig_genes 是 Ensembl ID，转换为 Entrez ID
if (all(grepl("^ENSG", sig_genes))) {
  library(org.Hs.eg.db)
  sig_genes <- mapIds(org.Hs.eg.db, keys = sig_genes, column = "ENTREZID", keytype = "ENSEMBL")
  sig_genes <- na.omit(sig_genes)
}

# 检查转换后的基因 ID
head(sig_genes)

# KEGG 富集分析
kegg_enrich <- enrichKEGG(gene = sig_genes, organism = "hsa", keyType = "kegg")
dotplot(kegg_enrich)


force = TRUE


# 安装并加载 Monocle
if (!require("monocle")) BiocManager::install("monocle")
packageVersion("monocle")
library(monocle)

library(monocle)
library(Seurat)

# 假设你的 Seurat 对象是 seurat_obj
# 将 Seurat 对象转换为 Monocle 对象
cds <- as.cell_data_set(seurat_obj)


exists("as.CellDataSet")

# 加载必要的包
library(monocle)
library(Seurat)

# 假设你的 Seurat 对象是 seurat_obj
# 将 Seurat 对象转换为 Monocle 对象
cds <- as.CellDataSet(seurat_obj)

# 检查转换后的 Monocle 对象
cds






# 将 Seurat 对象转换为 Monocle 对象
cds <- as.CellDataSet(seurat_obj)


seurat_obj <- NormalizeData(seurat_obj)

# 加载必要的包
library(monocle)
library(Seurat)

# 假设你的 Seurat 对象是 seurat_obj
# 标准化 Seurat 对象
seurat_obj <- NormalizeData(seurat_obj)


packageVersion("Seurat")
# 提取 Seurat 对象的表达矩阵
expression_matrix <- seurat_obj@assays$RNA$counts


cell_metadata <- seurat_obj@meta.data
gene_metadata <- data.frame(gene_short_name = rownames(seurat_obj))

rownames(expression_matrix)
rownames(gene_metadata)

rownames(gene_metadata) <- rownames(expression_matrix)
dim(expression_matrix)
dim(gene_metadata)
# 创建 Monocle 对象
cds <- newCellDataSet(
  expression_matrix,
  phenoData = AnnotatedDataFrame(cell_metadata),
  featureData = AnnotatedDataFrame(gene_metadata),
  expressionFamily = negbinomial.size()
)

# 设置大小因子
cds <- estimateSizeFactors(cds)

# 检查转换后的 Monocle 对象
cds

dim(cds)



# 运行轨迹分析
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds <- orderCells(cds)

# 绘制轨迹图
plot_cell_trajectory(cds, color_by = "State")

# 选择高变基因
cds <- detectGenes(cds)
disp_table <- dispersionTable(cds)
hvg_genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)
cds <- setOrderingFilter(cds, hvg_genes$gene_id)

# 使用 UMAP 降维
cds <- reduceDimension(cds, max_components = 2, method = 'UMAP')

# 运行轨迹分析
cds <- orderCells(cds)


# 可视化轨迹
plot_cell_trajectory(cds, color_by = "State")  # 按状态着色
plot_cell_trajectory(cds, color_by = "Pseudotime")  # 按伪时间着色

state_markers <- findMarkers(cds, group_col = "State")

# 加载必要的包
library(monocle)
library(Seurat)

# 假设 cds 是已经创建的单细胞数据集
# 检查细胞类型
head(pData(cds)$cell_type)

# 如果没有细胞注释，使用标记基因推断
plot_cell_trajectory(cds, markers = c("CD3D", "CD14", "CD19"), use_color_gradient = TRUE)



# 运行轨迹分析
cds <- orderCells(cds)

# 可视化轨迹
plot_cell_trajectory(cds, color_by = "State")  # 按状态着色
plot_cell_trajectory(cds, color_by = "Pseudotime")  # 按伪时间着色












if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("HSMMSingleCell")
options(BioC_mirror = "https://bioconductor.org")
n
BiocManager::install("HSMMSingleCell")
# 更换 Bioconductor 镜像
options(BioC_mirror = "https://bioconductor.org")
n
# 安装 HSMMSingleCell
BiocManager::install("HSMMSingleCell")

# 如果仍然失败，手动下载并安装
install.packages("D:/HSMMSingleCell_1.26.0.tar.gz", repos = NULL, type = "source")
library(monocle)
n

# 安装并加载 CellChat
if (!require("CellChat")) devtools::install_github("sqjin/CellChat")
library(CellChat)


packageVersion("Seurat")
# 创建 CellChat 对象
cellchat <- createCellChat(object = seurat_obj, group.by = "seurat_clusters")
library(CellChat)
exists("createCellChat")




# 加载必要的包
library(Seurat)
library(CellChat)

# 假设你的 Seurat 对象是 seurat_obj
# 创建 CellChat 对象
cellchat <- createCellChat(object = seurat_obj, group.by = "seurat_clusters")

# 检查创建的 CellChat 对象
cellchat

# 加载必要的包
library(Seurat)
library(CellChat)

# 假设你的 Seurat 对象是 seurat_obj
# 检查 seurat_clusters 标签
unique(seurat_obj$seurat_clusters)

# 修改 seurat_clusters 标签
seurat_obj$seurat_clusters <- paste0("Cluster", seurat_obj$seurat_clusters)

# 创建 CellChat 对象
cellchat <- createCellChat(object = seurat_obj, group.by = "seurat_clusters")

# 检查创建的 CellChat 对象
cellchat

# 加载必要的包
library(CellChat)


# 假设你的 Seurat 对象是 seurat_obj
# 创建 CellChat 对象
cellchat <- createCellChat(object = seurat_obj, group.by = "seurat_clusters")



print(cellchat)


str(cellchat)



# 运行数据预处理
cellchat <- subsetData(cellchat)
# 识别过表达的基因
cellchat <- identifyOverExpressedGenes(cellchat)

# 识别过表达的配体-受体对
cellchat <- identifyOverExpressedInteractions(cellchat)

library(CellChat)

# 创建 CellChat 对象
cellchat <- createCellChat(data = your_data, meta = your_meta)

# 加载数据库
cellchat@DB <- CellChatDB.human  # 或 CellChatDB.mouse

# 运行 subsetData
cellchat <- subsetData(cellchat)

# 检查 data.signaling
head(cellchat@data.signaling)

# 识别过表达的基因
cellchat <- identifyOverExpressedGenes(cellchat)

# 识别过表达的配体-受体对
cellchat <- identifyOverExpressedInteractions(cellchat)


# 检查 data.signaling 数据
head(cellchat@data.signaling)

# 运行细胞通讯分析
cellchat <- identifyOverExpressedGenes(cellchat)

# 检查分析结果
cellchat
















str(cellchat)
# 假设你的 CellChat 对象是 cellchat
# 运行数据预处理

# 检查 data.signaling 数据
head(cellchat@data.signaling)

# 运行细胞通讯分析
cellchat <- identifyOverExpressedGenes(cellchat)

# 检查分析结果
cellchat

# 运行细胞通讯分析
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

# 可视化细胞通讯网络
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize)

# 计算 groupSize
groupSize <- as.numeric(table(cellchat@idents))
print(groupSize)

# 可视化细胞通讯网络
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize)
cellchat@idents <- as.factor(cellchat@meta$seurat_clusters)

cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
print(groupSize)
print(cellchat@idents)
num_nodes <- nrow(cellchat@net$weight)
print(num_nodes)

print(head(cellchat@meta))
print(cellchat@idents)
# 假设元数据中有一列 celltype 表示细胞类型
cell_types <- cellchat@meta$celltype

# 将 cellchat@idents 映射到细胞类型
cluster_to_celltype <- data.frame(
  Cluster = levels(cellchat@idents),
  CellType = tapply(cell_types, cellchat@idents, function(x) names(sort(table(x), decreasing = TRUE))[1])
)

print(cluster_to_celltype)

print(length(cell_types))
print(length(cellchat@idents))

print(colnames(cellchat@meta))
cell_types <- cellchat@meta$seurat_clusters
print(head(cell_types))
cluster_to_celltype <- data.frame(
  Cluster = levels(cellchat@idents),
  CellType = tapply(cell_types, cellchat@idents, function(x) names(sort(table(x), decreasing = TRUE))[1])
)
print(cluster_to_celltype)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SingleR")
BiocManager::install("celldex")
n

library(SingleR)
library(celldex)

# 加载参考数据集
ref_data <- HumanPrimaryCellAtlasData()

# 进行注释
annotations <- SingleR(test = seurat_object@assays$RNA@data, ref = ref_data, labels = ref_data$label.main)
seurat_object$celltype <- annotations$labels
print(table(seurat_object$celltype))



library(Seurat)

# 加载数据
data <- Read10X(data.dir = "D:/GSE175453_RAW (1)")  # 替换为你的数据路径
seurat_object <- CreateSeuratObject(counts = data, project = "YourProject")

print(seurat_object)

print(head(rownames(seurat_object@assays$RNA@counts)))
# 数据标准化
seurat_object <- NormalizeData(seurat_object)

counts_data <- GetAssayData(seurat_object, assay = "RNA", slot = "counts")
print(head(rownames(counts_data)))




# 特征选择
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)

# 数据缩放
seurat_object <- ScaleData(seurat_object)

# 降维（PCA）
seurat_object <- RunPCA(seurat_object, npcs = 50)

# 聚类
seurat_object <- FindNeighbors(seurat_object, dims = 1:20)
seurat_object <- FindClusters(seurat_object, resolution = 0.5)