library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
setwd("E:/ScRNAseq_AdvanzaBio")
untar("GSE161529_RAW.tar", exdir = "GSE161529_raw")
list.files("GSE161529_raw")

list.files("GSE161529_raw", pattern = "features")
file.exists("GSE161529_features.tsv.gz")

#sample 1 ER+
counts_GSM4909296 = ReadMtx(
    mtx = "GSE161529_raw/GSM4909296_ER-MH0001-matrix.mtx.gz",
   cells = "GSE161529_raw/GSM4909296_ER-MH0001-barcodes.tsv.gz",
   features = "GSE161529_features.tsv.gz")
dim(counts_GSM4909296)
counts_GSM4909296[1:5, 1:5]

er1 = CreateSeuratObject(
   counts = counts_GSM4909296,
   project = "ER1",
   min.cells = 3,
   min.features = 200)
er1


#sample 2 ER+
counts_GSM4909297 = ReadMtx(
   mtx = "GSE161529_raw/GSM4909297_ER-MH0125-matrix.mtx.gz",
   cells = "GSE161529_raw/GSM4909297_ER-MH0125-barcodes.tsv.gz",
   features = "GSE161529_features.tsv.gz")
dim(counts_GSM4909297)

er2 = CreateSeuratObject(
   counts = counts_GSM4909297,
   project = "ER2",
   min.cells = 3,
   min.features = 200)


#sample 3 ER+
counts_GSM4909298 = ReadMtx(
   mtx = "GSE161529_raw/GSM4909298_ER-PM0360-matrix.mtx.gz",
   cells = "GSE161529_raw/GSM4909298_ER-PM0360-barcodes.tsv.gz",
   features = "GSE161529_features.tsv.gz")
dim(counts_GSM4909298)

er3 = CreateSeuratObject(
   counts = counts_GSM4909298,
   project = "ER3",
   min.cells = 3,
   min.features = 200)

#Sample 4 ER+
counts_GSM4909299 = ReadMtx(
   mtx = "GSE161529_raw/GSM4909299_ER-MH0114-T3-matrix.mtx.gz",
   cells = "GSE161529_raw/GSM4909299_ER-MH0114-T3-barcodes.tsv.gz",
   features = "GSE161529_features.tsv.gz")
dim(counts_GSM4909299)

er4 = CreateSeuratObject(
   counts = counts_GSM4909299,
   project = "ER4",
   min.cells = 3,
   min.features = 200)



#sample 1 control
counts_GSM4909253 = ReadMtx(
   mtx = "GSE161529_raw/GSM4909253_N-PM0092-Total-matrix.mtx.gz",
   cells = "GSE161529_raw/GSM4909253_N-PM0092-Total-barcodes.tsv.gz",
   features = "GSE161529_features.tsv.gz")
dim(counts_GSM4909253)


con1 = CreateSeuratObject(
   counts = counts_GSM4909253,
   project = "con1",
   min.cells = 3,
   min.features = 200)


#sample 2 control
counts_GSM4909254 = ReadMtx(
   mtx = "GSE161529_raw/GSM4909254_N-PM0019-Total-matrix.mtx.gz",
   cells = "GSE161529_raw/GSM4909254_N-PM0019-Total-barcodes.tsv.gz",
   features = "GSE161529_features.tsv.gz")
dim(counts_GSM4909254)

con2 = CreateSeuratObject(
   counts = counts_GSM4909254,
   project = "con2",
   min.cells = 3,
   min.features = 200)

#sample 3 Control
counts_GSM4909257 = ReadMtx(
   mtx = "GSE161529_raw/GSM4909257_N-PM0095-Total-matrix.mtx.gz",
   cells = "GSE161529_raw/GSM4909257_N-PM0095-Total-barcodes.tsv.gz",
   features = "GSE161529_features.tsv.gz")
dim(counts_GSM4909257)

con3 = CreateSeuratObject(
   counts = counts_GSM4909257,
   project = "con3",
   min.cells = 3,
   min.features = 200)



#sample 4 Control
counts_GSM4909261 = ReadMtx(
   mtx = "GSE161529_raw/GSM4909261_N-PM0230-Total-matrix.mtx.gz",
   cells = "GSE161529_raw/GSM4909261_N-PM0230-Total-barcodes.tsv.gz",
   features = "GSE161529_features.tsv.gz")
dim(counts_GSM4909261)


con4 = CreateSeuratObject(
   counts = counts_GSM4909261,
   project = "con4",
   min.cells = 3,
   min.features = 200)


#Make a list to make it easier instead of running the code 8 times
seurat_list = list(
   er1 = er1,
   er2 = er2,
   er3 = er3,
   er4 = er4,
   con1 = con1,
   con2 = con2,
   con3 = con3,
   con4 = con4)


#QC and seclecting cells for furthur analysis
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
seurat_list = lapply(seurat_list, function(x) {
   x[["percent.mt"]] = PercentageFeatureSet(x, pattern = "^MT-")
   x
 })


#visualization before filtration
# Visualize QC metrics as a violin plot
for (i in 1:length(seurat_list)) {
   print(
     VlnPlot(
       seurat_list[[i]],
       features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
       ncol = 3)) }

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
library(ggplot2)
for (i in 1:length(seurat_list)) {
   plot1 = FeatureScatter(
     seurat_list[[i]],
     feature1 = "nCount_RNA",
     feature2 = "percent.mt")
   
   plot2 = FeatureScatter(
     seurat_list[[i]],
     feature1 = "nCount_RNA",
     feature2 = "nFeature_RNA")
   
   print(plot1 + plot2 + plot_annotation(title = names(seurat_list)[i])) }


#Filtration
# QC thresholds for each sample
qc_thresholds = list(
  con1 = list(feature_min = 200, feature_max = 7000, mt_max = 20),
  con2 = list(feature_min = 200, feature_max = 7500, mt_max = 25),
  con3 = list(feature_min = 200, feature_max = 7000, mt_max = 30),
  con4 = list(feature_min = 200, feature_max = 8000, mt_max = 20),
  er1  = list(feature_min = 200, feature_max = 2500, mt_max = 20),
  er2  = list(feature_min = 200, feature_max = 6000, mt_max = 20),
  er3  = list(feature_min = 200, feature_max = 8000, mt_max = 25),
  er4  = list(feature_min = 200, feature_max = 6500, mt_max = 20))

# Filter cells using sample-specific QC thresholds
seurat_list_filtered = lapply(names(seurat_list), function(name) {
  
  thr = qc_thresholds[[name]]
  
  subset(
    seurat_list[[name]],
    subset = nFeature_RNA > thr$feature_min &
      nFeature_RNA < thr$feature_max &
      percent.mt < thr$mt_max) })

# Keep sample names
names(seurat_list_filtered) = names(seurat_list)

# Number of cells before filtering
sapply(seurat_list, ncol)

# Number of cells after filtering
sapply(seurat_list_filtered, ncol)


#Revisualize after filtration
for (i in 1:length(seurat_list_filtered)) { print(
    VlnPlot(
      seurat_list_filtered[[i]],
      features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
      ncol = 3)) }

#that only to check the number of cells on er3 beacuse the quality is not good (sensitivity check)
sum(seurat_list[["er3"]]$percent.mt < 15)
sum(seurat_list[["er3"]]$percent.mt < 20)
sum(seurat_list[["er3"]]$percent.mt < 25)
ncol(seurat_list[["er3"]])


#To edit the percentage of mt of er3 but keep the original filtration the same and make a new one
# Refine the QC threshold for ER3 after reviewing the QC distribution
qc_thresholds_refined = qc_thresholds

qc_thresholds_refined$er3$mt_max = 20


# Re-filter all samples using the refined QC thresholds
seurat_list_filtered_refined = lapply(names(seurat_list), function(name) {
  
  thr = qc_thresholds_refined[[name]]
  
  subset(
    seurat_list[[name]],
    subset = nFeature_RNA > thr$feature_min &
      nFeature_RNA < thr$feature_max &
      percent.mt < thr$mt_max) })

# Keep the original sample names
names(seurat_list_filtered_refined) = names(seurat_list)


# Compare the number of cells before and after refined filtering
before = sapply(seurat_list, ncol)
after_refined = sapply(seurat_list_filtered_refined, ncol)

data.frame(
  sample = names(before),
  before = before,
  after = after_refined,
  pct_retained = round(100 * after_refined / before, 1))
#To visualize the er3 after the final filtration
VlnPlot(
  seurat_list_filtered_refined[["er3"]],
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3)


#Normalization
# LogNormalize each sample separately
seurat_list_normalized = lapply(
  seurat_list_filtered_refined,
  function(x) {NormalizeData(
      x,
      normalization.method = "LogNormalize",
      scale.factor = 10000) })
saveRDS(seurat_list_normalized, "seurat_normalized.rds")
#Identification of highly variable features
seurat_list_hvf = lapply(
  seurat_list_normalized,
  function(x) {FindVariableFeatures(
      x,
      selection.method = "vst",
      nfeatures = 2000) })
saveRDS(seurat_list_hvf, "seurat_list_hvf.rds")
# Identify the 10 most highly variable genes for each sample
top10 = lapply(
  seurat_list_hvf,
  function(x) head(VariableFeatures(x), 10))

top10

# Plot variable features with  labels
# Plot variable features with sample names
for (i in 1:length(seurat_list_hvf)) {
  
  plot1 = VariableFeaturePlot(seurat_list_hvf[[i]]) +
    ggtitle(names(seurat_list_hvf)[i])
  
  plot2 = LabelPoints(
    plot = plot1,
    points = head(VariableFeatures(seurat_list_hvf[[i]]), 10),
    repel = TRUE)
  
  print(plot2) }

#######################################################

#for integration(selection the common features for integration)
integration_features = SelectIntegrationFeatures(
  object.list = seurat_list_hvf,
  nfeatures = 2000)

#find integration anchors
anchors = FindIntegrationAnchors(
  object.list = seurat_list_hvf,
  anchor.features = integration_features)
saveRDS(anchors, "anchors.rds")
#For the integration
anchors = readRDS("anchors.rds")
class(anchors)
seurat_integrated = IntegrateData(
  anchorset = anchors)
saveRDS(seurat_integrated, "seurat_integrated.rds")

#Scaling
seurat_integrated = ScaleData(
  seurat_integrated,
  features = rownames(seurat_integrated))


#perform linear dimensional reduction
seurat_integrated = RunPCA(
  seurat_integrated,
  features = integration_features)
saveRDS(seurat_integrated, "seurat_integrated_PCA.rds")
seurat_integrated = readRDS("seurat_integrated_PCA.rds")
#Elbowplot
ElbowPlot(seurat_integrated, ndims = 50)


#Clusters the cells
seurat_integrated = FindNeighbors(
  seurat_integrated,
  dims = 1:15)

seurat_integrated = FindClusters(
  seurat_integrated,
  resolution = 0.5)

# Run non-linear dimensional reduction (UMAP)
seurat_integrated = RunUMAP(
  seurat_integrated,
  dims = 1:15)

#UMAP
DimPlot(
  seurat_integrated,
  reduction = "umap",
  label = TRUE)

DimPlot(
  seurat_integrated,
  reduction = "umap",
  group.by = "orig.ident")

# find markers for every cluster compared to all remaining cells, report only the positive ones
seurat_markers = FindAllMarkers(
  seurat_integrated,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25)


library(dplyr)
top10_markers = seurat_markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 0.5) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10) %>%
  ungroup()

top10_markers
top10_markers %>%
  arrange(cluster, desc(avg_log2FC))
#heatmap for top 10 markers
DoHeatmap(
  seurat_integrated,
  features = top10_markers$gene) + NoLegend()


#top 5 markers
top5_markers = seurat_markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 0.5) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 5) %>%
  ungroup()

top5_markers %>%
  select(cluster, gene) %>%
  group_by(cluster) %>%
  summarise(
    markers = paste(gene, collapse = ", "))
library(ggplot2)
DoHeatmap(
  seurat_integrated,
  features = top5_markers$gene) +
  theme(
    axis.text.y = element_text(size = 4),
    axis.text.x = element_text(
      size = 5,
      angle = 0,
      hjust = 0.5,
      vjust = 0.5),
    axis.ticks.x = element_blank(),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8))
ggsave("top5_heatmap.png", width = 12, height = 10, dpi = 300)

