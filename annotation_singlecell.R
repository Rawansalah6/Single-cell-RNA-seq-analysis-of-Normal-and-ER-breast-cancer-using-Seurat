#Manual annotation vs Automated annotation
#Manual annotation
#define canonical marker genes for the main cell types
marker_list = list(
  
  Epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "MUC1"),
  
  Basal_Myoepithelial = c("KRT5", "KRT14", "KRT17", "TP63", "COL17A1",
                          "ACTA2", "TAGLN"),
  
  Luminal = c("KRT8", "KRT18", "KRT19", "EPCAM",
              "MUC1", "KRT7"),
  
  Luminal_Secretory = c("PIP", "SCGB1D2", "SCGB2A1",
                        "SCGB2A2", "SCGB3A1", "MUCL1"),
  
  Fibroblast = c("COL1A1", "COL1A2", "COL3A1",
                 "DCN", "LUM", "OGN", "PDGFRA"),
  
  Endothelial = c("PECAM1", "VWF", "KDR", "EMCN",
                  "SOX18", "ACKR1", "ENG"),
  
  T_cell = c("CD3D", "CD3E", "TRBC1", "TRBC2",
             "CD8A", "CD8B"),
  
  NK_cell = c("NKG7", "GNLY", "KLRD1",
              "FCGR3A", "XCL1", "XCL2"),
  
  Macrophage_Myeloid = c("LYZ", "CSF1R", "C1QC",
                         "LILRB1", "LILRB2", "LILRB4",
                         "TYROBP", "FCER1G"),
  
  B_cell = c("CD79A", "MS4A1", "CD37",
             "CD74", "HLA-DRA", "CD79B"),
  
  Plasma_cell = c("MZB1", "JCHAIN", "SDC1",
                  "XBP1", "CD38"),
  
  Cycling = c("MKI67", "TOP2A", "CENPA",
              "PBK", "ASPM", "NEK2"))
#visualize canonical marker expression across all clusters
DotPlot(
  seurat_integrated,
  features = unique(unlist(marker_list)),
  group.by = "seurat_clusters") +
  RotatedAxis() +
  theme(
    axis.text.x = element_text(size = 3),
    axis.text.y = element_text(size = 8),
    legend.title = element_text(
      size = 8,
      margin = margin(t = 15))) +
  ggtitle("Canonical marker expression across clusters")

#examine epithelial and epithelial subtype marker expression
epithelial_markers = c("EPCAM", "KRT8", "KRT18", "KRT19", "MUC1", "KRT7", "PIP",
                       "MUCL1", "SCGB1D2", "SCGB2A1", "SCGB2A2", "SCGB3A1", "PAEP", "KRT5",
                       "KRT14", "KRT17", "TP63", "COL17A1")

DotPlot(
  seurat_integrated,
  features = epithelial_markers,
  group.by = "seurat_clusters") +
  RotatedAxis() +
  theme(
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 8),
    legend.title = element_text(
      size = 8,
      margin = margin(t = 15))) +
  ggtitle("Epithelial marker expression across clusters")

#validate non-epithelial cell identities using lineage-specific markers
validation_markers = c("PECAM1", "VWF", "EMCN", "KDR", "ACKR1", "SOX18", "DLL4",
                       
                       "RGS5", "PDGFRB", "MCAM", "CSPG4", "ACTA2", "TAGLN",
                       
                       "CD3D", "CD3E", "TRBC1", "CD8A", "CD8B",
                       
                       "NKG7", "GNLY", "KLRD1", "LYZ", "TYROBP", "FCER1G", "CSF1R", "C1QC",
                       
                       "CD79A", "MS4A1", "CD74", "HLA-DRA", "MZB1", "JCHAIN", "SDC1")

DotPlot(
  seurat_integrated,
  features = validation_markers,
  group.by = "seurat_clusters") +
  RotatedAxis() +
  theme(
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 8),
    legend.title = element_text(
      size = 8,
      margin = margin(t = 15))) +
  ggtitle("Validation markers across clusters")
#additional validation of selected clusters
#validate the proliferative state and cell identity of cluster 12
s.genes = cc.genes$s.genes
g2m.genes = cc.genes$g2m.genes

seurat_integrated = CellCycleScoring(
  seurat_integrated,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = FALSE)

#examine S and G2/M cell cycle scores in cluster 12
VlnPlot(
  seurat_integrated,
  features = c("S.Score", "G2M.Score"),
  idents = "12")
#determine the cell cycle phase distribution in cluster 12
table(
  seurat_integrated$Phase[
    Idents(seurat_integrated) == "12" ])
#examine the top marker genes identified for cluster 12
cluster12.markers = seurat_markers %>%
  filter(cluster == 12)

head(cluster12.markers, 20)
#to validate the epithelial identity of cluster 12
FeaturePlot(
  seurat_integrated,
  features = c("EPCAM", "FDCSP", "LAMA3"),
  cells = WhichCells(seurat_integrated, idents = "12"))
#To validate the endothelial identity of cluster 13
endothelial_genes = c("PECAM1", "CLDN5", "VWF", "KDR", "EMCN", "ENG")

DotPlot(
  seurat_integrated,
  features = endothelial_genes,
  idents = "13") +
  RotatedAxis()
#define arterial endothelial markers
arterial_genes = c("GJA5", "HEY1", "EFNB2", "SOX17")
#define lymphatic endothelial markers
lymphatic_genes = c("PROX1", "LYVE1", "PDPN", "FLT4")
#compare arterial and lymphatic endothelial markers in cluster 13
DotPlot(
  seurat_integrated,
  features = c(arterial_genes, lymphatic_genes),
  idents = "13") +
  RotatedAxis()

#To check whether luminal clusters are affected by sample specific effects, cluster 0,1,4,6,9,10
DimPlot(
  subset(
    seurat_integrated,
    idents = c("0", "1", "4", "6", "9", "10")),
  group.by = "seurat_clusters",
  split.by = "orig.ident",
  ncol = 4)

table(
  seurat_integrated$seurat_clusters[
    seurat_integrated$seurat_clusters %in%
      c("0", "1", "4", "6", "9", "10") ],
  seurat_integrated$orig.ident[
    seurat_integrated$seurat_clusters %in%
      c("0", "1", "4", "6", "9", "10") ])
#define canonical markers for luminal epithelial cells
luminal_genes = c("EPCAM", "KRT8", "KRT18", "KRT19", "MUC1", "KRT7")
luminal_state_genes = c("ESR1", "PGR", "FOXA1", "GATA3", "ELF5", "MKI67")

DotPlot(
  seurat_integrated,
  features = c(
    luminal_genes,
    luminal_state_genes),
  idents = c("0", "1", "4", "6", "9", "10")) +
  RotatedAxis()

#validate cell identities using canonical marker expression
DotPlot(
  seurat_integrated,
  features = c("EPCAM","KRT8", "KRT18", "KRT19", "KRT14", "TP63",
               "COL17A1", "PECAM1", "CLDN5", "VWF", "COL1A1", "COL1A2", "DCN", "LST1",
               "CD68", "CD3D", "CD3E", "MS4A1", "CD79A", "NKG7", "GNLY")) +
  RotatedAxis()

#finally the manual annotation
#to save the original cluster numbers 
seurat_integrated$original_cluster = as.character(
  Idents(seurat_integrated))
colnames(seurat_integrated@meta.data)

annotation = c(
  "Luminal epithelial",
  "Luminal epithelial state",
  "Fibroblasts",
  "Basal/myoepithelial",
  "Mature luminal epithelial",
  "T/NK cells",
  "Luminal progenitor-like epithelial",
  "Pericytes/perivascular",
  "Endothelial",
  "Luminal secretory epithelial",
  "Luminal epithelial state",
  "Myeloid cells",
  "Proliferating epithelial",
  "Arterial/arteriolar endothelial",
  "B/plasma cells")

names(annotation) = as.character(0:14)

seurat_integrated$cell_type = unname(
  annotation[
    as.character(seurat_integrated$original_cluster) ])
"cell_type" %in% colnames(seurat_integrated@meta.data)
table(
  seurat_integrated$original_cluster,
  seurat_integrated$cell_type)

DimPlot(
  seurat_integrated,
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE,
  label.size = 1.8,
  pt.size = 0.2) +
  theme(
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7))


#########################################
#Automated annotation
#SingleR cell level
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SingleR")
BiocManager::install("BiocFileCache", ask = FALSE, update = FALSE)
BiocManager::install("celldex")
library(SingleR)
library(celldex)
library(Seurat)
DefaultAssay(seurat_integrated) = "RNA"
DefaultAssay(seurat_integrated)
#load the reference dataset
ref = HumanPrimaryCellAtlasData()

Layers(seurat_integrated[["RNA"]])
#Join RNA layers for singleR
seurat_singleR = seurat_integrated
seurat_singleR[["RNA"]] = JoinLayers(
  seurat_singleR[["RNA"]])

Layers(seurat_singleR[["RNA"]])
#Extract normalized RNA
test = GetAssayData(
  seurat_singleR,
  assay = "RNA",
  layer = "data")

dim(test)
#Run singleR at the cell level
pred = SingleR(
  test = test,
  ref = ref,
  labels = ref$label.main)
head(pred)
table(pred$labels)


seurat_singleR$SingleR_label = pred$labels

library(ggplot2)
DimPlot(
  seurat_singleR,
  group.by = "SingleR_label",
  label = TRUE,
  repel = TRUE,
  label.size = 1.8,
  pt.size = 0.2) +
  theme(
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 6))
#compare manual and automated annotations
table(
  seurat_singleR$cell_type,
  seurat_singleR$SingleR_label)

#SingleR cluster level annotation
BiocManager::install("scrapper", ask = FALSE, update = FALSE)
library(scrapper)

cluster_singleR = table(
  seurat_singleR$original_cluster,
  seurat_singleR$SingleR_label)

cluster_singleR
#Assign the dominant singleR label to each cluster
cluster_dominant = apply(
  cluster_singleR,
  1,
  function(x) names(which.max(x)))

cluster_dominant
#map cluster  level labels back to individual cells
cluster_labels_per_cell = cluster_dominant[
  as.character(seurat_singleR$original_cluster) ]

names(cluster_labels_per_cell) = colnames(seurat_singleR)

seurat_singleR$SingleR_cluster_label = cluster_labels_per_cell

DimPlot(
  seurat_singleR,
  group.by = "SingleR_cluster_label",
  label = TRUE,
  repel = TRUE,
  label.size = 1.8,
  pt.size = 0.2) +
  theme(
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 6))
