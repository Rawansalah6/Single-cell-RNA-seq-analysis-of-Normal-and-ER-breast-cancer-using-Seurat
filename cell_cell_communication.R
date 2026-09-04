install.packages("devtools")

pkgbuild::check_build_tools(debug = TRUE)
install.packages("NMF")
install.packages("circlize")
BiocManager::install("ComplexHeatmap")
devtools::install_github("jinworks/CellChat")
library(CellChat)

seurat_integrated = readRDS("seurat_integrated_PCA.rds")
seurat_integrated = readRDS("seurat_singleR")

# Prepare CellChat input
# Get normalized RNA expression from the joined object
data.input = GetAssayData(
  seurat_singleR,
  assay = "RNA",
  layer = "data")

# Get metadata from the manually annotated object
meta = seurat_integrated@meta.data

# Add manual cell-type annotation to metadata
meta$cell_type = seurat_integrated$cell_type

# Make sure metadata follows the expression matrix cell order
meta = meta[colnames(data.input), , drop = FALSE]

# Make sure cell_type is character
meta$cell_type = as.character(meta$cell_type)

# Check that everything matches

cat("Expression cells:", ncol(data.input), "\n")
cat("Metadata cells:", nrow(meta), "\n")

cat(
  "Cell names match:",
  identical(colnames(data.input), rownames(meta)),
  "\n")

cat(
  "Missing cell types:",
  sum(is.na(meta$cell_type)),
  "\n")

# Check cell-type distribution
print(table(meta$cell_type))

# Save CellChat input data

saveRDS(
  data.input,
  "CellChat_expression.rds",
  compress = FALSE)

saveRDS(
  meta,
  "CellChat_metadata.rds")

# Check that files were created
file.info(
  c("CellChat_expression.rds",
    "CellChat_metadata.rds"))[, c("size")]
data.input = readRDS("CellChat_expression.rds")
meta = readRDS("CellChat_metadata.rds")

# Create CellChat object
cellchat = createCellChat(
  object = data.input,
  meta = meta,
  group.by = "cell_type")

cellchat

# CellChat analysis - Human
# 1. Load CellChat database for human
CellChatDB = CellChatDB.human
CellChatDB

# 2. Use the full database
cellchat@DB = CellChatDB

# 3. Preprocess expression data
cellchat = subsetData(cellchat)

# 4. Identify overexpressed genes
devtools::install_github("immunogenomics/presto")
library(presto)
cellchat = identifyOverExpressedGenes(cellchat)

# 5. Identify overexpressed ligand-receptor interactions
cellchat = identifyOverExpressedInteractions(cellchat)
cellchat


# Compute cell-cell communication probability
cellchat = computeCommunProb(
  cellchat, type = "triMean")

# Remove interactions with very few cells
cellchat = filterCommunication(
  cellchat, min.cells = 10)

cellchat
saveRDS(cellchat, "CellChat_after_communication.rds")

cellchat = computeCommunProbPathway(cellchat)
cellchat = aggregateNet(cellchat)
saveRDS(cellchat, "CellChat_after_pathways.rds")

# Make the plotting area larger
par(mfrow = c(1, 1), mar = c(1, 1, 2, 1))

groupSize = as.numeric(table(cellchat@idents))

netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  vertex.label.cex = 0.6,
  title.name = "Number of interactions")

groupSize = as.numeric(table(cellchat@idents))

par(mfrow = c(1, 1), mar = c(1, 1, 2, 1))

netVisual_circle(
  cellchat@net$weight,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  vertex.label.cex = 0.6,
  title.name = "Interaction strength")

args(netVisual_heatmap)

netVisual_heatmap(cellchat,
  measure = "count",
  font.size = 5)

netVisual_heatmap(cellchat,
  measure = "weight",
  font.size = 5)

rankNet(cellchat,
  mode = "comparison",
  stacked = TRUE,
  do.stat = TRUE)

#############################################
#To compare between normal and ER+
# Create metadata for Normal cells
# Create metadata for Control/Normal cells
meta_normal = meta[
  meta$orig.ident %in% c("con1", "con2", "con3", "con4"),
  ,
  drop = FALSE]

# Create metadata for ER+ cells
meta_er = meta[
  meta$orig.ident %in% c("ER1", "ER2", "ER3", "ER4"),
  ,
  drop = FALSE]

nrow(meta_normal)
nrow(meta_er)
table(meta_normal$orig.ident)
table(meta_er$orig.ident)

# Extract normalized expression for Control/Normal cells
data.input_normal = data.input[
  ,
  rownames(meta_normal)]

# Extract normalized expression for ER+ cells
data.input_er = data.input[
  ,
  rownames(meta_er)]


identical(colnames(data.input_normal),
  rownames(meta_normal))

identical(colnames(data.input_er),
  rownames(meta_er))

# Create CellChat object for Control/Normal cells
cellchat_normal = createCellChat(
  object = data.input_normal,
  meta = meta_normal,
  group.by = "cell_type")
cellchat_normal


cellchat_normal@DB = CellChatDB.human
cellchat_normal = subsetData(cellchat_normal)
cellchat_normal = identifyOverExpressedGenes(cellchat_normal)
cellchat_normal = identifyOverExpressedInteractions(cellchat_normal)


cellchat_normal = computeCommunProb(
  cellchat_normal,
  type = "triMean")


cellchat_normal = filterCommunication(
  cellchat_normal,
  min.cells = 10)

cellchat_normal = computeCommunProbPathway(
  cellchat_normal)


cellchat_normal = aggregateNet(
  cellchat_normal)

###############
# Create CellChat object for ER+ cells
cellchat_er = createCellChat(
  object = data.input_er,
  meta = meta_er,
  group.by = "cell_type")
cellchat_er


cellchat_er@DB = CellChatDB.human


cellchat_er = subsetData(cellchat_er)


cellchat_er = identifyOverExpressedGenes(cellchat_er)


cellchat_er = identifyOverExpressedInteractions(cellchat_er)


cellchat_er = computeCommunProb(
  cellchat_er,
  type = "triMean")


cellchat_er = filterCommunication(
  cellchat_er,
  min.cells = 10)

cellchat_er = computeCommunProbPathway(cellchat_er)


cellchat_er = aggregateNet(cellchat_er)

dim(cellchat_er@net$count)
dim(cellchat_er@net$weight)



# Create a list containing the two CellChat objects
cellchat_list = list(
  Normal = cellchat_normal,
  ER = cellchat_er)

# Merge Normal and ER+ CellChat objects
cellchat_compare <- mergeCellChat(
  object.list = list(
    Normal = cellchat_normal,
    ER = cellchat_er),
  add.names = c("Normal", "ER"))

# Compare signaling pathways between Normal and ER+
rankNet(
  cellchat_compare,
  mode = "comparison",
  stacked = TRUE,
  do.stat = TRUE, font.size = 4)
args(rankNet)

# Compare number of interactions between Normal and ER+
netVisual_heatmap(
  cellchat_compare,
  comparison = c(1, 2),
  measure = "count",
  font.size = 5)

# Compare interaction strength between Normal and ER+
netVisual_heatmap(
  cellchat_compare,
  comparison = c(1, 2),
  measure = "weight",
  font.size = 5)
