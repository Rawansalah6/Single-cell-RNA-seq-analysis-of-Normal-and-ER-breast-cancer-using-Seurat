#Preparation of SC files
deconv_sc_celltype_table = data.frame(
  cell_id = rownames(pbulk_annotation_table),
  cell_type = pbulk_annotation_table$cell_type)
all(deconv_sc_celltype_table$cell_id == colnames(pbulk_raw_count_matrix))
#for bulk
bulk_expression_matrix <- read.delim(
  "Galaxy251-[Gene level summarization on dataset 250 and collection 238].tabular",
  header = TRUE, sep = "\t", check.names = FALSE)
dim(bulk_expression_matrix)

head(bulk_expression_matrix)

colnames(bulk_expression_matrix)

bulk_ensembl_ids = rownames(bulk_expression_matrix)

bulk_ensembl_ids_clean = sub("\\..*$", "",
  bulk_ensembl_ids)

head(bulk_ensembl_ids_clean)

library(rtracklayer)

gencode_gtf = import(
  "Galaxy250-[gencode.v50.annotation.gtf.gz].gtf.gz")
gencode_gtf

gene_mapping_table = unique(
  data.frame(gene_id = gencode_gtf$gene_id,
    gene_symbol = gencode_gtf$gene_name))

gene_mapping_table = gene_mapping_table[
  !is.na(gene_mapping_table$gene_id) &
    !is.na(gene_mapping_table$gene_symbol), ]

head(gene_mapping_table)

gene_mapping_table$gene_id_clean = sub(
  "\\..*$",
  "",
  gene_mapping_table$gene_id)

mapped_gene_count = sum(
  bulk_ensembl_ids_clean %in% gene_mapping_table$gene_id_clean)

total_bulk_gene_count = length(bulk_ensembl_ids_clean)

mapped_gene_count
total_bulk_gene_count

mapped_gene_percentage = 100 * mapped_gene_count / total_bulk_gene_count

mapped_gene_percentage

bulk_gene_mapping = data.frame(
  ensembl_id = bulk_ensembl_ids_clean,
  stringsAsFactors = FALSE)

bulk_gene_mapping$gene_symbol = gene_mapping_table$gene_symbol[
  match(bulk_gene_mapping$ensembl_id,
    gene_mapping_table$gene_id_clean) ]

sum(is.na(bulk_gene_mapping$gene_symbol))

head(bulk_gene_mapping)

sum(duplicated(bulk_gene_mapping$gene_symbol))

duplicate_symbols = unique(
  bulk_gene_mapping$gene_symbol[
    duplicated(bulk_gene_mapping$gene_symbol) ])
head(duplicate_symbols, 20)
bulk_gene_mapping[
  bulk_gene_mapping$gene_symbol == duplicate_symbols[1],]

duplicate_symbol_table = bulk_gene_mapping[
  duplicated(bulk_gene_mapping$gene_symbol) |
    duplicated(bulk_gene_mapping$gene_symbol, fromLast = TRUE), ]

head(duplicate_symbol_table, 30)
length(unique(duplicate_symbol_table$gene_symbol))
head(rownames(pbulk_raw_count_matrix))

# Prepare bulk matrix for deconvolution
# Keep protein-coding genes from the same GTF
protein_coding_mapping = gene_mapping_table[
  gene_mapping_table$gene_id_clean %in%
    unique(sub("\\..*$", "",
        gencode_gtf$gene_id[gencode_gtf$gene_type == "protein_coding"])), ]

# Keep only bulk genes that are protein-coding
bulk_protein_coding_mapping = bulk_gene_mapping[
  bulk_gene_mapping$ensembl_id %in%
    protein_coding_mapping$gene_id_clean, ]

# Remove duplicated Gene Symbols
# Keep one Ensembl gene ID per Gene Symbol
bulk_protein_coding_mapping = bulk_protein_coding_mapping[
  !duplicated(bulk_protein_coding_mapping$gene_symbol), ]

# Extract corresponding bulk expression values
bulk_protein_coding_matrix = bulk_expression_matrix[
  match(bulk_protein_coding_mapping$ensembl_id,
    bulk_ensembl_ids_clean),
  ,
  drop = FALSE ]

# Replace Ensembl row names with Gene Symbols
rownames(bulk_protein_coding_matrix) =
  bulk_protein_coding_mapping$gene_symbol
dim(bulk_protein_coding_matrix)
sum(duplicated(rownames(bulk_protein_coding_matrix)))

# Find genes shared between bulk and scRNA
common_deconv_genes = intersect(
  rownames(bulk_protein_coding_matrix),
  rownames(pbulk_raw_count_matrix))
length(common_deconv_genes)

# Keep only shared genes
bulk_deconv_matrix = bulk_protein_coding_matrix[
  common_deconv_genes,
  ,
  drop = FALSE ]
sc_deconv_matrix = pbulk_raw_count_matrix[
  common_deconv_genes,
  ,
  drop = FALSE]
dim(bulk_deconv_matrix)
dim(sc_deconv_matrix)
all(rownames(bulk_deconv_matrix) ==
      rownames(sc_deconv_matrix))
all(colnames(sc_deconv_matrix) ==
      deconv_sc_celltype_table$cell_id)
table(deconv_sc_celltype_table$cell_type)
all(
  colnames(sc_deconv_matrix) ==
    deconv_sc_celltype_table$cell_id)
length(unique(deconv_sc_celltype_table$cell_type))
table(deconv_sc_celltype_table$cell_type)


# Clean Environment for MuSiC
rm(list = setdiff(ls(),
    c("bulk_deconv_matrix",
      "sc_deconv_matrix",
      "deconv_sc_celltype_table")))
gc()
#MuSiC installation
BiocManager::install("SingleCellExperiment", ask = FALSE, update = FALSE)
BiocManager::install("TOAST", ask = FALSE, update = FALSE)
devtools::install_github("xuranw/MuSiC")
library(MuSiC)
#Another preparation after using MUsic
seurat_integrated = readRDS("seurat_integrated_PCA.rds")
deconv_sc_metadata = data.frame(
  cell_id = colnames(sc_deconv_matrix),
  cell_type = deconv_sc_celltype_table$cell_type,
  sample_id = seurat_integrated@meta.data[
    colnames(sc_deconv_matrix),
    "orig.ident"
  ],
  stringsAsFactors = FALSE,
  row.names = colnames(sc_deconv_matrix))
head(deconv_sc_metadata)
rm(seurat_integrated)
gc()
#create singlecellexperiment
library(SingleCellExperiment)
sc_deconv_sce = SingleCellExperiment(
  assays = list(
    counts = sc_deconv_matrix))
#  Add cell type information
colData(sc_deconv_sce)$cell_type =
  deconv_sc_metadata[
    colnames(sc_deconv_sce),
    "cell_type" ]
#  Add sample ID information
colData(sc_deconv_sce)$sample_id =
  deconv_sc_metadata[
    colnames(sc_deconv_sce),
    "sample_id" ]
sc_deconv_sce
table(is.na(colData(sc_deconv_sce)$cell_type))
table(is.na(colData(sc_deconv_sce)$sample_id))
table(colData(sc_deconv_sce)$cell_type)
table(colData(sc_deconv_sce)$sample_id)
###############################
# Separate ER bulk samples
bulk_case_matrix = bulk_deconv_matrix[
  ,
  grepl("_ER[0-9]+$", colnames(bulk_deconv_matrix)),
  drop = FALSE]
# Separate Normal bulk samples
bulk_control_matrix = bulk_deconv_matrix[
  ,
  grepl("_NOR[0-9]+$", colnames(bulk_deconv_matrix)),
  drop = FALSE ]
dim(bulk_case_matrix)
colnames(bulk_case_matrix)
dim(bulk_control_matrix)
colnames(bulk_control_matrix)

# Check that bulk and scRNA-seq before Music
all(rownames(bulk_case_matrix) == rownames(sc_deconv_sce))
all(rownames(bulk_control_matrix) == rownames(sc_deconv_sce))
length(unique(colData(sc_deconv_sce)$cell_type))
unique(colData(sc_deconv_sce)$cell_type)
length(unique(colData(sc_deconv_sce)$sample_id))
unique(colData(sc_deconv_sce)$sample_id)
class(assay(sc_deconv_sce, "counts"))
class(bulk_case_matrix)
class(bulk_control_matrix)


######################
# Load the packages required for MuSiC2 deconvolution
library(SingleCellExperiment)
library(TOAST)
library(MuSiC)

# Define the cell types to estimate
deconv_selected_celltypes = unique(
  colData(sc_deconv_sce)$cell_type)

# Define all genes shared between the scRNA-seq reference and bulk RNA-seq
deconv_music2_markers = intersect(
  rownames(sc_deconv_sce),
  rownames(bulk_deconv_matrix))

# Check the number of genes used as MuSiC2 markers
length(deconv_music2_markers)

# Run MuSiC2 on the Normal bulk samples as a diagnostic test
music2_normal_result = music2_prop(
  bulk.control.mtx = bulk_control_matrix,
  bulk.case.mtx = bulk_control_matrix,
  sc.sce = sc_deconv_sce,
  clusters = "cell_type",
  samples = "sample_id",
  select.ct = deconv_selected_celltypes,
  markers = deconv_music2_markers,
  expr_low = 5,
  n_resample = 10,
  maxiter = 50
)

# Build the MuSiC signature basis from the scRNA-seq reference
deconv_music_basis = music_basis(
  sc.sce = sc_deconv_sce,
  non.zero = TRUE,
  markers = deconv_music2_markers,
  clusters = "cell_type",
  samples = "sample_id",
  select.ct = deconv_selected_celltypes
)

# Check the dimensions of the MuSiC signature matrix
dim(deconv_music_basis$Disgn.mtx)

# Check how many genes are shared between the MuSiC signature
# matrix and the Normal bulk RNA-seq matrix
deconv_common_genes = intersect(
  rownames(deconv_music_basis$Disgn.mtx),
  rownames(bulk_control_matrix)
)

# Print the number of common genes
length(deconv_common_genes)

# Show the arguments accepted by music_basis
args(music_basis)
# Build the MuSiC signature basis from the scRNA-seq reference
deconv_music_basis = music_basis(
  x = sc_deconv_sce,
  non.zero = TRUE,
  markers = deconv_music2_markers,
  clusters = "cell_type",
  samples = "sample_id",
  select.ct = deconv_selected_celltypes
)

# Check the dimensions of the MuSiC signature matrix
dim(deconv_music_basis$Disgn.mtx)

# Find genes shared between the MuSiC signature matrix and the Normal bulk matrix
deconv_common_genes = intersect(
  rownames(deconv_music_basis$Disgn.mtx),
  rownames(bulk_control_matrix)
)

# Count the shared genes
length(deconv_common_genes)
# Check the components returned by music_basis
names(deconv_music_basis)
# Check the dimensions of the main MuSiC basis components
dim(deconv_music_basis$Disgn.mtx)
dim(deconv_music_basis$Sigma)
dim(deconv_music_basis$M.S)



# Temporarily patch music.iter so that gene names are restored when they are missing
trace(what = "music.iter",
  tracer = quote({
    if (is.null(names(Y)) && length(Y) == nrow(D)) {
      names(Y) = rownames(D) }
    }),
  at = 1,
  print = FALSE)
# Run MuSiC deconvolution on all Normal bulk samples
deconv_musicprop_normal = music_prop(
  bulk.mtx = bulk_control_matrix,
  sc.sce = sc_deconv_sce,
  clusters = "cell_type",
  samples = "sample_id",
  select.ct = deconv_selected_celltypes,
  markers = deconv_music2_markers,
  verbose = TRUE)

# Run MuSiC deconvolution on all ER+ bulk samples
deconv_musicprop_case = music_prop(
  bulk.mtx = bulk_case_matrix,
  sc.sce = sc_deconv_sce,
  clusters = "cell_type",
  samples = "sample_id",
  select.ct = deconv_selected_celltypes,
  markers = deconv_music2_markers,
  verbose = TRUE)

# Check the MuSiC deconvolution results for Normal samples
dim(deconv_musicprop_normal$Est.prop.weighted)
head(deconv_musicprop_normal$Est.prop.weighted)

# Check the MuSiC deconvolution results for ER+ samples
dim(deconv_musicprop_case$Est.prop.weighted)
head(deconv_musicprop_case$Est.prop.weighted)

# Extract weighted cell-type proportions from Normal samples
deconv_normal_proportions = as.data.frame(
  deconv_musicprop_normal$Est.prop.weighted)

# Extract weighted cell-type proportions from ER+ samples
deconv_case_proportions = as.data.frame(
  deconv_musicprop_case$Est.prop.weighted)

# Add condition labels to the two groups
deconv_normal_proportions$condition = "Normal"
deconv_case_proportions$condition = "ER+"

# Add sample IDs as an explicit column
deconv_normal_proportions$sample_id = rownames(deconv_normal_proportions)
deconv_case_proportions$sample_id = rownames(deconv_case_proportions)

# Combine Normal and ER+ deconvolution results
deconv_final_proportions = rbind(
  deconv_normal_proportions,
  deconv_case_proportions)

# Move sample ID and condition to the first columns
deconv_final_proportions = deconv_final_proportions[
  ,
  c("sample_id", "condition",
    setdiff(
      colnames(deconv_final_proportions),
      c("sample_id", "condition"))) ]
deconv_final_proportions

# Extract weighted cell-type proportions from Normal samples
deconv_normal_proportions = as.data.frame(
  deconv_musicprop_normal$Est.prop.weighted)

# Extract weighted cell-type proportions from ER+ samples
deconv_case_proportions = as.data.frame(
  deconv_musicprop_case$Est.prop.weighted)

# Add condition labels to the two groups
deconv_normal_proportions$condition = "Normal"
deconv_case_proportions$condition = "ER+"

# Add sample IDs as an explicit column
deconv_normal_proportions$sample_id = rownames(deconv_normal_proportions)
deconv_case_proportions$sample_id = rownames(deconv_case_proportions)

# Combine Normal and ER+ deconvolution results
deconv_final_proportions = rbind(
  deconv_normal_proportions,
  deconv_case_proportions)

# Move sample ID and condition to the first columns
deconv_final_proportions = deconv_final_proportions[
  ,
  c("sample_id", "condition",
    setdiff(colnames(deconv_final_proportions),
      c("sample_id", "condition"))) ]
deconv_final_proportionssaveRDS(
  deconv_final_proportions,
  "MuSiC_deconvolution_celltype_proportions.rds")


library(ggplot2)

# Convert the proportion table from wide format to long format
deconv_long_proportions = reshape2::melt(
  deconv_final_proportions,
  id.vars = c("sample_id", "condition"),
  variable.name = "cell_type",
  value.name = "proportion")

# Create a stacked bar plot of cell-type proportions
deconv_stacked_barplot = ggplot(
  deconv_long_proportions,
  aes(x = sample_id,
    y = proportion,
    fill = cell_type)) +
  geom_bar(
    stat = "identity",
    position = "fill") +
  labs(title = "MuSiC Estimated Cell-Type Proportions",
    x = "Bulk RNA-seq Sample",
    y = "Cell-Type Proportion") +
  theme_classic() + theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1),
    legend.title = element_blank())
deconv_stacked_barplot
#################################

library(ggplot2)
library(ggpubr)

# Convert the deconvolution table to long format
deconv_long_proportions = reshape2::melt(
  deconv_final_proportions,
  id.vars = c("sample_id", "condition"),
  variable.name = "cell_type",
  value.name = "proportion")

# Make sure condition is treated as a two-group factor
deconv_long_proportions$condition = factor(
  deconv_long_proportions$condition,
  levels = c("Normal", "ER+"))

# Create boxplots with individual sample points and statistical comparisons
deconv_boxplot = ggplot(
  deconv_long_proportions,
  aes(x = condition,
    y = proportion,
    fill = condition)) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA) +
  geom_jitter(
    width = 0.12,
    size = 2) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("Normal", "ER+")),
    label = "p.format",
    size = 4.5,
    label.y.npc = 0.95) +
  facet_wrap(
    ~ cell_type,
    scales = "free_y",
    ncol = 3) +
  labs(
    title = "MuSiC Cell-Type Proportions: ER+ vs Normal",
    x = NULL,
    y = "Estimated Cell-Type Proportion") +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 11),
    plot.title = element_text(size = 14))
deconv_boxplot
##############################################
# Prepare the cell-type proportion matrix for the heatmap
deconv_heatmap_matrix = as.matrix(
  deconv_final_proportions[
    ,
    !(colnames(deconv_final_proportions) %in% c("sample_id", "condition")) ])

# Use sample IDs as row names
rownames(deconv_heatmap_matrix) = deconv_final_proportions$sample_id

# Convert the matrix to long format
deconv_heatmap_long = reshape2::melt(
  deconv_heatmap_matrix,
  varnames = c("sample_id", "cell_type"),
  value.name = "proportion")

# Add the experimental condition
deconv_heatmap_long$condition = deconv_final_proportions$condition[
  match(deconv_heatmap_long$sample_id,
    deconv_final_proportions$sample_id) ]

# Create the cell-type proportion heatmap
deconv_heatmap = ggplot(
  deconv_heatmap_long,
  aes(x = cell_type,
    y = sample_id,
    fill = proportion)) +
  geom_tile() +
  labs(title = "MuSiC Cell-Type Proportions",
    x = "Cell Type",
    y = "Bulk RNA-seq Sample",
    fill = "Proportion") +
  theme_classic() +
  theme(axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 9),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(size = 14))


deconv_heatmap
