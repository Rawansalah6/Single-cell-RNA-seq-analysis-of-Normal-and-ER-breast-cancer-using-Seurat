#Preparation for pseusoBulk
#cteate condition labels

pbulk_condition_vector = ifelse(
  seurat_integrated$orig.ident %in% c("con1", "con2", "con3", "con4"),"Normal", "ER+")

identical(rownames(meta),
  colnames(seurat_integrated))
#extract raw RNA counts
pbulk_raw_count_matrix = LayerData(
  seurat_integrated,
  assay = "RNA",
  layer = "counts")
#create a compact metadata table for pseudobulk
pbulk_annotation_table = data.frame(
  sample_id = meta$orig.ident,
  condition = pbulk_condition_vector,
  cell_type = meta$cell_type,
  row.names = rownames(meta))

# Create a unique pseudobulk group for each sample and cell type
pbulk_sample_celltype_id = paste(
  pbulk_annotation_table$sample_id,
  pbulk_annotation_table$cell_type,
  sep = "__")

# Check cell order between metadata and raw count matrix
identical(rownames(pbulk_annotation_table),
  colnames(pbulk_raw_count_matrix))


library(Matrix)
# Create a sparse model matrix for sample-cell type groups
pbulk_group_matrix = sparse.model.matrix(
  ~ 0 + pbulk_sample_celltype_id)

# Aggregate raw counts while preserving sparse matrix format
pbulk_count_matrix = pbulk_raw_count_matrix %*% pbulk_group_matrix


dim(pbulk_count_matrix)
head(colnames(pbulk_count_matrix))
class(pbulk_count_matrix)

# Clean pseudobulk sample-cell type column names
colnames(pbulk_count_matrix) = sub(
  "^pbulk_sample_celltype_id",
  "",
  colnames(pbulk_count_matrix))
head(colnames(pbulk_count_matrix))


# Create metadata for the 112 pseudobulk profiles
pbulk_sample_metadata = unique(
  pbulk_annotation_table[, c("sample_id", "condition", "cell_type")])

# Create the same pseudobulk ID used in the count matrix
pbulk_sample_metadata$pbulk_id = paste(
  pbulk_sample_metadata$sample_id,
  pbulk_sample_metadata$cell_type,
  sep = "__")

# Use the pseudobulk ID as row names
rownames(pbulk_sample_metadata) = pbulk_sample_metadata$pbulk_id

# Reorder metadata to exactly match the count matrix columns
pbulk_sample_metadata = pbulk_sample_metadata[
  colnames(pbulk_count_matrix),
  ,
  drop = FALSE ]

identical(colnames(pbulk_count_matrix),
  rownames(pbulk_sample_metadata))

table(pbulk_sample_metadata$cell_type)

table(pbulk_sample_metadata$cell_type,
  pbulk_sample_metadata$condition)

# Get the 14 cell types available in the pseudobulk dataset
pbulk_celltype_list = unique(
  pbulk_sample_metadata$cell_type)
pbulk_celltype_list

# Create a list to store filtered pseudobulk count matrices
pbulk_filtered_counts_list = list()

# Filter low-count genes separately for each cell type
for (pbulk_current_celltype in pbulk_celltype_list) {
  pbulk_current_columns = rownames(
    pbulk_sample_metadata[
      pbulk_sample_metadata$cell_type == pbulk_current_celltype,
      ,
      drop = FALSE ] )
  pbulk_current_counts = pbulk_count_matrix[
    ,
    pbulk_current_columns,
    drop = FALSE ]
  pbulk_keep_genes = rowSums(
    pbulk_current_counts >= 10
    ) >= 4
  
  pbulk_filtered_counts_list[[pbulk_current_celltype]] <-
    pbulk_current_counts[pbulk_keep_genes, , drop = FALSE] }

# Check the number of genes remaining after filtering
pbulk_filtered_gene_numbers = sapply(
  pbulk_filtered_counts_list, nrow)
pbulk_filtered_gene_numbers

library(DESeq2)

# Select Fibroblast pseudobulk profiles
pbulk_fibroblast_ids = rownames(pbulk_sample_metadata[
    pbulk_sample_metadata$cell_type == "Fibroblasts",
    ,
    drop = FALSE ])

# Extract Fibroblast pseudobulk counts
pbulk_fibroblast_counts = pbulk_filtered_counts_list[["Fibroblasts"]][
  ,
  pbulk_fibroblast_ids,
  drop = FALSE]

# Extract matching metadata
pbulk_fibroblast_metadata = pbulk_sample_metadata[pbulk_fibroblast_ids,
  ,
  drop = FALSE]

# Make sure the condition is treated as a factor
pbulk_fibroblast_metadata$condition = factor(
  pbulk_fibroblast_metadata$condition,
  levels = c("Normal", "ER+"))
table(pbulk_fibroblast_metadata$condition)

# Create the DESeq2 dataset for Fibroblasts
pbulk_fibroblast_dds = DESeqDataSetFromMatrix(
  countData = pbulk_fibroblast_counts,
  colData = pbulk_fibroblast_metadata,
  design = ~ condition)

# Run differential expression analysis
pbulk_fibroblast_dds = DESeq(
  pbulk_fibroblast_dds)
# Extract differential expression results for ER+ versus Normal
pbulk_fibroblast_de_results = results(
  pbulk_fibroblast_dds,
  contrast = c("condition", "ER+", "Normal"))


pbulk_fibroblast_de_results = pbulk_fibroblast_de_results[
  order(pbulk_fibroblast_de_results$padj),]

head(as.data.frame(pbulk_fibroblast_de_results))
################
#for the 14 cell type
# Run DESeq2 separately for all 14 cell types
pbulk_deseq2_results_list = list()

for (pbulk_current_celltype in pbulk_celltype_list) {
  pbulk_current_ids = rownames(
    pbulk_sample_metadata[
      pbulk_sample_metadata$cell_type == pbulk_current_celltype,
      ,
      drop = FALSE ])
   pbulk_current_counts = pbulk_filtered_counts_list[[pbulk_current_celltype]][
    ,
    pbulk_current_ids,
    drop = FALSE ]
  pbulk_current_metadata = pbulk_sample_metadata[
    pbulk_current_ids,
    ,
    drop = FALSE ]
  pbulk_current_metadata$condition = factor(
    pbulk_current_metadata$condition,
    levels = c("Normal", "ER+"))
  pbulk_current_dds = DESeqDataSetFromMatrix(
    countData = pbulk_current_counts,
    colData = pbulk_current_metadata,
    design = ~ condition)
  pbulk_current_dds = DESeq(
    pbulk_current_dds,
    quiet = TRUE)
  
# Extract ER+ versus Normal results
  pbulk_current_results = results(
    pbulk_current_dds,
    contrast = c("condition", "ER+", "Normal"))
   pbulk_current_results = pbulk_current_results[
    order(pbulk_current_results$padj), ]
  pbulk_deseq2_results_list[[pbulk_current_celltype]] =
    pbulk_current_results }

# Create a summary table of significant DEGs for each cell type
pbulk_deg_summary = data.frame(
  cell_type = character(),
  significant_DEGs = integer(),
  up_in_ERplus = integer(),
  up_in_Normal = integer(),
  stringsAsFactors = FALSE)

for (pbulk_current_celltype in pbulk_celltype_list) {
  pbulk_current_deg_table = as.data.frame(
    pbulk_deseq2_results_list[[pbulk_current_celltype]])
  pbulk_current_deg_table = pbulk_current_deg_table[
    !is.na(pbulk_current_deg_table$padj),
    ,
    drop = FALSE ]
  pbulk_significant_genes = pbulk_current_deg_table[
    pbulk_current_deg_table$padj < 0.05 &
      abs(pbulk_current_deg_table$log2FoldChange) >= 1,
    ,
    drop = FALSE ]
  pbulk_up_erplus = sum(
    pbulk_significant_genes$log2FoldChange >= 1)
  
  pbulk_up_normal = sum(
    pbulk_significant_genes$log2FoldChange <= -1)
  
  pbulk_deg_summary = rbind(
    pbulk_deg_summary, data.frame(
      cell_type = pbulk_current_celltype,
      significant_DEGs = nrow(pbulk_significant_genes),
      up_in_ERplus = pbulk_up_erplus,
      up_in_Normal = pbulk_up_normal)) }

pbulk_deg_summary

# Extract the top 20 significant DEGs for every cell type
pbulk_top20_deg_list = list()

for (pbulk_current_celltype in pbulk_celltype_list) {
  pbulk_current_deg_table = as.data.frame(
    pbulk_deseq2_results_list[[pbulk_current_celltype]])
  pbulk_current_deg_table$gene = rownames(
    pbulk_current_deg_table)
  pbulk_current_deg_table = pbulk_current_deg_table[
    !is.na(pbulk_current_deg_table$padj) &
      pbulk_current_deg_table$padj < 0.05 &
      abs(pbulk_current_deg_table$log2FoldChange) >= 1,
    ,
    drop = FALSE ]
  pbulk_top20_deg_list[[pbulk_current_celltype]] =
    head(pbulk_current_deg_table[
        order(pbulk_current_deg_table$padj),
        c("gene", "log2FoldChange", "pvalue", "padj")], 20) }
# View the top 20 DEGs in Fibroblasts
pbulk_top20_deg_list[["Fibroblasts"]]
# View the top 20 DEGs in Endothelial cells
pbulk_top20_deg_list[["Endothelial"]]

# Check normalized expression of selected Fibroblast DEGs
pbulk_fibroblast_normalized_counts = counts(
  pbulk_fibroblast_dds,
  normalized = TRUE)

# Inspect selected genes across the 8 samples
pbulk_fibroblast_normalized_counts[
  c("TNFAIP6", "CXCL3", "COL8A1", "CXCL2", "SFRP2", "CTHRC1","MMP3",
    "CXCL1"), ]

# Perform variance stabilizing transformation on Fibroblast pseudobulk samples
pbulk_fibroblast_vst = vst(
  pbulk_fibroblast_dds,
  blind = FALSE)

# Plot PCA showing Normal versus ER+ samples
plotPCA(pbulk_fibroblast_vst,
  intgroup = "condition")

# Save DESeq2 results for all 14 cell types
saveRDS(pbulk_deseq2_results_list,
  "pseudobulk_DESeq2_results_all_celltypes.rds")

# Save the filtered pseudobulk counts
saveRDS(pbulk_filtered_counts_list,
  "pseudobulk_filtered_counts_all_celltypes.rds")

# Save the pseudobulk sample metadata
saveRDS(pbulk_sample_metadata,
  "pseudobulk_sample_metadata.rds")

# Create one combined table containing DEGs from all cell types
pbulk_all_deg_table = do.call(rbind,
  lapply(names(pbulk_deseq2_results_list),
    function(pbulk_current_celltype) {
      
      pbulk_current_table = as.data.frame(
        pbulk_deseq2_results_list[[pbulk_current_celltype]])
      pbulk_current_table$gene = rownames(
        pbulk_current_table)
      pbulk_current_table$cell_type = pbulk_current_celltype
      pbulk_current_table }))

# Keep only significant DEGs
pbulk_all_significant_deg_table = pbulk_all_deg_table[
  !is.na(pbulk_all_deg_table$padj) &
    pbulk_all_deg_table$padj < 0.05 &
    abs(pbulk_all_deg_table$log2FoldChange) >= 1,
  ,
  drop = FALSE]
dim(pbulk_all_significant_deg_table)


library(ggplot2)
# Create a folder for the 14 volcano plots
dir.create("pseudobulk_volcano_plots", showWarnings = FALSE)

# Generate and save a volcano plot for every cell type
for (pbulk_current_celltype in pbulk_celltype_list) {
  pbulk_current_volcano_table = as.data.frame(
    pbulk_deseq2_results_list[[pbulk_current_celltype]])
  pbulk_current_volcano_table$gene = rownames(
    pbulk_current_volcano_table)
  pbulk_current_volcano_table = pbulk_current_volcano_table[
    !is.na(pbulk_current_volcano_table$padj) &
      !is.na(pbulk_current_volcano_table$log2FoldChange),
    ,
    drop = FALSE]
  
  # Classify genes by significance and direction
  pbulk_current_volcano_table$DEG_status = "Not significant"
  
  pbulk_current_volcano_table$DEG_status[
    pbulk_current_volcano_table$padj < 0.05 &
      pbulk_current_volcano_table$log2FoldChange >= 1
  ] <- "Up in ER+"
  
  pbulk_current_volcano_table$DEG_status[
    pbulk_current_volcano_table$padj < 0.05 &
      pbulk_current_volcano_table$log2FoldChange <= -1
  ] <- "Up in Normal"
  
  # Create the volcano plot
  pbulk_current_volcano_plot = ggplot(
    pbulk_current_volcano_table,
    aes(x = log2FoldChange,
      y = -log10(padj), color = DEG_status)) +
    geom_point(alpha = 0.65, size = 1.6) +
    geom_vline(xintercept = c(-1, 1),
      linetype = "dashed", linewidth = 0.5) +
    geom_hline(yintercept = -log10(0.05),
      linetype = "dashed", linewidth = 0.5) +
    scale_color_manual(values = c(
        "Up in ER+" = "#D73027",
        "Up in Normal" = "#4575B4",
        "Not significant" = "grey75")) +
    labs(title = pbulk_current_celltype,
      subtitle = "ER+ vs Normal",
      x = "log2 Fold Change",
      y = "-log10 Adjusted P-value",
      color = NULL) +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(
        face = "bold",size = 15),
      plot.subtitle = element_text(
        size = 11),
      legend.position = "top",
      axis.title = element_text(
        face = "bold"))
  
  # Save the plot as a high-resolution PNG
  ggsave(filename = paste0(
      "pseudobulk_volcano_plots/",
      gsub("[^A-Za-z0-9]+",
        "_",  pbulk_current_celltype),
      "_volcano.png"),
    plot = pbulk_current_volcano_plot,
    width = 8, height = 6, dpi = 300) }

# Load packages for heatmap visualization
library(pheatmap)
library(grid)
library(patchwork)

# Create an empty list for the 14 heatmaps
pbulk_heatmap_plot_list = list()

# Generate a heatmap for each cell type
for (pbulk_current_celltype in pbulk_celltype_list) {
  pbulk_current_heatmap_genes =
    pbulk_top20_deg_list[[pbulk_current_celltype]]$gene
  pbulk_current_ids = rownames(
    pbulk_sample_metadata[
      pbulk_sample_metadata$cell_type == pbulk_current_celltype,
      ,
      drop = FALSE])
  pbulk_current_counts =
    pbulk_filtered_counts_list[[pbulk_current_celltype]][
      ,
      pbulk_current_ids,
      drop = FALSE]
   pbulk_current_metadata = pbulk_sample_metadata[
    pbulk_current_ids,
    ,
    drop = FALSE]
  pbulk_current_metadata$condition = factor(
    pbulk_current_metadata$condition,
    levels = c("Normal", "ER+"))
  pbulk_current_dds = DESeqDataSetFromMatrix(
    countData = pbulk_current_counts,
    colData = pbulk_current_metadata,
    design = ~ condition)
  pbulk_current_vst = varianceStabilizingTransformation(
    pbulk_current_dds,
    blind = FALSE)
  pbulk_current_heatmap_matrix = assay(
    pbulk_current_vst)[
    pbulk_current_heatmap_genes,
    ,
    drop = FALSE]
   pbulk_current_annotation = data.frame(
    Condition = pbulk_current_metadata$condition)
  
  rownames(pbulk_current_annotation) =
    rownames(pbulk_current_metadata)
  
 pbulk_current_heatmap_grob = grid.grabExpr(
    pheatmap(
      pbulk_current_heatmap_matrix,
      scale = "row",
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      annotation_col = pbulk_current_annotation,
      show_colnames = FALSE, show_rownames = TRUE,
      fontsize = 6, fontsize_row = 5,
      main = pbulk_current_celltype))
  pbulk_heatmap_plot_list[[pbulk_current_celltype]] =
    pbulk_current_heatmap_grob }
wrap_plots(lapply(
    pbulk_heatmap_plot_list,
    function(pbulk_current_grob) {
      wrap_elements(pbulk_current_grob) }), ncol = 4)
############################################
#GSEA
# Create ranked gene lists for GSEA from existing DESeq2 results
pbulk_gsea_ranked_genes_list = lapply(pbulk_deseq2_results_list,
  function(pbulk_current_results) {
    pbulk_current_results = pbulk_current_results[
      !is.na(pbulk_current_results$stat),]
    pbulk_current_ranked_genes = pbulk_current_results$stat
    names(pbulk_current_ranked_genes) =
      rownames(pbulk_current_results)
    sort(pbulk_current_ranked_genes, decreasing = TRUE) })

library(clusterProfiler)
install.packages("msigdbr")
library(msigdbr)
# Retrieve human Hallmark gene sets from MSigDB
pbulk_hallmark_sets = msigdbr(
  species = "Homo sapiens",
  collection = "H")

# Create the pathway-to-gene mapping required by GSEA
pbulk_hallmark_gene_list = split(
  pbulk_hallmark_sets$gene_symbol,
  pbulk_hallmark_sets$gs_name)

pbulk_gsea_results_list = list()

for (pbulk_current_celltype in pbulk_celltype_list) {
  pbulk_current_ranked_genes =
    pbulk_gsea_ranked_genes_list[[pbulk_current_celltype]]
  pbulk_current_gsea = GSEA(
    geneList = pbulk_current_ranked_genes,
    TERM2GENE = pbulk_hallmark_sets[, c("gs_name", "gene_symbol")],
    pvalueCutoff = 0.05,
    verbose = FALSE)
  pbulk_gsea_results_list[[pbulk_current_celltype]] =
    pbulk_current_gsea }


# Plot GSEA results for all 14 cell types
# Check the exact pathway names currently used in the GSEA plot
unique(pbulk_gsea_plot_table$Description)
# Calculate the best adjusted p-value for each pathway across all cell types
pbulk_top10_pathways = pbulk_gsea_plot_table %>%
  group_by(Description) %>%
  summarise(
    best_padj = min(p.adjust, na.rm = TRUE)
  ) %>%
  arrange(best_padj) %>%
  slice_head(n = 10) %>%
  pull(Description)

pbulk_gsea_top10_table = pbulk_gsea_plot_table %>%
  filter(Description %in% pbulk_top10_pathways)

pbulk_gsea_top10_table$Pathway_short = dplyr::recode(
  as.character(pbulk_gsea_top10_table$Description),
  
  "OXIDATIVE PHOSPHORYLATION" = "OxPhos",
  "TNFA SIGNALING VIA NFKB" = "TNFα–NFκB",
  "INFLAMMATORY RESPONSE" = "Inflammation",
  "EPITHELIAL MESENCHYMAL TRANSITION" = "EMT",
  "HYPOXIA" = "Hypoxia",
  "GLYCOLYSIS" = "Glycolysis",
  "MYC TARGETS V1" = "MYC V1",
  "MYC TARGETS V2" = "MYC V2",
  "MTORC1 SIGNALING" = "mTORC1",
  "UNFOLDED PROTEIN RESPONSE" = "UPR",
  "KRAS SIGNALING UP" = "KRAS ↑",
  "ALLOGRAFT REJECTION" = "Allograft",
  "KRAS SIGNALING DN" = "KRAS ↓",
  "ESTROGEN RESPONSE EARLY" = "Estrogen",
  "IL6 JAK STAT3 SIGNALING" = "IL6–JAK–STAT3",
  "IL2 STAT5 SIGNALING" = "IL2–STAT5",
  "INTERFERON ALPHA RESPONSE" = "IFN-α")

pbulk_gsea_top10_dotplot = ggplot(
  pbulk_gsea_top10_table,
  aes(x = cell_type,
    y = Pathway_short, size = setSize,
    color = NES)) +
  geom_point() + coord_flip() + theme_bw() + theme(
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    plot.title = element_text(size = 13)) +
  labs(title = "Top 10 Hallmark pathways across cell types",
    x = "Cell type", y = "Pathway",
    size = "Gene set size", color = "NES")
pbulk_gsea_top10_dotplot

# Prepare NES matrix for the top 10 GSEA pathways
# Replace missing NES values with NA-safe values for plotting
pbulk_gsea_nes_matrix_clean = pbulk_gsea_nes_matrix

# Remove pathways or cell types containing only missing values
pbulk_gsea_nes_matrix_clean <- pbulk_gsea_nes_matrix_clean[
  rowSums(!is.na(pbulk_gsea_nes_matrix_clean)) > 0,
  ,
  drop = FALSE ]
pbulk_gsea_nes_matrix_clean = pbulk_gsea_nes_matrix_clean[
  ,
  colSums(!is.na(pbulk_gsea_nes_matrix_clean)) > 0,
  drop = FALSE]

# Replace remaining missing values with zero for clustering/visualization
pbulk_gsea_nes_matrix_clean[
  is.na(pbulk_gsea_nes_matrix_clean) ] <- 0

# Plot the cleaned GSEA NES heatmap
pheatmap(pbulk_gsea_nes_matrix_clean,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  number_format = "%.2f",
  fontsize = 8,
  fontsize_row = 9,
  fontsize_col = 7,
  main = "Top 10 Hallmark pathways – GSEA NES")
##########################################
#ORA
# Prepare Hallmark gene sets for ORA
pbulk_ora_hallmark = msigdbr(
  species = "Homo sapiens",
  collection = "H")

pbulk_ora_term2gene = pbulk_ora_hallmark[
  ,
  c("gs_name", "gene_symbol") ]

pbulk_ora_directional_results = list()
pbulk_ora_celltypes = unique(
  pbulk_all_significant_deg_table$cell_type)

# Run ORA separately for ER+ upregulated and Normal upregulated genes
for (pbulk_current_celltype in pbulk_ora_celltypes) {
  pbulk_current_deg_table = pbulk_all_significant_deg_table[
    pbulk_all_significant_deg_table$cell_type ==
      pbulk_current_celltype,
    ,
    drop = FALSE ]
  pbulk_erplus_up_genes = unique(
    pbulk_current_deg_table$gene[
      pbulk_current_deg_table$log2FoldChange >= 1 ])
  pbulk_normal_up_genes = unique(
    pbulk_current_deg_table$gene[
      pbulk_current_deg_table$log2FoldChange <= -1 ])
  pbulk_current_universe = rownames(
    pbulk_filtered_counts_list[[pbulk_current_celltype]])
  
  # Run ORA for ER+ upregulated genes
  pbulk_erplus_ora = enricher(
    gene = pbulk_erplus_up_genes,
    universe = pbulk_current_universe,
    TERM2GENE = pbulk_ora_term2gene,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05)
  
  # Run ORA for Normal upregulated genes
  pbulk_normal_ora = enricher(
    gene = pbulk_normal_up_genes,
    universe = pbulk_current_universe,
    TERM2GENE = pbulk_ora_term2gene,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05)
  
  # Store both directional ORA results
  pbulk_ora_directional_results[[pbulk_current_celltype]] = list(
    ERplus = pbulk_erplus_ora,
    Normal = pbulk_normal_ora) }
length(pbulk_ora_directional_results)
names(pbulk_ora_directional_results)

pbulk_ora_table_list = list()

# Combine ER+ and Normal ORA results for each cell type
for (pbulk_current_celltype in names(pbulk_ora_directional_results)) {
  
  
  pbulk_erplus_table = as.data.frame(
    pbulk_ora_directional_results[[pbulk_current_celltype]]$ERplus)
  if (nrow(pbulk_erplus_table) > 0) {
    pbulk_erplus_table$direction <- "ER+ up"
    pbulk_erplus_table$cell_type <- pbulk_current_celltype
    pbulk_ora_table_list[[length(pbulk_ora_table_list) + 1]] =
      pbulk_erplus_table }
  pbulk_normal_table = as.data.frame(
    pbulk_ora_directional_results[[pbulk_current_celltype]]$Normal)
  if (nrow(pbulk_normal_table) > 0) {
    pbulk_normal_table$direction = "Normal up"
    pbulk_normal_table$cell_type = pbulk_current_celltype
    pbulk_ora_table_list[[length(pbulk_ora_table_list) + 1]] =
      pbulk_normal_table } }

# Combine all non-empty ORA results
pbulk_ora_combined_table = do.call(rbind,
  pbulk_ora_table_list)
rownames(pbulk_ora_combined_table) = NULL
dim(pbulk_ora_combined_table)
head(pbulk_ora_combined_table)

# Print all ORA pathway names without truncation
print(unique(pbulk_ora_combined_table$Description), max = 100)

# Create short biological names for ORA pathways
pbulk_ora_combined_table$Pathway_short = pbulk_ora_combined_table$Description

pbulk_ora_combined_table$Pathway_short = gsub("HALLMARK_", "",
  pbulk_ora_combined_table$Pathway_short)

pbulk_ora_combined_table$Pathway_short = gsub(
  "OXIDATIVE_PHOSPHORYLATION", "OxPhos",
  pbulk_ora_combined_table$Pathway_short)
pbulk_ora_combined_table$Pathway_short = gsub(
  "FATTY_ACID_METABOLISM", "FA Metabolism",
  pbulk_ora_combined_table$Pathway_short)

pbulk_ora_combined_table$Pathway_short = gsub(
  "TNFA_SIGNALING_VIA_NFKB", "TNFα–NFκB",
  pbulk_ora_combined_table$Pathway_short)

pbulk_ora_combined_table$Pathway_short = gsub(
  "INFLAMMATORY_RESPONSE", "Inflammation",
  pbulk_ora_combined_table$Pathway_short)

pbulk_ora_combined_table$Pathway_short = gsub(
  "EPITHELIAL_MESENCHYMAL_TRANSITION", "EMT",
  pbulk_ora_combined_table$Pathway_short)

pbulk_ora_combined_table$Pathway_short = gsub(
  "KRAS_SIGNALING_UP", "KRAS ↑",
  pbulk_ora_combined_table$Pathway_short)

pbulk_ora_combined_table$Pathway_short = gsub(
  "KRAS_SIGNALING_DN", "KRAS ↓",
  pbulk_ora_combined_table$Pathway_short)

pbulk_ora_combined_table$Pathway_short = gsub(
  "IL2_STAT5_SIGNALING", "IL2–STAT5",
  pbulk_ora_combined_table$Pathway_short)

pbulk_ora_combined_table$Pathway_short = gsub(
  "IL6_JAK_STAT3_SIGNALING", "IL6–JAK–STAT3",
  pbulk_ora_combined_table$Pathway_short)

pbulk_ora_combined_table$Pathway_short = gsub(
  "INTERFERON_ALPHA_RESPONSE", "IFN-α",
  pbulk_ora_combined_table$Pathway_short)

pbulk_ora_combined_table$Pathway_short = gsub(
  "INTERFERON_GAMMA_RESPONSE", "IFN-γ",
  pbulk_ora_combined_table$Pathway_short)

pbulk_ora_combined_table$Pathway_short = gsub(
  "MYC_TARGETS_V1", "MYC V1",
  pbulk_ora_combined_table$Pathway_short)

pbulk_ora_combined_table$Pathway_short = gsub(
  "MYC_TARGETS_V2", "MYC V2",
  pbulk_ora_combined_table$Pathway_short)

pbulk_ora_combined_table$Pathway_short = gsub(
  "MTORC1_SIGNALING", "mTORC1",
  pbulk_ora_combined_table$Pathway_short)

pbulk_ora_combined_table$Pathway_short = gsub(
  "UNFOLDED_PROTEIN_RESPONSE", "UPR",
  pbulk_ora_combined_table$Pathway_short)

pbulk_ora_combined_table$Pathway_short = gsub(
  "REACTIVE_OXYGEN_SPECIES_PATHWAY", "ROS",
  pbulk_ora_combined_table$Pathway_short)

pbulk_ora_combined_table$Pathway_short = gsub(
  "ESTROGEN_RESPONSE_EARLY", "Estrogen Early",
  pbulk_ora_combined_table$Pathway_short)

pbulk_ora_combined_table$Pathway_short = gsub(
  "ESTROGEN_RESPONSE_LATE", "Estrogen Late",
  pbulk_ora_combined_table$Pathway_short)

pbulk_ora_combined_table$Pathway_short = gsub(
  "G2M_CHECKPOINT", "G2/M",
  pbulk_ora_combined_table$Pathway_short)

pbulk_ora_combined_table$Pathway_short = gsub(
  "P53_PATHWAY", "p53",
  pbulk_ora_combined_table$Pathway_short)

pbulk_ora_combined_table$Pathway_short = gsub(
  "CHOLESTEROL_HOMEOSTASIS", "Cholesterol Homeostasis",
  pbulk_ora_combined_table$Pathway_short)

pbulk_ora_combined_table$Pathway_short = gsub(
  "_", " ",
  pbulk_ora_combined_table$Pathway_short)

# Select the 10 most significant ER+ ORA pathways overall
pbulk_ora_erplus_heatmap_table = pbulk_ora_combined_table %>%
  filter(direction == "ER+ up") %>%
  group_by(Pathway_short) %>%
  summarise(
    best_padj = min(p.adjust),
    .groups = "drop"
  ) %>%
  arrange(best_padj) %>%
  slice_head(n = 10) %>%
  pull(Pathway_short)

# Keep only the selected pathways
pbulk_ora_erplus_heatmap_data = pbulk_ora_combined_table %>%
  filter(
    direction == "ER+ up",
    Pathway_short %in% pbulk_ora_erplus_heatmap_table
  ) %>%
  mutate(
    score = -log10(p.adjust)
  ) %>%
  select(cell_type, Pathway_short, score) %>%
  distinct()

# Convert to wide format
pbulk_ora_erplus_heatmap_matrix = pbulk_ora_erplus_heatmap_data %>%
  tidyr::pivot_wider(
    names_from = cell_type,
    values_from = score,
    values_fill = 0
  ) %>%
  as.data.frame()

# Set unique pathway names as row names
rownames(pbulk_ora_erplus_heatmap_matrix) =
  pbulk_ora_erplus_heatmap_matrix$Pathway_short

pbulk_ora_erplus_heatmap_matrix$Pathway_short = NULL

# Plot the ER+ ORA heatmap
pheatmap(
  as.matrix(pbulk_ora_erplus_heatmap_matrix),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA,
  fontsize = 8,
  fontsize_row = 7,
  fontsize_col = 8,
  main = "ORA – ER+ Upregulated Genes")

###############
# Select the 10 most significant Normal ORA pathways overall
pbulk_ora_normal_heatmap_table = pbulk_ora_combined_table %>%
  filter(direction == "Normal up") %>%
  group_by(Pathway_short) %>%
  summarise(
    best_padj = min(p.adjust),
    .groups = "drop"
  ) %>%
  arrange(best_padj) %>%
  slice_head(n = 10) %>%
  pull(Pathway_short)

# Prepare the Normal ORA heatmap data
pbulk_ora_normal_heatmap_data = pbulk_ora_combined_table %>%
  filter(
    direction == "Normal up",
    Pathway_short %in% pbulk_ora_normal_heatmap_table
  ) %>%
  mutate(
    score = -log10(p.adjust)
  ) %>%
  select(cell_type, Pathway_short, score) %>%
  distinct()

# Convert to matrix format
pbulk_ora_normal_heatmap_matrix = pbulk_ora_normal_heatmap_data %>%
  tidyr::pivot_wider(
    names_from = cell_type,
    values_from = score,
    values_fill = 0
  ) %>%
  as.data.frame()

# Set pathway names as row names
rownames(pbulk_ora_normal_heatmap_matrix) =
  pbulk_ora_normal_heatmap_matrix$Pathway_short

pbulk_ora_normal_heatmap_matrix$Pathway_short = NULL

# Plot the Normal ORA heatmap
pheatmap(
  as.matrix(pbulk_ora_normal_heatmap_matrix),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA,
  fontsize = 8,
  fontsize_row = 7,
  fontsize_col = 8,
  main = "ORA – Normal Upregulated Genes")
