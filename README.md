# Single-cell-RNA-seq-analysis-of-Normal-and-ER-breast-cancer-using-Seurat
##Single-Cell and Bulk RNA-seq Integration in ER+ Breast Cancer

**Overview**

This project integrates single-cell RNA-seq (scRNA-seq) and bulk RNA-seq data to characterize the cellular and molecular landscape of ER-positive (ER+) breast cancer compared with normal breast tissue.

The scRNA-seq data provide cell-level resolution for identifying cellular populations and their molecular states, while bulk RNA-seq deconvolution enables estimation of cell-type composition across an independent bulk cohort.

**The analysis combines:**

- scRNA-seq preprocessing, QC, normalization, integration, and clustering
- Manual and automated cell-type annotation
- Cell-cell communication analysis
- Cell-type-specific pseudobulk differential expression
- GSEA and ORA pathway analysis
- Bulk RNA-seq cell-type deconvolution using the scRNA-seq reference

---

**Analysis Workflow**

The project follows a stepwise workflow from raw expression matrices to integrated biological interpretation.

scRNA-seq

Raw count matrices
↓
Quality Control & Filtering
↓
Normalization & HVG Selection
↓
Anchor-based Integration
↓
PCA & Clustering
↓
UMAP Visualization
↓
Cell-Type Annotation
↓
CellChat + Pseudobulk Analysis

Bulk RNA-seq

Raw sequencing reads
↓
Salmon Quantification
↓
Gene-level Count Matrix
↓
Gene Annotation & Filtering
↓
MuSiC Deconvolution
↓
Cell-Type Composition Comparison

Integrated Interpretation

Cellular composition + cell-cell communication + cell-type-specific expression + pathway activity

→ Characterization of ER+ breast cancer microenvironment

---

**Data**

Single-Cell RNA-seq

Source: GEO — "GSE161529" 

Eight samples were selected, including four ER+ tumors and four normal breast tissue samples.

Sample| GSM ID| Condition
ER1| GSM4909296| ER+
ER2| GSM4909297| ER+
ER3| GSM4909298| ER+
ER4| GSM4909299| ER+
con1| GSM4909253| Normal
con2| GSM4909254| Normal
con3| GSM4909257| Normal
con4| GSM4909261| Normal

Raw "matrix.mtx.gz", "barcodes.tsv.gz", and "features.tsv.gz" files were obtained from GEO and loaded into R using Seurat.

---

Bulk RNA-seq

The bulk RNA-seq cohort contains:

- 6 ER+ samples
- 5 Normal samples

Reads were quantified using Salmon through the Galaxy platform, generating a gene-level expression matrix.

Ensembl gene IDs were mapped to gene symbols using GENCODE v50. Protein-coding genes were retained, and duplicated gene symbols were collapsed before downstream analysis.

Bulk RNA-seq accession: To be added.

---

**Methods**

1. Quality Control & Preprocessing

Script: "01_preprocessing_QC.R"

Each scRNA-seq sample was processed independently in Seurat.

The following QC metrics were evaluated:

- "nFeature_RNA"
- "nCount_RNA"
- "percent.mt"

QC distributions were examined using violin and scatter plots before filtering.

Because the samples differed in sequencing depth and mitochondrial background, sample-specific QC thresholds were applied instead of using a single global cutoff.

Filtered datasets were re-evaluated to confirm the effect of QC, and cell numbers before and after filtering were recorded.

---

2. Normalization & Integration

Script: "02_integration_clustering.R"

Each sample was independently normalized using LogNormalize with a scale factor of 10,000.

The top 2,000 highly variable genes (HVGs) were identified using the "vst" method.

The eight samples were integrated using Seurat anchor-based CCA integration:

SelectIntegrationFeatures()
        ↓
FindIntegrationAnchors()
        ↓
IntegrateData()

The integrated dataset was subsequently scaled and used for dimensionality reduction and clustering.

---

3. Dimensionality Reduction & Clustering

PCA was performed on the integrated dataset, and an elbow plot was used to evaluate the informative principal components.

Graph-based clustering was performed using:

- 15 PCs
- Resolution = 0.5

The analysis identified 14 clusters (0–13).
UMAP was used to visualize the resulting cellular populations. UMAPs were examined both by cluster identity and sample of origin to assess the overall structure of the integrated dataset.

Cluster marker genes were identified using "FindAllMarkers()" and visualized using heatmaps and dot plots.

---

4. Cell-Type Annotation

Script: "03_annotation.R"

Cell identities were assigned using complementary manual and automated annotation approaches.

Manual Annotation

Clusters were annotated using canonical marker genes for major breast tissue populations, including:

- Luminal epithelial cells
- Luminal secretory cells
- Basal / myoepithelial cells
- Fibroblasts
- Endothelial cells
- Pericyte / perivascular cells
- T cells
- NK cells
- Myeloid / macrophage cells
- B / plasma cells
- Cycling cells

Marker expression was evaluated using cluster markers, DotPlots, FeaturePlots, and heatmaps.

Automated Annotation

SingleR was used as an independent annotation approach.

SingleR predictions were compared with marker-based annotations to support the final cell-type assignments.

---

Downstream Analyses

5. Cell-Cell Communication

Script: "04_CellChat.R"

Cell-cell communication was investigated using CellChat.

Communication networks were analyzed separately for Normal and ER+ samples and compared to identify changes in:

- Interaction number
- Interaction strength
- Signaling pathways
- Sender and receiver cell populations

This analysis was used to investigate how intercellular communication changes in ER+ breast cancer.

---

6. Pseudobulk Differential Expression

Script: "05_Pseudobulk_DESeq2.R"

To perform cell-type-specific differential expression while preserving biological replicates, single-cell expression data were aggregated into pseudobulk profiles by sample and cell type.

Pseudobulk expression profiles were compared between:

ER+ vs Normal

using DESeq2.

Significant transcriptional changes were subsequently investigated using:

- GSEA
- ORA

This enabled identification of biological pathways associated with ER+ breast cancer within specific cell populations.

---

7. Bulk RNA-seq Deconvolution

Script: "06_MuSiC_Deconvolution.R"

The bulk RNA-seq expression matrix was analyzed using MuSiC.

The annotated scRNA-seq dataset was used as a reference to estimate the relative abundance of cell types within the bulk RNA-seq samples.

Estimated cell-type proportions were compared between ER+ and Normal samples and visualized using stacked bar plots and boxplots.

This provides an independent bulk-level assessment of changes in cellular composition.

---

**Results & Visualizations**

The repository contains figures generated throughout the analysis, including:

scRNA-seq

- QC plots
- PCA and Elbow plots
- UMAP visualizations
- Cluster marker heatmaps
- Cell-type annotation plots
- DotPlots

CellChat

- Cell-cell interaction networks
- Interaction strength comparisons
- Signaling pathway analyses

Pseudobulk

- Differential expression results
- Volcano plots
- Heatmaps
- GSEA results
- ORA results

Bulk Deconvolution

- Estimated cell-type proportions
- Stacked bar plots
- Cell-type proportion boxplots

---


**Tools & Packages**

Analysis| Tool / Package
scRNA-seq analysis| Seurat
Cell-type annotation| SingleR
Cell-cell communication| CellChat
Differential expression| DESeq2
Pathway analysis| GSEA / ORA
Bulk deconvolution| MuSiC
Transcript quantification| Salmon
Gene annotation| GENCODE v50
Data processing| R / Galaxy
By combining cellular composition, cell-cell communication, cell-type-specific differential expression, pathway analysis, and bulk deconvolution, this project aims to connect cellular heterogeneity with molecular changes in the ER+ breast cancer microenvironment.
Project Objective

The overall objective is to integrate single-cell and bulk transcriptomic data to obtain a comprehensive view of ER+ breast cancer.
---
**Author**
**Rawan Salah**
Biotechnology graduate | Faculty of Science


