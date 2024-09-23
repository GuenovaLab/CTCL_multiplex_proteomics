# CTCL_multiplex_proteomics
Subcellular cell interactions reveal immune evasion in cutaneous T-cell lymphoma

# Description
This repository hosts the code for the paper "Subcellular cell interactions reveal immune evasion in cutaneous T-cell lymphoma," enabling the reproduction of figures and analyses presented in the study. The repository contains scripts in R, Python, and Bash for data processing, quality control, deep learning-based cell type classification, and spatial interaction analysis. The data used in the study can be accessed from the Dryad database at XXXXX.com, and an interactive tour of the results is available at YYYY.com.

# Step-by-step description

## 1. QC and segmentation
This folder contains scripts for preprocessing and quality control (QC) of raw imaging data. It includes helpers for marker quality assessment, input filtering, nuclear detection, pixel classification, and normalization of spatial data. Key steps involve creating the spatial experiment and segmenting cells based on pixel classification.

## 2. Deeplearning of cell type and PMD
This section covers cell type identification and positive marker detection (PMD) using deep learning. Scripts to guide the manual annotation of markers and cell types, training the CellSighter model for classification, and predicting cell types and marker positivity from images. QC tools for assessing cell phenotyping are also included.

## 3. Average intensity-based celltype and PMD (optional)
This folder provides optional scripts for an alternative method of cell type classification and PMD based on average cell intensity. It includes methods for clustering cells using Seurat and FlowSOM and average intensity-based marker expression analysis.

## 4. Region and Interactions
Scripts in this section allow for defining tissue regions and analyzing cell-to-cell interactions. This includes calculating distances between cells, identifying interacting cell borders (direct-contact), and performing region-based analysis of cellular interactions and spatial clustering.

## 5. Cellular Junction
This folder contains scripts to analyze cellular junctions, focusing on interactions at cell boundaries and spatial organization of neighboring cells.

## 6. Protein expression and celltype composition
This folder focuses on analyzing protein expression and cell type composition. It includes scripts to compare malignant versus benign cell types, protein co-expression analysis, and evaluating disease-specific and responder versus non-responder cases.

These folders together allow for a comprehensive analysis pipeline from raw image processing to spatial proteomics and cell type identification.

# How to cite
