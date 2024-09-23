# Subcellular cell interactions reveal immune evasion in cutaneous T-cell lymphoma

This repository hosts the code for the paper "Subcellular cell interactions reveal immune evasion in cutaneous T-cell lymphoma," enabling the reproduction of figures and analyses presented in the study. The repository contains scripts in R, Python, and Bash for data processing, quality control, deep learning-based cell type classification, and spatial interaction analysis. The data used in the study can be accessed from the Dryad database at XXXXX.com, and an interactive tour of the results is available at [bit.ly/skin-atlas]bit.ly/skin-atlas.

# Description

## 1. QC and segmentation
This folder contains scripts for preprocessing and quality control (QC) of raw imaging data. It includes helpers for marker quality assessment, input filtering, nuclear detection, pixel classification, and normalization of spatial data. Key steps involve creating the spatial experiment and segmenting cells based on pixel classification.

## 2. Deeplearning of cell type and PMD
This section covers cell type identification and positive marker detection (PMD) using deep learning. Scripts to guide the manual annotation of markers and cell types, training the CellSighter model for classification, and predicting cell types and marker positivity from images. QC tools for assessing cell phenotyping are also included.

## 3. Average intensity-based celltype and PMD (optional)
This folder provides optional scripts for an alternative method of cell type classification and PMD based on average cell intensity. It includes methods for clustering cells using Seurat and FlowSOM and average intensity-based marker expression analysis.

## 4. Regions and Interactions
Scripts in this section allow for defining tissue regions and analyzing cell-to-cell interactions. This includes calculating distances between cells, identifying interacting cell borders (direct-contact), and performing region-based analysis of cellular interactions and spatial clustering.

## 5. Cellular Junctions
This folder contains scripts to analyze cellular junctions, focusing on interactions at cell boundaries and spatial organization of neighboring cells.

## 6. Protein expression and cell type composition
This folder focuses on analyzing protein expression and cell type composition. It includes scripts to compare malignant versus benign cell types, protein co-expression analysis, and evaluating disease-specific and responder versus non-responder cases.

# How to cite

# R session info
```
R version 4.4.1 (2024-06-14)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 22.04.4 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0 
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0

locale:
 [1] LC_CTYPE=fr_CH.UTF-8       LC_NUMERIC=C               LC_TIME=fr_CH.UTF-8        LC_COLLATE=fr_CH.UTF-8     LC_MONETARY=fr_CH.UTF-8   
 [6] LC_MESSAGES=fr_CH.UTF-8    LC_PAPER=fr_CH.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=fr_CH.UTF-8 LC_IDENTIFICATION=C       

time zone: Europe/Zurich
tzcode source: system (glibc)

attached base packages:
[1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ComplexHeatmap_2.16.0       lubridate_1.9.2             forcats_1.0.0               stringr_1.5.0              
 [5] purrr_1.0.2                 readr_2.1.4                 tibble_3.2.1                tidyverse_2.0.0            
 [9] mclust_6.0.0                rstatix_0.7.2               ggbeeswarm_0.7.2            ggnewscale_0.4.9           
[13] ggpubr_0.6.0                ChromSCape_1.9.0            arrow_12.0.1                SpatialExperiment_1.10.0   
[17] SingleCellExperiment_1.22.0 SummarizedExperiment_1.30.2 Biobase_2.60.0              GenomicRanges_1.52.0       
[21] GenomeInfoDb_1.36.1         IRanges_2.34.1              S4Vectors_0.38.1            BiocGenerics_0.46.0        
[25] MatrixGenerics_1.12.2       matrixStats_1.0.0           Seurat_5.0.1.9004           SeuratObject_5.0.1         
[29] sp_2.0-0                    tidyr_1.3.0                 dplyr_1.1.2                 ggplot2_3.5.0              
[33] argparse_2.2.2             

loaded via a namespace (and not attached):
  [1] zoo_1.8-12                  doParallel_1.0.17           reticulate_1.30             pkgconfig_2.0.3            
  [5] pheatmap_1.0.12             later_1.3.1                 harmony_0.1.1               scales_1.3.0               
  [9] circlize_0.4.15             parallel_4.4.1              promises_1.2.0.1            spatstat.geom_3.2-1        
 [13] gtable_0.3.4                RcppHNSW_0.4.1              lattice_0.22-5              sctransform_0.4.1          
 [17] ggforce_0.4.1               dotCall64_1.0-2             systemfonts_1.0.5           RcppAnnoy_0.0.20           
 [21] glue_1.7.0                  S4Arrays_1.2.1              pbapply_1.7-2               KernSmooth_2.23-24         
 [25] bluster_1.10.0              shinydashboardPlus_2.0.3    munsell_0.5.0               ScaledMatrix_1.8.1         
 [29] data.table_1.14.8           shinycssloaders_1.0.0       bit_4.0.5                   restfulr_0.0.15            
 [33] Rhdf5lib_1.22.0             assertthat_0.2.1            pdist_1.2.1                 tidyselect_1.2.0           
 [37] memoise_2.0.1               scater_1.28.0               xml2_1.3.6                  processx_3.8.3             
 [41] sass_0.4.6                  shape_1.4.6                 vipor_0.4.5                 dqrng_0.3.0                
 [45] rlang_1.1.3                 shinyWidgets_0.7.6          parallelly_1.36.0           tzdb_0.4.0                 
 [49] openssl_2.0.6               generics_0.1.3              polyclip_1.10-4             qualV_0.3-4                
 [53] iterators_1.0.14            gridExtra_2.3               gggenes_0.5.0               jsonlite_1.8.5             
 [57] leiden_0.4.3                cluster_2.1.6               ResidualMatrix_1.10.0       viridisLite_0.4.2          
 [61] BiocParallel_1.34.2         FastKNN_0.0.1               shinyhelper_0.3.2           xtable_1.8-4               
 [65] RColorBrewer_1.1-3          DelayedArray_0.26.3         RcppParallel_5.1.7          rsvd_1.0.5                 
 [69] edgeR_3.42.4                beachmat_2.16.0             msigdbr_7.5.1               farver_2.1.1               
 [73] lmtest_0.9-40               rstudioapi_0.15.0           tiff_0.1-11                 RSpectra_0.16-1            
 [77] DT_0.28                     stringi_1.8.3               RProtoBufLib_2.12.0         shiny_1.7.4                
 [81] carData_3.0-5               pkgbuild_1.4.2              pillar_1.9.0                readxl_1.4.2               
 [85] HDF5Array_1.28.1            cachem_1.0.8                ggridges_0.5.4              car_3.1-2                  
 [89] zlibbioc_1.46.0             splines_4.4.1               Rsamtools_2.16.0            clue_0.3-64                
 [93] spatstat.utils_3.0-3        survival_3.6-4              yaml_2.3.8                  hms_1.1.3                  
 [97] uwot_0.1.15                 evaluate_0.21               codetools_0.2-19            BiocNeighbors_1.18.0       
[101] utf8_1.2.4                  ellipsis_0.3.2              knitr_1.43                  future_1.32.0              
[105] GetoptLong_1.0.5            viridis_0.6.3               dittoSeq_1.12.0             fitdistrplus_1.1-11        
[109] igraph_2.0.1.1              spam_2.9-1                  GenomicAlignments_1.36.0    beeswarm_0.4.0             
[113] magrittr_2.0.3              ggsignif_0.6.4              stringdist_0.9.10           batchelor_1.16.0           
[117] backports_1.4.1             babelgene_22.9              RANN_2.6.1                  progressr_0.13.0           
[121] ggfittext_0.10.0            ROCR_1.0-11                 cytolib_2.12.0              fs_1.6.3                   
[125] bslib_0.5.0                 crayon_1.5.2                deldir_1.0-9                cowplot_1.1.1              
[129] listenv_0.9.0               svglite_2.1.3               Rcpp_1.0.12                 askpass_1.2.0              
[133] labeling_0.4.3              abind_1.4-5                 prettyunits_1.2.0           reshape2_1.4.4             
[137] XVector_0.40.0              rlist_0.4.6.2               DelayedMatrixStats_1.22.1   scran_1.28.1               
[141] tensor_1.5                  rvest_1.0.3                 ica_1.0-3                   compiler_4.4.1             
[145] colorspace_2.1-0            patchwork_1.1.2             colorRamps_2.3.1            rhdf5filters_1.12.1        
[149] withr_3.0.0                 MASS_7.3-60                 digest_0.6.32               httpuv_1.6.11              
[153] spatstat.explore_3.2-1      png_0.1-8                   shinydashboard_0.7.2        matrixTests_0.2.2          
[157] ggsurvfit_1.0.0             Rtsne_0.16                  ggrepel_0.9.3               R.oo_1.25.0                
[161] httr_1.4.7                  nlme_3.1-165                jquerylib_0.1.4             spatstat.random_3.1-5      
[165] locfit_1.5-9.8              timechange_0.2.0            mime_0.12                   htmlwidgets_1.6.4          
[169] metapod_1.8.0               umap_0.2.10.0               flowCore_2.12.0             fastDummies_1.6.3          
[173] WriteXLS_6.4.0              RCurl_1.98-1.12             miniUI_0.1.1.1              shinyFiles_0.9.3           
[177] goftest_1.2-3               R6_2.5.1                    cli_3.6.2                   BiocSingular_1.16.0        
[181] globals_0.16.2              rhdf5_2.44.0                lazyeval_0.2.2              R.methodsS3_1.8.2          
[185] GlobalOptions_0.1.2         tweenr_2.0.2                lifecycle_1.0.4             foreach_1.5.2              
[189] kableExtra_1.3.4            rtracklayer_1.60.0          future.apply_1.11.0         limma_3.56.2               
[193] spatstat.sparse_3.0-2       GenomeInfoDbData_1.2.10     RApiSerialize_0.1.2         fansi_1.0.6                
[197] scuttle_1.10.1              coop_0.6-3                  magick_2.7.4                shinyjs_2.1.0              
[201] callr_3.7.3                 DropletUtils_1.20.0         Matrix_1.6-3                htmltools_0.5.7            
[205] R.utils_2.12.2              ps_1.7.6                    stringfish_0.15.8           Biostrings_2.68.1          
[209] bit64_4.0.5                 bitops_1.0-7                sparseMatrixStats_1.12.1    rjson_0.2.21               
[213] BiocIO_1.10.0               tools_4.4.1                 XML_3.99-0.14               irlba_2.3.5.1              
[217] spatstat.data_3.0-1         qs_0.25.5                   scattermore_1.2             statmod_1.5.0              
[221] fastmap_1.1.1               plotly_4.10.2               ConsensusClusterPlus_1.64.0 xfun_0.39                  
[225] colourpicker_1.2.0          vctrs_0.6.3                 FlowSOM_2.8.0               cellranger_1.1.0           
[229] plyr_1.8.8                  rmarkdown_2.22              broom_1.0.5                 webshot_0.5.5
```
