# Analyze Protein Co-Expression across samples
# authors: Pacome Prompsy
# contact: pacome.prompsy@chuv.ch
# Guenova Lab
# CHUV (Centre Hospitalier Universitaire Vaudois), Lausanne, Suisse

cat("Analyzing results from Protein Co-Expression across samples, condition and celltypes detection... \n")

# Loading packages --------------------------------------------------------
libraries = c("argparse",
              "ggplot2",
              "dplyr",
              "tidyr",
              "Seurat",
              "SpatialExperiment",
              "arrow",
              "ChromSCape",
              "ComplexHeatmap")
suppressPackageStartupMessages(invisible(lapply(libraries, require, character.only = TRUE)))
setwd("/mnt/RECHERCHE/GUENOVA_LAB/Project_Multiplex_Phenotyping/Miltenyi/Muliplex_Imaging_Pipeline/")

source("scripts/MultiplexImaging_utils.R")
source("scripts/GlobalVars.R")
library(ggplot2)
library(tidyverse)

# Reading in data  -------------------------------------------------------------
args = list(output = "output/")
cat("Output = ", args$output, "\n")

output_dir = file.path(args$output, "ProteinCoExpression")
if(!file.exists(output_dir)) dir.create(output_dir)

spe = qs::qread("output/SpatialExperiment.qs")

# Global Protein Co-expression  ------------------------------------------------
protein_cor = cor(t(spe@assays@data$normcounts))
# protein_cor = scale(protein_cor, center = T, scale = TRUE)


pdf(file.path(output_dir, paste0("Global_Protein_CoExpression_heatmap.pdf")),
    width = 15,
    height = 10)
protein_cor = cor(t(spe@assays@data$normcounts))
hc_cor = hclust(as.dist(1 - protein_cor), method = "ward.D2")
h = Heatmap( 
  protein_cor,
  cluster_rows = hc_cor,
  cluster_columns = hc_cor, 
  name = paste0("Global Protein Co-expression"),
  column_title = "Proteins",
  row_title = "Proteins",
  row_dend_side = "right",
  row_dend_width = unit(4, "cm"),
  column_dend_height = unit(5, "cm"),
  show_column_names = T,
  clustering_distance_columns ="pearson",
  clustering_distance_rows = "pearson",
  row_split = 10,
  column_split = 10,
  border = TRUE,
  row_names_gp = gpar(fontsize = 8),
  use_raster = FALSE
)
draw(h)
dev.off()

# Global Protein Co-expression - Binary ------------------------------------------------
jaccard. = function (x, y) {
  M.11 = sum(x == 1 & y == 1)
  M.10 = sum(x == 1 & y == 0)
  M.01 = sum(x == 0 & y == 1)
  return (M.11 / (M.11 + M.10 + M.01))
}

jaccard <- function(m){
  jac_mat = matrix(0,ncol = nrow(m), nrow= nrow(m),
                   dimnames = list(rownames(m),
                                   rownames(m))) 
  for(row in 1:nrow(m)){
    for(col in 1:nrow(m)){
      jac_mat[row,col] = jaccard.(m[row,],m[col,])
    }
  }
  return(jac_mat)
}

pdf(file.path(output_dir, paste0("Global_Protein_normcounts_correlation.pdf")),
    width = 15,
    height = 10)
protein_cor = spe@assays@data$normcounts
protein_cor = protein_cor[setdiff(rownames(protein_cor), markers_to_remove), setdiff(colnames(protein_cor), markers_to_remove)]
# protein_cor = jaccard(protein_cor)
protein_cor = cor(t(protein_cor))
corrplot::corrplot(protein_cor, type = "lower", order = "hclust", method = "color",
                   col = rev(corrplot::COL2('RdBu', 200)), col.lim = c(-1, 1),
                   diag = T, tl.cex = 0.8, tl.col = "black")
dev.off()

for(celltype in unique(spe$cell_type_CellSighter)){
  cells = c();
  for(samp in unique(spe$sample_id)){
    cells. = spe$cell_id[spe$sample_id == samp & spe$cell_type_CellSighter %in% "T_helper"]
    cells = c(cells, sample(cells., min(length(cells.), 500)))
  }
  
  spe. = spe[,match(cells, spe$cell_id)]
  spe. = reduce_dims_scExp(spe., n = 10, remove_PC = NULL, dimension_reductions = "PCA",
                           batch_correction = TRUE, batch_list = list(batch_1 = unique(spe.$sample_id[spe.$batch == "batch_1"]),
                                                                      batch_2 = unique(spe.$sample_id[spe.$batch == "batch_2"])
                                                                                      ))
  
  plot_reduced_dim_scExp(spe., color_by = "batch", reduced_dim = "PCA")
  
  # Cor based on PCA
  spe. = ChromSCape::correlation_and_hierarchical_clust_scExp(spe.)
  
  # Cor based on normcounts
  # mat <- spe.@assays@data$normcounts
  # cor_mat = coop::pcor(mat, inplace = TRUE)
  # reducedDim(spe., "Cor") = as(cor_mat,"dspMatrix")

  spe.. = spe.
  samp = intra_correlation_scExp(spe.., by = "condition", fullCor = TRUE)
  samp$condition = factor(samp$condition, levels = levels(spe..$condition))
  
  annot = SingleCellExperiment::colData(spe..)
  by = "condition"
  p <- ggplot(samp, aes(x = .data[[by]], y = .data$intra_corr, 
                        fill = .data[[by]])) + geom_violin(scale = "width") + theme_classic() + 
    theme(axis.text.x = element_text(angle = 90)) + ylab(paste0("Intra-", 
                                                                by, " correlation")) + xlab("")
  cols = unique(as.character(annot[, paste0(by, "_color")][match(rownames(samp), 
                                                                 annot$cell_id)]))
  names(cols) = unique(as.character(annot[, by][match(rownames(samp), 
                                                      annot$cell_id)]))
  p <- p + scale_fill_manual(values = cols)
  
  
  library(ggpubr)
  png(file.path(output_dir, "intracorrelation_violin_plot_normcounts_based_on_pca_batch_corrected_T_helper.png"), 
      width = 2000, height = 1500, res = 450)
  p + stat_compare_means(aes(label = after_stat(p.signif)),
                                method = "t.test", ref.group = "HD") + ylim(c(0.25, 1))
  dev.off()
  
  library(ggpubr)
  png(file.path(output_dir, "correlation_heatmap_normcounts_based_on_pca_batch_corrected_T_helper.png"), 
      width = 2000, height = 2000, res = 300)
  ChromSCape::plot_heatmap_scExp(spe.., color_by = c("condition", "sample_id", "batch") )
  dev.off()
  
  
}

# In Tumor
pdf(file.path(output_dir, paste0("Global_Protein_Jaccard_CoExpression_CellSighter_corrplot_tumor.pdf")),
    width = 15,
    height = 10)
spe. = spe[,spe$condition %in% c("MF", "SS")]
protein_cor = spe.@assays@data$CellSighter_marker_mat
protein_cor = protein_cor[setdiff(rownames(protein_cor), markers_to_remove), setdiff(colnames(protein_cor), markers_to_remove)]
protein_cor = jaccard(protein_cor)
corrplot::corrplot(protein_cor, type = "upper", order = "hclust",
                   col = rev(corrplot::COL2("RdBu", 200)), col.lim = c(-0.4, 1),
                   diag = F, tl.cex = 0.8, tl.col = "black")
dev.off()

# Non Tumor
pdf(file.path(output_dir, paste0("Global_Protein_Jaccard_CoExpression_CellSighter_corrplot_non_tumor.pdf")),
    width = 15,
    height = 10)
spe. = spe[,!spe$condition %in% c("MF", "SS")]
protein_cor = spe.@assays@data$CellSighter_marker_mat
protein_cor = protein_cor[setdiff(rownames(protein_cor), markers_to_remove), setdiff(colnames(protein_cor), markers_to_remove)]
protein_cor = jaccard(protein_cor)
corrplot::corrplot(protein_cor, type = "upper", order = "hclust",
                   col = rev(corrplot::COL2("RdBu", 200)), col.lim = c(-0.4, 1),
                   diag = F, tl.cex = 0.8, tl.col = "black")
dev.off()

# TMA1
pdf(file.path(output_dir, paste0("Global_Protein_Jaccard_CoExpression_CellSighter_corrplot_TMA1.pdf")),
    width = 15,
    height = 10)
spe. = spe[,spe$batch %in% c("batch_1")]
protein_cor = spe.@assays@data$CellSighter_marker_mat
protein_cor = protein_cor[setdiff(rownames(protein_cor), markers_to_remove), setdiff(colnames(protein_cor), markers_to_remove)]
protein_cor = jaccard(protein_cor)
corrplot::corrplot(protein_cor, type = "upper", order = "hclust",
                   col = rev(corrplot::COL2("RdBu", 200)), col.lim = c(-0.4, 1),
                   diag = F, tl.cex = 0.8, tl.col = "black")
dev.off()

# TMA2
pdf(file.path(output_dir, paste0("Global_Protein_Jaccard_CoExpression_CellSighter_corrplot_TMA2.pdf")),
    width = 15,
    height = 10)
spe. =  spe[,spe$batch %in% c("batch_2")]
protein_cor = spe.@assays@data$CellSighter_marker_mat
protein_cor = protein_cor[setdiff(rownames(protein_cor), markers_to_remove), setdiff(colnames(protein_cor), markers_to_remove)]
protein_cor = jaccard(protein_cor)
corrplot::corrplot(protein_cor, type = "upper", order = "hclust",
                   col = rev(corrplot::COL2("RdBu", 200)), col.lim = c(-0.4, 1),
                   diag = F, tl.cex = 0.8, tl.col = "black")
dev.off()

# Per condition
pdf(file.path(output_dir, paste0("Global_Protein_Jaccard_CoExpression_CellSighter_corrplot_per_condition.pdf")),
    width = 15,
    height = 10)

for(i in unique(spe$condition)){
  print(i)
  spe. =  spe[,spe$condition == i]
  protein_cor = spe.@assays@data$CellSighter_marker_mat
  protein_cor = protein_cor[setdiff(rownames(protein_cor), markers_to_remove), setdiff(colnames(protein_cor), markers_to_remove)]
  protein_cor = jaccard(protein_cor)
  print(corrplot::corrplot(protein_cor, type = "upper", order = "hclust",
                     col = rev(corrplot::COL2("RdBu", 200)), col.lim = c(-0.4, 1),
                     diag = F, tl.cex = 0.8, tl.col = "black", title = i))
}
dev.off()

# Sample By Sample Protein Co-expression  ------------------------------------------------
pdf(file.path(output_dir, paste0("Sample_Protein_CoExpression_heatmap.pdf")),
    width = 15,
    height = 10)
for(samp in unique(spe$sample_id)){
  spe. = spe[,spe$sample_id == samp]
  protein_cor = cor(t(spe.@assays@data$normcounts))
  hc_cor = hclust(as.dist(1 - protein_cor), method = "ward.D2")
  h = Heatmap( 
    protein_cor,
    cluster_rows = hc_cor,
    cluster_columns = hc_cor, 
    name = paste0("Protein Co-expression: ", samp),
    column_title = "Proteins",
    row_title = "Proteins",
    row_dend_side = "right",
    row_dend_width = unit(4, "cm"),
    column_dend_height = unit(5, "cm"),
    show_column_names = T,
    clustering_distance_columns ="pearson",
    clustering_distance_rows = "pearson",
    row_split = 10,
    column_split = 10,
    border = TRUE,
    row_names_gp = gpar(fontsize = 8),
    use_raster = FALSE
  )
  draw(h)
}
dev.off()


# Condition By Condition Protein Co-expression  ------------------------------------------------
pdf(file.path(output_dir, paste0("Condition_Protein_CoExpression_heatmap.pdf")),
    width = 15,
    height = 10)
for(condition in unique(spe$condition)){
  spe. = spe[,spe$condition == condition]
  protein_cor = cor(t(spe.@assays@data$normcounts))
  hc_cor = hclust(as.dist(1 - protein_cor), method = "ward.D2")
  h = Heatmap( 
    protein_cor,
    cluster_rows = hc_cor,
    cluster_columns = hc_cor, 
    name = paste0("Protein Co-expression: ", condition),
    column_title = "Proteins",
    row_title = "Proteins",
    row_dend_side = "right",
    row_dend_width = unit(4, "cm"),
    column_dend_height = unit(5, "cm"),
    show_column_names = T,
    clustering_distance_columns ="pearson",
    clustering_distance_rows = "pearson",
    row_split = 10,
    column_split = 10,
    border = TRUE,
    row_names_gp = gpar(fontsize = 8),
    use_raster = FALSE
  )
  draw(h)
}
dev.off()


# CellType By CellType Protein Co-expression  ------------------------------------------------
pdf(file.path(output_dir, paste0("CellType_Protein_CoExpression_heatmap.pdf")),
    width = 15,
    height = 10)
for(celltype in unique(spe$cell_type_CellSighter)){
  spe. = spe[,spe$cell_type_CellSighter == celltype]
  protein_cor = cor(t(spe.@assays@data$normcounts))
  hc_cor = hclust(as.dist(1 - protein_cor), method = "ward.D2")
  h = Heatmap( 
    protein_cor,
    cluster_rows = hc_cor,
    cluster_columns = hc_cor, 
    name = paste0("Protein Co-expression: ", celltype),
    column_title = "Proteins",
    row_title = "Proteins",
    row_dend_side = "right",
    row_dend_width = unit(4, "cm"),
    column_dend_height = unit(5, "cm"),
    show_column_names = T,
    clustering_distance_columns ="pearson",
    clustering_distance_rows = "pearson",
    row_split = 10,
    column_split = 10,
    border = TRUE,
    row_names_gp = gpar(fontsize = 8),
    use_raster = FALSE
  )
  draw(h)
}
dev.off()

# Condition by CellType Protein Co-expression  ------------------------------------------------
pdf(file.path(output_dir, paste0("Condition_by_CellType_Protein_CoExpression_heatmap.pdf")),
    width = 15,
    height = 10)
for(condition in unique(spe$condition)){
  spe. = spe[,spe$condition == condition]
  for(celltype in unique(spe.$cell_type_CellSighter)){
    spe.. = spe.[,spe.$cell_type_CellSighter == celltype]
    if(ncol(spe..) > 50){
      protein_cor = cor(t(spe..@assays@data$normcounts))
      hc_cor = hclust(as.dist(1 - protein_cor), method = "ward.D2")
      h = Heatmap( 
        protein_cor,
        cluster_rows = hc_cor,
        cluster_columns = hc_cor, 
        name = paste0(condition, " x ", celltype),
        column_title = "Proteins",
        row_title = "Proteins",
        row_dend_side = "right",
        row_dend_width = unit(4, "cm"),
        column_dend_height = unit(5, "cm"),
        show_column_names = T,
        clustering_distance_columns ="pearson",
        clustering_distance_rows = "pearson",
        row_split = 10,
        column_split = 10,
        border = TRUE,
        row_names_gp = gpar(fontsize = 8),
        use_raster = FALSE
      )
      draw(h)
    }
  }
}
dev.off()


# ------------------------------------------------
# Protein - Condition - Marker Heatmap - T cells
# ------------------------------------------------
spe. = spe[,spe$cell_type_CellSighter %in% T_cells & !spe$condition %in% c("THY", "TON", "LN")]


df = as.data.frame(t(spe.@assays@data$CellSighter_marker_mat[setdiff(rownames(spe.), markers_to_remove ),]))
df$condition = spe.$condition
df = df %>% group_by(condition) %>%
  dplyr::summarise(across(all_of(colnames(df)[-length(colnames(df))]), mean))
mat = as.matrix(df[,-1])
rownames(mat) = df$condition


pdf(file.path(output_dir, paste0("Protein_x_Condition_in_Tcells_PMD.pdf")),
    width = 15,
    height = 10)
ha = rowAnnotation(df = data.frame(condition = df$condition),
                   col = list("condition" = setNames(unique(spe.$condition_color), unique(spe.$condition))))

colsum = colSums(mat)
rowsum = rowSums(mat)

mat = mat / colsum
mat = mat / rowsum

h = Heatmap( 
  mat,
  name = paste0("Condition x Markers"),
  column_title = "Proteins",
  row_title = "Conditions",
  row_dend_side = "right",
  row_dend_width = unit(4, "cm"),
  column_dend_height = unit(5, "cm"),
  show_column_names = T,
  clustering_distance_columns ="pearson",
  clustering_distance_rows = "pearson",
  row_split = 10,
  column_split = 10,
  border = TRUE,
  row_names_gp = gpar(fontsize = 8),
  use_raster = FALSE,
  right_annotation = ha
)
draw(h)
dev.off()

dev.off()




