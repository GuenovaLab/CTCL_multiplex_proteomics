# Finds cell clusters using cell average intensities
# Seurat or FlowSOM
# authors: Pacome Prompsy
# contact: pacome.prompsy@chuv.ch
# Guenova Lab
# CHUV (Centre Hospitalier Universitaire Vaudois), Lausanne, Suisse

cat("Finds cell clusters using cell average intensities (Seurat, FlowSOM)... \n")

# Loading packages --------------------------------------------------------
libraries = c("argparse",
              "ggplot2",
              "dplyr",
              "tidyr",
              "ChromSCape",
              "SpatialExperiment",
              "FlowSOM",
              "Seurat")
suppressPackageStartupMessages(invisible(lapply(libraries, require, character.only = TRUE)))
setwd("/mnt/RECHERCHE/GUENOVA_LAB/Project_Multiplex_Phenotyping/Miltenyi/Muliplex_Imaging_Pipeline/")

source("scripts/MultiplexImaging_utils.R")
source("scripts/GlobalVars.R")

# Directorires --------------------------------------------
output = "./output/" 

cat("Output = ", output, "\n")


output_dir = file.path(output, "cell_based_clustering")
if (!dir.exists(output_dir))
  dir.create(output_dir)

# Loading the dataset ----------------------------------------------------------
spe = qs::qread(file.path(output, "SpatialExperiment.qs"))

# Dimensionality Reduction -----------------------------------------------------
spe = reduce_dims_scExp(spe, n = 10, remove_PC = NULL)
plot_correlation_PCA_scExp(spe,
                           correlation_var = "log10_area_normalized_intensity",
                           color_by = "sample_id", topPC = 10)

plot_correlation_PCA_scExp(spe,
                           correlation_var = "log10_area_normalized_intensity", topPC = 10)
pdf(
  file.path(output_dir, "DimensionalityReduction.pdf"),
  width =  10,
  height = 6
)
plot_reduced_dim_scExp(spe, color_by = "sample_id", reduced_dim = "PCA", downsample = 1e6, size = 0.25, min_quantile = 0, max_quantile = 1)
plot_reduced_dim_scExp(spe, color_by = "log10_area_normalized_intensity", reduced_dim = "PCA", downsample = 1e6, size = 0.25, min_quantile = 0, max_quantile = 1)
plot_reduced_dim_scExp(spe, color_by = "sample_id", reduced_dim = "UMAP", downsample = 1e6, size = 0.25, min_quantile = 0, max_quantile = 1)
plot_reduced_dim_scExp(spe, color_by = "log10_area_normalized_intensity", reduced_dim = "UMAP", downsample = 1e6, size = 0.25, min_quantile = 0, max_quantile = 1)
plot_reduced_dim_scExp(spe, color_by = "condition", reduced_dim = "UMAP", downsample = 1e6, size = 1, min_quantile = 0, max_quantile = 1)
dev.off()

pdf(file.path(output_dir, "UMAP_normcounts_sample_highlight.pdf"))
for(i in unique(spe$sample_id)){
  spe$highlight = spe$sample_id == i 
  print(scater::plotReducedDim(spe, dimred = "UMAP", colour_by = "highlight", order_by = "highlight") + ggtitle(i) +
          scale_color_manual(name = i, values = c("grey75", "darkblue"))
  )
}
dev.off()


pdf(
  file.path(output_dir, "DimensionalityReduction-CellType-Marker-Based.pdf"),
  width =  7,
  height = 6
)
plot_reduced_dim_scExp(spe, color_by = "cell_type_positive_marker", reduced_dim = "UMAP", size = 0.25,  downsample = 1e6, min_quantile = 0, max_quantile = 1)
dev.off()

################################################################################
# Dimensionality Reduction with Seurat (using all proteins)
################################################################################
run = TRUE
if(run == TRUE){
  library(Seurat)
  library(harmony)
  library(ggrastr)
  logcounts(spe) = normcounts(spe)
  counts(spe) = logcounts(spe)
  
  Seu = as.Seurat(spe)
  for(i in names(Seu@reductions)) Seu@reductions[[i]] = NULL
  
  Seu = ScaleData(Seu, features = rownames(Seu))
  Seu = RunPCA(Seu, npcs = 10, features = rownames(Seu))
  Seu = RunHarmony(Seu, group.by.vars = "sample_id", dims.use = 1:10, reduction = "pca")
  source('../../../Pacome/GitLab/FIt-SNE-master/fast_tsne.R', chdir=T)
  tsne = fftRtsne(Seu@reductions$harmony@cell.embeddings, perplexity = 20, nthreads = 8)
  Seu@reductions$tsne = SeuratObject::CreateDimReducObject(tsne, key = "tsne_")
  # Seu = RunTSNE(Seu,tsne.method = "FIt-SNE", reduction = "harmony", dims = 1:10)
  Seu = RunUMAP(Seu, dims = 1:10,  reduction = "harmony")
  
  png(
    file.path(output_dir, paste0("TSNE_full_sample_id_QN.png")),
    width =  2500,
    height = 1800,
    res = 300
  )
  df = as.data.frame(Seu@reductions$tsne@cell.embeddings[,1:2])
  df$sample_id = Seu$sample_id
  p = df %>% ggplot() + geom_point(aes(x = tsne_1, y = tsne_2, color = sample_id),
                                   pch = 20, size = 0.05, alpha = 0.5) + 
    theme_bw() + scale_color_manual(values = unique(Seu$sample_id_color))
  print(p)
  dev.off()
  
  
  png(
    file.path(output_dir, paste0("TSNE_full_celltype_QN.png")),
    width =  2500,
    height = 1800,
    res = 300
  )
  df = as.data.frame(Seu@reductions$tsne@cell.embeddings[,1:2])
  df$cell_type_CellSighter = Seu$cell_type_CellSighter
  p = df %>% ggplot() + geom_point(aes(x = tsne_1, y = tsne_2, color = cell_type_CellSighter),
                                   pch = 20, size = 0.05, alpha = 0.5) + 
    theme_bw() + scale_color_manual(values = setNames(unique(Seu$cell_type_CellSighter_color), unique(Seu$cell_type_CellSighter)))
  print(p)
  dev.off()
  
  png(
    file.path(output_dir, paste0("TSNE_full_condition_QN.png")),
    width =  2500,
    height = 1800,
    res = 300
  )
  df = as.data.frame(Seu@reductions$tsne@cell.embeddings[,1:2])
  df$condition = Seu$condition
  p = df %>% ggplot() + geom_point(aes(x = tsne_1, y = tsne_2, color = condition),
                                   pch = 20, size = 0.05, alpha = 0.5) + 
    theme_bw() + scale_color_manual(values = setNames(unique(Seu$condition_color), unique(Seu$condition)))
  print(p)
  dev.off()
  
  png(
    file.path(output_dir, paste0("UMAP_full_sample_id_QN.png")),
    width =  2500,
    height = 1800,
    res = 300
  )
  df = as.data.frame(Seu@reductions$umap@cell.embeddings[,1:2])
  df$sample_id = Seu$sample_id
  p = df %>% ggplot() + geom_point(aes(x = UMAP_1, y = UMAP_2, color = sample_id),
                                   pch = 20, size = 0.05, alpha = 0.5) + 
    theme_bw() + scale_color_manual(values = unique(Seu$sample_id_color))
  print(p)
  dev.off()
  
  
  png(
    file.path(output_dir, paste0("UMAP_full_celltype_QN.png")),
    width =  2500,
    height = 1800,
    res = 300
  )
  df = as.data.frame(Seu@reductions$umap@cell.embeddings[,1:2])
  df$cell_type_CellSighter = Seu$cell_type_CellSighter
  p = df %>% ggplot() + geom_point(aes(x = UMAP_1, y = UMAP_2, color = cell_type_CellSighter),
                                   pch = 20, size = 0.05, alpha = 0.5) + 
    theme_bw() + scale_color_manual(values = setNames(unique(Seu$cell_type_CellSighter_color), unique(Seu$cell_type_CellSighter)))
  print(p)
  dev.off()
  
  png(
    file.path(output_dir, paste0("UMAP_full_condition_QN.png")),
    width =  2500,
    height = 1800,
    res = 300
  )
  df = as.data.frame(Seu@reductions$umap@cell.embeddings[,1:2])
  df$condition = Seu$condition
  p = df %>% ggplot() + geom_point(aes(x = UMAP_1, y = UMAP_2, color = condition),
                                   pch = 20, size = 0.05, alpha = 0.5) + 
    theme_bw() + scale_color_manual(values = setNames(unique(Seu$condition_color), unique(Seu$condition)))
  print(p)
  dev.off()
  qs::qsave(Seu, "output/Seu_all_QN.qs")

  
}

# Save in SPE
spe = qs::qread("output/SpatialExperiment.qs")

tsne_harmony = DataFrame(Seu@reductions$tsne@cell.embeddings)
colnames(tsne_harmony) = c("Component_1", "Component_2")
spe@int_colData$reducedDims$TSNE_harmony = tsne_harmony

umap_harmony = DataFrame(Seu@reductions$umap@cell.embeddings)
colnames(umap_harmony) = c("Component_1", "Component_2")
spe@int_colData$reducedDims$UMAP_harmony = umap_harmony

pca_harmony = DataFrame(Seu@reductions$harmony@cell.embeddings)
colnames(pca_harmony) = paste0("Component_", 1:10)
spe@int_colData$reducedDims$PCA_harmony = pca_harmony

qs::qsave(spe, "output/SpatialExperiment.qs")

################################################################################
# Dimensionality Reduction with Seurat (using phenotyping markers only)
################################################################################
if(run == TRUE){
  library(Seurat)
  library(harmony)
  library(ggrastr)
  logcounts(spe) = limma::normalizeQuantiles(counts(spe))
  normcounts(spe) = logcounts(spe)
  counts(spe) = logcounts(spe)
  spe. = spe[unique(cell_type_markers$marker), !spe$condition %in% c("THY", "LN", "TON")]
  Seu = as.Seurat(spe.)
  for(i in names(Seu@reductions)) Seu@reductions[[i]] = NULL
  
  Seu = ScaleData(Seu, features = rownames(Seu))
  Seu = RunPCA(Seu, npcs = 10, features = rownames(Seu))
  Seu = RunHarmony(Seu, group.by.vars = "sample_id", dims.use = 1:10, reduction = "pca")
  source('../../../Pacome/GitLab/FIt-SNE-master/fast_tsne.R', chdir=T)
  tsne = fftRtsne(Seu@reductions$harmony@cell.embeddings, perplexity = 20, nthreads = 8)
  Seu@reductions$tsne = SeuratObject::CreateDimReducObject(tsne, key = "tsne_")
  # Seu = RunTSNE(Seu,tsne.method = "FIt-SNE", reduction = "harmony", dims = 1:10)
  Seu = RunUMAP(Seu, dims = 1:10,  reduction = "harmony")
  
  png(
    file.path(output_dir, paste0("TSNE_cellmarkers_sample_id_QN.png")),
    width =  2500,
    height = 1800,
    res = 300
  )
  df = as.data.frame(Seu@reductions$tsne@cell.embeddings[,1:2])
  df$sample_id = Seu$sample_id
  p = df %>% ggplot() + geom_point(aes(x = tsne_1, y = tsne_2, color = sample_id),
                                   pch = 20, size = 0.05, alpha = 0.5) + 
    theme_bw() + scale_color_manual(values = unique(Seu$sample_id_color))
  print(p)
  dev.off()
  
  
  png(
    file.path(output_dir, paste0("TSNE_cellmarkers_celltype_QN.png")),
    width =  2500,
    height = 1800,
    res = 300
  )
  df = as.data.frame(Seu@reductions$tsne@cell.embeddings[,1:2])
  df$cell_type_CellSighter = Seu$cell_type_CellSighter
  p = df %>% ggplot() + geom_point(aes(x = tsne_1, y = tsne_2, color = cell_type_CellSighter),
                                   pch = 20, size = 0.05, alpha = 0.5) + 
    theme_bw() + scale_color_manual(values = setNames(unique(Seu$cell_type_CellSighter_color), unique(Seu$cell_type_CellSighter)))
  print(p)
  dev.off()
  
  png(
    file.path(output_dir, paste0("TSNE_cellmarkersl_condition_QN.png")),
    width =  2500,
    height = 1800,
    res = 300
  )
  df = as.data.frame(Seu@reductions$tsne@cell.embeddings[,1:2])
  df$condition = Seu$condition
  p = df %>% ggplot() + geom_point(aes(x = tsne_1, y = tsne_2, color = condition),
                                   pch = 20, size = 0.05, alpha = 0.5) + 
    theme_bw() + scale_color_manual(values = setNames(unique(Seu$condition_color), unique(Seu$condition)))
  print(p)
  dev.off()
  
  png(
    file.path(output_dir, paste0("UMAP_cellmarkers_sample_id_QN.png")),
    width =  2500,
    height = 1800,
    res = 300
  )
  df = as.data.frame(Seu@reductions$umap@cell.embeddings[,1:2])
  df$sample_id = Seu$sample_id
  p = df %>% ggplot() + geom_point(aes(x = UMAP_1, y = UMAP_2, color = sample_id),
                                   pch = 20, size = 0.05, alpha = 0.5) + 
    theme_bw() + scale_color_manual(values = unique(Seu$sample_id_color))
  print(p)
  dev.off()
  
  
  png(
    file.path(output_dir, paste0("UMAP_cellmarkers_celltype_QN.png")),
    width =  2500,
    height = 1800,
    res = 300
  )
  df = as.data.frame(Seu@reductions$umap@cell.embeddings[,1:2])
  df$cell_type_CellSighter = Seu$cell_type_CellSighter
  p = df %>% ggplot() + geom_point(aes(x = UMAP_1, y = UMAP_2, color = cell_type_CellSighter),
                                   pch = 20, size = 0.05, alpha = 0.5) + 
    theme_bw() + scale_color_manual(values = setNames(unique(Seu$cell_type_CellSighter_color), unique(Seu$cell_type_CellSighter)))
  print(p)
  dev.off()
  
  png(
    file.path(output_dir, paste0("UMAP_cellmarkers_condition_QN.png")),
    width =  2500,
    height = 1800,
    res = 300
  )
  df = as.data.frame(Seu@reductions$umap@cell.embeddings[,1:2])
  df$condition = Seu$condition
  p = df %>% ggplot() + geom_point(aes(x = UMAP_1, y = UMAP_2, color = condition),
                                   pch = 20, size = 0.05, alpha = 0.5) + 
    theme_bw() + scale_color_manual(values = setNames(unique(Seu$condition_color), unique(Seu$condition)))
  print(p)
  dev.off()
  qs::qsave(Seu, "output/Seu_cellmarkers_QN.qs")
  
  Seu = qs::qread("output/Seu_cellmarkers_QN.qs")
  
  # Clustering and Heatmap
  Seu = FindNeighbors(Seu, reduction = "harmony", dims = 1:10)
  Seu = FindClusters(Seu, res = 0.25)
  
  spe$seurat_cluster = Seu$seurat_clusters[match(colnames(spe), colnames(Seu))]
  qs::qsave(spe, "output/SpatialExperiment.qs")
  

  png(
    file.path(output_dir, paste0("UMAP_cellmarkers_cluster_QN.png")),
    width =  2500,
    height = 1800,
    res = 300
  )
  df = as.data.frame(Seu@reductions$umap@cell.embeddings[,1:2])
  df$seurat_clusters = Seu$seurat_clusters
  p = df %>% ggplot() + geom_point(aes(x = UMAP_1, y = UMAP_2, color = seurat_clusters),
                                   pch = 20, size = 0.05, alpha = 0.5) + 
    theme_bw() + scale_color_manual(values = discrette_colors_50[1:length(unique(Seu$seurat_clusters))])
  print(p)
  dev.off()
  
  
  library(ComplexHeatmap)
  
  cells = c()
  spe. = spe[,which(!is.na(spe$seurat_cluster))]
  for(i in unique(spe.$seurat_cluster)){
    cells = c(cells, sample(spe.$cell_id[which(spe.$seurat_cluster == i)], min(500, length(which(spe.$seurat_cluster == i))) ))
  }
  spe. = spe.[,match(cells, spe.$cell_id)]
  
  png(
    file.path(output_dir, paste0("Heatmap_cellmarkers_cluster_QN_sampled.png")),
    width =  2500,
    height = 1800,
    res = 300
  )
  
  set.seed(47)
  mat = normcounts(spe.)
  mat = t(mat)[order(spe.$seurat_cluster), c(golden_markers, "CollagenI")]
  
  h = Heatmap(
    scale(mat),
    name = "Quantile Normalized Counts",
    column_title = "Markers",
    row_title = "Seurat cluster",
    row_dend_side = "right",
    row_dend_width = unit(4, "cm"),
    column_dend_height = unit(5, "cm"), 
    cluster_columns = F,
    cluster_rows = FALSE,
    show_column_names = TRUE,
    show_row_names = FALSE,
    border = TRUE,
    row_names_gp = gpar(fontsize = 8),
    use_raster = FALSE,
    right_annotation = rowAnnotation(
      "Seurat cluster" = spe.$seurat_cluster[order(spe.$seurat_cluster)],
      col = list("Seurat cluster" = setNames(
        unique(discrette_colors_50[1:length(unique(spe.$seurat_cluster))]),
        unique(spe.$seurat_cluster)
      )))
  )
  draw(h)
  dev.off()                   
  

  
  
}

cat("Finished running cell clustering using Seurat ! \n")
cat("\n")

################################################################################
# FlowSOM clustering
################################################################################
library(FlowSOM)
library(flowCore)
library(RColorBrewer)

# Set the number of clusters (cell types) you expect
num_clusters <- 17
spe. = spe[,!spe$condition %in% c("THY", "TON", "LN")]

# Matrix
matrix = as.data.frame(t(normcounts(spe.)))
matrix$sample_id = spe.$sample_id
matrix$condition = spe.$condition
matrix$batch = spe.$batch
matrix$cell_type_CellSighter = spe.$cell_type_CellSighter

matrix = matrix[,unique(cell_type_markers$marker)]
matrix = as.matrix(matrix)

# Run FlowSOM clustering
flowsom_model <- FlowSOM(input = matrix, nClus = num_clusters, transform = FALSE, scale = FALSE, seed = 47)
qs::qsave(flowsom_model, file.path(output_dir, "flowsom_model_counts_qn.qs"))

# Assign cluster labels to each cell
cluster_labels <- FlowSOM::GetMetaclusters(flowsom_model)
table(cluster_labels)

# Create a data frame combining the cluster labels and the scaled intensity matrix
clustered_data <- data.frame(cluster = factor(cluster_labels), matrix)

# Sort the data frame by the cluster labels for better visualization in the heatmap
clustered_data <- clustered_data[order(clustered_data$cluster), ]

heatmap_col <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)  # Define heatmap colors

library(ComplexHeatmap)

spe$FlowSom_clusters = clustered_data$cluster[match(colnames(spe), rownames(clustered_data))]
qs::qsave(spe, "output/SpatialExperiment.qs")


cells = c()
spe. = spe[,which(!is.na(spe$FlowSom_clusters))]
for(i in unique(spe.$FlowSom_clusters)){
  cells = c(cells, sample(spe.$cell_id[which(spe.$FlowSom_clusters == i)], min(500, length(which(spe.$FlowSom_clusters == i))) ))
}
spe. = spe.[,match(cells, spe.$cell_id)]

png(file.path(output_dir,paste0("Heatmap_FlowSom_Clusters_counts_QN_sampled.png")),
    width = 2500,
    height = 1800, res = 300)

mat = t(normcounts(spe.))
mat = mat[order(spe.$FlowSom_clusters), c(golden_markers, "CollagenI")]

h = Heatmap(
  scale(mat),
  name = "Quantile Normalized Counts",
  column_title = "Markers",
  row_title = "FlowSOM cluster",
  row_dend_side = "right",
  row_dend_width = unit(4, "cm"),
  column_dend_height = unit(5, "cm"),
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  show_column_names = TRUE,
  show_row_names = FALSE,
  border = TRUE,
  row_names_gp = gpar(fontsize = 8),
  use_raster = FALSE,
  right_annotation = rowAnnotation(
    "FlowSOM cluster" = spe.$FlowSom_clusters,
    col = list("FlowSOM cluster" = setNames(
      unique(discrette_colors_50[1:length(unique(spe.$FlowSom_clusters))]),
      unique(spe.$FlowSom_clusters)
    )))
)
draw(h)
dev.off()                   


Seu$FlowSom_clusters = clustered_data$cluster[match(colnames(Seu), rownames(clustered_data))]

png(file.path(output_dir,paste0("Heamtap_Correspondance_FlowSom_Clusters_counts.png")),
    width = 1500,
    height = 1500, res = 250)
confusion_df = as.data.frame((table(cluster_labels, spe.$cell_type_CellSighter)))
colnames(confusion_df) <- c("FlowSOM", "CellSighter", "Count")

p = (ggplot(data = confusion_df, aes(x = FlowSOM, y = CellSighter, fill = log10(Count+1), label = Count)) +
       geom_tile() +
       scale_fill_gradient2(midpoint = 0.5, low = "royalblue4", mid = "royalblue4",high = "gold") +
       labs(title = "Celltype confusion matrix",
            x = "True Class", y = "Predicted Class", fill = "log10(Cell Overlap + 1)") +
       geom_text(color = "white", size = 3.5)  + theme(axis.text.x = element_text(angle = 90)) +theme_classic()
) + ggtitle("CellSighter CellType identification ") + theme(axis.text.x = element_text(angle=90))
p
dev.off()



cat("Finished running Positive Marker Detection ! \n")
cat("\n")



