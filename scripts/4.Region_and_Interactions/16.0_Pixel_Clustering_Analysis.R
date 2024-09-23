# Analyze results from Pixie pixel clustering
# authors: Pacome Prompsy
# contact: pacome.prompsy@chuv.ch
# Guenova Lab
# CHUV (Centre Hospitalier Universitaire Vaudois), Lausanne, Suisse

cat("Running Positive Marker Detection based on average cell expression distribution... \n")

# Loading packages --------------------------------------------------------
libraries = c("argparse",
              "ggplot2",
              "dplyr",
              "tidyr",
              "Seurat",
              "SpatialExperiment",
              "arrow",
              "ChromSCape")
suppressPackageStartupMessages(invisible(lapply(libraries, require, character.only = TRUE)))
setwd("/mnt/RECHERCHE/GUENOVA_LAB/Project_Multiplex_Phenotyping/Miltenyi/Muliplex_Imaging_Pipeline/")

source("scripts/MultiplexImaging_utils.R")
source("scripts/GlobalVars.R")
# Reading in data  -------------------------------------------------------------
args = list(output = "output/")
cat("Output = ", args$output, "\n")

cell_output_dir = file.path(args$output, "pixie", "region_cell_output_dir")
pixel_output_dir = file.path(args$output, "pixie", "region_pixel_output_dir")
spe = qs::qread("output/SpatialExperiment.qs")

#################################################################################
# Plotting pixel heatmap 
#################################################################################

samples = unique(spe$sample_id)

pixel_list = list()
for(samp in samples){
  print(samp)
  pixel_mat. = arrow::read_feather(file.path(pixel_output_dir, "pixel_mat_data", paste0(samp, ".feather"))) 
  print(barplot(table(pixel_mat.$pixel_meta_cluster), col = discrette_colors_50[1:20]))
  
  set.seed(47); 
  # for(clust in unique(pixel_mat.$pixel_meta_cluster)){
  #   print(clust)
  #   pixel_mat.. = pixel_mat.[pixel_mat.$pixel_meta_cluster == clust,]
  #   idx = sample(nrow(pixel_mat..), 100, replace = T)
  #   pixel_mat.. = pixel_mat..[idx,]
  #   pixel_mat..$condition = spe$condition[match(pixel_mat..$fov, spe$sample_id)]
  #   pixel_list[[paste0(samp,"_",clust)]] = pixel_mat..
  # }
  gc()
}



pixel_mat = do.call("rbind", pixel_list)
qs::qsave(pixel_mat, file.path(pixel_output_dir, "pixel_mat.qs"))
###

pixel_mat = qs::qread(file.path(pixel_output_dir, "pixel_mat.qs"))
pixel_mat = read.csv(file.path(pixel_output_dir, "pixel_channel_avg_meta_cluster.csv"))
pixel_mat = read.csv(file.path(cell_output_dir, "cell_meta_cluster_channel_avg.csv"))

mat = as.matrix(pixel_mat[,7:17])
mat = scale(mat, center = TRUE, scale = TRUE)

pixel_cluster_colors = setNames(
    discrette_colors_50[1:length(unique(pixel_mat$cell_meta_cluster))], unique(pixel_mat$cell_meta_cluster)
)

library(ComplexHeatmap)
#mat = mat[ !(rownames(mat)  %in% "CD79a"), ]
ht_opt$raster_temp_image_max_width = 1000000
ht_opt$raster_temp_image_max_height = 1000000

h = Heatmap(
  t(mat),
  name = paste0("Pixel Clustering ", samp),
  column_title = "Pixels",
  row_title = "Proteins",
  row_dend_side = "right",
  row_dend_width = unit(4, "cm"),
  column_dend_height = unit(5, "cm"),
  show_column_names = FALSE,
  cluster_columns = FALSE,
  column_order = order(pixel_mat$cell_meta_cluster),
  clustering_distance_columns = "pearson",
  clustering_distance_rows =  "pearson",
  row_split = 10,
  # column_split = length(unique(pixel_mat$pixel_meta_cluster)),
  border = TRUE,
  row_names_gp = gpar(fontsize = 8),
  use_raster = FALSE,
  top_annotation = HeatmapAnnotation(
    fov = pixel_mat$fov,
    condition = pixel_mat$condition,
    metacluster = as.factor(pixel_mat$cell_meta_cluster),
    col = list(metacluster = pixel_cluster_colors,
    fov = setNames(
      unique(spe$sample_id_color), unique(spe$sample_id)
    ),
    condition = setNames(
      unique(spe$condition_color), unique(spe$condition)
    )
    
    ))
)

pdf(file.path(pixel_output_dir, paste0("Pixel_heatmap.pdf")),
    width = 20,
    height = 15)
draw(h)
dev.off()


png(file.path(pixel_output_dir, paste0("Pixel_heatmap.png")),
    width = 1500,
    height = 1000)
draw(h)
dev.off()


# Mapping 
mapping = setNames(c(6,7,10,21,23,24,25,26), c(
  "Epidermis Outer",
  "Blood Vessels Outer",
  "Lymphatic Vessels",
  "Blood Vessels Inner",
  "Dermis",
  "Epidermis Inner",
  "Background",
  "Immune"
))

colors_region = data.frame(
  region = 
    c(  "Epidermis Outer",
        "Epidermis Inner",
        "Blood Vessels Outer",
        "Lymphatic Vessels",
        "Blood Vessels Inner",
        "Dermis",
        "Immune",
        "Background"),
  region_color =
    c(
      "#930805ff",
      "#4a0100ff",
      "#C23EB5",
      "#fec911ff",
      "#78108F",
      "grey85",
      "#0d522aff",
      "white")
)

# Plot PDF
for(samp in unique(spe$sample_id)[1:10]){
  celltype_img_mat = getImageAsMatrix(file.path(pixel_output_dir, "pixel_masks", paste0(samp, "_pixel_mask.tiff")))
  cell_cluster_color_df = data.frame(cell_type = names(mapping),
                                     cell_type_color = colors_region$region_color[match(names(mapping), colors_region$region)])
  
  pdf(file.path(pixel_output_dir, "pixel_masks", paste0(samp,"_type.pdf")), height = 5, width = 6)
  rgb_color_cell_raster(img_mat = celltype_img_mat, color_df = cell_cluster_color_df, main = samp) 
  dev.off()
}

#################################################################################
#################################################################################

cluster_counts_size_norm = arrow::read_feather(file.path(cell_output_dir, "cluster_counts_size_norm.feather"))


mat = as.matrix(cluster_counts_size_norm[,3:22])
rownames(mat)

cluster_to_cell_type = setNames(c(1:25), c("T_helper", 
                                           "Macrophage",
                                           "Monocyte",
                                           "Basophil",
                                           "Macrophage2",
                                           "Epithelium-Corneum",
                                           "T_regulatory",
                                           "Mixed2",
                                           "T_helper2",
                                           "Epithelium-Spinosum",
                                           "Endothelial",
                                           "T_helper_Bcl2+",
                                           "pDC",
                                           "Endothelial2",
                                           "Epithelium-Basale",
                                           "T_cytotoxic",
                                           "pDC2",
                                           "B_cell",
                                           "Fibroblast",
                                           "Monocyte2",
                                           "CD123+_CD134+",
                                           "NK",
                                           "Lymphatic",
                                           "Fibroblast2"
                                           ))



cluster_counts_size_norm$cell_type = names(cluster_to_cell_type)[match(cluster_counts_size_norm$cell_meta_cluster_rename, cluster_to_cell_type)]
cluster_counts_size_norm$cell_id = paste0(cluster_counts_size_norm$fov, "-", cluster_counts_size_norm$segmentation_label)

spe = qs::qread("output/SpatialExperiment.qs")
spe$cell_type_pixel = cluster_counts_size_norm$cell_type[match(spe$cell_id,
                                                               cluster_counts_size_norm$cell_id)]

# spe = spe[,-which(is.na(spe$cell_type_pixel))]
set.seed(48)
cell_cluster_color_df = data.frame("cell_type_pixel" = sort(unique(spe$cell_type_pixel)),
                                   "cell_type_pixel_color" = sample(discrette_colors_50[1:25], length(sort(unique(spe$cell_type_pixel)))))
spe = colors_scExp(spe, annotCol = "cell_type_pixel", color_by = "cell_type_pixel", color_df = cell_cluster_color_df)
spe. = spe
spe.$cell_cluster = spe.$cell_type_pixel
plot_reduced_dim_scExp(spe., color_by = "cell_type_pixel", reduced_dim = "UMAP", size = 0.25,  downsample = 1e6, min_quantile = 0, max_quantile = 1, annotate_clusters = T)

par(mai = c(1,2,1,1))
tab = table(spe$cell_type_pixel)
barplot(tab,  col = cell_cluster_color_df$cell_type_pixel_color[match(names(tab), cell_cluster_color_df$cell_type_pixel)], horiz = T, las = 1)

for(i in unique(spe$sample_id)){
  cell_overlay_mat = getImageAsMatrix(spe, i, segment = "whole_cell")
  celltype_img_mat = get_metadata_image_overlay(spe, cell_overlay_mat, i, metadata = "cell_type_pixel")
  
  pdf(file.path(output_dir, "cell_masks", paste0(i, "_cell_type.pdf")), width = 6, height = 5)
  print(rgb_color_cell_raster(celltype_img_mat, cell_cluster_color_df, i, cex.main = 0.6, line = 0))
  dev.off()
}

# Cell Overlay of positive markers
# Plot 
for (i in unique(spe$sample_id)) {
  cell_overlay_mat = getImageAsMatrix(spe, i, segment = "whole_cell")
  celltype_img_mat = get_metadata_image_overlay(spe, cell_overlay_mat, i, metadata = "cell_type")
  tiff::writeTIFF(as.matrix(celltype_img_mat),
                  file.path(output_dir, "cell_masks",  paste0(i,"_cell_type.tiff")),
                  bits.per.sample = 32)
  
}

par(mai = c(1,2,1,4))
heatmap(table(spe$cell_type, spe$cell_type_pixel), scale = "row")



spe. = spe[,spe$sample_id=="ROI-15"]
table(spe.$cell_type, (labels_Christoph[match(spe.$cell_id, names(labels_Christoph))]))

confusion_matrix <- table(as.factor(spe.$cell_type), 
                          as.factor((labels_Christoph[match(spe.$cell_id, names(labels_Christoph))])))


# convert confusion matrix to data frame
confusion_df <- as.data.frame(confusion_matrix)

# rename columns
colnames(confusion_df) <- c("PositiveMarker.Class", "ManualAnnot.Class", "Count")

# create heatmap with text labels
png(file.path(output_dir, "Manual_vs_PositiveMarker.png"), height = 1200, width = 1200, res = 300)
ggplot(data = confusion_df, aes(x = PositiveMarker.Class, y = ManualAnnot.Class, fill = log10(Count+1), label = Count)) + 
  geom_tile() + 
  scale_fill_viridis_c(begin = 0.4, direction = -1) +
  labs(title = "Confusion Matrix Heatmap",
       x = "PositiveMarker Class", y = " ManualAnnot Class", fill = "Count") +
  geom_text(color = "white", size = 5)  + theme(axis.text.x = element_text(angle = 90))
dev.off()

# Composition analysis ---------------------------------------------------------
library(dittoSeq)
pdf(file.path(output_dir, "Composition_analysis.pdf"))
dittoBarPlot(spe, "condition", group.by = "sample_id", color.panel = unique(spe$condition_color[order(spe$condition)]), main = "Cell Count", scale = "count")  + Seurat::DarkTheme()
dittoBarPlot(spe, "cell_type_pixel", group.by = "condition", color.panel =  unique(spe$cell_type_pixel_color[order(spe$cell_type_pixel)]), main = "Condition x Cell Type")  + Seurat::DarkTheme()
dittoBarPlot(spe, "condition", group.by = "condition", scale = "count",  color.panel = unique(spe$condition_color[order(spe$condition)]))
dev.off()

png(file.path(output_dir, "Sample x Clusters x Condition.png"), width = 3500, height = 1600, res = 300)
dittoBarPlot(spe, "cell_cluster", group.by = "sample_id", color.panel = unique(spe$sample_id_color[order(spe$sample_id)]), main = "Sample x Clusters", 
             split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1) + Seurat::DarkTheme()
dev.off()

d = dittoBarPlot(spe, "cell_type_pixel", group.by = "sample_id", color.panel = unique(spe$cell_type_pixel_color)[order(spe$cell_type_pixel)],  main = "Sample x Clusters",
                 split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)

png(file.path(output_dir, "Sample x CellType x Condition.png"), width = 3500, height = 1600, res = 300)
colors =  setNames(unique(spe$cell_type_pixel_color), unique(spe$cell_type_pixel))
dittoBarPlot(spe, "cell_type_pixel", group.by = "sample_id", color.panel = colors[levels(d$data$label)],  main = "Sample x Clusters",
             split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)
dev.off()

spe. = spe[,!spe$cell_type_pixel %in% c("Epidermis_Mixed", "Epidermis", "Endothelial", "Fibroblast")]
d = dittoBarPlot(spe., "cell_type_pixel", group.by = "sample_id", color.panel = unique(spe.$cell_type_pixel_color[order(spe.$cell_type_pixel)]), main = "Sample x Clusters", 
                 split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)

png(file.path(output_dir, "Sample x CellType x Condition ImmuneCompartment.png"), width = 3500, height = 1600, res = 300)
colors =  setNames(unique(spe.$cell_type_pixel_color), unique(spe.$cell_type_pixel))
dittoBarPlot(spe., "cell_type_pixel", group.by = "sample_id", color.panel = colors[levels(d$data$label)],  main = "Sample x Clusters",
             split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)

dev.off()

spe. = spe[,spe$cell_type_pixel %in% c("T_regulatory", "T_helper", "T_cytotoxic")]
d = dittoBarPlot(spe., "cell_type_pixel", group.by = "sample_id", color.panel = unique(spe.$cell_type_pixel_color[order(spe.$cell_type_pixel)]), main = "Sample x Clusters", 
                 split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)

png(file.path(output_dir, "Sample x CellType x Condition ImmuneCompartment.png"), width = 3500, height = 1600, res = 300)
colors =  setNames(unique(spe.$cell_type_pixel_color), unique(spe.$cell_type_pixel))
dittoBarPlot(spe., "cell_type_pixel", group.by = "sample_id", color.panel = colors[levels(d$data$label)],  main = "Sample x Clusters",
             split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)

dev.off()
