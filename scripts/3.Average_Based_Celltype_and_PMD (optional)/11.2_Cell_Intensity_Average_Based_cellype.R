# Finds cell type using marker positivity
# authors: Pacome Prompsy
# contact: pacome.prompsy@chuv.ch
# Guenova Lab
# CHUV (Centre Hospitalier Universitaire Vaudois), Lausanne, Suisse

cat("Running celltype classification using average-expression based PMD... \n")

# Loading packages --------------------------------------------------------
libraries = c("argparse",
              "ggplot2",
              "dplyr",
              "tidyr",
              "Seurat",
              "SpatialExperiment")
suppressPackageStartupMessages(invisible(lapply(libraries, require, character.only = TRUE)))
setwd("/mnt/RECHERCHE/GUENOVA_LAB/Project_Multiplex_Phenotyping/Miltenyi/Muliplex_Imaging_Pipeline/")

source("scripts/MultiplexImaging_utils.R")
source("scripts/GlobalVars.R")
devtools::load_all("annotation/silvermantest/")
devtools::load_all("/media/localadmin//T7/InstitutCurie/Documents/GitLab/ChromSCape/")

# Directorires --------------------------------------------
output = "./output/" 

output_dir = file.path(output, "markers")
if (!dir.exists(output_dir))
  dir.create(output_dir)

# Loading the dataset -----------------------------------------------------
spe = qs::qread(file.path(output, "SpatialExperiment.qs"))

# Reading the markers ----------------------------------------------
cell_markers = read.csv("annotation/cell_markers.csv")
cell_markers = cell_markers[-which(cell_markers$marker %in% c("CollagenI", "DAPI")),]
cell_markers = cell_markers %>% filter(important == TRUE)
cell_markers = cell_markers[nrow(cell_markers):1,]

cell_type = rep("Fibroblast", ncol(spe))
pos_mark = spe@assays@data$CellSighter_marker_mat

for(type in unique(cell_markers$cell_type)){
  df = cell_markers[cell_markers$cell_type == type,]
  vec_cell_type = setNames(rep(TRUE, ncol(spe)), colnames(spe))
  for(i in 1:nrow(df)){
    if(df$positive[i]) vec_cell_type = vec_cell_type & (pos_mark[df$marker[i],] == 1)
    else vec_cell_type = vec_cell_type & (pos_mark[df$marker[i],] == 0)
  }
  cell_type[vec_cell_type] = type
  print(table(cell_type))
}

spe$cell_type_positive_marker = cell_type
colnames(cell_type_color_df) = c("cell_type_positive_marker", "cell_type_positive_marker_color")
spe = colors_scExp(spe, annotCol = "cell_type_positive_marker", color_by = "cell_type_positive_marker", color_df = cell_type_color_df)

par(mai = c(1,2,1,1))
tab = table(spe$cell_type_positive_marker)
tab = tab[c("Fibroblast","Leukocyte", "B_cell", "T_helper", "T_cytotoxic", "T_regulatory", "pDC", "Basophil",  "NKT", "APC", "Monocytes", "Monocytic_Lineage", "Neutrophils", "Macrophages", "Keratinocyte", "Endothelial", "Lymphatic")]
barplot(tab,  col = cell_type_color_df$cell_type_positive_marker_color[match(names(tab), cell_type_color_df$cell_type_positive_marker)], horiz = T, las = 1)

qs::qsave(spe, "output/SpatialExperiment.qs")


# Cell Overlay of positive markers
# Plot 
for (samp in unique(spe$sample_id)) {
  cell_overlay_mat = getImageAsMatrix(file.path(spe@int_metadata$segmentation_dir, paste0(samp,"_whole_cell.tiff")))
  celltype_img_mat = get_metadata_image_overlay(spe, cell_overlay_mat, samp, metadata = "cell_type_positive_marker")
  tiff::writeTIFF(as.matrix(celltype_img_mat),
                  file.path(output_dir, paste0(samp,"_cell_type.tiff")),
                  bits.per.sample = 32)
  
}

# Plot PNG
for (samp in unique(spe$sample_id)) {
  cell_overlay_mat = getImageAsMatrix(file.path(spe@int_metadata$segmentation_dir, paste0(samp,"_whole_cell.tiff")))
  celltype_img_mat = get_metadata_image_overlay(spe, cell_overlay_mat, samp, metadata = "cell_type_positive_marker")
  
  pdf(file.path(output_dir, paste0(samp,"_cell_type_positive_marker.pdf")), height = 5, width = 6)
  rgb_color_cell_raster(celltype_img_mat, cell_type_color_df, samp) 
  dev.off()
}

# Composition analysis ---------------------------------------------------------
library(dittoSeq)
pdf(file.path(output_dir, "Composition_analysis.pdf"))
dittoBarPlot(spe, "condition", group.by = "sample_id", color.panel = unique(spe$condition_color[order(spe$condition)]), main = "Cell Count", scale = "count")  + Seurat::DarkTheme()
dittoBarPlot(spe, "cell_type_positive_marker", group.by = "condition", color.panel =  unique(spe$cell_type_positive_marker_color[order(spe$cell_type_positive_marker)]), main = "Condition x Cell Type")  + Seurat::DarkTheme()
dittoBarPlot(spe, "condition", group.by = "condition", scale = "count",  color.panel = unique(spe$condition_color[order(spe$condition)]))
dev.off()

d = dittoBarPlot(spe, "cell_type_positive_marker", group.by = "sample_id", color.panel = unique(spe$cell_type_positive_marker_color)[order(spe$cell_type_positive_marker)],  main = "Sample x Clusters",
                 split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)

png(file.path(output_dir, "Sample x CellType x Condition.png"), width = 3500, height = 1600, res = 300)
colors =  setNames(unique(spe$cell_type_positive_marker_color), unique(spe$cell_type_positive_marker))
dittoBarPlot(spe, "cell_type_positive_marker", group.by = "sample_id", color.panel = colors[levels(d$data$label)],  main = "Sample x Clusters",
             split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)
dev.off()

spe. = spe[,!spe$cell_type_positive_marker %in% c("Keratinocyte", "Endothelial", "Fibroblast", "Lymphatic")]
d = dittoBarPlot(spe., "cell_type_positive_marker", group.by = "sample_id", color.panel = unique(spe.$cell_type_positive_marker_color[order(spe.$cell_type_positive_marker)]), main = "Sample x Clusters", 
             split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)

png(file.path(output_dir, "Sample x CellType x Condition ImmuneCompartment.png"), width = 3500, height = 1600, res = 300)
colors =  setNames(unique(spe.$cell_type_positive_marker_color), unique(spe.$cell_type_positive_marker))
dittoBarPlot(spe., "cell_type_positive_marker", group.by = "sample_id", color.panel = colors[levels(d$data$label)],  main = "Sample x Clusters",
             split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)

dev.off()

spe. = spe[,spe$cell_type_positive_marker %in% c("T_regulatory", "T_helper", "T_cytotoxic")]
d = dittoBarPlot(spe., "cell_type_positive_marker", group.by = "sample_id", color.panel = unique(spe.$cell_type_positive_marker_color[order(spe.$cell_type_positive_marker)]), main = "Sample x Clusters", 
                 split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)

png(file.path(output_dir, "Sample x CellType x Condition ImmuneCompartment.png"), width = 3500, height = 1600, res = 300)
colors =  setNames(unique(spe.$cell_type_positive_marker_color), unique(spe.$cell_type_positive_marker))
dittoBarPlot(spe., "cell_type_positive_marker", group.by = "sample_id", color.panel = colors[levels(d$data$label)],  main = "Sample x Clusters",
             split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)

dev.off()


intensities = as.data.frame(t(spe@assays@data$positive_marker))
intensities$sample_id = spe$sample_id
intensities$condition = spe$condition
intensities$cell_type_positive_marker = spe$cell_type_positive_marker
intensities$cell_id = rownames(intensities)

intensities = intensities %>% tidyr::pivot_longer(-c(sample_id,condition, cell_type_positive_marker, cell_id), names_to = c("marker"))

for(celltype in c("T_helper", "T_cytotoxic", "T_regulatory", "Keratinocyte")){
  colors = setNames(unique(spe$condition_color), unique(spe$condition))
  colors_samp = setNames(unique(spe$sample_id_color), unique(spe$sample_id))
  pdf(file.path(output_dir, paste0("Functional_Marker_in_cell_type_positive_marker_",celltype,".pdf")), height = 20, width = 17)
   p =  intensities %>% filter(cell_type_positive_marker %in% c(celltype)) %>%
    group_by(sample_id, condition, marker) %>% 
    summarise(Percent = 100 * sum(value) / n()) %>%
    ggplot(aes(y = Percent, fill = condition, x = condition)) + geom_boxplot(outlier.alpha = 0) +
     geom_jitter(aes(color = sample_id)) + 
     theme_classic() +
    theme(axis.text.x = element_text(angle = 90)) +
     scale_color_manual(values = colors_samp) +
     scale_fill_manual(values = colors) +
    facet_wrap(~marker, ncol = 5) + ylab(paste0("% of ", celltype , " cells positive for marker")) + xlab("") +
     ggtitle(celltype)
  print(p)
  dev.off()
}

