# Finds positive markers in cell detections by distribution of the average intensities
# authors: Pacome Prompsy
# contact: pacome.prompsy@chuv.ch
# Guenova Lab
# CHUV (Centre Hospitalier Universitaire Vaudois), Lausanne, Suisse

cat("Running Positive Marker Detection based on average cell expression distribution... \n")

# Loading packages -------------------------------------------------------------
libraries = c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "ggspavis",
  "SpatialExperiment",
  "ComplexHeatmap",
  "seriation",
  "EBImage",
  "ComplexHeatmap"
)
suppressPackageStartupMessages(invisible(lapply(libraries, require, character.only = TRUE)))

setwd("/mnt/RECHERCHE/GUENOVA_LAB/Project_Multiplex_Phenotyping/Miltenyi/Muliplex_Imaging_Pipeline/")

source("scripts/MultiplexImaging_utils.R")
source("scripts/GlobalVars.R")
devtools::load_all("annotation/silvermantest/")
marker_metadata = read.csv("annotation/marker_metadata.csv")

# Directorires --------------------------------------------
output = "./output/" 

output_dir = file.path(output, "markers")
if (!dir.exists(output_dir))
  dir.create(output_dir)


# Loading the dataset ----------------------------------------------------------
spe = qs::qread(file.path(output, "SpatialExperiment.qs"))

# Finding positive markers -----------------------------------------------------
cat("Finding positive markers...\n")

pdf(file.path(output_dir, "Histogram_positive_marker_findCutoffs.pdf"))
marker_cutoffs = cell_marker_identification(spe)
dev.off()

write.csv(marker_cutoffs, file.path(output_dir, "marker_cutoffs.csv"), row.names = T)


#Cyto per sample
positive_marker_mat_cyto = normcounts(spe)
for (marker in names(marker_cutoffs) ) {
  positive_marker_mat_cyto[marker, ] =
    ifelse(positive_marker_mat_cyto[marker, ] > marker_cutoffs[marker], 1, 0)
}
rownames(positive_marker_mat_cyto) = rownames(spe)
colnames(positive_marker_mat_cyto) = colnames(spe)
spe@assays@data$positive_marker = positive_marker_mat_cyto

# Cell Overlay of positive markers
spe@int_metadata$segmentation_dir = "output/segmentation/"

pdf(file.path(output_dir, "postivive_marker_overlay.pdf"))
for(i in rownames(spe)){
  print(i)
  print(plotSPE(
    spe,
    feature = i,
    assay = "positive_marker", 
    x_coord = "centroid.0",
    y_coord = "centroid.1",
    sample_id = "sample_id", size = 0.25
  )
  )
}
dev.off()

# Plot Tiffs
for (samp in "ROI-19") {
  print(samp)
  sampdir = file.path(output_dir, samp)
  if (!dir.exists(sampdir))
    dir.create(sampdir)
  cell_overlay_mat = getImageAsMatrix(file.path(spe@int_metadata$segmentation_dir, paste0(samp,"_whole_cell.tiff")))
  for (i in rownames(spe)) {
    print(i)
    positive_marker_img_mat = get_metadata_image_overlay(spe = spe, img_mat = cell_overlay_mat,
                                                         sample = samp, metadata = i, levels = c(1, 2))
    tiff::writeTIFF(positive_marker_img_mat,
                    file.path(sampdir, paste0(i, ".tiff")),
                    bits.per.sample = 8)
  }
}


# Plotting Heatmaps -----------------------------------------------------------

cat("Plotting heatmaps All Markers...\n")
set.seed(48)
idx = sample(ncol(spe), size = 6000, replace = F)
annot = as.data.frame(colData(spe)[idx, ])

for (assay in c("normcounts", "counts", "positive_marker")) {
  if (assay == "normcounts")
    mat <- (normcounts(spe))
  if (assay == "counts")
    mat <- (counts(spe))
  if (assay == "positive_marker")
    mat <- (spe@assays@data$positive_marker)
  
  mat = mat[, idx]
  if (assay != "positive_marker")
    mat = scale(mat, center = TRUE, scale = TRUE)
  #mat = mat[ !(rownames(mat)  %in% "CD79a"), ]
  pdf(file.path(output_dir, paste0("Marker_heatmap_", assay, ".pdf")),
      width = 15,
      height = 10)
  h = Heatmap(
    mat,
    name = "Markers",
    column_title = "Cells",
    row_title = "Proteins",
    row_dend_side = "right",
    row_dend_width = unit(4, "cm"),
    column_dend_height = unit(5, "cm"),
    show_column_names = FALSE,
    clustering_distance_columns = ifelse(assay == "positive_marker", "binary", "pearson"),
    clustering_distance_rows = ifelse(assay == "positive_marker", "binary", "pearson"),
    row_split = 5,
    column_split = 5,
    border = TRUE,
    row_names_gp = gpar(fontsize = 8),
    use_raster = FALSE,
    top_annotation = HeatmapAnnotation(
      sample_id = annot$sample_id,
      col = list(sample_id = setNames(
        unique(spe$sample_id_color), unique(spe$sample_id)
      ))
    )
  )
  draw(h)
  dev.off()
}

cat("Plotting heatmaps Cell Markers...\n")
set.seed(48)
idx = sample(ncol(spe), size = 4000, replace = F)
annot = as.data.frame(colData(spe)[idx, ])
markers = read.csv("annotation/cell_markers.csv")
markers = unique(markers$marker)

for (assay in c("normcounts", "counts", "positive_marker")) {
  if (assay == "normcounts")
    mat <- (normcounts(spe))
  if (assay == "counts")
    mat <- (counts(spe))
  if (assay == "positive_marker")
    mat <- (spe@assays@data$positive_marker)
  
  mat = mat[, idx]
  mat = mat[(rownames(mat)  %in% markers), ]
  
  if (assay != "positive_marker")
    mat = scale(mat, center = TRUE, scale = TRUE)
  
  pdf(file.path(output_dir, paste0("Marker_heatmap_cell_markers_", assay, ".pdf")),
      width = 15,
      height = 10)
  h = Heatmap(
    mat,
    name = "Markers",
    column_title = "Cells",
    row_title = "Proteins",
    row_dend_side = "right",
    row_dend_width = unit(4, "cm"),
    column_dend_height = unit(5, "cm"),
    show_column_names = FALSE,
    clustering_distance_columns = ifelse(assay == "positive_marker", "binary", "pearson"),
    clustering_distance_rows = ifelse(assay == "positive_marker", "binary", "pearson"),
    row_split = 5,
    column_split = 5,
    border = TRUE,
    row_names_gp = gpar(fontsize = 8),
    use_raster = F,
    top_annotation = HeatmapAnnotation(
      sample_id = annot$sample_id,
      col = list(sample_id = setNames(
        unique(spe$sample_id_color), unique(spe$sample_id)
      ))
    )
  )
  draw(h)
  dev.off()
}
  



