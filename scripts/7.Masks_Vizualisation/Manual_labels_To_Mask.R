# Finds cell clusters using cell average intensities
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
              "ChromSCape",
              "SpatialExperiment",
              "Matrix",
              "Seurat")
suppressPackageStartupMessages(invisible(lapply(libraries, require, character.only = TRUE)))
setwd("/media/localadmin/T7/CHUV/Documents/Muliplex_Imaging_Pipeline/")

source("scripts/MultiplexImaging_utils.R")
source("scripts/GlobalVars.R")
spe = qs::qread("output/SpatialExperiment.qs")

markers = rownames(spe)

output_dir = file.path("output/CellSighter/marker/marker_classification/")
mask_dir = file.path(output_dir, "Masks")
if(!dir.exists(mask_dir)) dir.create(mask_dir)

library(reticulate)
np = reticulate::import("numpy")


for(marker in rownames(spe)){
  print(marker)
  for(image in c("ROI-01-Pacome","ROI-02-Pacome","ROI-10-Pacome","ROI-11-Pacome","ROI-20-Pacome","ROI-21-Pacome")){
    print(image)
    
    cell2label_dir = file.path(output_dir, marker, "CellTypes/cells2labels/")
    labels <- np$load(file.path(cell2label_dir, paste0(image, ".npz")), allow_pickle = T)["data"]
    image = gsub("-Pacome|-Christoph|Ionoss","",image)
    names(labels) = paste0(image, "-", 0:(length(labels)-1))
    labels = labels[names(labels) %in% spe$cell_id]
    
    cell_overlay_mat = getImageAsMatrix(file.path("output", "segmentation", paste0(image,"_whole_cell.tiff")))
    
    spe. = spe[,grep(image, spe$sample_id)]
    colData(spe.)[marker] = 0
    colData(spe.)[marker] = labels[match(spe.$cell_id, names(labels))]
    
    celltype_img_mat = get_metadata_image_overlay(spe., cell_overlay_mat, image, metadata = marker)
    
    if(!dir.exists(file.path(mask_dir, marker))) dir.create(file.path(mask_dir, marker))
    tiff::writeTIFF(as.matrix(celltype_img_mat),
                    file.path(mask_dir, marker, paste0(image,"_",marker,"_mask.tiff")),
                    bits.per.sample = 32)
    
  }
}





