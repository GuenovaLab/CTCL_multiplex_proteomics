# Finds cell clusters using cell average intensities
# authors: Pacome Prompsy
# contact: pacome.prompsy@chuv.ch
# Guenova Lab
# CHUV (Centre Hospitalier Universitaire Vaudois), Lausanne, Suisse

cat("Running Predicted label to masks... \n")
setwd("/mnt/RECHERCHE/GUENOVA_LAB/Project_Multiplex_Phenotyping/Miltenyi/Muliplex_Imaging_Pipeline/")

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

source("scripts/MultiplexImaging_utils.R")
source("scripts/GlobalVars.R")
spe = qs::qread("output/SpatialExperiment.qs")

################################################################################
# Positive Marker Detection (DeepLearning)
################################################################################

output_dir = file.path("output/CellSighter/marker/Predictions/")
mask_dir = file.path(output_dir, "Masks")
if(!dir.exists(mask_dir)) dir.create(mask_dir)


cell_markers = read.csv("annotation/cell_markers.csv")
marker_metadata = read.csv("annotation/marker_metadata.csv")
markers = marker_metadata$Marker[which(marker_metadata$PassOverallQuality == "True")]
sample_location = read.csv("annotation/Sample_metadata.csv", sep =";")
samples = sample_location$ROI[which(sample_location$KeepForAnalysis == TRUE)]

markers = setdiff(sort(markers), c("DAPI", "CXCR4", "CXCL12", "Ki67", "TRBC1"))

for(marker in markers) { 
  print(marker)
  
  for(image in samples){
    print(image)
    pred_file = file.path(output_dir, marker, paste0(image, "_val_results.csv"))
    pred_file_batch = file.path(output_dir, marker, ifelse(as.numeric(gsub(".*-", "",image)) <= 28, "batch_1", "batch_2"), paste0(image, "_val_results.csv"))
    
    if(file.exists(pred_file)){
      predictions = read.csv(pred_file)
      predictions = predictions[grep(image, predictions$image_id),]
      
      cell_overlay_mat = getImageAsMatrix(file.path("output", "segmentation", paste0(image,"_whole_cell.tiff")))
      boundary_overlay_mat = getImageAsMatrix(file.path("output", "segmentation", paste0(image,"_whole_cell_segmentation_borders.tiff")))
      
      spe. = spe[,grep(image, spe$sample_id)]
      colData(spe.)[marker] = 0
      colData(spe.)[marker] = predictions$pred[match(spe.$cell_id, paste0(image, "-", predictions$cell_id))]
      
      celltype_img_mat = get_metadata_image_overlay(spe = spe., img_mat = cell_overlay_mat,
                                                    sample = image, metadata = marker, levels = c(0, 1, 2))
      
      celltype_img_mat[celltype_img_mat == 1] = 0
      celltype_img_mat[boundary_overlay_mat == 0] = 0
      
      if(!dir.exists(file.path(mask_dir, image))) dir.create(file.path(mask_dir, image))
      
      tiff::writeTIFF(as.matrix(celltype_img_mat),
                      file.path(mask_dir, image, paste0(image,"_",marker,"_mask.tiff")),
                      bits.per.sample = 8)
      
    } else if(file.exists(pred_file_batch)){
      predictions = read.csv(pred_file_batch)
      predictions = predictions[grep(image, predictions$image_id),]
      
      cell_overlay_mat = getImageAsMatrix(file.path("output", "segmentation", paste0(image,"_whole_cell.tiff")))
      boundary_overlay_mat = getImageAsMatrix(file.path("output", "segmentation", paste0(image,"_whole_cell_segmentation_borders.tiff")))
      
      spe. = spe[,grep(image, spe$sample_id)]
      colData(spe.)[marker] = 0
      colData(spe.)[marker] = predictions$pred[match(spe.$cell_id, paste0(image, "-", predictions$cell_id))]
      
      celltype_img_mat = get_metadata_image_overlay(spe = spe., img_mat = cell_overlay_mat,
                                                    sample = image, metadata = marker, levels = c(0, 1, 2))
      
      celltype_img_mat[celltype_img_mat == 1] = 0
      celltype_img_mat[boundary_overlay_mat == 0] = 0
      
      if(!dir.exists(file.path(mask_dir, image))) dir.create(file.path(mask_dir, image))
      tiff::writeTIFF(as.matrix(celltype_img_mat),
                      file.path(mask_dir, image, paste0(image,"_",marker,"_mask.tiff")),
                      bits.per.sample = 32)
    }
  }
}


################################################################################
# Positive Marker Detection (Distribution-based)
################################################################################

output_dir = file.path("output/markers/")
if(!dir.exists(output_dir)) dir.create(output_dir)

mask_dir = file.path(output_dir, "Masks")
if(!dir.exists(mask_dir)) dir.create(mask_dir)

cell_markers = read.csv("annotation/cell_markers.csv")
marker_metadata = read.csv("annotation/marker_metadata.csv")
markers = marker_metadata$Marker[which(marker_metadata$PassOverallQuality == "True")]
sample_location = read.csv("annotation/Sample_metadata.csv", sep =";")
samples = sample_location$ROI[which(sample_location$KeepForAnalysis == TRUE)]

markers = setdiff(sort(markers), c("DAPI", "CXCR4", "CXCL12", "Ki67", "TRBC1"))

# markers =  marker_metadata$Marker[marker_metadata$BackgroundInKeratynocyteNotExpected == "True"]
samples = c("ROI-11", "ROI-42")

for(marker in markers) { 
  print(marker)
  for(image in samples){
    
    cell_overlay_mat = getImageAsMatrix(file.path("output", "segmentation", paste0(image,"_whole_cell.tiff")))
    boundary_overlay_mat = getImageAsMatrix(file.path("output", "segmentation", paste0(image,"_whole_cell_segmentation_borders.tiff")))
    
    
    
    spe. = spe[,grep(image, spe$sample_id)]
    colData(spe.)[marker] = 0
    colData(spe.)[marker] = spe.@assays@data$positive_marker[marker,]
    
    celltype_img_mat = get_metadata_image_overlay(spe = spe., img_mat = cell_overlay_mat,
                                                  sample = image, metadata = marker, levels = c(0, 1, 2))
    
    celltype_img_mat[celltype_img_mat == 1] = 0
    celltype_img_mat[boundary_overlay_mat == 0] = 0
    
    if(!dir.exists(file.path(mask_dir, image))) dir.create(file.path(mask_dir, image))
    
    tiff::writeTIFF(as.matrix(celltype_img_mat),
                    file.path(mask_dir, image, paste0(image,"_",marker,"_mask.tiff")),
                    bits.per.sample = 8)
  }
}



################################################################################
# Cell Phenotyping (DeepLearning)
################################################################################
output_dir = file.path("output/CellSighter/celltype/Predictions/")
mask_dir = file.path(output_dir, "Masks")
if(!dir.exists(mask_dir)) dir.create(mask_dir)

images = c("ROI-06","ROI-54", "ROI-17")
images = unique(spe$sample_id)
images = unique(spe$sample_id[spe$condition == "AD"])

for(image in images){
  print(image)
  
  cell_overlay_mat = getImageAsMatrix(file.path("output", "segmentation", paste0(image,"_whole_cell.tiff")))
  boundary_overlay_mat = getImageAsMatrix(file.path("output", "segmentation", paste0(image,"_whole_cell_segmentation_borders.tiff")))
  
  spe. = spe[,grep(image, spe$sample_id)]
  
  colData(spe.)["celltype"] = 0
  colData(spe.)[, "celltype"] = as.numeric(factor(spe.$cell_type_CellSighter, levels = celltype_levels))
  colData(spe.)[, "celltype"] = floor(changeRange(colData(spe.)[, "celltype"], newmin = 1, newmax = 255))
  
  celltype_img_mat = get_metadata_image_overlay(spe = spe., img_mat = cell_overlay_mat,
                                                sample = image, metadata = "celltype", levels = seq(0,255))
  celltype_img_mat@x = changeRange(celltype_img_mat@x, newmin = 0, newmax = 1)
  
  # celltype_img_mat[celltype_img_mat == 1] = 0
  # celltype_img_mat[boundary_overlay_mat == 0] = 0
  
  if(!dir.exists(file.path(mask_dir, image))) dir.create(file.path(mask_dir, image))
  
  tiff::writeTIFF(as.matrix(celltype_img_mat),
                  file.path(mask_dir, image, paste0(image,"_celltype_mask.tiff")),
                  bits.per.sample = 8)
  
  # color_df = data.frame(cell_type = c("Rest","T_cytotoxic", "NKT"), cell_type_color = c("grey", "#D72638", "#8A6C0B"))
  # pdf(file.path(mask_dir, image, paste0(image,"_", "cytotoxic","_mask.pdf")))
  # rgb_color_cell_raster(img_mat = celltype_img_mat, color_df = color_df,
  #                       main = "Cytotoxic")
  # dev.off()
}



################################################################################
# Unique Cell Snapshots
################################################################################
output_dir = file.path("output/Snapshots/CellSnapshots/")
if(!dir.exists(output_dir)) dir.create(output_dir)

cells=c("ROI-06-4625", "ROI-06-4162", "ROI-06-4677")

library(stringr)
for(cell in cells) { 
  
  image = stringr::str_replace(cell, "(?<=ROI-..).*", "")
  cell_overlay_mat = getImageAsMatrix(file.path("output", "segmentation", paste0(image,"_whole_cell.tiff")))
  boundary_overlay_mat = getImageAsMatrix(file.path("output", "segmentation", paste0(image,"_whole_cell_segmentation_borders.tiff")))
  
  spe. = spe[,grep(image, spe$sample_id)]
  colData(spe.)["cell"] = 0
  colData(spe.)[which(spe.$cell_id == cell), "cell"] = 1
  
  celltype_img_mat = get_metadata_image_overlay(spe = spe., img_mat = cell_overlay_mat,
                                                sample = image, metadata = "cell", levels = c(0, 1))
  
  celltype_img_mat[celltype_img_mat == 1] = 0
  celltype_img_mat[boundary_overlay_mat == 0] = 0
  
  tiff::writeTIFF(as.matrix(celltype_img_mat),
                  file.path(output_dir, paste0(cell,".tiff")),
                  bits.per.sample = 8)
  
}





################################################################################
# Cell Interactions 
################################################################################
output_dir = file.path("output/cell_interactions/enriched_interactions/")
mask_dir = file.path(output_dir, "Masks")
if(!dir.exists(mask_dir)) dir.create(mask_dir)

images = c("ROI-06","ROI-54", "ROI-17")
images = unique(spe$sample_id)
images = unique(spe$sample_id[spe$condition == "AD"])

# for(celltype in celltypes) { 
#   print(celltype)
#   celltype = c("NKT", "T_cytotoxic")
# 
celltype1 = "T_helper" 
celltype2 =  "T_cytotoxic"
mark =  "CLA"
matrix_list = qs::qread(file.path(output_dir, "matrix_list.qs"))


for(image in images){
  print(image)
  
  cell_overlay_mat = getImageAsMatrix(file.path("output", "segmentation", paste0(image,"_whole_cell.tiff")))
  boundary_overlay_mat = getImageAsMatrix(file.path("output", "segmentation", paste0(image,"_whole_cell_segmentation_borders.tiff")))
  
  spe. = spe[,grep(image, spe$sample_id)]
  
  # Enriched interaction mat
  mat = matrix_list[[mark]]
  mat = mat[match(spe.$cell_id, rownames(mat)), match(spe.$cell_id, rownames(mat))]
  mat = mat[match(spe.$cell_id[spe.$cell_type_CellSighter == celltype1], rownames(mat)),
            match(spe.$cell_id[spe.$cell_type_CellSighter == celltype2],  colnames(mat))]
  
  # Celltype1 
  cells_1 = rownames(mat)[which(rowSums(mat) > 0)]
  
  colData(spe.)["enriched"] = 0
  colData(spe.)[match(cells_1, spe.$cell_id), "enriched"] = 1
  celltype_img_mat = get_metadata_image_overlay(spe = spe., img_mat = cell_overlay_mat,
                                                sample = image, metadata = "enriched", levels = c(0, 1, 2))
  
  celltype_img_mat[celltype_img_mat == 1] = 0
  celltype_img_mat[boundary_overlay_mat == 0] = 0
  
  if(!dir.exists(file.path(mask_dir, image))) dir.create(file.path(mask_dir, image))
  tiff::writeTIFF(as.matrix(celltype_img_mat),
                  file.path(mask_dir, image, paste0(image,"_", celltype1,"_mask.tiff")),
                  bits.per.sample = 8)
  
  # Celltype2 
  cells_2 = colnames(mat)[which(colSums(mat) > 0)]
  
  colData(spe.)["enriched"] = 0
  colData(spe.)[match(cells_2, spe.$cell_id), "enriched"] = 1
  celltype_img_mat = get_metadata_image_overlay(spe = spe., img_mat = cell_overlay_mat,
                                                sample = image, metadata = "enriched", levels = c(0, 1, 2))
  
  celltype_img_mat[celltype_img_mat == 1] = 0
  celltype_img_mat[boundary_overlay_mat == 0] = 0
  
  if(!dir.exists(file.path(mask_dir, image))) dir.create(file.path(mask_dir, image))
  tiff::writeTIFF(as.matrix(celltype_img_mat),
                  file.path(mask_dir, image, paste0(image,"_", celltype2,"_mask.tiff")),
                  bits.per.sample = 8)
  
  
  
  # color_df = data.frame(cell_type = c("Rest","T_cytotoxic", "NKT"), cell_type_color = c("grey", "#D72638", "#8A6C0B"))
  # pdf(file.path(mask_dir, image, paste0(image,"_", "cytotoxic","_mask.pdf")))
  # rgb_color_cell_raster(img_mat = celltype_img_mat, color_df = color_df,
  #                       main = "Cytotoxic")
  # dev.off()
  
}




################################################################################
# CD3+ Th vs CD3- Th
################################################################################
output_dir = file.path("output/Snapshots/CD3+_vs_CD3-_Th/")
if(!dir.exists(output_dir)) dir.create(output_dir)

for(image in unique(spe$sample_id[spe$disease != "Malign"])) { 
  cell_overlay_mat = getImageAsMatrix(file.path("output", "segmentation", paste0(image,"_whole_cell.tiff")))
  
  spe. = spe[,grep(image, spe$sample_id)]
  spe. = spe.[,spe.$cell_type_CellSighter == "T_helper"]
  
  colData(spe.)["CD3"] = spe.@assays@data$CellSighter_marker_mat["CD3",] 
  
  celltype_img_mat = get_metadata_image_overlay(spe = spe., img_mat = cell_overlay_mat,
                                                sample = image, metadata = "CD3", levels = c(0, 1))
  tiff::writeTIFF(as.matrix(celltype_img_mat),
                  file.path(output_dir, paste0(image,".tiff")),
                  bits.per.sample = 8)
  
  
  
}
