# Finds cell clusters using cell average intensities
# authors: Pacome Prompsy
# contact: pacome.prompsy@chuv.ch
# Guenova Lab
# CHUV (Centre Hospitalier Universitaire Vaudois), Lausanne, Suisse

cat("Running Numerical Values to to masks... \n")
setwd("/mnt/RECHERCHE/GUENOVA_LAB/Project_Multiplex_Phenotyping/Miltenyi/Muliplex_Imaging_Pipeline/")

# Loading packages --------------------------------------------------------
libraries = c("argparse",
              "ggplot2",
              "dplyr",
              "tidyr",
              "ChromSCape",
              "SpatialExperiment",
              "Matrix",
              "Seurat",
              "detrendr")
suppressPackageStartupMessages(invisible(lapply(libraries, require, character.only = TRUE)))

source("scripts/MultiplexImaging_utils.R")
source("scripts/GlobalVars.R")
spe = qs::qread("output/SpatialExperiment.qs")

output_dir = file.path("output", "Snapshots", "RawIntensity_Masks")
dir.create(output_dir)

# Values to mask
markers = setdiff(rownames(spe), c("DAPI", "CXCR4", "CXCL12", "Ki67", "TRBC1"))
samples = c("ROI-16")

for(marker in markers) { 
  print(marker)
  for(image in samples){
    
    cell_overlay_mat = getImageAsMatrix(file.path("output", "segmentation", paste0(image,"_whole_cell.tiff")))
    
    spe. = spe[,grep(image, spe$sample_id)]
    colData(spe.)[marker] = 0
    # colData(spe.)[marker] = round(changeRange(spe.@assays@data$counts[marker,], , 0(255^2 -1)))
    colData(spe.)[marker] = round(spe.@assays@data$counts[marker,])
    
    celltype_img_mat = get_numeric_image_overlay(spe = spe., img_mat = cell_overlay_mat,
                                                  sample = image, metadata = marker)
    
    if(!dir.exists(file.path(output_dir, image))) dir.create(file.path(output_dir, image))
    
    # tiff::writeTIFF(as.matrix(celltype_img_mat),
    #                 compression = "none", reduce = FALSE,
    #                 file.path(output_dir, image, paste0(image,"_",marker,"_mask.tiff")),
    #                 bits.per.sample = 16)
    
    ijtiff::write_tif(img = as.matrix(celltype_img_mat),
                      bits_per_sample = 16, compression = "none",
                      overwrite = TRUE, 
                      path = file.path(output_dir, image, paste0(image,"_",marker,"_mask.tiff"))
    )
    
    hist(celltype_img_mat@x, breaks = 150)
  }
}


