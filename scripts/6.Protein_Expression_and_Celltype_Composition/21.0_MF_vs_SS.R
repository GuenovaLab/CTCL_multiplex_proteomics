cat("Analyzing MF vs SS skins... \n")

# Loading packages --------------------------------------------------------
libraries = c("argparse",
              "ggplot2",
              "dplyr",
              "tidyr",
              "Seurat",
              "SpatialExperiment",
              "arrow",
              "ChromSCape", 
              "ggpubr",
              "ggnewscale")
suppressPackageStartupMessages(invisible(lapply(libraries, require, character.only = TRUE)))
setwd("/mnt/RECHERCHE/GUENOVA_LAB/Project_Multiplex_Phenotyping/Miltenyi/Muliplex_Imaging_Pipeline/")

source("scripts/MultiplexImaging_utils.R")
source("scripts/GlobalVars.R")
library(ggplot2)
library(tidyverse)
library(ComplexHeatmap)
library(ggbeeswarm)

# Reading in data  -------------------------------------------------------------
output = file.path("output", "MF_vs_SS")
dir.create(file.path(output))
spe = qs::qread("output/SpatialExperiment.qs")


# Select MF and SS
spe_tumor = spe[, spe$condition %in% c("MF", "SS")]

################################################################################
# Composition
################################################################################
dir.create(file.path(output, "Composition"))
for(i in unique(spe_tumor$cell_type_CellSighter)){
  spe_tumor. = spe_tumor
  name = "all"
  if(!(i %in% struct_celltype)){
    spe_tumor. = spe_tumor[, !spe_tumor$cell_type_CellSighter %in% struct_celltype]
    name = "immune"
  }
  spe_tumor.$cell_type_CellSighter[spe_tumor.$cell_type_CellSighter != i] = "Aa"
  meta = as.data.frame(colData(spe_tumor.))
  
  meta = meta %>% 
    dplyr::group_by(sample_id, condition, batch) %>% 
    dplyr::summarise(percent_CN = 100 * length(
      which(cell_type_CellSighter == i)) / length(cell_type_CellSighter))
  
  png(file.path(output, "Composition", paste0("Boxplot_composition_",i,".png")), width = 1200, height = 1150, res = 300)
  p = grouped_dotplot(meta, y = "percent_CN", categ1 = "condition", categ2 = NULL,
                      ref.group = "MF", add_violin = F,
                      color_categ1 = setNames(unique(spe$condition_color), unique(spe$condition)))
  p = p + xlab("") + ylab(paste0("Percentage of ", gsub("_"," ",i), " (% of immune cells)")) +
    guides(fill="none", color = "none")
  print(p)
  dev.off()
  
}

################################################################################
# Plotting Functional Markers in cell types
################################################################################
intensities = as.data.frame(t(spe_tumor@assays@data$CellSighter_marker_mat))
intensities$sample_id = spe_tumor$sample_id
intensities$condition = spe_tumor$condition
intensities$batch = spe_tumor$batch
intensities$cell_type_CellSighter = spe_tumor$cell_type_CellSighter
intensities = intensities %>% group_by(sample_id, condition, batch, cell_type_CellSighter) %>% 
  summarise(across(Actin:Ki67, mean))

dir.create(file.path(output, "Markers"))
for(celltype in unique(spe_tumor$cell_type_CellSighter)){
  print(celltype)
  
  pdf(file.path(output,  "Markers", paste0("Functional_Marker_in_cell_type_CellSighter_",celltype,".pdf")), height = 3, width = 3)
  for(mark in colnames(intensities)[5:65]){
    intensities. = intensities %>% filter(cell_type_CellSighter %in% celltype)
    intensities.[[mark]] = 100 * intensities.[[mark]]
    
    p = grouped_dotplot(intensities., y =  mark, categ1 = "condition", categ2 = NULL,
                        ref.group = "MF", add_violin = T,
                        color_categ1 = setNames(unique(spe$condition_color), unique(spe$condition)))
    p = p + xlab("") + ylab(paste0("", celltype, " expressing ", mark, " (%)")) +
      guides(fill="none", color = "none")
    print(p)
  }
  dev.off()
}


