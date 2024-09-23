# Analyze Protein Expression across samples
# authors: Pacome Prompsy
# contact: pacome.prompsy@chuv.ch
# Guenova Lab
# CHUV (Centre Hospitalier Universitaire Vaudois), Lausanne, Suisse

cat("Analyzing protein expression across samples, condition and celltypes... \n")

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
setwd("/media/localadmin/T7/CHUV/Documents/Muliplex_Imaging_Pipeline/")

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

pdf(file.path(output_dir, paste0("Celltype_Protein_CoExpression_heatmap.pdf")),
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