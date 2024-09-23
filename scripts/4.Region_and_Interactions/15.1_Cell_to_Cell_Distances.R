# Calculate distances between each and every cell
# authors: Pacome Prompsy
# contact: pacome.prompsy@chuv.ch
# Guenova Lab
# CHUV (Centre Hospitalier Universitaire Vaudois), Lausanne, Suisse

# Loading packages --------------------------------------------------------
libraries = c("argparse",
              "ggplot2",
              "plyr",
              "dplyr",
              "tidyr",
              "Seurat",
              "SpatialExperiment")
suppressPackageStartupMessages(invisible(lapply(libraries, require, character.only = TRUE)))
setwd("/mnt/RECHERCHE/GUENOVA_LAB/Project_Multiplex_Phenotyping/Miltenyi/Muliplex_Imaging_Pipeline/")

source("scripts/MultiplexImaging_utils.R")
source("scripts/GlobalVars.R")
devtools::load_all("annotation/silvermantest/")

# Change when running from R:
args = list(output = "output/")

cat("Output = ", args$output, "\n")

output_dir = file.path(args$output, "cell_neighborhood")
if (!dir.exists(output_dir))
  dir.create(output_dir)

spe = qs::qread("output/SpatialExperiment.qs")

###############################################################################
# Distances between cells & calculate K-Nearest Neighbors clusters
###############################################################################
library(FastKNN)
graph_list = list()
all_clusters = c()

library(igraph)

list_distances = list()

pdf(file.path(output_dir, "Cell_networks.pdf"))
for(samp in unique(spe$sample_id)){
  print(samp)
  spe. = spe[,spe$sample_id == samp]
  df = cbind(colData(spe.), spe.@int_colData$spatialCoords)
  m = as.matrix(df[,c("centroid.0","centroid.1")])
  d = FastKNN::Distance_for_KNN_test(m,m)
  max_dist = 400
  
  rownames(d) = colnames(d) = df$cell_id
  d = apply(d, 1, function(i){
    v = i[i <= max_dist]
    if(length(v) > 0){
      if(all(v == 0)) {
        i[which.min(i)] = max_dist
      }
    }
    i
  })
  d[d>max_dist] = Inf
  list_distances[[samp]] = d
  
  graph = d %>% as.data.frame %>% rownames_to_column %>% pivot_longer(-rowname) %>%
    filter(value !=Inf) %>% graph_from_data_frame(directed = FALSE) %>%  simplify(edge.attr.comb = "mean")
  clusters = igraph::cluster_louvain(graph, weights = E(graph)$value/ max_dist, resolution = 2)$membership
  
  spe.$cell_neighborhood = (clusters)
  spe.$cell_neighborhood = paste0(samp,"_CN", spe.$cell_neighborhood)
  spe. = colors_scExp(spe., annotCol = "cell_neighborhood")
  
  graph_list[[samp]] = graph
  
  # Plot the graph
  layout = layout_with_fr(graph, weights = E(graph)$value/ max_dist)
  
  par(bg = "white", mar = c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  vertex.color = spe.$cell_type_positive_marker_color[match(V(graph)$name, spe.$cell_id)]
  print(plot(graph,  vertex.color = vertex.color,
             layout = layout, vertex.size = 3, vertex.label = "", main = samp))
  print(legend(x = "topright", inset = c(-0.2, 0.1),  cex = 0.6, title = "Legend", ncol = 1,
               legend = cell_type_color_df[,1], col = cell_type_color_df[,2],
               fill = cell_type_color_df[,2]))
  
  vertex.color = spe.$cell_neighborhood_color[match(V(graph)$name, spe.$cell_id)]
  print(plot(graph,  vertex.color = vertex.color,
             layout = layout, vertex.size = 3, vertex.label = "", main = samp))
  print(legend(x = "topright", inset = c(-0.2, 0.1),  cex = 0.6, title = "Legend", ncol = 1,
               legend = unique(spe.$cell_neighborhood), col = unique(spe.$cell_neighborhood_color),
               fill = unique(spe.$cell_neighborhood_color)))
  
  print(plotSPE(spe., assay = "metadata", feature = "cell_neighborhood", size = 2) +
          scale_color_manual(values = unique(spe.$cell_neighborhood_color[order(spe.$cell_neighborhood)])) + ggtitle(samp))
  
  all_clusters = c(all_clusters, setNames(spe.$cell_neighborhood, spe.$cell_id))
}
dev.off()
qs::qsave(list_distances, file.path(output_dir, "list_distances.qs"))
