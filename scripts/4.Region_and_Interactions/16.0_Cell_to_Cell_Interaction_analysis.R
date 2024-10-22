# Analysis direct interaction between cells
# Find Immune Pathological Infiltrate in samples
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
              "SpatialExperiment",
              "ggpubr",
              "ggbeeswarm",
              "rstatix")
suppressPackageStartupMessages(invisible(lapply(libraries, require, character.only = TRUE)))
setwd("/mnt/RECHERCHE/GUENOVA_LAB/Project_Multiplex_Phenotyping/Miltenyi/Muliplex_Imaging_Pipeline/")

source("scripts/MultiplexImaging_utils.R")
source("scripts/GlobalVars.R")

# Change when running from R:
args = list(output = "output/")

cat("Output = ", args$output, "\n")

output_dir = file.path(args$output, "cell_interactions")
if (!dir.exists(output_dir))
  dir.create(output_dir)

spe = qs::qread("output/SpatialExperiment.qs")

list_distances = qs::qread("output/cell_neighborhood/list_distances.qs")

################################################################################
# Finding Interacting cell communities - IPI 
################################################################################

library(Matrix)
interaction_list = list()
for(i in unique(spe$sample_id)){
  interaction_list[[i]] = read.csv(file.path("output/cell_neighborhood/interactions/", paste0(i, ".csv.gz")), header = TRUE, row.names = 1)
  interaction_list[[i]] = as(as.matrix(interaction_list[[i]]), "sparseMatrix")
}

for(i in unique(spe$sample_id)){
  colnames(interaction_list[[i]]) = paste0(i, "-", gsub("X","", colnames(interaction_list[[i]])))
  rownames(interaction_list[[i]]) = paste0(i, "-", gsub("X","", rownames(interaction_list[[i]])))
}

interaction_mat = Matrix(data = 0, nrow = ncol(spe), ncol = ncol(spe), dimnames = list(colnames(spe), colnames(spe)))
interaction_mat = as(interaction_mat, "sparseMatrix")

for(i in unique(spe$sample_id)){
  print(i)
  mat = interaction_list[[i]]
  mat = mat[intersect(colnames(mat), colnames(spe)),intersect(colnames(mat), colnames(spe))]
  mat@x = mat@x
  interaction_mat[colnames(mat), colnames(mat)] = mat 
}

spe@assays@data$interaction = interaction_mat

# Find list 
for(i in sample(spe$cell_id,10)){
  spe. = spe
  spe.$touching = FALSE
  spe.$touching[spe@assays@data$interaction[i,]>0] = TRUE
  samp = spe$sample_id[spe$cell_id == i]
  
  # Plotting
  cell_overlay_mat = getImageAsMatrix(file.path("output", "segmentation", paste0(samp,"_whole_cell.tiff")))
  celltype_img_mat = get_metadata_image_overlay(spe., cell_overlay_mat, sample = samp, metadata = "touching", levels = c(FALSE, TRUE))
  
  tiff::writeTIFF(as.matrix(celltype_img_mat),
                  file.path(output_dir, paste0(i,"_touching.tiff")))
}

################################################################################
# Number of interacting celltype1 vs celltype2 - immune cells - by disease
################################################################################
spe. = spe[,!spe$condition %in% c("LN", "THY", "TON") & spe$cell_type_CellSighter != "Keratinocyte"]

dir.create(file.path(output_dir, "interacting_cells_by_disease"))
list_DA = list()
for(celltype1 in c("T_helper", "T_cytotoxic","APC","Macrophages", "Monocytic_Lineage", "Monocytes", "pDC", "B_cell", "T_regulatory")){
  print(celltype1)
  dir.create(file.path(output_dir, "interacting_cells_by_disease", celltype1))
  
  for(celltype2 in unique(spe.$cell_type_CellSighter)){
    print(celltype2)
    interaction_list = list()
    interaction_list_permuted = list()
    
    for(samp in unique(spe.$sample_id)){
      spe.. = spe.[,spe.$sample_id == samp]
      mat = spe..@metadata$interaction[spe..$cell_id,spe..$cell_id]
      n_interactions  = sum(apply(mat, 1, function(i) length(which(i > 0))))
      
      cells1 = spe..$cell_id[spe..$cell_type_CellSighter == celltype1]
      cells2 = spe..$cell_id[spe..$cell_type_CellSighter == celltype2]
      
      if(length(cells1) > 1 & length(cells2)>1){
        
        # Random permutations - same cell composition - same network
        for(i in 1:10){
          spe_permuted = spe..
          set.seed(i + 47)
          spe_permuted$cell_type_CellSighter = sample(spe_permuted$cell_type_CellSighter)
          cells1_p = spe_permuted$cell_id[spe_permuted$cell_type_CellSighter == celltype1]
          cells2_p = spe_permuted$cell_id[spe_permuted$cell_type_CellSighter == celltype2]
          mat = spe..@metadata$interaction[cells1_p,cells2_p]
          
          interaction_list_permuted[[paste0(samp, "_",i)]] =
            data.frame(
              sample_id = samp,
              interacting = length(which(apply(mat, 1, function(i) length(which(i>0))) > 0)),
              total1 = length(cells1_p),
              total2 = length(cells2_p),
              iteration = paste0("random_",i),
              n_total_interactions = n_interactions,
              condition = unique(spe..$condition))
        }
        
        mat = spe..@metadata$interaction[cells1,cells2]
        
        interaction_list[[samp]] = data.frame(
          sample_id = samp,
          interacting = length(which(apply(mat, 1, function(i) length(which(i>0))) > 0)),
          total1 = length(cells1),
          total2 = length(cells2),
          iteration = "real",
          n_total_interactions = n_interactions,
          condition = unique(spe..$condition)
        )
      }
    }
    
    df = do.call("rbind", interaction_list)
    df_p = do.call("rbind", interaction_list_permuted)
    df_p = df_p %>% group_by(sample_id) %>% 
      summarise(
        interacting = mean(interacting),
        total1 = mean(total1),
        total2 = mean(total2),
        iteration = "random",
        n_total_interactions = mean(n_total_interactions),
        condition= head(condition,1)
      )
    df = rbind(df, df_p)
    df = df %>% mutate(percent_in_contact = 100 * interacting / n_total_interactions)
    df$batch = spe$batch[match(df$sample_id, spe$sample_id)]
    df$disease = spe$disease[match(df$sample_id, spe$sample_id)]
    
    png(file.path(output_dir, "interacting_cells_by_disease", celltype1, paste0("Percent_in_contact_", celltype1,"_with_",celltype2,".png")), width = 1800, height = 1500, res = 300)
    p = df %>%
      ggplot(aes(x = disease, y = percent_in_contact,
                 fill = iteration)) +
      scale_fill_manual(values = c("#593E3CFE", "darkgreen")) +
      theme_classic() +
      theme(
        axis.text.x = element_text(angle = 90), strip.background = element_blank(),
        panel.background = element_blank(),
        panel.ontop = F,
        panel.spacing = unit(0.75, "lines"), strip.text.x = element_text(size = 12.5))
    p = p + geom_boxplot(outlier.colour = "white")
    p = p + geom_jitter(position = position_jitterdodge(),
                        colour="black",pch=21, size=2)
    p = p + scale_shape_manual(values = setNames(c(21, 24), c("TMA1", "TMA2"))) +
      guides(shape = guide_legend("TMA"))
    p = p + stat_compare_means(aes(label = after_stat(p.format)),
                               method = "t.test",  paired = TRUE)
    p = p + facet_grid(formula(paste0("~disease")), scales = "free_x", axes = "all_x",
                       axis.labels = "margins", space = "free_x", switch="both")
    
    p = p + xlab("") + ylab(paste0("Interactions ",celltype1,"<-->",celltype2," (% of all interactions)"))
    print(p)
    dev.off()
    
    list_DA[[paste0(celltype1,"_",celltype2)]] = df
  }
}

for(i in names(list_DA)){
  list_DA[[i]]$celltype_1 = gsub(paste(paste0("_", unique(spe$cell_type_CellSighter)), collapse = "|"), "", i)
  list_DA[[i]]$celltype_2 = gsub(paste(paste0(unique(spe$cell_type_CellSighter), "_"), collapse = "|"), "", i)
}

DA = do.call("rbind", list_DA)
WriteXLS::WriteXLS(DA, file.path(output_dir, "interacting_cells", "Real_vs_Permuted_all_by_disease.xlsx"))

################################################################################
# Graph based Communities based on interactions: detection of the IPI
################################################################################
library(igraph)
library(ChromSCape)

all_cluster = c()
all_cells_in_infiltrate = c()
for(samp in unique(spe$sample_id)){
  
  spe. = spe[,spe$sample_id == samp & !(spe$cell_type_CellSighter %in% struct_celltype)]
  mat = spe.@assays@data$interaction
  mat = mat[colnames(spe.), colnames(spe.)]
  
  g = igraph::graph_from_adjacency_matrix(mat)
  g = simplify(g)
  g = as.undirected(g)
  
  clusters = igraph::cluster_louvain(g, resolution = 0.001)$membership
  tab = table(clusters)
  clusters[ clusters %in% as.numeric(names(tab)[tab <= 30])] = 0
  clusters = as.character(clusters)
  n = 0
  clusters. = clusters
  for(i in sort(unique(clusters))){
    if(i == "0")
      clusters.[clusters == i] = paste0("isolated_cells")
    else
      clusters.[clusters == i] = paste0("infiltrate")
    n = n+1
  }
  spe.$immune_interaction_community = clusters.
  color_df = data.frame("immune_interaction_community" =  c("isolated_cells", "infiltrate"),
                        "immune_interaction_community_color" = c("#d1d1d1ff", "#681FC2"))
  spe. = colors_scExp(spe., annotCol = "immune_interaction_community", color_df =  color_df, color_by = "immune_interaction_community")
  
  all_cluster = c(all_cluster, setNames(spe.$immune_interaction_community, spe.$cell_id) )
  
  not_isolated = which(spe.$immune_interaction_community != "isolated_cells")
  
  if(length(not_isolated) ==  0 ){
    not_isolated = colnames(spe.)
  }
  spe.. = spe.[,not_isolated]
  g = subgraph(g,not_isolated)
  
  # Plot the graph
  png(file.path(output_dir, paste0(samp, "_interaction_graph_by_celltype.png")), width = 1000, height = 1000, res = 200)
  par(bg = "white", mar = c(1.1, 1.1, 1.1, 1.1), xpd=TRUE)
  set.seed(49)
  layout = layout_nicely(g)
  vertex.color = spe.$cell_type_CellSighter_color[match(V(g)$name, spe..$cell_id)]
  print(plot(g,  layout = layout, vertex.size = 1, vertex.label = ""))
  dev.off()
  
  # Find the regions around the immune infitlrate, also containing non immune cells
  spe_samp = spe[,spe$sample_id == samp]
  
  for(clust in setdiff(unique(clusters), "0")){
    cells = spe.$cell_id[clusters == clust]
    polygon = get_convex_hull(spe., cells = cells)
    coords = spe_samp@int_colData$spatialCoords
    cells_in_IPI_region = point.in.polygon(coords[,1], coords[,2], polygon[,1], polygon[,2])
    cells_in_IPI_region = spe_samp$cell_id[cells_in_IPI_region == 1]
    all_cells_in_infiltrate = c(all_cells_in_infiltrate, cells_in_IPI_region)
  }
  
  spe_samp$immune_interaction_region = "isolated_cells"
  spe_samp$immune_interaction_region[spe_samp$cell_id %in% all_cells_in_infiltrate] = "infiltrate"
  
  png(file.path(output_dir, paste0(samp, "_interaction_graph_by_community_all.png")), width = 1000, height = 1000, res = 200)
  print(plotSPE(spe_samp, assay = "metadata", feature = "immune_interaction_region", size = 1) +
          scale_color_manual(values = unique(spe.$immune_interaction_community_color[order(spe.$immune_interaction_community)])) +
          ggtitle(samp))
  dev.off()
  
  # Plotting
  cell_overlay_mat = getImageAsMatrix(file.path("output", "segmentation", paste0(samp,"_whole_cell.tiff")))
  celltype_img_mat = get_metadata_image_overlay(spe_samp, cell_overlay_mat, sample = samp,
                                                metadata = "immune_interaction_region",
                                                levels = c("isolated_cells", "infiltrate"))
  
  tiff::writeTIFF(as.matrix(celltype_img_mat),
                  file.path(output_dir, paste0(samp,"_infiltrate.tiff")))
}

spe$immune_interaction_community = NA
spe$immune_interaction_community[match(names(all_cluster), spe$cell_id)] = all_cluster

spe$immune_interaction_region = "isolated_cells"
spe$immune_interaction_region[spe$cell_id %in% all_cells_in_infiltrate] = "infiltrate"

qs::qsave(spe, "output/SpatialExperiment.qs")

# Boxplot of % infiltrate / all immune per condition with stats
meta = as.data.frame(colData(spe)) %>% filter(!condition %in% c("THY", "TON", "LN"))
meta = meta %>% filter(!cell_type_CellSighter %in% struct_celltype) %>%
  dplyr::group_by(sample_id, disease, batch, condition) %>% dplyr::summarise(ratio = length(which(immune_interaction_community == "infiltrate")) /
                                                                               length(which(!is.na(immune_interaction_community))))

library(ggpubr)
png(file.path(output_dir, paste0("Ratio_Immune_Infiltrate.png")), width = 1300, height = 1100, res = 250)
p = grouped_dotplot(meta, y = "ratio", categ1 = "condition", categ2 = "disease", 
                    ref.group = "HD",
                    color_categ1 = setNames(c("grey", colCondition$condition_color), c("Rest", colCondition$condition)))
p = p + xlab("") + ylab(paste0("Immune infiltrate (% of immune cells)")) +
  guides(fill="none", color = "none")
print(p)
dev.off()

  meta %>% ungroup %>%
  t_test(formula = ratio ~ disease, paired = FALSE, p.adjust.method = "none")

  meta %>% ungroup %>% mutate(condition = as.character(condition)) %>% 
  t_test(formula = ratio ~ condition, ref.group = "HD", paired = FALSE, p.adjust.method = "none")

################################################################################
# Grouped by IPI group
################################################################################
library(ggpubr)
png(file.path(output_dir, paste0("Ratio_Immune_Infiltrate_grouped.png")),
    width = 1200, height = 1200, res = 300)
meta = meta %>% filter(!condition %in% c("THY", "TON", "LN"))
meta$condition = factor(meta$condition, c("HD", "AD", "PS", "DAR", "LP", "MF", "SS"))
p = grouped_dotplot(meta, y = "ratio", categ1 = "condition", categ2 = "disease", 
                    ref.group = "HD",
                    color_categ1 = setNames(c("grey", colCondition$condition_color), c("Rest", colCondition$condition)))
p = p + xlab("") + ylab(paste0("Immune infiltrate (% of immune cells)")) +
  guides(fill="none", color = "none")
print(p)
dev.off()

mean(meta$ratio[meta$condition %in% c("AD", "PS", "DAR")])
sd(meta$ratio[meta$condition %in% c("AD", "PS", "DAR")])
mean(meta$ratio[meta$condition %in% c("LP", "MF", "SS")])
sd(meta$ratio[meta$condition %in% c("LP", "MF", "SS")])
test = t.test(meta$ratio[meta$condition %in% c("AD", "PS", "DAR")],
              meta$ratio[meta$condition %in% c("LP", "MF", "SS")])
cat("P-value = ", test$p.value)

################################################################################
# Variance of the infiltrate 
################################################################################
library(ggpubr)
library(car)
library(caret)

meta$disease = "NonCTCL"
meta$disease[meta$condition %in% c("MF", "SS")] = "CTCL"
levene_test_result <- leveneTest(ratio ~ condition, data = meta %>% filter(condition != "HD"))
levene_test_result

meta. = meta %>% dplyr::group_by(condition) %>% dplyr::summarise(variance = var(100*ratio))


png(file.path(output_dir, paste0("Variance_Immune_Infiltrate_grouped.png")),
    width = 1200, height = 1200, res = 300)
meta. = meta. %>% filter(!condition %in% c("THY", "TON", "LN"))
meta.$condition = factor(meta.$condition, c("HD", "AD", "PS", "DAR", "LP", "MF", "SS"))
p = meta. %>% ggplot(aes(x = condition, y = variance, fill = condition)) +
  geom_bar(stat =  "identity") +
  scale_fill_manual(values = unique(spe$condition_color[match(levels(meta.$condition), spe$condition)])) +
  theme_classic() + xlab("") + ylab(paste0("Variance")) +
  NoLegend()  + 
  theme(axis.text.x = element_text(angle = 90))
print(p)
dev.off()

################################################################################
# Barplot of cell composition in infiltrate vs not  immune per condition 
################################################################################
library(dittoSeq)
spe. = spe[,!spe$cell_type_CellSighter %in% c("Leukocyte", struct_celltype)]
spe. = spe.[,!spe.$condition %in% c("HD","THY", "TON", "LN")]
spe.$cell_type_simplified = spe.$cell_type_CellSighter
spe.$cell_type_simplified[spe.$cell_type_CellSighter %in% c("T_cytotoxic","T_regulatory", "T_helper", "NKT")] = "Lymphoid"
spe.$cell_type_simplified[!spe.$cell_type_CellSighter %in% c("T_cytotoxic","T_regulatory", "T_helper", "NKT")] = "Myeloid"

png(file.path(output_dir, paste0("Composition_Immune_Infiltrate.png")), width = 1800, height = 1600, res = 300)
colors =  setNames(unique(spe.$cell_type_CellSighter_color), unique(spe.$cell_type_CellSighter))
dittoBarPlot(spe., "cell_type_CellSighter", group.by = "condition",  main = "Composition of Immune Infiltrate",
             color.panel = colors,
             split.adjust = list(scales = 'free'), retain.factor.levels = F, 
             split.by = "immune_interaction_community", split.nrow = 1)
dev.off()

# Fold Change
meta = as.data.frame(colData(spe.))[,c("condition", "immune_interaction_community", "cell_type_simplified")]
meta = meta %>% group_by(condition, immune_interaction_community) %>%
  dplyr::summarise(n_tot = n(), n_lymphoid = length(which(cell_type_simplified == "Lymphoid"))) %>%
  mutate(percent = 100 * n_lymphoid / n_tot )
t.test(meta$percent[meta$immune_interaction_community == "infiltrate"], meta$percent[meta$immune_interaction_community != "infiltrate"])

spe.$condition = factor(spe.$condition, levels = c("AD", "PS", "DAR", "LP",  "MF", "SS"))
png(file.path(output_dir, paste0("Composition_Immune_Infiltrate_simplified_grouped_IPI.png")), width = 1400, height = 1200, res = 300)
dittoBarPlot(spe., "cell_type_simplified", group.by = "condition",  main = "Composition of Immune Infiltrate",
             split.adjust = list(scales = 'free'), retain.factor.levels = T, 
             split.by = "immune_interaction_community", split.nrow = 1, ylab = "Fraction of cells", 
             color.panel = cell_type_simplified_color_df[c("Myeloid", "Lymphocyte")])
dev.off()


################################################################################
# Compare markers from IPI vs non-IPI
################################################################################
spe. = spe[,!is.na(spe$immune_interaction_community)]
spe. = spe.[,!spe.$condition %in% c("HD","THY", "TON", "LN")]
spe.$cell_type_simplified = spe.$cell_type_CellSighter
spe.$cell_type_simplified[spe.$cell_type_CellSighter %in% c("T_cytotoxic","T_regulatory", "T_helper", "NKT")] = "Lymphoid"
spe.$cell_type_simplified[!spe.$cell_type_CellSighter %in% c("T_cytotoxic","T_regulatory", "T_helper", "NKT")] = "Myeloid"

table(spe.$cell_type_simplified, spe.$condition)
cells = c()
for(cond in unique(spe.$condition)){
  spe_cond = spe.[,spe.$condition == cond]
  myeloid = sample(spe_cond$cell_id[spe_cond$cell_type_simplified == "Myeloid"], 3500, replace = F)
  lymphoid = sample(spe_cond$cell_id[spe_cond$cell_type_simplified == "Lymphoid"], 3500, replace = F)
  
  cells = c(cells, myeloid, lymphoid)
}
spe. = spe.[,match(cells, spe.$cell_id)]

mat = spe.@assays@data$CellSighter_marker_mat
mat = t(mat)
mat. = mat[,setdiff(colnames(mat), c("Ki67","CD209", "CD270", "CXCR4"))]

intensities. = as.data.frame(mat.) 
intensities.$sample_id = spe.$sample_id[match(rownames(intensities.), spe.$cell_id)] 
intensities.$condition = spe.$condition[match(rownames(intensities.), spe.$cell_id)] 
intensities.$immune_interaction_community = spe.$immune_interaction_community[match(rownames(intensities.), spe.$cell_id)] 
intensities.$batch = spe.$batch[match(rownames(intensities.), spe.$cell_id)] 

intensities.. = intensities. %>% group_by(condition, immune_interaction_community) %>%
  summarise_at(1:57, mean) %>% arrange(immune_interaction_community, condition)

intensities.. = intensities.. %>% filter(!condition %in% c("HD","THY", "TON", "LN"))

mat. = as.matrix(intensities..[,setdiff(colnames(mat.), unique(cell_type_markers$marker))]) #c("HLAABC", "CD2", "CD4", "CD8a", "CD3", "CD5", "HLADR")]


FC = apply(mat., 2, function(i) {
  FC = mean(i[intensities..$immune_interaction_community == "infiltrate"]) /
    mean(i[intensities..$immune_interaction_community != "infiltrate"])
  FC
} )

pvals = sapply(colnames(intensities..)[4:57], function(i) t.test(intensities..[intensities..$immune_interaction_community == "infiltrate",i],
                                                                 intensities..[intensities..$immune_interaction_community != "infiltrate",i])$p.value)


cormat = cor(mat.)
cormat[is.na(cormat)] = 0
hc_cor_marker = hclust(as.dist(1 - cormat), method = "ward.D2")

color_df = data.frame(
  "condition" = intensities..$condition,
  # "sample_id" = intensities..$sample_id,
  "disease_type" = intensities..$immune_interaction_community
)

col = list(
  condition = setNames(unique(spe.$condition_color), unique(spe.$condition)),
  # sample_id = setNames(unique(spe.$sample_id_color), unique(spe.$sample_id)),
  disease_type = setNames(c("#d1d1d1ff", "#681FC2"), unique(spe.$immune_interaction_community))
  
)

mat. = mat.[,order(FC, decreasing = TRUE)]
# Select top and bottom 5 non-celltype markers
mat_scaled = scaleMat(mat., scale = "column")[,c("CD2","TRBC1", "CD5", "Bcl2", "ZAP70", "CollagenIII", "CD324", "CD204", "TIM3", "CD13")]

library(circlize)
library("RColorBrewer")
library("ComplexHeatmap")
# c("royalblue4", "white", "gold")
col_fun = colorRamp2(c(-2:2), rev(brewer.pal(5,"RdBu")))

h = Heatmap( 
  mat_scaled,
  cluster_columns = F,
  cluster_rows = F,
  col = col_fun, 
  name = "Positive Marker",
  column_title = "Samples",
  row_title = "Markers",
  row_dend_side = "right",
  row_dend_width = unit(4, "cm"),
  column_dend_height = unit(5, "cm"),
  show_column_names = T,
  show_row_names = F,
  clustering_distance_columns ="pearson",
  clustering_distance_rows = "pearson",
  column_split = c(rep(1,5), rep(2,5)),
  row_split = c(rep(1,6), rep(2,6)),
  border = TRUE,
  use_raster = FALSE,
  right_annotation = rowAnnotation(df = color_df,
                                   col = col)
)
draw(h)

signature_pos = c("CD2","TRBC1", "CD5", "Bcl2", "ZAP70")
signature_neg = c("CollagenIII", "CD324", "CD204", "TIM3", "CD13")

intensities.. = intensities. %>% group_by(sample_id, condition, immune_interaction_community) %>%
  summarise_at(1:57, mean) %>% arrange(immune_interaction_community, condition)

pvals = sapply(colnames(mat_scaled), function(i) t.test(intensities..[intensities..$immune_interaction_community == "infiltrate",i],
                                                        intensities..[intensities..$immune_interaction_community != "infiltrate",i])$p.value)

pvals[pvals < 0.0001] # ****
pvals[pvals < 0.001] # ***
pvals[pvals < 0.01] # **
pvals[pvals < 0.05] # *

png(file.path(output_dir, paste0("Heatmap_markers_infiltrate_vs_not_infiltrate_inflammatory_scaled_col.png")),
    width = 1400,
    height = 1200, res = 200)
draw(h)
dev.off()


################################################################################
# Grouped
# Boxplots of marker infiltrate vs rest
################################################################################

pdf(file.path(output_dir, paste0("Boxplot_markes_infiltrate_vs_rest.pdf")), width = 6, height = 5)
for(marker in setdiff(rownames(spe.),c("Ki67","CD209","CD279","CD270", "CXCR4"))){
  
  p =  intensities.. %>%
    mutate(immune_interaction_community = factor(immune_interaction_community, levels = c("isolated_cells", "infiltrate"))) %>%
    ggplot(aes(x = immune_interaction_community, y = .data[[marker]], fill = immune_interaction_community)) +
    geom_boxplot(outlier.colour =  "white") +
    geom_line(aes(group=sample_id), linewidth = 0.2) +
    scale_fill_manual(values =  setNames(c("#d1d1d1ff", "#681fc2"), c("isolated_cells", "infiltrate"))) +
    stat_compare_means(aes(label = after_stat(p.signif)),
                       method = "t.test", paired = TRUE, ref.group = "isolated_cells") +
    ggnewscale::new_scale("fill") +
    geom_jitter(aes(fill = condition),   colour="black",pch=21,  size=2, width = 0) +
    scale_fill_manual(values = setNames(unique(spe.$condition_color), unique(spe.$condition))) +
    theme_classic() + xlab("") + ylab(paste0("Cells expressing ",marker, " (%)")) + 
    ggtitle(marker) +
    theme(axis.text.x = element_text(angle = 90))  +
    geom_hline(yintercept = median(intensities..[[marker]]), lwd = 0.35, lty = 2, col = "grey20")
  print(p)
}
dev.off()

#################################################################################
# Comparing infiltrate vs isolated immune cells 
################################################################################# 
spe. = spe[setdiff(rownames(spe), markers_to_remove), !spe$condition %in% c("THY", "TON", "LN", "HD")]

list_DA_IPI_vs_isolated = list()

pdf(file.path(output_dir, paste0("list_DA_IPI_vs_isolated_per_condition.pdf")), width = 6, height = 5)
for(cond in setdiff(spe.$condition, "HD")){
    
    spe_isolated = spe.[,which(spe.$condition == cond & spe.$immune_interaction_community == "isolated_cells")]
    isolated_cells = spe_isolated$cell_id
  
    
    spe_disease = spe.[, which(spe.$condition %in% cond  & spe.$immune_interaction_community == "infiltrate")]
    disease_cells = spe_disease$cell_id

    if(length(isolated_cells) > 10 & length(disease_cells) > 10){
      
      mat1 = spe.@assays@data$CellSighter_marker_mat[,isolated_cells]
      mat2 = spe.@assays@data$CellSighter_marker_mat[,disease_cells]
      
      values1 = rowSums(mat1) 
      values2 = rowSums(mat2) 
      
      percent1 = 100 * rowSums(mat1) / ncol(mat1)
      percent2 = 100 * rowSums(mat2) / ncol(mat2)
      
      p.values = sapply(1:nrow(mat1), function(i){
        contingency_table = data.frame(c(values1[i], values2[i]),
                                       c(ncol(mat1) - values1[i], ncol(mat2) - values2[i]))
        fisher.test(contingency_table)$p.value
      })
      p.values = p.adjust(p.values, method = "bonferroni")
      
      log2FC = log2(percent2 / percent1)
      effect = abs(percent2 - percent1)
      
      df = data.frame(condition = cond, gene = names(values1), value.HD = values1, 
                      value.condition = values2, percent.HD = percent1,
                      percent.cond = percent2, log2FC = log2FC, 
                      p.value.adjusted = p.values, effect = effect)
      
      df = df %>% mutate(categ = ifelse( p.values < 0.01, ifelse(log2FC > 1, "signifplus", ifelse(log2FC < -1, "signifminus", "unsig")), "unsig"))
      
      p = df  %>%
        ggplot(aes(x = log2FC, y = -log10(p.values), color = categ, label = gene)) +
        geom_point() +
        scale_color_manual(values = c("darkgreen","red","black")) +
        theme_classic() + xlab("Percent - log2FC") + ylab("-log10(adjusted p-value)") + 
        ggtitle(paste0(cond, " - volcano plot - Infiltrate vs HD")) + 
        geom_hline(yintercept = 2, lwd = 0.35, lty = 2, col = "grey20") +
        geom_vline(xintercept = -1, lwd = 0.35, lty = 2, col = "grey20") +
        geom_vline(xintercept = 1,  lwd = 0.35, lty = 2, col = "grey20") +
        ggrepel::geom_text_repel()
      print(p)
      
    }
    list_DA_IPI_vs_isolated[[paste0(cond)]] = df
}
dev.off()

DA_IPI_vs_isolated = do.call("rbind", list_DA_IPI_vs_isolated)
WriteXLS::WriteXLS(DA_IPI_vs_isolated, file.path(output_dir, "DA_IPI_vs_isolated_per_condition.xlsx"))

pdf(file.path(output_dir, paste0("top_7_over_genes.pdf")), width = 6, height = 5)
DA_IPI_vs_isolated %>% filter(categ == "signifplus", !gene %in% unique(cell_type_markers$marker)) %>% group_by(gene) %>% summarise(n_diff = n()) %>%
  arrange(desc(n_diff)) %>% head(7) %>% ggplot() + 
  geom_bar(aes(x = forcats::fct_reorder(gene, n_diff) ,y = n_diff), stat = "identity", fill = "red") +
  xlab("") + ylab("# Diseases overrexpressed in") + theme_pubr()
dev.off()

pdf(file.path(output_dir, paste0("top_7_under_genes.pdf")), width = 6, height = 5)
DA_IPI_vs_isolated %>%  filter(categ == "signifminus" ,!gene %in% unique(cell_type_markers$marker)) %>% group_by(gene) %>% summarise(n_diff = n()) %>%
  arrange(desc(n_diff)) %>% head(7) %>% ggplot() + 
  geom_bar(aes(x = forcats::fct_reorder(gene, desc(n_diff)) ,y = n_diff), stat = "identity", fill = "forestgreen") +
  xlab("") + ylab("# Diseases underexpressed in") + theme_pubr()
dev.off()


#################################################################################
# Vascularization of blood vessels
################################################################################# 
run_vasc_blood = TRUE

if(run_vasc_blood == TRUE){
  library(MASS)
  library(cluster)
  source("scripts/MultiplexImaging_utils.R")
  
  list_vascular_df = list()
  spe. = spe[,!spe$condition %in% c("THY", "TON", "LN")]
  
  for(samp in unique(spe$sample_id)){
    print(samp)
    
    spe.. = spe.[,spe.$sample_id == samp]
    if(ncol(spe..) > 2){
      
      mat = spe..@metadata$interaction
      mat = mat[colnames(spe..), colnames(spe..)]
      
      vascular_cells = spe..$cell_id[spe..$cell_type_CellSighter == "Endothelial"]
      mat_vascular = mat[vascular_cells, vascular_cells]
      
      g = igraph::graph_from_adjacency_matrix(mat_vascular)
      g = igraph::simplify(g)
      g = igraph::as.undirected(g)
      
      clusters = igraph::cluster_louvain(g, resolution = 0.001)$membership
      names(clusters) = vascular_cells
      
      print(plot(g, vertex.size = 0.5, vertex.label = NA, main = samp))
      tab = table(clusters)
      vessels = tab[tab > 2]
      clusters = clusters[clusters %in% names(vessels)]
      
      vascular_cells = names(clusters)
      mat_vascular = mat[vascular_cells, vascular_cells]
      
      infiltrate_cells = spe..$cell_id[which(spe..$immune_interaction_community == "infiltrate")]
      
      if(length(clusters) > 0) {
        
        normalized_diameter = rep(0, length(unique(clusters)))
        names(normalized_diameter) = paste0("vessel_", unique(clusters))
        vessels_in_infiltrate = c()
        
        for(vessel in unique(clusters)){
          cells = names(clusters[clusters == vessel])
          
          ell = get_ellipse(spe = spe.., cells = cells)
          
          # Adjust the margins to ensure the entire ellipse is visible
          xlim_min <- min(ell$ellipse_points[,1])
          xlim_max <- max(ell$ellipse_points[,1])
          ylim_min <- min(ell$ellipse_points[,2])
          ylim_max <- max(ell$ellipse_points[,2])
          
          points_matrix = as.matrix(spe@int_colData$spatialCoords[cells,])
          
          print(plot(points_matrix[,1], points_matrix[,2], pch = 19, col = "blue", xlab = "X", ylab = "Y", main = "Smallest Ellipse",
                     xlim = c(xlim_min, xlim_max), ylim = c(ylim_min, ylim_max)))
          
          # Plot the ellipse
          print(lines(ell$ellipse_points[,1], ell$ellipse_points[,2], col = "red", lwd = 2))
          
          normalized_diameter[ paste0("vessel_", vessel)] = ell$normalized_diameter
          
          if(any(cells %in% spe..$cell_id[spe..$immune_interaction_region == "infiltrate"])) 
            vessels_in_infiltrate = c(vessels_in_infiltrate, vessel)
        }
        
        list_vascular_df[[samp]] = data.frame(sample_id = samp,
                                              n_total = ncol(spe..),
                                              n_infiltrate = length(infiltrate_cells),
                                              n_vessels = length(vessels),
                                              n_vessels_in_infiltrate = length(vessels_in_infiltrate),
                                              n_vessels_in_isolated =  length(vessels) - length(vessels_in_infiltrate),
                                              mean_normalized_diameter_in_infiltrate = mean(normalized_diameter[paste0("vessel_", vessels_in_infiltrate)]),
                                              mean_normalized_diameter_in_isolate = mean(normalized_diameter[setdiff(paste0("vessel_", names(vessels)),
                                                                                                                     paste0("vessel_", vessels_in_infiltrate))])
        )
        
      }
    }
  }
  
  vascular_df = do.call("rbind", list_vascular_df)
  vascular_df$condition = spe$condition[match(vascular_df$sample_id, spe$sample_id)]
  vascular_df$fraction_infiltrate = vascular_df$n_infiltrate / vascular_df$n_total
  vascular_df = vascular_df %>% mutate(vascularization_score_in_infiltrate = (mean_normalized_diameter_in_infiltrate * (n_vessels_in_infiltrate / n_vessels))  / fraction_infiltrate)
  vascular_df = vascular_df %>% mutate(vascularization_score_in_isolated = (mean_normalized_diameter_in_isolate * (n_vessels_in_isolated / n_vessels))  / (1-fraction_infiltrate) )
  vascular_df = vascular_df %>% filter(
    !is.nan(vascularization_score_in_infiltrate), !is.na(vascularization_score_in_infiltrate), !is.infinite(vascularization_score_in_infiltrate),
    !is.nan(vascularization_score_in_isolated), !is.na(vascularization_score_in_isolated), !is.infinite(vascularization_score_in_isolated))
  
  
  # Compare infiltrate vascularization vs non infiltrate vascularization
  vascular_df. = vascular_df %>% pivot_longer(cols = c(vascularization_score_in_infiltrate, vascularization_score_in_isolated),
                                              names_to = "infiltrate", values_to = "vascularization_score", names_prefix =  "vascularization_score_")
  vascular_df.$infiltrate = factor(vascular_df.$infiltrate, levels =c("in_isolated", "in_infiltrate") )
  library(ggpubr)
  
  png(file.path(output_dir, paste0("Vascularization_Score_of_infiltrate.png")), width = 1400, height = 1200, res = 200)
  p = vascular_df.  %>% filter(!condition %in% c("HD","TON", "THY", "LN") ) %>%
    ggplot(aes(x = infiltrate, y = vascularization_score, fill = infiltrate)) +
    geom_boxplot(outlier.colour =  "white") +
    scale_fill_manual(values =  setNames(c("#d1d1d1ff", "#681fc2"), c("in_isolated", "in_infiltrate"))) +
    geom_jitter(aes(color = sample_id),  colour="black",pch=21, size=2, width = 0.3) +
    theme_classic() + xlab("") + ylab(paste0(" Vascularization Score")) + 
    ggtitle(" Vascularization score of Infiltrates") +
    theme(axis.text.x = element_text(angle = 90))  + ylim(c(0,120)) + 
    stat_compare_means(aes(label = after_stat(p.signif)), label.y.npc = 0.75,
                       method = "t.test", ref.group = "in_isolated") +
    facet_wrap(~condition)
  print(p)
  dev.off()
  
  png(file.path(output_dir, paste0("Vascularization_Score_of_infiltrate_grouped.png")), width = 1300, height = 1400, res = 300)
  p = vascular_df.  %>% filter(!condition %in% c("HD","TON", "THY", "LN") ) %>%
    ggplot(aes(x = infiltrate, y = vascularization_score, fill = infiltrate)) +
    geom_boxplot(outlier.colour =  "white") +
    scale_fill_manual(values =  setNames(c("#d1d1d1ff", "#681fc2"), c("in_isolated", "in_infiltrate"))) +
    geom_jitter(aes(color = sample_id),  colour="black",pch=21, size=2, width = 0.3) +
    theme_classic() + xlab("") + ylab(paste0(" Vascularization Score")) + 
    ggtitle(" Vascularization score of Infiltrates") +
    theme(axis.text.x = element_text(angle = 90))  + ylim(c(0,120)) + 
    stat_compare_means( label.y.npc = 1,
                       method = "t.test", ref.group = "in_isolated")
  print(p)
  dev.off()
  
  
  png(file.path(output_dir, paste0("Relative_vessel_number.png")), width = 1400, height = 1200, res = 200)
  p = vascular_df  %>% filter(!condition %in% c("HD","TON", "THY", "LN") ) %>%
    ggplot(aes(x = condition, y = n_vessels/n_total, fill = condition)) +
    geom_boxplot(outlier.colour =  "white") +
    scale_fill_manual(values =  setNames(unique(spe$condition_color), unique(spe$condition))) +
    geom_jitter(aes(color = sample_id),  colour="black",pch=21, size=2, width = 0.3) +
    theme_classic() + xlab("") + ylab(paste0(" Vascularization Score")) + 
    ggtitle(" Vascularization score of Infiltrates") +
    theme(axis.text.x = element_text(angle = 90))  + 
    stat_compare_means(aes(label = after_stat(p.signif)), label.y.npc = 0.75,
                       method = "t.test", ref.group = ".all.")  
  print(p)
  dev.off()
  
  write.csv(vascular_df, file.path(output_dir, "vascular_df.csv"))
  
  # Correlation of vascularization of infiltrate with cell composition
  spe. = spe[,!is.na(spe$immune_interaction_community)]
  spe.$cell_type_simplified = spe.$cell_type_CellSighter
  spe.$cell_type_simplified[spe.$cell_type_CellSighter %in% c("T_cytotoxic","T_regulatory", "T_helper", "NKT")] = "Lymphoid"
  spe.$cell_type_simplified[!spe.$cell_type_CellSighter %in% c("T_cytotoxic","T_regulatory", "T_helper", "NKT")] = "Myeloid"
  
  meta = as.data.frame(colData(spe.))
  meta = meta %>% group_by(sample_id, immune_interaction_community, cell_type_simplified) %>%
    summarise(celltype = n())
  meta = meta %>% group_by(sample_id, immune_interaction_community) %>% 
    summarise(n_Lymphoid = celltype[which(cell_type_simplified == "Lymphoid")],
              n_Myeloid = celltype[which(cell_type_simplified == "Myeloid")],
              percent_Lymphoid = 100 * n_Lymphoid  / (n_Lymphoid + n_Myeloid)
    )
  meta$infiltrate = ifelse(meta$immune_interaction_community == "infiltrate", "in_infiltrate", "in_isolated")
  vascular_df.. = left_join(vascular_df., meta, by = c("sample_id", "infiltrate"))
  vascular_df..$percent_Lymphoid
  
  png(file.path(output_dir, paste0("DotPlot_percent_Lymphoid_vascularization_in_infiltrate.png")), width = 1400, height = 1200, res = 200)
  vascular_df.. = vascular_df..  %>% filter(!condition %in% c("HD","TON", "THY", "LN"),
                                            immune_interaction_community == "infiltrate")
  p = vascular_df.. %>%
    ggplot(aes(x = (vascularization_score), y = percent_Lymphoid, fill = condition)) +
    geom_point(color = "black", pch = 21) +
    scale_fill_manual(values =  setNames(unique(spe$condition_color), unique(spe$condition))) +
    theme_classic()  + ylab(paste0(" Percent Lymphoid")) + 
    ggtitle("Percent lymphoid vs Vascularization in infiltrate ")  + xlim(c(0,150))
  
  print(p)
  dev.off()
  
  # Compare infiltrate Mean_Normalized_Diameter vs non infiltrate vascularization
  vascular_df. = vascular_df %>% pivot_longer(cols = c(mean_normalized_diameter_in_infiltrate, mean_normalized_diameter_in_isolate),
                                              names_to = "infiltrate", values_to = "mean_normalized_diameter", names_prefix =  "mean_normalized_diameter_")
  vascular_df.$infiltrate = factor(vascular_df.$infiltrate, levels =c("in_isolate", "in_infiltrate") )
  library(ggpubr)
  
  png(file.path(output_dir, paste0("Mean_Normalized_Diameter_of_infiltrate.png")), width = 1400, height = 1200, res = 200)
  p = vascular_df.  %>% filter(!condition %in% c("HD","TON", "THY", "LN") ) %>%
    ggplot(aes(x = infiltrate, y = mean_normalized_diameter, fill = infiltrate)) +
    geom_boxplot(outlier.colour =  "white") +
    scale_fill_manual(values =  setNames(c("#d1d1d1ff", "#681fc2"), c("in_isolate", "in_infiltrate"))) +
    geom_jitter(aes(color = sample_id),  colour="black",pch=21, size=2, width = 0.3) +
    theme_classic() + xlab("") + ylab(paste0(" Mean Normalized Diameter")) + 
    ggtitle(" Mean Normalized Diameter of Infiltrates") +
    theme(axis.text.x = element_text(angle = 90))  + 
    stat_compare_means(aes(label = after_stat(p.signif)), label.y.npc = 0.75,
                       method = "t.test", ref.group = "in_isolate", paired = TRUE) +
    facet_wrap(~condition)
  print(p)
  dev.off()
  
  # Vascularization
  vascular_df. = vascular_df %>% pivot_longer(cols = c(n_vessels_in_infiltrate, n_vessels_in_isolated),
                                              names_to = "infiltrate", values_to = "vessels", names_prefix =  "n_vessels_")
  vascular_df.$infiltrate = factor(vascular_df.$infiltrate, levels =c("in_isolated", "in_infiltrate") )
  
  png(file.path(output_dir, paste0("Number_of_vessels_in_infiltrate.png")), width = 1400, height = 1200, res = 200)
  p = vascular_df.  %>% filter(!condition %in% c("HD","TON", "THY", "LN") ) %>%
    ggplot(aes(x = infiltrate, y = vessels, fill = infiltrate)) +
    geom_boxplot(outlier.colour =  "white") +
    scale_fill_manual(values =  setNames(c("#d1d1d1ff", "#681fc2"), c("in_isolated", "in_infiltrate"))) +
    geom_jitter(aes(color = sample_id),  colour="black",pch=21, size=2, width = 0.3) +
    theme_classic() + xlab("") + ylab(paste0(" Number of vessels")) + 
    ggtitle(" Number of vessels in infiltrate") +
    theme(axis.text.x = element_text(angle = 90))  + 
    stat_compare_means(aes(label = after_stat(p.signif)), label.y.npc = 0.75,
                       method = "t.test", ref.group = "in_isolated", paired = TRUE) +
    facet_wrap(~condition)
  print(p)
  dev.off()

  qs::qsave(vascular_df, file.path(output_dir, "vascular_df.qs"))
}


