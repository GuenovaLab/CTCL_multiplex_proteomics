# Celltype and protein expression analysis using deep-learning based approaches
# authors: Pacome Prompsy
# contact: pacome.prompsy@chuv.ch
# Guenova Lab
# CHUV (Centre Hospitalier Universitaire Vaudois), Lausanne, Suisse
cat("Analyzing results from deep-lerarning based cell type and marker detection... \n")

# Loading packages --------------------------------------------------------
libraries = c("argparse",
              "ggplot2",
              "dplyr",
              "tidyr",
              "Seurat",
              "SpatialExperiment",
              "arrow",
              "ChromSCape",
              "ggbeeswarm",
              "rstatix")
suppressPackageStartupMessages(invisible(lapply(libraries, require, character.only = TRUE)))
setwd("/mnt/RECHERCHE/GUENOVA_LAB/Project_Multiplex_Phenotyping/Miltenyi/Muliplex_Imaging_Pipeline/")

source("scripts/MultiplexImaging_utils.R")
source("scripts/GlobalVars.R")
library(ggplot2)
library(tidyverse)
library(ComplexHeatmap)

# Reading in data  -------------------------------------------------------------
output = file.path("output", "CellSighter")
output_marker = file.path(output, "marker")
output_celltype = file.path(output, "celltype")

config_marker = jsonlite::parse_json(file(file.path(output_marker, "Predictions", "config.json")))
corresponding_marker = unlist(config_marker$hierarchy_match)

spe = qs::qread("output/SpatialExperiment.qs")

# Make simplified cell types ---------------------------------------------------
spe$cell_type_simplified = spe$cell_type_CellSighter
spe$cell_type_simplified[spe$cell_type_CellSighter %in% c("T_cytotoxic","T_regulatory", "T_helper", "NKT", "Leukocyte","B_cell")] = "Lymphocyte"
spe$cell_type_simplified[spe$cell_type_CellSighter %in% c("Macrophages","Monocytic_Lineage","Basophil", "APC", "pDC", "Neutrophils", "Monocytes")] = "Myeloid"
spe$cell_type_simplified[spe$cell_type_CellSighter %in% c("Keratinocyte")] = "Epidermis"
spe$cell_type_simplified[spe$cell_type_CellSighter %in% c("Unknown-Stroma")] = "Dermis"
spe$cell_type_simplified[spe$cell_type_CellSighter %in% c("Endothelial", "Lymphatic")] = "Vessels"
spe$disease = "Benign"
spe$disease[spe$condition %in% c("MF", "SS")] = "Malign"
spe$disease[spe$condition %in% c("HD", "TON", "LN", "THY")] = "Healthy"
spe$disease = factor(spe$disease, levels = c("Healthy", "Benign", "Malign"))
spe$disease_color = setNames(c(c("#7699D4","#B3A632FA", "#b81f53ff")), c("Healthy", "Benign", "Malign"))[match(spe$disease,c("Healthy", "Benign", "Malign"))]
spe$batch[spe$batch == "batch_1"] = "TMA1"
spe$batch[spe$batch == "batch_2"] = "TMA2"

qs::qsave(spe, "output/SpatialExperiment.qs")

# Composition analysis ---------------------------------------------------------
library(dittoSeq)

pdf(file.path(output, "Plots", "Composition_analysis.pdf"))
dittoBarPlot(spe, "cell_type_CellSighter", group.by = "condition",
             color.panel =  unique(spe$cell_type_CellSighter_color[order(spe$cell_type_CellSighter)]),
             main = "Condition x Cell Type")  + Seurat::DarkTheme()
dev.off()


d = dittoBarPlot(spe, "cell_type_CellSighter", group.by = "sample_id",
                 color.panel = unique(spe$cell_type_CellSighter_color)[order(spe$cell_type_CellSighter)],  main = "Sample x Clusters",
                 split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)

png(file.path(output, "Plots", "Sample x CellType x Condition.png"), width = 3500, height = 1600, res = 300)
colors =  setNames(unique(spe$cell_type_CellSighter_color), unique(spe$cell_type_CellSighter))
dittoBarPlot(spe, "cell_type_CellSighter", group.by = "sample_id", color.panel = colors[levels(d$data$label)],  main = "Sample x Clusters",
             split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)
dev.off()

spe. = spe[, !spe$condition %in% c("LN", "TON", "THY")]
spe. = spe.[, !spe.$cell_type_CellSighter %in% c("Unknown-Stroma")]
celltype_levels = c("Keratinocyte", "Endothelial", "Lymphatic", "Leukocyte",
                    "APC", "Monocytic_Lineage",  "Monocytes", "Macrophages", "Neutrophils",
                    "Basophil", "pDC", "B_cell", "NKT", "T_regulatory", "T_cytotoxic", "T_helper")
spe.$cell_type_CellSighter = factor(spe.$cell_type_CellSighter,
                                    levels = celltype_levels)
png(file.path(output, "Plots", "Sample x CellType x Condition Filtered.png"), width = 3500, height = 1600, res = 300)
colors =  setNames(unique(spe.$cell_type_CellSighter_color), unique(spe.$cell_type_CellSighter))
dittoBarPlot(spe., "cell_type_CellSighter", group.by = "sample_id", color.panel = colors[levels(spe.$cell_type_CellSighter)],  main = "Sample x Clusters",
             split.adjust = list(scales = 'free'), retain.factor.levels = TRUE, split.by = "condition", split.nrow = 1)
dev.off()

meta = as.data.frame(colData(spe))

meta %>% group_by(condition) %>% mutate(total = n()) %>%
  group_by(condition,cell_type_CellSighter) %>%
  summarise(percent = 100*n() / mean(total)) %>% filter(cell_type_CellSighter == "Unknown-Stroma")

for(i in unique(spe$cell_type_CellSighter)){
  png(file.path(output, "Plots", paste0("Sample x Condition x ",i,".png")), width = 3500, height = 1600, res = 300)
  spe. = spe
  spe.$cell_type_CellSighter[spe.$cell_type_CellSighter != i] = "Aa"
  print(dittoBarPlot(spe., "cell_type_CellSighter", group.by = "sample_id", color.panel = c(colors[i],"grey"),  main = "Sample x Clusters",
                     split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)
  )
  dev.off()
}

spe. = spe[,!spe$cell_type_CellSighter %in% c("Keratinocyte","Endothelial", "Unknown-Stroma", "Lymphatic")]
d = dittoBarPlot(spe., "cell_type_CellSighter", group.by = "sample_id", color.panel = unique(spe.$cell_type_CellSighter_color[order(spe.$cell_type_CellSighter)]), main = "Sample x Clusters", 
                 split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)

for(i in unique(spe.$cell_type_CellSighter)){
  png(file.path(output, "Plots", paste0("Sample x Condition x ",i," (Immune).png")), width = 3500, height = 1600, res = 300)
  spe.. = spe.
  spe..$cell_type_CellSighter[spe..$cell_type_CellSighter != i] = "Aa"
  print(dittoBarPlot(spe.., "cell_type_CellSighter", group.by = "sample_id", color.panel = c(colors[i],"grey"),  main = "Sample x Clusters",
                     split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)
  )
  dev.off()
}


spe. = spe
spe.$type = "NonImmune"
spe.$type[!spe.$cell_type_CellSighter %in% struct_celltype] = "Immune"
spe.$type = factor(spe.$type, levels = c( "NonImmune", "Immune"))
png(file.path(output, "Plots", "Sample x CellType x Condition - Immune vs Non Immune.png"), width = 3500, height = 1600, res = 300)
dittoBarPlot(spe., "type", group.by = "sample_id",  main = "Sample x Clusters", retain.factor.levels = T,
             split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)

dev.off()

png(file.path(output, "Plots", "Sample x CellType x Condition ImmuneCompartment.png"), width = 3500, height = 1600, res = 300)
colors =  setNames(unique(spe.$cell_type_CellSighter_color), unique(spe.$cell_type_CellSighter))
dittoBarPlot(spe., "cell_type_CellSighter", group.by = "sample_id", color.panel = colors[levels(d$data$label)],  main = "Sample x Clusters",
             split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)

dev.off()

spe. = spe[,!spe$cell_type_CellSighter %in% struct_celltype]
d = dittoBarPlot(spe., "cell_type_CellSighter", group.by = "sample_id", color.panel = unique(spe.$cell_type_CellSighter_color[order(spe.$cell_type_CellSighter)]), main = "Sample x Clusters", 
                 split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)

png(file.path(output, "Plots", "Sample x CellType x Condition ImmuneCompartment.png"), width = 3500, height = 1600, res = 300)
colors =  setNames(unique(spe.$cell_type_CellSighter_color), unique(spe.$cell_type_CellSighter))
dittoBarPlot(spe., "cell_type_CellSighter", group.by = "sample_id", color.panel = colors[levels(d$data$label)],  main = "Sample x Clusters",
             split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)

dev.off()

library(ggpubr)
for(i in unique(spe$cell_type_CellSighter)){
  spe. = spe[, !spe$condition %in% c("LN", "TON", "THY")]
  spe.$cell_type_CellSighter[spe.$cell_type_CellSighter != i] = "Aa"
  meta = as.data.frame(colData(spe.)) %>% dplyr::group_by(sample_id, condition, disease, batch) %>% 
    dplyr::summarise(percent_CN = 100 * length(
      which(cell_type_CellSighter == i)) / length(cell_type_CellSighter))
  
  png(file.path(output, "Plots", paste0("Boxplot_composition_",i,"_no_temoin_w_stats.png")), width = 1800, height = 1600, res = 300)
  p = grouped_dotplot(meta, y = "percent_CN", categ1 = "condition", categ2 = "disease",
                      ref.group = "HD", add_violin = F,
                      color_categ1 = setNames(unique(spe$condition_color), unique(spe$condition))
  )
  p = p + xlab("") + ylab(paste0("Percentage of ", gsub("_"," ",i), " (%)")) +
    guides(fill="none", color = "none")
  print(p)
  dev.off()
  
  # Immune cells
  spe. = spe[, !spe$cell_type_CellSighter %in% c(
    "Keratinocyte",
    "Lymphatic",
    "Endothelial",
    "Unknown-Stroma"
  )]
  spe.$cell_type_CellSighter[spe.$cell_type_CellSighter != i] = "Aa"
  
  meta = as.data.frame(colData(spe.)) %>% filter(!condition %in% c("LN", "TON", "THY")) %>%
    dplyr::group_by(sample_id, condition, disease, batch) %>% 
    dplyr::summarise(percent_CN = 100 * length(
      which(cell_type_CellSighter == i)) / length(cell_type_CellSighter))
  
  png(file.path(output, "Plots", paste0("Boxplot_composition_",i,"_no_temoin_w_stats_immune.png")), width = 1800, height = 1600, res = 300)
  p = grouped_dotplot(meta, y = "percent_CN", categ1 = "condition", categ2 = "disease",
                      ref.group = "HD", add_violin = F,
                      color_categ1 = setNames(unique(spe$condition_color), unique(spe$condition))
  )
  p = p + xlab("") + ylab(paste0("Percentage of ", gsub("_"," ",i), " (% immune cells)")) +
    guides(fill="none", color = "none")
  print(p)
  dev.off()
  
}

# Grouping Myeloid cells
spe. = spe. = spe[, !spe$cell_type_CellSighter %in% c(
  "Keratinocyte",
  "Lymphatic",
  "Endothelial",
  "Unknown-Stroma"
)]
myeloid = c("pDC", "Basophil", "Macrophages", "Monocytes", "Monocytic_Lineage", "Neutrophils")
spe.$cell_type_CellSighter[!spe.$cell_type_CellSighter %in% myeloid] = "Aa"
meta = as.data.frame(colData(spe.)) %>% filter(!condition %in% c("LN", "TON", "THY")) %>%
  dplyr::group_by(sample_id, condition, disease, batch) %>% 
  dplyr::summarise(percent_CN = 100 * length(
    which(cell_type_CellSighter %in% myeloid)) / length(cell_type_CellSighter))

png(file.path(output, "Plots", paste0("Boxplot_composition_Myeloid_no_temoin_w_stats_immune.png")), width = 1800, height = 1600, res = 300)
p = grouped_dotplot(meta, y = "percent_CN", categ1 = "condition", categ2 = "disease",
                    ref.group = "HD", add_violin = F,
                    color_categ1 = setNames(unique(spe$condition_color), unique(spe$condition))
)
p = p + xlab("") + ylab(paste0("Percentage of Myeloid cells (% immune cells)")) +
  guides(fill="none", color = "none")
print(p)
dev.off()

# Grouping Myeloid cells
spe. = spe. = spe[, !spe$cell_type_CellSighter %in% c(
  "Keratinocyte",
  "Lymphatic",
  "Endothelial",
  "Unknown-Stroma"
)]
lymphocyte = c("T_helper", "T_regulatory", "NKT", "T_cytotoxic")
spe.$cell_type_CellSighter[!spe.$cell_type_CellSighter %in% lymphocyte] = "Aa"

meta = as.data.frame(colData(spe.)) %>% filter(!condition %in% c("LN", "TON", "THY")) %>%
  dplyr::group_by(sample_id, condition, disease, batch) %>% 
  dplyr::summarise(percent_CN = 100 * length(
    which(cell_type_CellSighter %in% lymphocyte)) / length(cell_type_CellSighter))

png(file.path(output, "Plots", paste0("Boxplot_composition_Lymphoid_no_temoin_w_stats_immune.png")), width = 1800, height = 1600, res = 300)
p = grouped_dotplot(meta, y = "percent_CN", categ1 = "condition", categ2 = "disease",
                    ref.group = "HD", add_violin = F,
                    color_categ1 = setNames(unique(spe$condition_color), unique(spe$condition))
)
p = p + xlab("") + ylab(paste0("Percentage of Lymphocyte cells (% immune cells)")) +
  guides(fill="none", color = "none")
print(p)
dev.off() 

################################################################################
# Condition vs all else 
################################################################################
library(ggpubr)

# Composition
dir.create(file.path(output, "Plots", "One_vs_rest"))
for(cond in unique(spe$condition)){
  dir.create(file.path(output, "Plots", "One_vs_rest", cond))
  
  for(i in unique(spe$cell_type_CellSighter)){
    spe. =  spe
    name = "all"
    if(!(i %in% struct_celltype)){
      spe. = spe.[, !spe.$cell_type_CellSighter %in% struct_celltype]
      name = "immune"
    }
    spe.$cell_type_CellSighter[spe.$cell_type_CellSighter != i] = "Aa"
    
    meta = as.data.frame(colData(spe.)) %>%
      dplyr::group_by(sample_id, condition, batch) %>% 
      dplyr::summarise(percent_CN = 100 * length(
        which(cell_type_CellSighter == i)) / length(cell_type_CellSighter))  %>%
      mutate(condition = ifelse(condition == cond, cond, "Rest"))
    meta$condition = factor(meta$condition, levels = c("Rest", cond))
    png(file.path(output, "Plots", "One_vs_rest", cond, paste0("Boxplot_composition_",i,".png")), width = 1200, height = 1150, res = 300)
    p = grouped_dotplot(meta, y = "percent_CN", categ1 = "condition", categ2 = NULL,
                        ref.group = "Rest", add_violin = F,
                        color_categ1 = setNames(c("grey", unique(spe$condition_color)),
                                                c("Rest", as.character(unique(spe$condition)))))
    p = p + xlab("") + ylab(paste0("", gsub("_"," ",i), " (%)")) +
      guides(fill="none", color = "none")
    print(p)
    dev.off()
  }
}

##############################################################################
# Markers
##############################################################################
intensities = as.data.frame(t(spe@assays@data$CellSighter_marker_mat))
intensities$sample_id = spe$sample_id
intensities$condition = spe$condition
intensities$batch = spe$batch
intensities$cell_type_CellSighter = spe$cell_type_CellSighter
intensities = intensities %>% dplyr::group_by(sample_id, condition, batch, cell_type_CellSighter) %>% 
  dplyr::summarise(across(Actin:Ki67, mean))

colors_disease = setNames(unique(spe$condition_color), unique(spe$condition))
colors_samp = setNames(unique(spe$sample_id_color), unique(spe$sample_id))

for(cond in c("AD", "LP", "PS", "DAR", "MF", "SS")){
  print(cond)
  intensities. = intensities  %>%  mutate(condition = ifelse(condition == cond, cond, "Rest"))
  intensities.$condition = factor(intensities.$condition, levels = c("Rest", cond))
  
  for(celltype in c("T_helper", "T_cytotoxic", "T_regulatory", "B_cell", "Keratinocyte", "Macrophages")){
    print(celltype)
    
    pdf(file.path(output, "Plots", "One_vs_rest", cond, paste0("Functional_Marker_in_cell_type_CellSighter_",celltype,".pdf")), height = 4, width = 4)
    for(mark in colnames(intensities.)[5:65]){
      intensities.. = intensities. %>% filter(cell_type_CellSighter %in% celltype)
      intensities..[[mark]] = 100 * intensities..[[mark]]
      
      p = grouped_dotplot(intensities.., y =  mark, categ1 = "condition", categ2 = NULL,
                          ref.group = "Rest", add_violin = T,
                          color_categ1 = setNames(c("grey", unique(spe$condition_color)),
                                                  c("Rest", as.character(unique(spe$condition))))
      )
      p = p + xlab("") + ylab(paste0("Percentage of ", celltype, " expressing ", mark, " (%)")) +
        guides(fill="none", color = "none")
      print(p)
    }
    dev.off()
  }
  
}


##############################################################################
# Markers - normcounts
##############################################################################
intensities = as.data.frame(t(spe@assays@data$normcounts))
intensities$sample_id = spe$sample_id
intensities$condition = spe$condition
intensities$batch = spe$batch
intensities$cell_type_CellSighter = spe$cell_type_CellSighter
intensities = intensities %>% dplyr::group_by(sample_id, condition, batch, cell_type_CellSighter) %>% 
  dplyr::summarise(across(Actin:Ki67, mean))

colors_disease = setNames(unique(spe$condition_color), unique(spe$condition))
colors_samp = setNames(unique(spe$sample_id_color), unique(spe$sample_id))
dir.create(file.path(output, "Plots", "One_vs_rest_normcounts"))

for(cond in c("AD", "LP", "PS", "DAR", "MF", "SS")){
  print(cond)
  intensities. = intensities  %>%  mutate(condition = ifelse(condition == cond, cond, "Rest"))
  intensities.$condition = factor(intensities.$condition, levels = c("Rest", cond))
  
  for(celltype in c("T_helper", "T_cytotoxic", "T_regulatory", "B_cell", "Keratinocyte", "Macrophages")){
    print(celltype)
    dir.create(file.path(output, "Plots", "One_vs_rest_normcounts", cond))
    pdf(file.path(output, "Plots", "One_vs_rest_normcounts", cond, paste0("Functional_Marker_in_cell_type_CellSighter_",celltype,".pdf")), height = 4, width = 4)
    for(mark in colnames(intensities.)[5:65]){
      intensities.. = intensities. %>% filter(cell_type_CellSighter %in% celltype)
      intensities..[[mark]] = 100 * intensities..[[mark]]
      
      p = grouped_dotplot(intensities.., y =  mark, categ1 = "condition", categ2 = NULL,
                          ref.group = "Rest", add_violin = T,
                          color_categ1 = setNames(c("grey", unique(spe$condition_color)),
                                                  c("Rest", as.character(unique(spe$condition))))
      )
      p = p + xlab("") + ylab(paste0("Percentage of ", celltype, " expressing ", mark, " (%)")) +
        guides(fill="none", color = "none")
      print(p)
    }
    dev.off()
  }
  
}

##############################################################################
# Doghnut charts of the main cell types found 
##############################################################################


spe. = spe

df = as.data.frame(colData(spe.)) %>% group_by(condition, cell_type_simplified) %>%
  summarise(n_celltype =  n())

for(i in unique(spe.$condition)){
  data = df %>% filter(condition == i)
  
  # Compute percentages
  data$fraction <- data$n_celltype / sum(data$n_celltype)
  
  # Compute the cumulative percentages (top of each rectangle)
  data$ymax <- cumsum(data$fraction)
  
  # Compute the bottom of each rectangle
  data$ymin <- c(0, head(data$ymax, n=-1))
  
  # Compute label position
  data$labelPosition <- (data$ymax + data$ymin) / 2
  
  # Compute a good label
  # data$label <- paste0(data$cell_type_simplified, "\n (", round(100*data$fraction,1) ,"%)")
  data$label <- paste0(round(100*data$fraction,1), "%")
  
  pdf(file.path(output, "Plots", paste0("Doughnut_",i,"_nolabel.pdf")), width = 5, height =5)
  # Make the plot
  print(ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=2, fill=cell_type_simplified)) +
          geom_rect() +
          geom_text(x=5, aes(y=labelPosition, label=label), color="black", size=10) + # x here controls label position (inner / outer)
          scale_fill_manual(values = cell_type_simplified_color_df) +
          scale_color_manual(values = cell_type_simplified_color_df) +
          annotate("text", x = -1, y = 0, label =paste0(i), size = 14) +
          coord_polar(theta="y") +
          xlim(c(-1, 4)) +
          theme_void() +
          theme(legend.position = "none"))
  dev.off()
}

# Calculate p-values :
df = as.data.frame(colData(spe.)) %>% group_by(condition, sample_id, disease, batch, cell_type_simplified) %>%
  summarise(n_celltype =  n())
df = df %>% group_by(sample_id) %>% mutate(total = sum(n_celltype), fraction = n_celltype/total) %>%
  filter( !condition %in% c("LN", "THY", "TON"))

for(i in unique(df$cell_type_simplified)){
  library(ggpubr)
  
  pdf(file.path(output, "Plots", paste0("Composition_simplified_",i,".pdf")), width = 5, height =5)
  p = grouped_dotplot(df %>% filter(cell_type_simplified == i), y = "fraction", categ1 = "condition", categ2 = "disease",
                      ref.group = "HD", add_violin = F,
                      color_categ1 = setNames(unique(spe$condition_color), unique(spe$condition))
  )
  p = p + xlab("") + ylab(paste0("Percentage of ", gsub("_", " ", i), " cells (% immune cells)")) +
    guides(fill="none", color = "none")
  print(p)
  dev.off()
}


##############################################################################
# Normcounts x Celltype Heatmap
##############################################################################
spe. = spe[,!spe$condition %in% c("THY", "LN", "TON")]

cells = c()
for(i in unique(spe.$cell_type_CellSighter)){
  cells = c(cells, sample(spe.$cell_id[which(spe.$cell_type_CellSighter == i)], min(500, length(which(spe.$cell_type_CellSighter == i))) ))
}
spe. = spe.[,match(cells, spe.$cell_id)]

png(file.path("output", "CellSighter", "Plots", paste0("Heatmap_Celltype_Clusters_counts_QN.png")),
    width = 2500,
    height = 1800, res = 300)
set.seed(47)
mat = normcounts(spe.)
celltype_order = c( "Endothelial", "Lymphatic",
             "B_cell", "NKT",  "pDC", "Basophil",
             "Macrophages", "Monocytes", "Neutrophils",
                "T_cytotoxic", "T_regulatory", "T_helper",
             "Keratinocyte",  "Monocytic_Lineage",
                "APC", "Leukocyte", "Unknown-Stroma")
cells = spe.$cell_id[order(match(spe.$cell_type_CellSighter, celltype_order))]

mat = t(mat)[cells, c(golden_markers, "CollagenI")]

h = Heatmap(
  scale(mat),
  name = "Quantile Normalized Counts",
  column_title = "Markers", 
  row_title = "Celltype",
  row_dend_side = "right",
  row_dend_width = unit(4, "cm"),
  column_dend_height = unit(5, "cm"),
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  show_column_names = TRUE,
  show_row_names = FALSE,
  border = TRUE,
  row_names_gp = gpar(fontsize = 8),
  use_raster = FALSE,
  right_annotation = rowAnnotation(
    "Celltype" = spe.$cell_type_CellSighter[match(cells, spe.$cell_id)],
    col = list("Celltype" = setNames(
      unique(spe.$cell_type_CellSighter_color),
      unique(spe.$cell_type_CellSighter)
    )))
)
draw(h)
dev.off()  

##############################################################################
# Correlation heatmap based on celltype density correlations  - Immune
##############################################################################
library(dittoSeq)
density = read.csv("output/pixie/region_pixel_output_dir/areas_per_sample.csv")
meta = as.data.frame(colData(spe))
meta$area = density$tissue_area[match(meta$sample_id, density$sample_id)]

meta. = meta  %>%
  dplyr::group_by(sample_id, cell_type_CellSighter, condition, batch, Th_group) %>% 
  dplyr::summarise(density = n() / mean(area)) 

prop = meta. %>% filter(!condition %in% c("THY", "TON", "LN"))
prop = pivot_wider(prop, names_from = cell_type_CellSighter, values_from = density, values_fill = 0)

mat = as.matrix(as.data.frame(prop[,5:ncol(prop)]))
rownames(mat) = prop$sample_id
cormat = cor(t(mat), method = "pearson")

hc = hclust(1-dist(cormat), method = "ward.D2")

df = data.frame("sample_id" = prop$sample_id,
                "condition" = prop$condition,
                "batch" = prop$batch,
                "Th_group" = prop$Th_group
)
list_colors = list(
  "sample_id" = setNames(unique(spe$sample_id_color), unique(spe$sample_id)),
  "condition" = setNames(unique(spe$condition_color), unique(spe$condition)),
  "batch" = setNames(unique(spe$batch_color), unique(spe$batch)),
  "Th_group" = setNames(unique(spe$Th_group_color), unique(spe$Th_group))
)

pdf(file.path("output", "CellSighter", "Plots", paste0("Heatmap_correlations_of_density.pdf")),
    width = 13,
    height = 12)
h = Heatmap(
  cormat[hc$order, hc$order],
  name = "Correlations of density",
  column_title = "Celltype",
  row_title = "Samples",
  cluster_rows = F,
  cluster_columns = F,
  row_dend_side = "right",
  # col =c("royalblue4", "gold"),
  show_column_names = TRUE, 
  show_row_names = FALSE, 
  border = TRUE,
  column_names_gp = gpar(fontsize = 6),
  use_raster = FALSE,
  right_annotation = rowAnnotation(df = df[hc$order, ],
                                   col = list_colors
  )
)
draw(h)
dev.off()

# Dittoplot, grouped by correlation of celltype proportions
spe. = spe[, !spe$condition %in% c("LN", "TON", "THY")]
spe. = spe.[, !spe.$cell_type_CellSighter %in% c(struct_celltype)]
spe.$cell_type_CellSighter = factor(spe.$cell_type_CellSighter,
                                    levels = celltype_levels)
spe.$sample_id = factor(spe.$sample_id, levels = hc$labels[hc$order])
spe.$area = density$tissue_area[match(spe.$sample_id, density$sample_id)]

png(file.path(output, "Plots", "Sample x CellType x Condition Filtered - grouped by density correlations.png"), width = 3500, height = 1600, res = 300)
colors =  setNames(unique(spe.$cell_type_CellSighter_color), unique(spe.$cell_type_CellSighter))
dittoBarPlot(spe., "cell_type_CellSighter", group.by = "sample_id", color.panel = colors[levels(spe.$cell_type_CellSighter)],  main = "Sample x Clusters",
             split.adjust = list(scales = 'free'), retain.factor.levels = TRUE)
dev.off()

png(file.path(output, "Plots", 
              "Sample x CellType x Condition Filtered - density - bar.png"), width = 3500, height = 800, res = 300)
barplot(rep(1, length(hc$labels)), border = NA, col = spe$condition_color[match(hc$labels[hc$order], spe$sample_id)])
dev.off()

pdf(file.path(output, "Plots", 
              "Sample x CellType x Condition Filtered - grouped by density correlations - hc.pdf"), width = 10, height = 3)
plot(hc, cex = 0.5, main = "", ylab = "", xlab = "")
dev.off()


##############################################################################
# Correlation heatmap  based on celltype proportions correlations - Immune
##############################################################################
spe. = spe[, !spe$condition %in% c("LN", "TON", "THY")]
spe. = spe.[, !spe.$cell_type_CellSighter %in% c(struct_celltype, "Leukocyte")]
# spe.$cell_type_CellSighter[which(spe.$cell_type_CellSighter %in%
#                                    c("Monocytes", "Monocytic_Lineage"))] =
#   "Mono"
spe.$cell_type_CellSighter = factor(spe.$cell_type_CellSighter, levels = celltype_levels)

meta. = as.data.frame(colData(spe.))  %>%
    dplyr::group_by(sample_id, condition, batch, Th_group) %>% 
  mutate(n_all = n()) %>%
  dplyr::group_by(sample_id, cell_type_CellSighter, condition, batch, Th_group) %>% 
    dplyr::summarise(percent_CN = 100 * n() / mean(n_all))
prop = meta. # %>% filter(!condition %in% c("THY", "TON", "LN"))
prop = pivot_wider(prop, names_from = cell_type_CellSighter, values_from = percent_CN, values_fill = 0)

mat = as.matrix(as.data.frame(prop[,5:ncol(prop)]))
rownames(mat) = prop$sample_id
cormat = cor(t(mat), method = "pearson")

hc = hclust(1-dist(cormat), method = "ward.D2")

df = data.frame("sample_id" = prop$sample_id,
                "condition" = prop$condition,
                "batch" = prop$batch,
                "Th_group" = prop$Th_group
)
list_colors = list(
  "sample_id" = setNames(unique(spe$sample_id_color), unique(spe$sample_id)),
  "condition" = setNames(unique(spe$condition_color), unique(spe$condition)),
  "batch" = setNames(unique(spe$batch_color), unique(spe$batch)),
  "Th_group" = setNames(unique(spe$Th_group_color), unique(spe$Th_group))
)

pdf(file.path("output", "CellSighter", "Plots", paste0("Heatmap_correlations_of_cell_proportions.pdf")),
    width = 13,
    height = 12)
h = Heatmap(
  cormat,
  name = "Correlations of proportions",
  column_title = "Celltype",
  row_title = "Samples",
  cluster_rows = hc,
  cluster_columns = hc,
  row_dend_side = "right",
  # col =c("royalblue4", "gold"),
  show_column_names = TRUE, 
  show_row_names = FALSE, 
  border = TRUE,
  column_names_gp = gpar(fontsize = 6),
  use_raster = FALSE,
  right_annotation = rowAnnotation(df = df,
                                   col = list_colors
  )
)
draw(h)
dev.off()

# Dittoplot, grouped by correlation of celltype proportions
spe.$sample_id = factor(spe.$sample_id, levels = hc$labels[hc$order])

png(file.path(output, "Plots", "Sample x CellType x Condition Filtered - grouped by correlations.png"), width = 3500, height = 1600, res = 300)
colors =  setNames(unique(spe.$cell_type_CellSighter_color), unique(spe.$cell_type_CellSighter))
dittoBarPlot(spe., "cell_type_CellSighter", group.by = "sample_id",
             color.panel = colors[levels(spe.$cell_type_CellSighter)],  main = "Sample x Clusters",
             split.adjust = list(scales = 'free'), retain.factor.levels = TRUE)
dev.off()

png(file.path(output, "Plots", 
              "Sample x CellType x Condition Filtered - bar.png"), width = 3500, height = 800, res = 300)
barplot(rep(1, length(hc$labels)), border = NA, col = spe$condition_color[match(hc$labels[hc$order], spe$sample_id)])
dev.off()

clust = cutree(hc, h = 2)
set.seed(47)
cols = setNames(discrette_colors_50[1:length(unique(clust))], 1:6)

png(file.path(output, "Plots", 
              "Sample x CellType x Condition Filtered - bar - cluster.png"), width = 3500, height = 800, res = 300)
barplot(rep(1, length(hc$labels)), border = NA,
        col = cols[match(clust[hc$order], names(cols))])
dev.off()

pdf(file.path(output, "Plots", 
              "Sample x CellType x Condition Filtered - grouped by correlations - hc.pdf"), width = 10, height = 3)
plot(hc, cex = 0.5, main = "", ylab = "", xlab = "")
dev.off()



##############################################################################
# Correlation heatmap  based on celltype proportions correlations - ALL
##############################################################################
spe. = spe[, !spe$condition %in% c("LN", "TON", "THY")]
spe.$cell_type_CellSighter = factor(spe.$cell_type_CellSighter, levels = celltype_levels)

meta. = as.data.frame(colData(spe.))  %>%
  dplyr::group_by(sample_id, condition, batch, Th_group) %>% 
  mutate(n_all = n()) %>%
  dplyr::group_by(sample_id, cell_type_CellSighter, condition, batch, Th_group) %>% 
  dplyr::summarise(percent_CN = 100 * n() / mean(n_all))
prop = meta. 
prop = pivot_wider(prop, names_from = cell_type_CellSighter, values_from = percent_CN, values_fill = 0)

mat = as.matrix(as.data.frame(prop[,5:ncol(prop)]))
rownames(mat) = prop$sample_id
cormat = cor(t(mat), method = "pearson")

hc = hclust(1-dist(cormat), method = "ward.D2")

df = data.frame("sample_id" = prop$sample_id,
                "condition" = prop$condition,
                "batch" = prop$batch,
                "Th_group" = prop$Th_group
)
list_colors = list(
  "sample_id" = setNames(unique(spe$sample_id_color), unique(spe$sample_id)),
  "condition" = setNames(unique(spe$condition_color), unique(spe$condition)),
  "batch" = setNames(unique(spe$batch_color), unique(spe$batch)),
  "Th_group" = setNames(unique(spe$Th_group_color), unique(spe$Th_group))
)

pdf(file.path("output", "CellSighter", "Plots", paste0("Heatmap_correlations_of_cell_proportions_ALL.pdf")),
    width = 13,
    height = 12)
h = Heatmap(
  cormat,
  name = "Correlations of proportions",
  column_title = "Celltype",
  row_title = "Samples",
  cluster_rows = hc,
  cluster_columns = hc,
  row_dend_side = "right",
  # col =c("royalblue4", "gold"),
  show_column_names = TRUE, 
  show_row_names = FALSE, 
  border = TRUE,
  column_names_gp = gpar(fontsize = 6),
  use_raster = FALSE,
  right_annotation = rowAnnotation(df = df,
                                   col = list_colors
  )
)
draw(h)
dev.off()

# Dittoplot, grouped by correlation of celltype proportions
spe.$sample_id = factor(spe.$sample_id, levels = hc$labels[hc$order])

png(file.path(output, "Plots", "Sample x CellType x Condition Filtered - grouped by correlations - ALL.png"), width = 3500, height = 1600, res = 300)
colors =  setNames(unique(spe.$cell_type_CellSighter_color), unique(spe.$cell_type_CellSighter))
dittoBarPlot(spe., "cell_type_CellSighter", group.by = "sample_id",
             color.panel = colors[levels(spe.$cell_type_CellSighter)],  main = "Sample x Clusters",
             split.adjust = list(scales = 'free'), retain.factor.levels = TRUE)
dev.off()

png(file.path(output, "Plots", 
              "Sample x CellType x Condition Filtered - bar - ALL.png"), width = 3500, height = 800, res = 300)
barplot(rep(1, length(hc$labels)), border = NA, col = spe$condition_color[match(hc$labels[hc$order], spe$sample_id)])
dev.off()

pdf(file.path(output, "Plots", 
              "Sample x CellType x Condition Filtered - grouped by correlations - hc - ALL.pdf"), width = 10, height = 3)
plot(hc, cex = 0.5, main = "", ylab = "", xlab = "")
dev.off()

# Early vs late
meta = as.data.frame(colData(spe[,which(spe$sample_id %in% hc$labels)]))

colors_early_late = setNames(c("#EDCC5F", "#941D1B"), c("Early", "Late"))

meta$stage_color = "grey"
meta$stage_color[meta$stage == "Early"] = colors_early_late[1]
meta$stage_color[meta$stage == "Late"] = colors_early_late[2]

meta = meta %>% group_by(sample_id, condition, condition_color, stage, stage_color) %>% 
  summarise(total = n())

labels = hc$labels[hc$order]

png(file.path(output, "Plots", "stage_bar_for_cellcomposition_cor.png"), height = 200, width = 1600)
barplot(rep(1, nrow(meta)), col = meta$stage_color[match(labels, meta$sample_id)])
dev.off()

#################################################################################
# Marker expression
#################################################################################

spe. = spe[, !spe$condition %in% c("LN", "TON", "THY")]
intensities = as.data.frame(t(spe.@assays@data$CellSighter_marker_mat))
intensities$sample_id = spe.$sample_id
intensities$condition = spe.$condition
intensities$disease = spe.$disease
intensities$batch = spe.$batch
intensities$cell_type_CellSighter = spe.$cell_type_CellSighter
intensities = intensities %>% tidyr::pivot_longer(-c(sample_id,condition, cell_type_CellSighter, disease, batch), names_to = c("marker"))
intensities$condition = factor(intensities$condition, levels = colCondition$condition)

for(celltype in unique(spe$cell_type_CellSighter)){
  colors_cond = setNames(unique(spe$condition_color), unique(spe$condition))
  colors_samp = setNames(unique(spe$sample_id_color), unique(spe$sample_id))
  
  pdf(file.path(output, "Plots", paste0("Functional_Marker_in_cell_type_CellSighter_",celltype,".pdf")), height = 20, width = 17)
  intensities. = intensities %>% filter(cell_type_CellSighter %in% c(celltype)) %>%
    dplyr::group_by(sample_id, condition, disease, batch, marker) %>% 
    dplyr::summarise(Percent = 100 * sum(value) / dplyr::n())
  
  p = grouped_dotplot(intensities., y = "Percent", categ1 = "condition",,
                      ref.group = "HD", add_violin = F,
                      color_categ1 = setNames(unique(spe$condition_color), unique(spe$condition))
  )
  p = p + facet_wrap(~marker, ncol = 5) +  xlab("") +
    ylab(paste0("% of ", celltype , " cells positive for marker")) +
    guides(fill="none", color = "none") + ggtitle(celltype) + ylim(c(0,120))
  print(p)
  
  dev.off()
}



#################################################################################
# Heatmap of distribution-based markers vs celltype predictions (celltype markers)
#################################################################################

intensities = as.data.frame(t(spe@assays@data$positive_marker))
intensities$sample_id = spe$sample_id
intensities$condition = spe$condition
intensities$cell_type_CellSighter = spe$cell_type_CellSighter
intensities$cell_id = rownames(intensities)
intensities = intensities %>% arrange(factor(cell_type_CellSighter, levels = gsub("Fibroblast","Unknown-Stroma",unique(cell_type_markers$cell_type))))

intensities. = data.frame()
set.seed(47)
for(celltype in unique(cell_type_markers$cell_type) ){
  df = intensities[which(intensities$cell_type_CellSighter == celltype),]
  intensities. = rbind(intensities., df[sample(nrow(df), min(1000, nrow(df))),])
}

mat = as.matrix(intensities.[, colnames(intensities.) %in% golden_markers])
mat = mat[, golden_markers]
mat[which(is.na(mat))] = 0
mat = mat[-which(rowSums(mat) == 0),]

df = data.frame("Celltype" = factor(intensities.$cell_type_CellSighter[match(rownames(mat), intensities.$cell_id)], levels = gsub("Fibroblast","Unknown-Stroma",unique(cell_type_markers$cell_type))))

pdf(file.path("output", "CellSighter", "Plots", paste0("Heatmap_celltype_vs_markers_Classic.pdf")),
    width = 13,
    height = 12)
h = Heatmap(
  mat,
  name = "CellSighter Celltype vs Markers",
  column_title = "Celltype (DeepLearning)",
  row_title = "Positive Markers (Histogram-based)",
  cluster_rows = F,
  column_order = 1:ncol(mat),
  row_dend_side = "right",
  col =c("royalblue4", "gold"),
  show_column_names = TRUE, 
  show_row_names = FALSE, 
  border = TRUE,
  column_names_gp = gpar(fontsize = 6),
  use_raster = TRUE,
  right_annotation = rowAnnotation(df = df,
                                   col = list("Celltype" = setNames(unique(spe$cell_type_CellSighter_color[order(spe$cell_type_CellSighter)]),
                                                                    unique(spe$cell_type_CellSighter[order(spe$cell_type_CellSighter)])
                                   )))
)
draw(h)
dev.off()

#################################################################################
# Heatmap of distribution-based markers vs celltype predictions (all markers)
#################################################################################

mat = as.matrix(intensities.[, colnames(intensities.) %in% rownames(spe)])
mat = mat[, setdiff(colnames(mat), markers_to_remove)]
mat[which(is.na(mat))] = 0
mat = mat[-which(rowSums(mat) == 0),]

df = data.frame("Celltype" = factor(intensities.$cell_type_CellSighter[match(rownames(mat), intensities.$cell_id)], levels = gsub("Fibroblast","Unknown-Stroma",unique(cell_type_markers$cell_type))))

pdf(file.path("output", "CellSighter", "Plots", paste0("Heatmap_celltype_vs_markers_Classic_all.pdf")),
    width = 16,
    height = 12)
h = Heatmap(
  mat,
  name = "CellSighter Celltype vs Markers",
  column_title = "Celltype (DeepLearning)",
  row_title = "Positive Markers (Histogram-based)",
  cluster_rows = F,
  cluster_columns = T,
  row_dend_side = "right",
  col =c("royalblue4", "gold"),
  show_column_names = TRUE, 
  show_row_names = FALSE, 
  border = TRUE,
  column_names_gp = gpar(fontsize = 6),
  use_raster = TRUE,
  right_annotation = rowAnnotation(df = df,
                                   col = list("Celltype" = setNames(unique(spe$cell_type_CellSighter_color[order(spe$cell_type_CellSighter)]),
                                                                    unique(spe$cell_type_CellSighter[order(spe$cell_type_CellSighter)])
                                   )))
)
draw(h)
dev.off()


#################################################################################
# Heatmap of deeplearning-based markers vs celltype predictions (celltype markers)
#################################################################################

cell_type_markers = read.csv("annotation/cell_markers.csv")
intensities = as.data.frame(t(spe@assays@data$CellSighter_marker_mat))
intensities$sample_id = spe$sample_id
intensities$condition = spe$condition
intensities$cell_type_CellSighter = spe$cell_type_CellSighter
intensities$cell_id = rownames(intensities)
intensities = intensities %>% arrange(factor(cell_type_CellSighter))
intensities$cell_type_CellSighter_color = spe$cell_type_CellSighter

intensities. = data.frame()
set.seed(47)
for(celltype in unique(cell_type_markers$cell_type)){
  df = intensities[which(intensities$cell_type_CellSighter == celltype),]
  intensities. = rbind(intensities., df[sample(nrow(df), min(500, nrow(df))),])
}

# Plotting heatmaps CellSighter Cell Type vs CellSighter Marker
mat = as.matrix(intensities.[, colnames(intensities.) %in% cell_type_markers$marker])
mat = mat[, cell_type_markers$marker]
mat[which(is.na(mat))] = 0
mat = mat[-which(rowSums(mat) == 0),]

mat = mat[,golden_markers]

pdf(file.path("output", "CellSighter", "Plots", paste0("Heatmap_celltype_vs_markers_CellSighter.pdf")),
    width = 13,
    height = 12)
color_df = data.frame("Celltype" = (intensities$cell_type_CellSighter[match(rownames(mat), rownames(intensities))]))

h = Heatmap(
  mat,
  name = "Positive marker",
  column_title = "Positive Markers (DeepLearning)",
  row_title = "Celltype (DeepLearning)",
  cluster_rows = F,
  column_order = 1:ncol(mat),
  row_dend_side = "right",
  col =c("royalblue4",  "royalblue4", "gold"),
  show_column_names = TRUE, 
  show_row_names = FALSE, 
  border = TRUE,
  column_names_gp = gpar(fontsize = 9),
  use_raster = FALSE,
  right_annotation = rowAnnotation(df = color_df, 
                                   col = list(Celltype = setNames(unique(spe$cell_type_CellSighter_color), 
                                                                  unique(spe$cell_type_CellSighter))))
)
draw(h)
dev.off()


#################################################################################
# Heatmap of deeplearning-based markers vs celltype predictions (all markers)
#################################################################################
mat = as.matrix(intensities.[, colnames(intensities.) %in% rownames(spe)])
mat = mat[, setdiff(colnames(mat), markers_to_remove)]
mat[which(is.na(mat))] = 0
mat = mat[-which(rowSums(mat) == 0),]

pdf(file.path("output", "CellSighter", "Plots", paste0("Heatmap_celltype_vs_markers_CellSighter_all.pdf")),
    width = 16,
    height = 12)
color_df = data.frame(
  "Celltype" = (intensities.$cell_type_CellSighter[match(rownames(mat),rownames(intensities.))])
)

h = Heatmap(
  mat,
  name = "Positive marker",
  column_title = "Positive Markers (DeepLearning)",
  row_title = "Celltype (DeepLearning)",
  cluster_rows = F,
  cluster_columns = T,
  row_dend_side = "right",
  col =c("royalblue4",  "royalblue4", "gold"),
  show_column_names = TRUE, 
  show_row_names = FALSE, 
  border = TRUE,
  column_names_gp = gpar(fontsize = 9),
  use_raster = TRUE,
  right_annotation = rowAnnotation(df = color_df, 
                                   col = list(Celltype = setNames(unique(spe$cell_type_CellSighter_color), 
                                                                  unique(spe$cell_type_CellSighter))))
)
draw(h)
dev.off()

#################################################################################
# Heatmap of normalized intensitiy of markers vs celltype predictions (all markers)
#################################################################################
intensities = as.data.frame(t(spe@assays@data$normcounts))
intensities$sample_id = spe$sample_id
intensities$condition = spe$condition
intensities$cell_type_CellSighter = spe$cell_type_CellSighter
intensities$cell_id = rownames(intensities)
intensities = intensities %>% arrange(factor(cell_type_CellSighter))
intensities$cell_type_CellSighter_color = spe$cell_type_CellSighter_color

intensities. = intensities %>% group_by(cell_type_CellSighter) %>%
  summarise_at(rownames(spe), mean)

mat = as.matrix(intensities.[,-1])
rownames(mat) = intensities.$cell_type_CellSighter
mat = mat[unique(spe$cell_type_CellSighter),]
mat = scale(mat)

pdf(file.path(output, "Plots", paste0("Heatmap_celltype_vs_markers_normcounts_all.pdf")),
    width = 16,
    height = 12)
hc_cor_marker = hclust(as.dist(1 - cor(mat)), method = "ward.D2")

color_df = data.frame("Celltype" = rownames(mat))

h = Heatmap( 
  mat,
  cluster_columns = hc_cor_marker,
  cluster_rows = F,
  name = paste0("Protein Normalized Intensity"),
  column_title = "Celltype (DeepLearning)",
  row_title = "Proteins",
  row_dend_side = "right",
  row_dend_width = unit(4, "cm"),
  column_dend_height = unit(5, "cm"),
  show_column_names = T,
  show_row_names = T,
  clustering_distance_columns ="pearson",
  clustering_distance_rows = "pearson",
  row_split = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),
  column_split = 10,
  border = TRUE,
  use_raster = FALSE,
  right_annotation = rowAnnotation(df = color_df, 
                                   col = list(Celltype = setNames(unique(spe$cell_type_CellSighter_color), 
                                                                  unique(spe$cell_type_CellSighter))))
)
draw(h)
dev.off()

#################################################################################
# Heatmap of normalized intensitiy of markers vs celltype predictions (celltype markers)
#################################################################################
mat = as.matrix(intensities.[,golden_markers])
rownames(mat) = intensities.$cell_type_CellSighter
mat = mat[ gsub("Fibroblast","Unknown-Stroma",unique(cell_type_markers$cell_type)),]
mat = scale(mat)

pdf(file.path(output, "Plots", paste0("Heatmap_celltype_vs_markers_normcounts.pdf")),
    width = 13,
    height = 12)
h = Heatmap( 
  mat,
  cluster_columns = F,
  cluster_rows = F,
  name = paste0(" Protein Normalized Intensity"),
  column_title = "Celltype (DeepLearning)",
  row_title = "Proteins",
  row_dend_side = "right",
  row_dend_width = unit(4, "cm"),
  column_dend_height = unit(5, "cm"),
  show_column_names = T,
  show_row_names = T,
  clustering_distance_columns ="pearson",
  clustering_distance_rows = "pearson",
  row_split = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),
  column_split = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),
  border = TRUE,
  use_raster = FALSE
)
draw(h)
dev.off()

#################################################################################
# Heatmap of normalized intensitiy of markers vs celltype predictions 
# per condition (all markers)
#################################################################################
pdf(file.path(output, "Plots", paste0("Heatmap_celltype_vs_markers_normcounts_per_condition.pdf")),
    width = 16,
    height = 10)

intensities. = intensities %>%
  group_by(condition) %>%
  summarise_at(rownames(spe), mean)
intensities. = intensities. %>% filter(!condition %in% c("THY", "TON", "LN") )

mat = as.matrix(intensities.[,-1])
rownames(mat) = intensities.$condition
mat = mat[,setdiff(colnames(mat), markers_to_remove)]
mat = scale(mat)
hc_cor_marker = hclust(as.dist(1 - cor(mat)), method = "ward.D2")
hc_cor_sample = hclust(as.dist(1 - cor(t(mat))), method = "ward.D2")

color_df = data.frame("condition" = rownames(mat))
h = Heatmap( 
  mat,
  cluster_columns = hc_cor_marker,
  cluster_rows = hc_cor_sample,
  name = paste0("Condition"),
  column_title = "Proteins",
  row_title = "Conditions",
  row_dend_side = "right",
  row_dend_width = unit(4, "cm"),
  column_dend_height = unit(5, "cm"),
  show_column_names = T,
  show_row_names = T,
  clustering_distance_columns ="pearson",
  clustering_distance_rows = "pearson",
  column_split = 10,
  border = TRUE,
  use_raster = FALSE,
  right_annotation = rowAnnotation(df = color_df, 
                                   col = list(condition = setNames(unique(spe$condition_color), 
                                                                   unique(spe$condition))))
  
)
draw(h)
dev.off()

#################################################################################
# Heatmap of normalized intensitiy of markers vs celltype predictions 
# per condition per celltype (all markers)
#################################################################################

pdf(file.path(output, "Plots", paste0("Heatmap_celltype_vs_markers_normcounts_per_condition_per_celltype.pdf")),
    width = 16,
    height = 12)

# Global hcor
intensities. = intensities %>%
  group_by(condition) %>%
  summarise_at(rownames(spe), mean)
intensities. = intensities. %>% filter(!condition %in% c("THY", "TON", "LN") )
mat = spe.@assays@data$normcounts
mat = t(mat)

for(celltype in unique(spe$cell_type_CellSighter)){
  
  mat. = mat[spe.$cell_id[spe.$cell_type_CellSighter == celltype],]
  mat. = mat.[,setdiff(colnames(mat), markers_to_remove)]
  
  
  intensities. = as.data.frame(mat.) 
  intensities.$condition = spe$condition[match(rownames(intensities.), spe$cell_id)] 
  
  intensities. = intensities. %>%
    group_by(condition) %>%
    summarise_at(colnames(mat.), mean)
  intensities. = intensities. %>% filter(!condition %in% c("THY", "TON", "LN") )
  
  hc_cor_markers = hclust(as.dist(1 - cor(mat.)), method = "ward.D2")
  hc_cor_sample = hclust(as.dist(1 - cor(t(as.matrix(intensities.[,-1])))), method = "ward.D2")
  
  
  color_df = data.frame("condition" = intensities.$condition)
  mat.. = as.matrix(intensities.[,-1])
  mat.. = scale(mat..)
  rownames(mat..) = intensities.$condition
  h = Heatmap( 
    mat..,
    cluster_columns = hc_cor_markers,
    cluster_rows = hc_cor_sample,
    name = paste0(celltype),
    column_title = celltype,
    row_title = "Conditions",
    row_dend_side = "right",
    row_dend_width = unit(4, "cm"),
    column_dend_height = unit(5, "cm"),
    show_column_names = TRUE,
    show_row_names = TRUE,
    clustering_distance_columns ="pearson",
    clustering_distance_rows = "pearson",
    column_split = 10,
    border = TRUE,
    use_raster = FALSE,
    right_annotation = rowAnnotation(df = color_df, 
                                     col = list(condition = setNames(unique(spe.$condition_color), 
                                                                     unique(spe.$condition))))
    
  )
  draw(h)
}
dev.off()

#################################################################################
# Heatmap of fraction of positive deeplearning-based marker vs celltype predictions 
# per condition per celltype (all markers)
#################################################################################
intensities = as.data.frame(t(spe@assays@data$CellSighter_marker_mat))
intensities$sample_id = spe$sample_id
intensities$condition = spe$condition
intensities$cell_type_CellSighter = spe$cell_type_CellSighter
intensities$cell_id = rownames(intensities)
intensities = intensities %>% arrange(factor(cell_type_CellSighter))
intensities$cell_type_CellSighter_color = spe$cell_type_CellSighter_color


pdf(file.path(output, "Plots", paste0("Heatmap_celltype_vs_markers_CellSighter_per_condition_per_celltype.pdf")),
    width = 16,
    height = 12)

for(celltype in setdiff(unique(spe$cell_type_CellSighter), "Leukocyte")){
  
  intensities. = intensities %>% filter(cell_type_CellSighter == celltype)
  all_cells = c()
  for(i in unique(intensities.$sample_id)){
    int_samp = intensities.[intensities.$sample_id == i,]
    all_cells = c(all_cells, int_samp$cell_id[sample(1:nrow(int_samp), min(100, nrow(int_samp)))])
  }
  intensities. = intensities.[match(all_cells, intensities.$cell_id), ]
  
  intensities. = intensities.  %>%
    group_by(condition) %>%
    summarise_at(rownames(spe), mean)
  intensities. = intensities. %>% filter(!condition %in% c("THY", "TON", "LN") )
  
  mat = as.matrix(intensities.[,-1])
  rownames(mat) = intensities.$condition
  mat = mat[,setdiff(colnames(mat), markers_to_remove)]
  if(celltype == "Neutrophils")
    mat = mat[,setdiff(colnames(mat), c("Cytokeratin"))]
  
  cormat = cor(mat)
  cormat[is.na(cormat)] = 0
  hc_cor_marker = hclust(as.dist(1 - cormat), method = "ward.D2")
  cormat = cor(t(mat))
  cormat[is.na(cormat)] = 0
  hc_cor_sample = hclust(as.dist(1 - cormat), method = "ward.D2")
  
  color_df = data.frame("condition" = rownames(mat))
  h = Heatmap( 
    mat,
    cluster_columns = hc_cor_marker,
    cluster_rows = hc_cor_sample,
    name = "Fraction of cells",
    column_title = paste0(celltype, " - n= ", length(all_cells)),
    row_title = "Conditions",
    row_dend_side = "right",
    row_dend_width = unit(4, "cm"),
    column_dend_height = unit(5, "cm"),
    show_column_names = T,
    show_row_names = T,
    clustering_distance_columns ="pearson",
    clustering_distance_rows = "pearson",
    column_split = 10,
    border = TRUE,
    use_raster = FALSE, 
    col = c(viridis::magma(100)[100], viridis::magma(100)[50], viridis::magma(100)[1]),
    right_annotation = rowAnnotation(df = color_df, 
                                     col = list(condition = setNames(unique(spe$condition_color), 
                                                                     unique(spe$condition))))
    
  )
  draw(h)
}
dev.off()

#################################################################################
# T helper / T cyto ratio per condition
#################################################################################

meta = as.data.frame(colData(spe))
meta = meta  %>% filter(!condition %in% c("THY", "TON", "LN")) %>%
  filter(cell_type_CellSighter %in% c("T_helper", "T_cytotoxic")) %>%
  dplyr::group_by(sample_id, condition, disease, batch) %>% 
  dplyr::summarise(ratio = length(which(cell_type_CellSighter == "T_helper")) /
                     length(which(cell_type_CellSighter == "T_cytotoxic")))

png(file.path(output, "Plots", paste0("Ratio_T_helper_T_cytotoxic_cells_per_condition.png")), width = 1800, height = 1600, res = 300)
p = grouped_dotplot(meta, y = "ratio", categ1 = "condition", categ2 = "disease",
                          ref.group = "HD", add_violin = F,
                          color_categ1 = setNames(unique(spe$condition_color), unique(spe$condition))
)
p = p + xlab("") + ylab(paste0("Ratio T helper/T_cytotoxic")) +
  guides(fill="none", color = "none")

print(p)
dev.off()




#################################################################################
# PCA & UMAPS
#################################################################################
spe. = spe
mat = spe.@assays@data$CellSighter_marker_mat
dir.create(file.path(output, "Plots", "UMAP"))

# !! 
# mat_all = read.csv("output/cell_table/cell_table_arcsinh_transformed.csv.gz")
# mat = mat_all
# rownames(mat) = paste0( mat$fov, "-", mat$label)
# mat = mat[match(spe$cell_id, rownames(mat)),]
# mat = mat[,match(setdiff(rownames(spe), c("CD209","Ki67", "CD270", "CXCR4")) , colnames(mat))]
# mat = mat[match(spe.$cell_id, rownames(mat)),]
# mat = t(mat)
# !!

intensities = as.data.frame(t(mat))
intensities$sample_id = spe.$sample_id
intensities$cell_type_CellSighter = spe.$cell_type_CellSighter

intensities = intensities %>% group_by(sample_id, cell_type_CellSighter) %>%
  summarise_at(rownames(mat), mean) 
mat. = as.matrix(intensities[,-c(1:2)])
pca = prcomp(mat., center = F, scale. = F)
var = 100 * (pca$sdev) / sum(pca$sdev)

harmony  = harmony::HarmonyMatrix(data_mat = as.matrix(pca$x[,1:30]), meta_data = spe.$batch[match(intensities$sample_id, spe.$sample_id)], do_pca = F)

config <- umap::umap.defaults
config$metric <- "cosine"
umap <- umap::umap(harmony, config = config, method = "naive")
umap <- as.data.frame(umap$layout)
colnames(umap) <- c("Component_1", "Component_2")
umap$sample_id = intensities$sample_id
umap$cell_type_CellSighter = intensities$cell_type_CellSighter
umap$condition = spe.$condition[match(umap$sample_id, spe.$sample_id)]
umap$batch = spe.$batch[match(umap$sample_id, spe.$sample_id)]
umap$Th_group = spe.$Th_group[match(umap$sample_id, spe.$sample_id)]

umap_no_correction <- umap::umap(pca$x[,1:30], config = config, method = "naive")
umap_no_correction <- as.data.frame(umap_no_correction$layout)
colnames(umap_no_correction) <- c("Component_1", "Component_2")
umap_no_correction$sample_id = intensities$sample_id
umap_no_correction$cell_type_CellSighter = intensities$cell_type_CellSighter
umap_no_correction$condition = spe.$condition[match(umap_no_correction$sample_id, spe.$sample_id)]
umap_no_correction$batch = spe.$batch[match(umap_no_correction$sample_id, spe.$sample_id)]
umap_no_correction$Th_group = spe.$Th_group[match(umap_no_correction$sample_id, spe.$sample_id)]

umap = umap %>% filter(!condition %in% c("THY", "TON", "LN"))
umap_no_correction = umap_no_correction %>% filter(!condition %in% c("THY", "TON", "LN"))

png(file.path(output, "Plots", "UMAP",  paste0("UMAP_by_condition.png")), width = 1600, height = 1600, res = 300)
umap %>% ggplot(aes(x = Component_1, y = Component_2, color = condition)) +
  geom_point() +
  scale_color_manual(values = unique(spe$condition_color[match(levels(umap$condition), spe$condition)])) +
  theme_classic() + ggtitle("UMAP")
dev.off()

png(file.path(output, "Plots", "UMAP",  paste0("UMAP_by_batch.png")), width = 1600, height = 1600, res = 300)
umap %>% ggplot(aes(x = Component_1, y = Component_2, color = batch)) +
  geom_point() +
  theme_classic() + ggtitle("UMAP")
ylab(paste0("Component_2 (",round(var[2],2),"%)"))
dev.off()

png(file.path(output, "Plots", "UMAP",  paste0("UMAP_by_Thgroup.png")), width = 1600, height = 1600, res = 300)
umap %>% ggplot(aes(x = Component_1, y = Component_2, color = Th_group)) +
  geom_point() +
  theme_classic() + ggtitle("UMAP")
dev.off()

png(file.path(output, "Plots", "UMAP",  paste0("UMAP_by_celltype.png")), width = 2000, height = 1600, res = 300)
umap %>% ggplot(aes(x = Component_1, y = Component_2, color = cell_type_CellSighter)) +
  geom_point() +
  scale_color_manual(values = setNames(spe$cell_type_CellSighter_color, spe$cell_type_CellSighter)) +
  theme_classic() + ggtitle("UMAP") 
dev.off()

png(file.path(output, "Plots", "UMAP",  paste0("UMAP_by_sample_id.png")), width = 2500, height = 1600, res = 300)
umap %>% ggplot(aes(x = Component_1, y = Component_2, color = sample_id)) +
  geom_point() +
  scale_color_manual(values = setNames(spe$sample_id_color, spe$sample_id)) +
  theme_classic() + ggtitle("UMAP")
dev.off()

# UMAP_no_correction no correction ------------------------------------------------------------
png(file.path(output, "Plots", "UMAP",  paste0("UMAP_no_correction_by_condition.png")), width = 1600, height = 1600, res = 300)
umap_no_correction %>% ggplot(aes(x = Component_1, y = Component_2, color = condition)) +
  geom_point() +
  scale_color_manual(values = unique(spe$condition_color[match(levels(umap_no_correction$condition), spe$condition)])) +
  theme_classic() + ggtitle("UMAP_no_correction")
dev.off()

png(file.path(output, "Plots", "UMAP",  paste0("UMAP_no_correction_by_batch.png")), width = 1600, height = 1600, res = 300)
umap_no_correction %>% ggplot(aes(x = Component_1, y = Component_2, color = batch)) +
  geom_point() +
  theme_classic() + ggtitle("UMAP_no_correction")
dev.off()

png(file.path(output, "Plots", "UMAP",  paste0("UMAP_no_correction_by_Thgroup.png")), width = 1600, height = 1600, res = 300)
umap_no_correction %>% ggplot(aes(x = Component_1, y = Component_2, color = Th_group)) +
  geom_point() +
  theme_classic() + ggtitle("UMAP_no_correction")
dev.off()

png(file.path(output, "Plots", "UMAP",  paste0("UMAP_no_correction_by_celltype.png")), width = 2000, height = 1600, res = 300)
umap_no_correction %>% ggplot(aes(x = Component_1, y = Component_2, color = cell_type_CellSighter)) +
  geom_point() +
  scale_color_manual(values = setNames(spe$cell_type_CellSighter_color, spe$cell_type_CellSighter)) +
  theme_classic() + ggtitle("UMAP_no_correction") 
dev.off()

png(file.path(output, "Plots", "UMAP",  paste0("UMAP_no_correction_by_sample_id.png")), width = 2000, height = 1600, res = 300)
umap_no_correction %>% ggplot(aes(x = Component_1, y = Component_2, color = sample_id)) +
  geom_point() +
  scale_color_manual(values = setNames(spe$sample_id_color, spe$sample_id)) +
  theme_classic() + ggtitle("UMAP")
dev.off()

#################################################################################
# Morphology analyis
#################################################################################
library(ggbeeswarm)
library(rstatix)

# Comparing T helpers sizes
for( celltype in unique(spe$cell_type_CellSighter)){
  meta = as.data.frame(colData(spe))
  meta = meta[which(meta$area_nuclear != 0),]
  meta = meta %>% mutate(ratio_nuclear_area = area / area_nuclear)
  
  for(type in c("area", "equivalent_diameter", "num_concavities", "convex_hull_resid","ratio_nuclear_area", "log10_area_normalized_intensity")){
    meta. = meta %>% filter(!condition %in% c("TON", "THY", "LN")) %>%
      filter(cell_type_CellSighter == celltype) %>%
      dplyr::group_by(sample_id, condition, disease, batch) %>%
      dplyr::summarise(type = mean(.data[[type]]))
    
    png(file.path(output, "Plots", "Morphology", paste0("Cell_", type ,"_",celltype,".png")), width = 1800, height = 1600, res = 300)
    p = grouped_dotplot(meta., y = "type", categ1 = "condition", categ2 = "disease",
                        ref.group = "HD", add_violin = F,
                        color_categ1 = setNames(unique(spe$condition_color), unique(spe$condition))
    )
    p = p + xlab("") + ylab(paste0(type," of ", celltype)) +
      guides(fill="none", color = "none")
    print(p)
    dev.off()
  }
  
}


# Double positive CD4 / CD8
meta = as.data.frame(t(spe@assays@data$CellSighter_marker_mat))
meta$sample_id = spe$sample_id
meta$condition = spe$condition
meta$cell_type_CellSighter = spe$cell_type_CellSighter
meta$batch = spe$batch
meta$disease = spe$disease

meta. = meta  %>%
  dplyr::group_by(sample_id, condition, disease, batch) %>% 
  filter(cell_type_CellSighter %in% c("T_helper", "T_cytotoxic", "T_regulatory", "NKT")) %>%
  filter(!condition %in% c("TON", "THY", "LN")) %>%
  dplyr::summarise(ratio = 100 * length(which(CD4 == 1 & CD8a == 1)) /
                     length(which(CD4 > -1)))

png(file.path(output, "Plots", paste0("Double_positive_CD4_CD8.png")), width = 1800, height = 1600, res = 300)
 p = grouped_dotplot(meta., y = "ratio", categ1 = "condition", categ2 = "disease",
                         ref.group = "HD", add_violin = F,
                         color_categ1 = setNames(unique(spe$condition_color), unique(spe$condition))
)
p = p + xlab("") + ylab(paste0("Double positive CD4+ CD8a+ (% of all T cells)")) +
  guides(fill="none", color = "none")
print(p)
dev.off()


#################################################################################
# Correlation between number of macrophages and lymphatic cell
#################################################################################
meta = as.data.frame(colData(spe)) %>% filter(!condition %in% c("THY", "TON", "LN"))

pdf(file.path(output,"Plots", paste0("Correlation_celltype.pdf")), height = 5, width = 5)
for(type2 in unique(spe$cell_type_CellSighter)){
  for(type1 in unique(spe$cell_type_CellSighter)){
    
    if(type1 != type2){
      meta. = meta %>% group_by(sample_id, condition) %>%
        summarise(p_1 = length(which(cell_type_CellSighter == type1)) / n(), 
                  p_2 = length(which(cell_type_CellSighter == type2)) / n()
        )
      cor = round(cor(meta.$p_1, meta.$p_2), 3)
      pval = round(cor.test(meta.$p_1, meta.$p_2)$p.value, 5)
      
      if(abs(cor) > 0.4 & pval < 0.01 & !("Leukocyte" %in% c(type1, type2))){
        print(meta. %>% ggplot(aes(x = 100 * p_1, y = 100 * p_2)) + 
                geom_point(aes(color = condition))  + geom_smooth(lty = 3, se = F, method = "lm", show.legend = F, color = "grey50") +
                scale_color_manual(values = setNames(unique(spe$condition_color), unique(spe$condition))) + 
                theme_classic() +
                ggtitle(paste0("Cor = ", cor, "; p.val = ", pval)) + 
                xlab(paste0(type1, " (%)")) + ylab(paste0(type2, " (%)"))
        )
      }
    }
  }
}
dev.off()

#################################################################################
# Plotting TRBC1 levels in T cells / T helper in batch 1 # OR CD4+ / CD8+
#################################################################################
spe. = spe[,which(!spe$condition %in% c("LN", "THY","TON"))]
meta = as.data.frame(colData(spe.))
meta$TRBC1 = spe.@assays@data$CellSighter_marker_mat["TRBC1",]
meta$CD4 = spe.@assays@data$CellSighter_marker_mat["CD4",]
meta$CD8a = spe.@assays@data$CellSighter_marker_mat["CD8a",]
meta$CD3 = spe.@assays@data$CellSighter_marker_mat["CD3",]
meta = meta %>% filter(batch == "TMA1")

png(file.path(output, "Plots", paste0("TRBC1_TMA1_T_cells.png")), width = 1800, height = 1600, res = 300)
meta. = meta %>% filter(cell_type_CellSighter %in% T_cells) %>%
  group_by(sample_id, condition, disease, batch) %>%
  summarise(TRBC1 = mean(TRBC1)) 

p = grouped_dotplot(meta., y = "TRBC1", categ1 = "condition", categ2 = "disease",
                    ref.group = "HD", add_violin = F, shape_by = "batch",
                    color_categ1 = setNames(unique(spe$condition_color), unique(spe$condition))
)
p = p + xlab("") + ylab(paste0("TRBC1+ (% of all T cells)"))  +
  guides(fill="none", color = "none")
print(p)
dev.off()

png(file.path(output, "Plots", paste0("TRBC1_TMA1_T_helper_clonality.png")), width = 1800, height = 1600, res = 300)
meta. = meta %>% filter(condition != "HD") %>%
  filter(cell_type_CellSighter %in% T_cells) %>%
  group_by(sample_id, condition, disease, batch) %>%
  summarise(TRBC1 = mean(TRBC1)) %>% 
  mutate(Clonality = abs(0.5 - TRBC1))

p = grouped_dotplot(meta., y = "Clonality", categ1 = "condition", categ2 = "disease",
                    ref.group = "HD", add_violin = F, shape_by = "batch",
                    color_categ1 = setNames(unique(spe$condition_color), unique(spe$condition))
)
p = p + xlab("") + ylab(paste0("Clonality  in T cells")) +
  guides(fill="none", color = "none")
print(p)
dev.off()


png(file.path(output, "Plots", paste0("TRBC1_TMA1_T_helper.png")), width = 1800, height = 1600, res = 300)
meta. = meta %>% filter(cell_type_CellSighter %in% "T_helper") %>%
  group_by(sample_id, condition, batch, disease) %>%
  summarise(TRBC1 = mean(TRBC1)) 

p = grouped_dotplot(meta., y = "TRBC1", categ1 = "condition", categ2 = "disease",
                    ref.group = "HD", add_violin = F, shape_by = "batch",
                    color_categ1 = setNames(unique(spe$condition_color), unique(spe$condition))
)
p = p + xlab("") + ylab(paste0("TRBC1+ (% of T helper)"))  +
  guides(fill="none", color = "none")
print(p)
dev.off()


png(file.path(output, "Plots", paste0("TRBC1_TMA1_CD4+_cells.png")), width = 1800, height = 1600, res = 300)
meta. = meta %>% filter(CD4 > 0) %>% group_by(sample_id, condition, batch, disease) %>%
  summarise(TRBC1 = mean(TRBC1)) 
p = grouped_dotplot(meta., y = "TRBC1", categ1 = "condition", categ2 = "disease",
                    ref.group = "HD", add_violin = F, shape_by = "batch",
                    color_categ1 = setNames(unique(spe$condition_color), unique(spe$condition))
)
p = p + xlab("") + ylab(paste0("TRBC1+ (% of CD4+)"))  +
  guides(fill="none", color = "none")
print(p)
dev.off()

png(file.path(output, "Plots", paste0("TRBC1_TMA1_CD3+_CD4+_cells.png")), width = 1800, height = 1600, res = 300)
meta. = meta %>% filter(CD4 > 0 & CD3 > 0) %>% group_by(sample_id, condition, batch, disease) %>%
  summarise(TRBC1 = mean(TRBC1)) 

p = grouped_dotplot(meta., y = "TRBC1", categ1 = "condition", categ2 = "disease",
                    ref.group = "HD", add_violin = F, shape_by = "batch",
                    color_categ1 = setNames(unique(spe$condition_color), unique(spe$condition))
)
p = p + xlab("") + ylab(paste0("TRBC1+ (% of CD3+CD4+)"))  +
  guides(fill="none", color = "none")
print(p)
dev.off()