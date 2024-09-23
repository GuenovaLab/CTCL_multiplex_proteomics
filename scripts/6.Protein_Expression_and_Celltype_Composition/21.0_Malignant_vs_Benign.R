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
              "ggnewscale",
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
output = file.path("output", "Tumor_vs_Rest")
spe = qs::qread("output/SpatialExperiment.qs")

# Select MF and SS
spe_tumor = spe[, spe$condition %in% c("AD","PS","LP","DAR","MF", "SS")]

#################################################################################
# Cell Composition
#################################################################################
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
  
  meta = meta %>% filter(!disease %in% c("LN", "TON", "THY")) %>%
    dplyr::group_by(sample_id, disease, batch) %>% 
    dplyr::summarise(percent_CN = 100 * length(
      which(cell_type_CellSighter == i)) / length(cell_type_CellSighter))
  
  png(file.path(output, "Composition", paste0("Boxplot_composition_",i,".png")), width = 1200, height = 1150, res = 300)
  p = grouped_dotplot(meta, y = "percent_CN", categ1 = "disease", categ2 = NULL,
                      ref.group = "Benign", add_violin = F,
                      color_categ1 = setNames(unique(spe$disease_color), unique(spe$disease)))
  p = p + xlab("") + ylab(paste0("Percentage of ", gsub("_"," ",i), " (%)")) +
    guides(fill="none", color = "none")
  print(p)
  dev.off()
}

#################################################################################
# Cell Composition - T cells
#################################################################################
dir.create(file.path(output, "Composition_T_cells"))

for(i in c("T_helper", "T_regulatory", "T_cytotoxic", "NKT")){
  spe_tumor. = spe_tumor[, spe_tumor$cell_type_CellSighter %in% c("T_helper", "T_regulatory", "T_cytotoxic", "NKT")]
  
  spe_tumor.$cell_type_CellSighter[spe_tumor.$cell_type_CellSighter != i] = "Aa"
  meta = as.data.frame(colData(spe_tumor.))
  
  meta = meta %>% filter(!disease %in% c("LN", "TON", "THY")) %>%
    dplyr::group_by(sample_id, disease, batch) %>% 
    dplyr::summarise(percent_CN = 100 * length(
      which(cell_type_CellSighter == i)) / length(cell_type_CellSighter))
  
  png(file.path(output, "Composition_T_cells", paste0("Boxplot_composition_",i,".png")), width = 1200, height = 1150, res = 300)
  p = grouped_dotplot(meta, y = "percent_CN", categ1 = "disease", categ2 = NULL,
                      ref.group = "Benign", add_violin = F,
                      color_categ1 = setNames(unique(spe$disease_color), unique(spe$disease)))
  p = p + xlab("") + ylab(paste0("Percentage of ", gsub("_"," ",i), " (% of all T cells)")) +
    guides(fill="none", color = "none")
  print(p)
  dev.off()
}


#################################################################################
# Cell Composition - simplified
#################################################################################
dir.create(file.path(output, "Composition_simplified"))
for(i in unique(spe_tumor$cell_type_simplified)){
  spe_tumor. = spe_tumor
  name = "all"
  
  spe_tumor.$cell_type_simplified[spe_tumor.$cell_type_simplified != i] = "Aa"
  meta = as.data.frame(colData(spe_tumor.))
  
  meta = meta %>% filter(!disease %in% c("LN", "TON", "THY")) %>%
    dplyr::group_by(sample_id, disease, batch) %>% 
    dplyr::summarise(percent_CN = 100 * length(
      which(cell_type_simplified == i)) / length(cell_type_simplified))
  
  png(file.path(output, "Composition_simplified", paste0("Boxplot_composition_",i,".png")), width = 1200, height = 1150, res = 300)
  p = grouped_dotplot(meta, y = "percent_CN", categ1 = "disease", categ2 = NULL,
                      ref.group = "Benign", add_violin = F,
                      color_categ1 = setNames(unique(spe$disease_color), unique(spe$disease)))
  p = p + xlab("") + ylab(paste0("Percentage of ", gsub("_"," ",i), " (%)")) +
    guides(fill="none", color = "none")
  print(p)
  dev.off()
}

###############################################################################
# Cell Composition of IPIs
###############################################################################
spe_tumor_IPI = spe_tumor[,which(spe_tumor$immune_interaction_community == "infiltrate")]
dir.create(file.path(output, "Composition_IPI"))

for(i in unique(spe_tumor_IPI$cell_type_CellSighter)){
  spe_tumor_IPI. = spe_tumor_IPI
  name = "all"
  if(!(i %in% struct_celltype)){
    spe_tumor_IPI. = spe_tumor_IPI[, !spe_tumor_IPI$cell_type_CellSighter %in% struct_celltype]
    name = "immune"
  }
  spe_tumor_IPI.$cell_type_CellSighter[spe_tumor_IPI.$cell_type_CellSighter != i] = "Aa"
  meta = as.data.frame(colData(spe_tumor_IPI.))
  
  meta = meta %>% filter(!disease %in% c("LN", "TON", "THY")) %>%
    dplyr::group_by(sample_id, disease, batch) %>% 
    dplyr::summarise(percent_CN = 100 * length(
      which(cell_type_CellSighter == i)) / length(cell_type_CellSighter))
  
  png(file.path(output, "Composition_IPI", paste0("Boxplot_composition_",i,".png")), width = 1200, height = 1150, res = 300)
  p = grouped_dotplot(meta, y = "percent_CN", categ1 = "disease", categ2 = NULL,
                      ref.group = "Benign", add_violin = F,
                      color_categ1 = setNames(unique(spe$disease_color), unique(spe$disease)))
  p = p + xlab("") + ylab(paste0("Percentage of ", gsub("_"," ",i), " (%)")) +
    guides(fill="none", color = "none")
  print(p)
  dev.off()
}

#################################################################################
# Functional Markers in cell types
#################################################################################
intensities = as.data.frame(t(spe_tumor@assays@data$CellSighter_marker_mat))
intensities$sample_id = spe_tumor$sample_id
intensities$condition = spe_tumor$condition
intensities$disease = spe_tumor$disease
intensities$batch = spe_tumor$batch
intensities$cell_type_CellSighter = spe_tumor$cell_type_CellSighter
intensities = intensities %>% dplyr::group_by(sample_id, disease, condition, batch, cell_type_CellSighter) %>% 
  dplyr::summarise(across(Actin:Ki67, mean))

colors_disease = setNames(unique(spe_tumor$disease_color), unique(spe_tumor$disease))
colors_samp = setNames(unique(spe_tumor$sample_id_color), unique(spe_tumor$sample_id))

dir.create(file.path(output, "Markers_by_cond"))
for(celltype in unique(spe$cell_type_CellSighter)){
  
  pdf(file.path(output,  "Markers_by_cond", paste0("Functional_Marker_in_cell_type_CellSighter_",celltype,".pdf")), height = 4, width = 4)
  for(mark in colnames(intensities)[6:66]){
    intensities. = intensities %>% filter(cell_type_CellSighter %in% celltype)
    intensities.[[mark]] = 100 * intensities.[[mark]]
    
    p = grouped_dotplot(intensities., y =  mark, categ1 = "condition", categ2 = "disease",
                        ref.group = "MF", add_violin = T, stats = F,
                        color_categ1 = setNames(unique(spe$condition_color), unique(spe$condition)))
    p = p + xlab("") + ylab(paste0("Percentage of ", celltype, " expressing ", mark, " (%)")) +
      guides(fill="none", color = "none")
    print(p)
  }
  dev.off()
}

#################################################################################
# Functional Markers in cell types - normcounts
#################################################################################
intensities = as.data.frame(t(spe_tumor@assays@data$normcounts))
intensities$sample_id = spe_tumor$sample_id
intensities$condition = spe_tumor$condition
intensities$disease = spe_tumor$disease
intensities$batch = spe_tumor$batch
intensities$cell_type_CellSighter = spe_tumor$cell_type_CellSighter
intensities = intensities %>% dplyr::group_by(sample_id, disease, condition, batch, cell_type_CellSighter) %>% 
  dplyr::summarise(across(Actin:Ki67, mean))

colors_disease = setNames(unique(spe_tumor$disease_color), unique(spe_tumor$disease))
colors_samp = setNames(unique(spe_tumor$sample_id_color), unique(spe_tumor$sample_id))

dir.create(file.path(output, "Markers_normcounts"))
for(celltype in unique(spe$cell_type_CellSighter)){
  print(celltype)
  
  pdf(file.path(output,  "Markers_normcounts", paste0("Functional_Marker_in_cell_type_CellSighter_",celltype,".pdf")), height = 4, width = 4)
  for(mark in colnames(intensities)[6:65]){
    intensities. = intensities %>% filter(cell_type_CellSighter %in% celltype)
    
    p = grouped_dotplot(intensities., y =  mark, categ1 = "disease", categ2 = NULL,
                        ref.group = "Benign", add_violin = T,
                        color_categ1 = setNames(unique(spe$disease_color), unique(spe$disease)))
    p = p + xlab("") + ylab(paste0("Expression of ", mark, " in ", celltype)) +
      guides(fill="none", color = "none")
    print(p + geom_text(aes(label = sample_id), color = "black",position = position_jitterdodge(jitter.width = 0.5, jitter.height = 0.02)))
    
  }
  dev.off()
}


#################################################################################
# Functional Markers in cell types - quantNorm
#################################################################################
intensities = as.data.frame(t(limma::normalizeQuantiles(counts(spe))))
intensities$sample_id = spe$sample_id
intensities$condition = spe$condition
intensities$disease = spe$disease
intensities$batch = spe$batch
intensities$cell_type_CellSighter = spe$cell_type_CellSighter

intensities = intensities %>% dplyr::group_by(sample_id, disease, condition, batch, cell_type_CellSighter) %>% 
  dplyr::summarise(across(Actin:Ki67, mean))

colors_disease = setNames(unique(spe_tumor$disease_color), unique(spe_tumor$disease))
colors_samp = setNames(unique(spe_tumor$sample_id_color), unique(spe_tumor$sample_id))

dir.create(file.path(output, "Markers_quantNorm"))
for(celltype in c("T_helper", "T_cytotoxic", "NKT")){
  print(celltype)
  
  pdf(file.path(output,  "Markers_quantNorm", paste0("Functional_Marker_in_cell_type_CellSighter_",celltype,".pdf")), height = 4, width = 4)
  for(mark in colnames(intensities)[6:65]){
    intensities. = intensities %>% filter(cell_type_CellSighter %in% celltype, disease != "Healthy")
    
    p = grouped_dotplot(intensities., y =  mark, categ1 = "disease", categ2 = NULL,
                        ref.group = "Benign", add_violin = T,
                        color_categ1 = setNames(unique(spe$disease_color), unique(spe$disease)))
    p = p + xlab("") + ylab(paste0("quantNorm expression of ", mark, " in ", celltype, "")) +
      guides(fill="none", color = "none")
    print(p)
    
  }
  dev.off()
}

#################################################################################
# Functional Markers in cell types - IPI or isolated
#################################################################################

for(type in c("IPI", "isolated")){
  dir.create(file.path(output, paste0("Markers_", type)))
  spe_tumor. = spe_tumor[, !is.na(spe_tumor$immune_interaction_community)]
  
  if(type == "IPI"){
    spe_tumor. = spe_tumor.[,spe_tumor.$immune_interaction_community == "infiltrate"]
  } else{
    spe_tumor. = spe_tumor.[,spe_tumor.$immune_interaction_community == "isolated_cells"]
  }
  
  intensities_IPI = as.data.frame(t(spe_tumor.@assays@data$CellSighter_marker_mat))
  intensities_IPI$sample_id = spe_tumor.$sample_id
  intensities_IPI$condition = spe_tumor.$condition
  intensities_IPI$disease = spe_tumor.$disease
  intensities_IPI$batch = spe_tumor.$batch
  intensities_IPI$cell_type_CellSighter = spe_tumor.$cell_type_CellSighter
  intensities_IPI = intensities_IPI %>% dplyr::group_by(sample_id, disease, condition, batch, cell_type_CellSighter) %>% 
    dplyr::summarise(across(Actin:Ki67, mean))
  
  colors_disease = setNames(unique(spe_tumor.$disease_color), unique(spe_tumor.$disease))
  colors_samp = setNames(unique(spe_tumor.$sample_id_color), unique(spe_tumor.$sample_id))
  
  for(celltype in "all"){ #unique(spe_tumor.$cell_type_CellSighter)
    print(celltype)
    
    pdf(file.path(output,  paste0("Markers_", type), paste0("Functional_Marker_in_cell_type_CellSighter_",celltype,".pdf")), height = 4, width = 4)
    for(mark in colnames(intensities_IPI)[6:65]){
      if(celltype == "all") {
        intensities_IPI. = intensities_IPI %>% group_by(sample_id, disease, condition, batch)  %>% 
          dplyr::summarise(across(Actin:Ki67, mean))
      } else{
        intensities_IPI. = intensities_IPI %>% filter(cell_type_CellSighter %in% celltype)
      }
      intensities_IPI.[[mark]] = 100 * intensities_IPI.[[mark]]
      
      p = grouped_dotplot(intensities_IPI., y =  mark, categ1 = "disease", categ2 = NULL,
                          ref.group = "Benign", add_violin = T,
                          color_categ1 = setNames(unique(spe$disease_color), unique(spe$disease)))
      p = p + xlab("") + ylab(paste0("Percentage of ", celltype, " expressing ", mark, " (%)")) +
        guides(fill="none", color = "none")
      print(p)
    }
    dev.off()
  }
}



#################################################################################
# Functional Markers in cell types - IPI vs isolated
#################################################################################

dir.create(file.path(output, paste0("Markers_IPI_vs_isolated")))
spe_tumor. = spe_tumor[, !is.na(spe_tumor$immune_interaction_community)]
spe_tumor. = spe_tumor.[, spe_tumor.$disease %in% c("Benign", "Malign")]

intensities_IPI = as.data.frame(t(spe_tumor.@assays@data$CellSighter_marker_mat))
intensities_IPI$sample_id = spe_tumor.$sample_id
intensities_IPI$condition = spe_tumor.$condition
intensities_IPI$disease = spe_tumor.$disease
intensities_IPI$immune_interaction_community = factor(spe_tumor.$immune_interaction_community, levels =  c("isolated_cells", "infiltrate"))
intensities_IPI$batch = spe_tumor.$batch
intensities_IPI$cell_type_CellSighter = spe_tumor.$cell_type_CellSighter
intensities_IPI = intensities_IPI %>% dplyr::group_by(sample_id, immune_interaction_community, condition, disease, batch, cell_type_CellSighter) %>% 
  dplyr::summarise(across(Actin:Ki67, mean))

colors_immune_interaction_community = setNames( c("#d1d1d1ff", "#681FC2"),c("isolated_cells", "infiltrate"))

for(celltype in c("all", unique(spe_tumor.$cell_type_CellSighter))){
  print(celltype)
  
  pdf(file.path(output,  "Markers_IPI_vs_isolated", paste0("Functional_Marker_in_cell_type_CellSighter_",celltype,".pdf")), height = 4, width = 5)
  for(mark in colnames(intensities_IPI)[7:67]){
    if(celltype == "all") {
      intensities_IPI. = intensities_IPI %>% group_by(sample_id, immune_interaction_community, condition, batch, disease)  %>% 
        dplyr::summarise(across(Actin:Ki67, mean))
    } else{
      intensities_IPI. = intensities_IPI %>% filter(cell_type_CellSighter %in% celltype)
    }
    intensities_IPI.[[mark]] = 100 * intensities_IPI.[[mark]]
    intensities_IPI. = intensities_IPI.[which(intensities_IPI.$sample_id != "ROI-17"),]
    p = grouped_dotplot(intensities_IPI., y =  mark, categ1 = "immune_interaction_community", categ2 = "disease",
                        ref.group = "isolated_cells", add_violin = T,
                        color_categ1 = colors_immune_interaction_community)
    p = p + stat_compare_means(aes(label = after_stat(p.format)),
                               method = "t.test", paired = T)
    p = p + xlab("") + ylab(paste0("Percentage of ", celltype, " expressing ", mark, " (%)")) +
      guides(fill="none", color = "none")
    print(p)
  }
  dev.off()
}


#################################################################################
# Functional Markers in cell types - IPI vs isolated - simplified
#################################################################################

dir.create(file.path(output, paste0("Markers_IPI_vs_isolated")))
spe_tumor. = spe_tumor[, !is.na(spe_tumor$immune_interaction_community)]
spe_tumor. = spe_tumor.[, spe_tumor.$disease %in% c("Benign", "Malign")]

intensities_IPI_simplified = as.data.frame(t(spe_tumor.@assays@data$CellSighter_marker_mat))
intensities_IPI_simplified$sample_id = spe_tumor.$sample_id
intensities_IPI_simplified$condition = spe_tumor.$condition
intensities_IPI_simplified$disease = spe_tumor.$disease
intensities_IPI_simplified$immune_interaction_community = factor(spe_tumor.$immune_interaction_community, levels =  c("isolated_cells", "infiltrate"))
intensities_IPI_simplified$batch = spe_tumor.$batch
intensities_IPI_simplified$cell_type_simplified = spe_tumor.$cell_type_simplified
intensities_IPI_simplified = intensities_IPI_simplified %>% dplyr::group_by(sample_id, immune_interaction_community, condition, disease, batch, cell_type_simplified) %>% 
  dplyr::summarise(across(Actin:Ki67, mean))

colors_immune_interaction_community = setNames( c("#d1d1d1ff", "#681FC2"),c("isolated_cells", "infiltrate"))

for(celltype in c(unique(spe_tumor.$cell_type_simplified))){
  print(celltype)
  
  pdf(file.path(output,  "Markers_IPI_vs_isolated", paste0("Functional_Marker_in_cell_type_simplified_",celltype,".pdf")), height = 4, width = 5)
  for(mark in colnames(intensities_IPI_simplified)[7:67]){
    if(celltype == "all") {
      intensities_IPI_simplified. = intensities_IPI_simplified %>% group_by(sample_id, immune_interaction_community, condition, batch, disease)  %>% 
        dplyr::summarise(across(Actin:Ki67, mean))
    } else{
      intensities_IPI_simplified. = intensities_IPI_simplified %>% filter(cell_type_simplified %in% celltype)
    }
    intensities_IPI_simplified.[[mark]] = 100 * intensities_IPI_simplified.[[mark]]
    intensities_IPI_simplified. = intensities_IPI_simplified.[which(intensities_IPI_simplified.$sample_id != "ROI-17"),]
    p = grouped_dotplot(intensities_IPI_simplified., y =  mark, categ1 = "immune_interaction_community", categ2 = "disease",
                        ref.group = "isolated_cells", add_violin = T,
                        color_categ1 = colors_immune_interaction_community)
    p = p + stat_compare_means(aes(label = after_stat(p.format)),
                               method = "t.test", paired = T)
    p = p + xlab("") + ylab(paste0("Percentage of ", celltype, " expressing ", mark, " (%)")) +
      guides(fill="none", color = "none")
    print(p)
  }
  dev.off()
}

#################################################################################
# Functional Markers in cell types simplified 
#################################################################################

intensitites_simplified = as.data.frame(t(spe_tumor@assays@data$CellSighter_marker_mat))
intensitites_simplified$sample_id = spe_tumor$sample_id
intensitites_simplified$disease = spe_tumor$disease
intensitites_simplified$batch = spe_tumor$batch
intensitites_simplified$condition = spe_tumor$condition
intensitites_simplified$cell_type_simplified = spe_tumor$cell_type_simplified
intensitites_simplified = intensitites_simplified %>% group_by(sample_id, disease, condition, batch, cell_type_simplified) %>% 
  summarise(across(Actin:Ki67, mean))

dir.create(file.path(output, "Markers_Celltype_Simplififed"))
for(celltype in c("Myeloid", "Lymphocyte")){
  print(celltype)
  
  if( celltype != "all") p = intensitites_simplified %>% filter(cell_type_simplified %in% c(celltype)) else p = intensitites_simplified
  
  pdf(file.path(output,  "Markers_Celltype_Simplififed", paste0("Functional_Marker_in_cell_type_CellSighter_",celltype,".pdf")), height = 4, width = 4)
  for(mark in colnames(intensitites_simplified)[6:66]){
    intensitites_simplified. = intensitites_simplified %>% filter(cell_type_simplified == celltype)
    intensitites_simplified.[[mark]] = 100 * intensitites_simplified.[[mark]]
    
    p = grouped_dotplot(intensitites_simplified., y =  mark, categ1 = "disease", categ2 = NULL,
                        ref.group = "Benign", add_violin = T, shape_by = "batch",
                        color_categ1 = setNames(unique(spe$disease_color), unique(spe$disease)))
    p = p + xlab("") + ylab(paste0("Percentage of ", celltype, " expressing ", mark, " (%)")) +
      guides(fill="none", color = "none")
    print(p)
  }
  dev.off()
}

#################################################################################
# Correlations of markers
#################################################################################
intensities$condition = spe$condition[match(intensities$sample_id, spe$sample_id)]

dir.create(file.path(output, paste0("Marker_Cor")))

markers = c("Bcl2", "CD5", "CD3")
# markers = setdiff(rownames(spe), markers_to_remove)

pdf(file.path(output, paste0("Marker_Cor"), "Marker_Correlations_by_cond.pdf"), height = 4, width = 4)
for(celltype1 in c("T_helper")){
  print(celltype1)
  for(celltype2 in c("T_cytotoxic")){
    
    print(celltype2)
    
    for(mark1 in markers){
      for(mark2 in markers){
        
        if( !(mark1 == mark2)){
          
          x =  intensities %>% filter(cell_type_CellSighter == celltype1) %>% select(mark1)
          y = intensities %>% filter(cell_type_CellSighter == celltype2) %>% select(mark2)
          
          df = cbind(x, y[,5])
          
          
          COR = cor.test(df[[mark1]], df[[mark2]])
          if(abs(COR$estimate) > 0. & COR$p.value < 0.5){
            # R = summary(lm(df[[mark1]] ~ df[[celltype2]]))$r.squared
            p = df %>%
              ggplot(aes(x = 100 * .data[[mark1]], y = 100 * .data[[mark2]])) +
              geom_point(aes(fill = condition, shape = batch), size = 3) + 
              scale_fill_manual(values = setNames(colCondition$condition_color, colCondition$condition)) +
              scale_shape_manual(values = setNames(c(21, 24), c("TMA1", "TMA2"))) +
              guides(shape = guide_legend("TMA")) +
              geom_smooth(method = "lm", color = "black") + 
              xlab(paste0(mark1, " in ", gsub("_", " ",celltype1), " (%)")) +
              ylab(paste0(mark2, " in ", gsub("_", " ",celltype2), " (%)")) +
              theme_classic() +
              theme(
                axis.text.x = element_text(angle = 90), strip.background = element_blank(),
                panel.background = element_blank(),
                panel.ontop = F,
                panel.spacing = unit(0.75, "lines"), strip.text.x = element_text(size = 12.5))
            p = p +  annotate("text", x = 75, y = 75, size = 5, 
                              label = paste0("Pearson Cor Test, p = ", round(COR$p.value,4))) +
              annotate("text",  x = 80, y = 65, size = 5, 
                       label = paste0("Pearson Cor = ", round(COR$estimate,2)))
            
            print(p + NoLegend())
          }
        }
      }
    }
    
    
  }
}
dev.off()



#################################################################################
# Odd ratios
#################################################################################

pairs = list(
  "Keratinocyte" = c("Slug", "TOM22", "ZAP70", "betaActin"),
  "T_cytotoxic" = c("Bcl2", "CD195", "CD196", "CD5", "CD3", "betaActin", "p16INK4a"),
  "T_helper" = c("Bcl2", "CD3", "CD5", "p16INK4a", "p16INK4a", "betaActin"),
  "Monocytes" = c("p16INK4a"),
  "NKT" = c("Bcl2", "CD195", "betaActin"),
  "T_regulatory" = c("Bcl2", "CD195", "betaActin", "CD134", "TOM22", "p16INK4a"),
  "B_cell" = c("CD181", "CD195", "CD196", "CD204", "CD1c", "CD88", "TOM22", "HLAABC", "HLADR", "CD45RA"),
  "Macrophages" = c("betaActin", "CD195", "CD204", "p16INK4a", "CD88", "ZAP70"),
  "pDC" = c("CD195", "CD4", "CD2")
  
)

set.seed(47)
cells = c()
for(i in unique(spe_tumor$sample_id)){
  cells = c(cells, sample(spe_tumor$cell_id[spe_tumor$sample_id == i], min(1500, length(which(spe_tumor$sample_id == i)))))
}
spe_tumor_subsampled = spe_tumor[,match(cells, spe_tumor$cell_id)]
intensities_subsampled = as.data.frame(t(spe_tumor_subsampled@assays@data$CellSighter_marker_mat))
intensities_subsampled$sample_id = spe_tumor_subsampled$sample_id
intensities_subsampled$disease = spe_tumor_subsampled$disease
intensities_subsampled$batch = spe_tumor_subsampled$batch
intensities_subsampled$cell_type_CellSighter = spe_tumor_subsampled$cell_type_CellSighter
intensities_subsampled = intensities_subsampled %>%
  tidyr::pivot_longer(-c(sample_id,disease, cell_type_CellSighter, batch), names_to = c("marker"))

library(ggpattern)
library(ggtext)
library(glue)
list_OR_df = list()

for(celltype in names(pairs)){
  for(mark in setdiff(rownames(spe_tumor),c("Ki67", "CXCR4", "CD270", "CD209")) ){ # pairs[[celltype]]
    
    p =  intensities %>% filter(cell_type_CellSighter %in% c(celltype)) %>%
      dplyr::group_by(sample_id, disease, batch, marker) %>%
      filter(marker == mark) %>%
      dplyr::mutate(count = n())  %>% filter(count > 10) %>%
      dplyr::summarise(Percent = 100 * sum(value) / dplyr::n()) %>%
      ggplot(aes(y = Percent, fill = disease, x = disease)) + geom_boxplot(outlier.colour = "white") +
      scale_fill_manual(values = setNames(unique(spe_tumor$disease_color), unique(spe_tumor$disease))) +
      ggnewscale::new_scale("fill") +
      geom_jitter(aes(fill = batch),  colour="black",pch=21, size=2, width = 0.3) +
      scale_fill_manual(values = unique(spe_tumor$batch_color[match(meta$batch, spe_tumor$batch)])) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90)) +
      stat_compare_means(aes(label = after_stat(p.signif)),
                         method = "t.test", ref.group = "Inflammatory")  + ylab(paste0("% of ", celltype , " cells positive for marker")) + xlab("") +
      ggtitle(celltype)
    
    df = intensities_subsampled %>% filter(marker == mark)
    value = (sum(df$value[df$disease != "Inflammatory" & df$cell_type_CellSighter == celltype]) / length(which(df$disease != "Inflammatory" & df$cell_type_CellSighter == celltype))) /
      (sum(df$value[df$disease == "Inflammatory" & df$cell_type_CellSighter == celltype]) / length(which(df$disease == "Inflammatory" & df$cell_type_CellSighter == celltype)))
    alternative = ifelse( value > 1,  "greater","less")
    
    odd_ratios = lapply(unique(spe_tumor$cell_type_CellSighter), function(celltype.){
      contingency = data.frame(
        "Positive"  = c(sum(df$value[df$disease != "Inflammatory" & df$cell_type_CellSighter == celltype.]),
                        sum(df$value[df$disease == "Inflammatory" & df$cell_type_CellSighter == celltype.])
        ),
        "Negative"  = c(length(which(df$value[df$disease != "Inflammatory" & df$cell_type_CellSighter == celltype.] == 0)),
                        length(which(df$value[df$disease  == "Inflammatory" & df$cell_type_CellSighter == celltype.] == 0))
        ), row.names = c("Inflammatory", "SS")
      )
      df = data.frame("log2OddRatio" = log2((contingency[1,1] / contingency[1,2])  /  (contingency[2,1] / contingency[2,2])),
                      "Pval" = (fisher.test(contingency, alternative = alternative)$p.value),
                      "Difference" = 100 * abs((contingency[2,1] / (contingency[2,1] + contingency[2,2])) -
                                                 (contingency[1,1] / (contingency[1,1] + contingency[1,2])))
      )
      if(sum(contingency) < 100)
        df = data.frame("log2OddRatio" = 0, "Pval" = 1, "Difference" = 0)
      return(df)
    })
    
    odd_ratios = as.data.frame(do.call("rbind", odd_ratios) )
    odd_ratios$celltype = unique(spe_tumor$cell_type_CellSighter)
    odd_ratios$Significant = ifelse(odd_ratios$Pval < 0.01, "TRUE", "FALSE")
    odd_ratios$Styled = ifelse(odd_ratios$celltype == celltype, glue("<b>{celltype}</b>"), odd_ratios$celltype )
    odd_ratios$Marker = mark
    
    p2 = odd_ratios %>% ggplot() +
      geom_bar_pattern(aes(x = Styled, y = log2OddRatio, fill = Difference, pattern = Significant),
                       stat = "identity",
                       width = 0.5,
                       position = position_dodge(preserve = "single"),
                       color = "black",
                       pattern_fill = "white",
                       pattern_angle = 45,
                       pattern_density = 0.1,
                       pattern_spacing = 0.025,
                       pattern_key_scale_factor = 0.6) +
      scale_fill_gradientn(colours=c("grey", rev(viridis::magma(20)))) +
      scale_pattern_manual(values = c("FALSE" = "stripe", "TRUE" = "none")) +
      xlab("") + ylab("Log2OddRatio") + theme_classic() +
      theme(legend.title = element_text(size = 8),
            axis.text.x = ggtext::element_markdown(angle = 90, size=8, hjust = 1, vjust = 0.5),
            axis.title.x = element_blank(),
            legend.text = element_text(size = 6), legend.key.size = unit(0.2, 'cm'),
      ) +
      ggtitle("Specificity") +
      guides(pattern = guide_legend(override.aes = list(fill = "white")))
    
    png(file.path(output, paste0(celltype,"-",mark,"_specificity.png")), height = 1200, width = 1400, res = 300)
    print(p2)
    dev.off()
    
    # png(file.path(output, paste0(celltype,"-",mark,".png")), height = 2000, width = 1600, res = 300)
    # print(ggarrange(p, p2, heights = c(2, 1.3),
    #                 ncol = 1, nrow = 2))
    # dev.off()
    
    list_OR_df[[mark]] = odd_ratios
    
  }
  
}

write.csv(do.call("rbind",list_OR_df), file.path(output, "odd_ratios_celltype_marker.csv"))


#################################################################################
# Interacting cell types
#################################################################################

spe_tumor. = spe_tumor
spe_tumor.$cell_type_CellSighter[spe_tumor.$basal] = "Basal"

for(celltype1 in unique(spe_tumor.$cell_type_CellSighter)){ # c("Macrophages", "T_helper", "Basal", "T_cytotoxic", "NKT", "Lymphatic", "T_regulatory")
  for(celltype2 in unique(spe_tumor.$cell_type_CellSighter)){
    print(celltype2)
    interaction_list = list()
    
    for(samp in unique(spe_tumor.$sample_id)){
      spe_tumor.. = spe_tumor.[,spe_tumor.$sample_id == samp]
      
      cells1 = spe_tumor..$cell_id[spe_tumor..$cell_type_CellSighter == celltype1]
      cells2 = spe_tumor..$cell_id[spe_tumor..$cell_type_CellSighter == celltype2]
      
      if(length(cells1) > 1 & length(cells2)>1){
        mat = spe_tumor..@metadata$interaction[cells1,cells2]
        
        interaction_list[[samp]] = data.frame(interacting = length(which(apply(mat, 1, function(i) length(which(i>0))) > 0)),
                                              total1 = length(cells1),
                                              total2 = length(cells2),
                                              sample_id = samp,
                                              disease = unique(spe_tumor..$disease)
        )
      }
    }
    
    df = do.call("rbind", interaction_list)
    df = df %>% mutate(score_in_contact = 10000 * interacting / (total1 * total2))
    df$batch = spe_tumor$batch[match(df$sample_id, spe_tumor$sample_id)]
    
    png(file.path(output, "Score_in_contact", paste0("Score_in_contact_", celltype1,"_with_",celltype2,".png")), width = 1100, height = 1200, res = 300)
    p = df %>% 
      ggplot(aes(x = disease, y = score_in_contact, fill = disease)) +
      geom_boxplot(outlier.color = "white") +
      scale_fill_manual(values = setNames(unique(spe_tumor$disease_color), unique(spe_tumor$disease))) +
      ggnewscale::new_scale("fill") +
      geom_jitter(aes(fill = batch),  colour="black",pch=21, size=2, width = 0.3) +
      scale_fill_manual(values = setNames(unique(spe_tumor$batch_color), unique(spe_tumor$batch))) +
      theme_classic() + xlab("") + ylab(paste0("% of ", celltype1, " in direct contact with ", celltype2)) +
      theme(axis.text.x = element_text(angle = 90))  +
      stat_compare_means(aes(label = after_stat(p.signif)),
                         method = "t.test", ref.group = "Inflammatory")
    print(p)
    dev.off()
  }
}


###############################################################################
# Differential Analysis of IPI
###############################################################################
spe_tumor. = spe_tumor[setdiff(rownames(spe_tumor), markers_to_remove),
                       which(spe_tumor$immune_interaction_community == "infiltrate")]

list_DA = list()

for(celltype in c("all", unique(spe_tumor.$cell_type_CellSighter))){
  
  set.seed(47)
  spe_tumor_inflammatory = spe_tumor.[,spe_tumor.$disease == "Inflammatory"]
  inflammatory_cells = c()
  
  list_samp = list()
  
  # subsample 
  for(inflammatory in unique(spe_tumor_inflammatory$sample_id)){
    
    if(celltype == "all")
      cells = spe_tumor_inflammatory$cell_id[spe_tumor_inflammatory$sample_id == inflammatory]
    else
      cells = spe_tumor_inflammatory$cell_id[spe_tumor_inflammatory$sample_id == inflammatory & spe_tumor_inflammatory$cell_type_CellSighter == celltype]
    other_cells = spe_tumor$cell_id[spe_tumor$sample_id == inflammatory & spe_tumor$cell_type_CellSighter %in% struct_celltype]
    
    mat = spe_tumor_inflammatory@assays@data$CellSighter_marker_mat[,cells, drop = F]
    values = rowSums(mat) 
    percent = 100 * rowSums(mat) / ncol(mat)
    
    mat_other = spe_tumor@assays@data$CellSighter_marker_mat[rownames(spe_tumor.),other_cells, drop = F]
    values_other = rowSums(mat_other) 
    percent_other = 100 * rowSums(mat_other) / ncol(mat_other)
    
    if(length(cells) > 50){
      
      df = data.frame(
        celltype = celltype,
        sample_id = inflammatory,
        gene = names(values),
        group = "Inflammatory",
        n_positive = values,
        n_total = length(cells),
        percent = percent,
        percent_other = percent_other
      )
      cat(inflammatory, " TRBC1 = ", percent["TRBC1"],"%\n")
      list_samp[[inflammatory]] = df
      
    }
  }
  
  set.seed(47)
  spe_tumor_CTCL = spe_tumor.[, spe_tumor.$disease == "Lymphoma"  ]
  
  for(CTCL in unique(spe_tumor_CTCL$sample_id)){
    if(celltype == "all")
      cells = spe_tumor_CTCL$cell_id[spe_tumor_CTCL$sample_id == CTCL]
    else
      cells = spe_tumor_CTCL$cell_id[spe_tumor_CTCL$sample_id == CTCL & spe_tumor_CTCL$cell_type_CellSighter == celltype]
    other_cells = spe_tumor$cell_id[spe_tumor$sample_id == CTCL & spe_tumor$cell_type_CellSighter %in% struct_celltype]
    
    mat = spe_tumor_CTCL@assays@data$CellSighter_marker_mat[,cells, drop = F]
    values = rowSums(mat) 
    percent = 100 * rowSums(mat) / ncol(mat)
    
    mat_other = spe_tumor@assays@data$CellSighter_marker_mat[rownames(spe_tumor_CTCL),other_cells, drop = F]
    values_other = rowSums(mat_other) 
    percent_other = 100 * rowSums(mat_other) / ncol(mat_other)
    
    if(length(cells) > 50){
      df = data.frame(
        celltype = celltype,
        sample_id = CTCL,
        gene = names(values),
        group = "Lymphoma",
        n_positive = values,
        n_total = length(cells),
        percent = percent,
        percent_other = percent_other
      )
      cat(CTCL, " TRBC1 = ", percent["TRBC1"],"%\n")
      
      list_samp[[CTCL]] = df
    }
  }
  
  if(length(list_samp) > 0){
    df = do.call("rbind", list_samp)
    df_grouped = df %>% group_by(celltype, group, gene) %>% 
      summarise(avg_percent = mean(percent), avg_percent_other = mean(percent_other), n_positive = sum(n_positive), n_total = sum(n_total)) 
    
    df_grouped = df_grouped %>% ungroup %>% select(-celltype) %>%
      pivot_wider(names_from = group, values_from = c(avg_percent, avg_percent_other, n_positive, n_total))
    
    if(length(unique(df$sample_id[df$group == "Lymphoma"])) > 2 & length(unique(df$sample_id[df$group == "Inflammatory"])) > 2){
      
      p.values = sapply(1:nrow(df_grouped), function(i){
        contingency_table = data.frame(c(df_grouped$n_positive_Inflammatory[i], df_grouped$avg_percent_Lymphoma[i]),
                                       c(df_grouped$n_total_Inflammatory[i] - df_grouped$n_positive_Inflammatory[i], df_grouped$n_positive_Lymphoma[i] - df_grouped$n_positive_Lymphoma[i]))
        fisher.test(contingency_table)$p.value
      })
      df_grouped$p.values = p.adjust(p.values, method = "bonferroni")
      
      df_grouped = df_grouped %>% mutate(
        log2FC = log2(avg_percent_Lymphoma / avg_percent_Inflammatory),
        log2FC_other = log2(avg_percent_other_Lymphoma / avg_percent_other_Inflammatory),
        effect_other = avg_percent_other_Lymphoma- avg_percent_other_Inflammatory,
        effect = avg_percent_Lymphoma - avg_percent_Inflammatory,
        specificity =  abs(log2FC/log2FC_other) * abs(effect - effect_other)
      )
      
      df_grouped = df_grouped %>%
        mutate(categ = ifelse(p.values < 0.01,
                              ifelse(log2FC > 1, "signifplus",
                                     ifelse(log2FC < -1, "signifminus", "unsig")), "unsig"))
      
      changeRange = function (v, newmin = 0, newmax = 1) 
      {
        oldmin <- min(v, na.rm = TRUE)
        oldmax <- max(v, na.rm = TRUE)
        return(newmin + ((newmax - newmin) * (v - oldmin)/(oldmax - oldmin)))
      }
      
      df_grouped$specificity = changeRange(log2(df_grouped$specificity), 0, 1)
      
      
      list_DA[[celltype]] = df_grouped
      
    }    
  }
}


pdf(file.path(output, paste0("DA_per_celltype.pdf")), width = 6, height = 5)
for(celltype in names(list_DA)){
  df_grouped = list_DA[[celltype]]
  
  p = df_grouped  %>%
    ggplot(aes(x = log2FC, y = specificity, color = categ, label = gene, size = abs(effect))) +
    geom_point() +
    scale_color_manual(values = c("darkgreen","red","black")) +
    theme_classic() + xlab("Percent - log2FC") + ylab("Specificity") + 
    ggtitle(paste0(celltype, " - volcano plot - Lymphoma vs Inflammatory")) + 
    geom_hline(yintercept = 0.6, lwd = 0.35, lty = 2, col = "grey20") +
    geom_vline(xintercept = -1, lwd = 0.35, lty = 2, col = "grey20") +
    geom_vline(xintercept = 1,  lwd = 0.35, lty = 2, col = "grey20") +
    ggrepel::geom_text_repel() 
  print(p)
}
dev.off()


#################################################################################
# Interaction maps per celltype 
#################################################################################
dir.create(file.path(output, "Interaction_Maps"))

for(celltype1 in unique(spe_tumor$cell_type_CellSighter)){
  print(celltype1)
  
  dir.create(file.path(output, "Interaction_Maps", "absolute"))
  
  pdf(file.path(output, "Interaction_Maps",  "absolute", paste0("Interaction_Map_", celltype1, "_by_disease.pdf")))
  interaction_doughnut_plot(spe = spe_tumor, celltype = celltype1, 
                            group_by = "cell_type_CellSighter", enriched_protein = NULL,
                            plot_by = "disease", weighted = "absolute",
                            colors = setNames(unique(spe_tumor$cell_type_CellSighter_color), unique(spe_tumor$cell_type_CellSighter)))
  
  dev.off()
  
}


#################################################################################
# Ratio T-Helper / T-Cytotoxic
#################################################################################

meta = as.data.frame(colData(spe_tumor))
meta = meta  %>% filter(!disease %in% c("THY", "TON", "LN")) %>%
  filter(cell_type_CellSighter %in% c("T_helper", "T_cytotoxic")) %>%
  dplyr::group_by(sample_id, disease) %>% dplyr::summarise(ratio = length(which(cell_type_CellSighter == "T_helper")) /
                                                             length(which(cell_type_CellSighter == "T_cytotoxic")))

png(file.path(output, paste0("Ratio_T_helper_T_cytotoxic_cells_per_disease.png")), width = 1800, height = 1600, res = 300)
p = meta %>% ggplot(aes(x = disease, y = ratio, fill = disease)) +
  geom_boxplot(outlier.colour =  "white") +
  geom_jitter(aes(color = sample_id),  colour="black",pch=21, size=2, width = 0.3) +
  scale_fill_manual(values = setNames(unique(spe_tumor$disease_color), unique(spe_tumor$disease))) +
  theme_classic() + xlab("") + ylab(paste0("Ratio T helper/T_cytotoxic")) +
  NoLegend()  + ggtitle("Ratio T helper/T_cytotoxic") + theme(axis.text.x = element_text(angle = 90))  +
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "t.test", ref.group = "Inflammatory")  +
  geom_hline(yintercept = 1, lwd = 0.35, lty = 2, col = "grey20") 
print(p)
dev.off()


#################################################################################
# Molecular profile of T cells interacting with macrophages 
#################################################################################

mat_raw = read.csv("output/cell_table/cell_table_arcsinh_transformed.csv.gz")
for(celltype1 in c("Macrophages", "T_helper", "Basal", "T_cytotoxic", "NKT", "Lymphatic", "T_regulatory")){
  for(celltype2 in unique(spe_tumor.$cell_type_CellSighter)){
    print(celltype2)
    interacting = list()
    non_interacting = list()
    
    for(samp in unique(spe_tumor.$sample_id)){
      spe_tumor.. = spe_tumor.[,spe_tumor.$sample_id == samp]
      
      cells1 = spe_tumor..$cell_id[spe_tumor..$cell_type_CellSighter == celltype1]
      cells2 = spe_tumor..$cell_id[spe_tumor..$cell_type_CellSighter == celltype2]
      
      if(length(cells1) > 1 & length(cells2)>1){
        mat_interacting = spe_tumor..@assays@data$interaction[cells1,cells2]
        
        cells_interacting = rownames(mat_interacting)[which(apply(mat_interacting, 1, function(i) length(which(i>0))) > 0)]
        cells_non_interacting = setdiff(rownames(mat_interacting), cells_interacting)
        
        interacting[[samp]] = cells_interacting
        non_interacting[[samp]] = cells_non_interacting
        
      }
    }
    
    interacting = interacting[sapply(interacting, length) > 20]
    interacting = unlist(lapply(interacting, function(i) sample(i, 50, replace = T)))
    non_interacting = non_interacting[sapply(non_interacting, length) > 20]
    non_interacting = unlist(lapply(non_interacting, function(i) sample(i, 50, replace = T)))
    
    
    mat = spe_tumor@assays@data$CellSighter_marker_mat[setdiff(rownames(spe_tumor),c("CD270","Ki67", "CXCR4", "CD279", "CD209")),
                                                       c(interacting, non_interacting)]
    
    hist(rowSums(mat), breaks = 150)
    mat = mat[rowSums(mat) > 100,]
    
    png(file.path(output, paste0("Heatmap_", celltype1,"_interacting_with_",celltype2,".png")), width = 2200, height = 2000, res = 300)
    df = as.data.frame(colData(spe_tumor[,match(colnames(mat), spe_tumor$cell_id)]))[,c("sample_id", "disease")]
    df$interacting = FALSE
    df$interacting[match(interacting, rownames(df))] = TRUE
    
    h = Heatmap(
      t(mat),
      name = paste0(celltype1),
      column_title =  "Markers",
      row_title =paste0(celltype1),
      cluster_rows = T,
      cluster_columns = T,
      row_dend_side = "right",
      col =c("royalblue4", "gold"),
      show_column_names = TRUE, 
      show_row_names = FALSE, 
      border = TRUE,
      column_names_gp = gpar(fontsize = 6),
      use_raster = TRUE,
      right_annotation = rowAnnotation(df = df)
    )
    
    # right_annotation = rowAnnotation(df = df,
    #                                  col = list("Celltype" = setNames(unique(spe$cell_type_CellSighter_color[order(spe$cell_type_CellSighter)]),
    #                                                                   unique(spe$cell_type_CellSighter[order(spe$cell_type_CellSighter)])
    #                                  )))
    draw(h)
    dev.off()
  }
}


#################################################################################
# Functional Markers in cell types - Dual Mark
#################################################################################
intensities = as.data.frame(t(spe_tumor@assays@data$CellSighter_marker_mat))
intensities$sample_id = spe_tumor$sample_id
intensities$condition = spe_tumor$condition
intensities$disease = spe_tumor$disease
intensities$batch = spe_tumor$batch
intensities$cell_type_CellSighter = spe_tumor$cell_type_CellSighter


dir.create(file.path(output, "Markers_Dual"))

for(celltype in unique(spe$cell_type_CellSighter)){
  print(celltype)
  
  pdf(file.path(output,  "Markers_Dual", paste0("Functional_Marker_in_cell_type_CellSighter_",celltype,".pdf")), height = 4, width = 4)
  for(mark1 in c("CD5", "CD2", "CD3", "ZAP70", "Bcl2", "STING")){
    for(mark2 in c("CD5", "CD2", "CD3", "ZAP70", "Bcl2", "STING")){
      if(mark1 != mark2){
        intensities. = intensities %>% filter(cell_type_CellSighter %in% celltype) %>%
          select(any_of(mark1), any_of(mark2), sample_id, disease, condition, batch, cell_type_CellSighter) %>%
          mutate(double_pos = (.data[[mark1]] + .data[[mark2]] == 2) )
        intensities. = intensities. %>% dplyr::group_by(sample_id, disease, condition, batch, cell_type_CellSighter) %>% 
          dplyr::summarise(double_pos = 100 * mean(double_pos))
        
        p = grouped_dotplot(intensities., y =  "double_pos", categ1 = "disease", categ2 = NULL,
                            ref.group = "Benign", add_violin = T,
                            color_categ1 = setNames(unique(spe$disease_color), unique(spe$disease)))
        p = p + xlab("") + ylab(paste0("", celltype, " expressing both ", mark1," & ", mark2, " (%)")) +
          guides(fill="none", color = "none")
        print(p)
      }
      
    }
  }
  dev.off()
}


#################################################################################
# Functional Markers in cell types - Tumor CD3 vs Non CD3 T helper
#################################################################################

dir.create(file.path(output, "T_helper_CD3"))

celltype = "T_helper"

for(type in c("Benign", "Malign")){
  pdf(file.path(output,  "T_helper_CD3", paste0(type,"_Th_CD3-_vs_Rest.pdf")), height = 4, width = 4)
  for(mark in setdiff(rownames(spe), c("CD3"))){
    
    
    intensities. = intensities %>% filter(cell_type_CellSighter %in% celltype, disease == type) %>%
      # filter((CD3 == 0 & CD5 == 0) | (CD3 == 1 & CD5 == 1)) %>%
      dplyr::mutate(CD3 = ifelse(CD3 == 0, "CD3m", "CD3p")) %>% 
      dplyr::group_by(CD3, sample_id, disease, condition, batch, cell_type_CellSighter) %>%
      dplyr::summarise(marker = 100 * mean(.data[[mark]])) 
    if(mark == "TRBC1") intensities. = intensities. %>% filter(batch == "TMA1")
    if(mark == "Ki67") intensities. = intensities. %>% filter(batch == "TMA2")
    if(mark == "CXCR4") intensities. = intensities. %>% filter(batch == "TMA1")
    
    intensities.$CD3 = factor(intensities.$CD3, levels = c("CD3p", "CD3m"))
    p = grouped_dotplot(intensities., y =  "marker", categ1 = "CD3", categ2 = NULL,
                        ref.group = "CD3p", add_violin = T, stats = TRUE,  paired = TRUE,
                        color_categ1 =  c("CD3p" = "#74C276", "CD3m" = "#B54A48")) 
    
    p = p + xlab("") + ylab(paste0("", celltype, " expressing ", mark, "(%")) +
      guides(fill="none", color = "none") + 
      theme(axis.text.x = element_text(size=12))+
      scale_x_discrete(labels=c("CD3p" = "CD3+", "CD3m" = "CD3-"), )
    print(p)
  }
  dev.off()
}

#################################################################################
# Functional Markers in cell types - Tumor CD3 vs Non CD3 T helper
#################################################################################
mark = "ZAP70"

intensities. = intensities %>% filter(cell_type_CellSighter %in% celltype, disease == type) %>%
  # filter((CD3 == 0 & CD5 == 0) | (CD3 == 1 & CD5 == 1)) %>%
  dplyr::mutate(CD3 = ifelse(CD3 == 0, "CD3m", "CD3p")) %>% 
  dplyr::group_by(CD3, sample_id, disease, condition, batch, cell_type_CellSighter) %>%
  dplyr::summarise(marker = 100 * mean(.data[[mark]])) 
if(mark == "TRBC1") intensities. = intensities. %>% filter(batch == "TMA1")
if(mark == "Ki67") intensities. = intensities. %>% filter(batch == "TMA2")
if(mark == "CXCR4") intensities. = intensities. %>% filter(batch == "TMA1")

intensities.$CD3 = factor(intensities.$CD3, levels = c("CD3p", "CD3m"))
p = grouped_dotplot(intensities., y =  "marker", categ1 = "CD3", categ2 = NULL,
                    ref.group = "CD3p", add_violin = T, stats = TRUE,  paired = TRUE,
                    color_categ1 =  c("CD3p" = "#74C276", "CD3m" = "#B54A48")) 

p = p + xlab("") + ylab(paste0("", celltype, " expressing ", mark, "(%")) +
  guides(fill="none", color = "none") + 
  theme(axis.text.x = element_text(size=12))+
  scale_x_discrete(labels=c("CD3p" = "CD3+", "CD3m" = "CD3-"), )
print(p)


plot_reduced_dim_scExp(spe_tumor., color_by = "condition")

