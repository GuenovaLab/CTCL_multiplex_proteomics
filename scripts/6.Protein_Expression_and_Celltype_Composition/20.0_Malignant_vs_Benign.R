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

#################################################################################
# Functional Markers in cell types - Tumor CD3 vs Non CD3 T helper
#################################################################################

dir.create(file.path(output, "T_helper_CD3"))

celltype = "T_helper"

for(type in c("Benign", "Malign")){
  pdf(file.path(output,  "T_helper_CD3", paste0(type,"_Th_CD3-_vs_Rest.pdf")), height = 4, width = 4)
  for(mark in setdiff(rownames(spe), c("CD3"))){
    
    
    intensities. = intensities %>% filter(cell_type_CellSighter %in% celltype, disease == type) %>%
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
