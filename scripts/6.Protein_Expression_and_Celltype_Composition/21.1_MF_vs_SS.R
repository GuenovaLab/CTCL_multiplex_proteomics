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



################################################################################
# Number of directly interacting celltype 1 vs other celltypes
################################################################################
spe_tumor. = spe_tumor
spe_tumor.$cell_type_CellSighter[spe_tumor.$basal] = "Basal"

for(celltype1 in c("Macrophages", "T_helper", "Basal", "T_cytotoxic", "NKT", "Lymphatic", "T_regulatory")){
  for(celltype2 in unique(spe_tumor.$cell_type_CellSighter)){
    print(celltype2)
    interaction_list = list()
    
    for(samp in unique(spe_tumor.$sample_id)){
      spe_tumor.. = spe_tumor.[,spe_tumor.$sample_id == samp]
      
      cells1 = spe_tumor..$cell_id[spe_tumor..$cell_type_CellSighter == celltype1]
      cells2 = spe_tumor..$cell_id[spe_tumor..$cell_type_CellSighter == celltype2]
      
      if(length(cells1) > 1 & length(cells2)>1){
        mat = spe_tumor..@assays@data$interaction[cells1,cells2]
        
        interaction_list[[samp]] = data.frame(interacting = length(which(apply(mat, 1, function(i) length(which(i>0))) > 0)),
                                              total1 = length(cells1),
                                              total2 = length(cells2),
                                              sample_id = samp,
                                              condition = unique(spe_tumor..$condition)
        )
      }
    }
    
    df = do.call("rbind", interaction_list)
    df = df %>% mutate(percent_in_contact = 100 * interacting / total1)
    df$batch = spe_tumor$batch[match(df$sample_id, spe_tumor$sample_id)]
    png(file.path(output, "Percent_in_contact", paste0("Percent_in_contact_", celltype1,"_with_",celltype2,".png")), width = 1100, height = 1200, res = 300)
    p = df %>% 
      ggplot(aes(x = condition, y = percent_in_contact, fill = condition)) +
      geom_boxplot(outlier.color = "white") +
      scale_fill_manual(values = setNames(unique(spe_tumor$condition_color), unique(spe_tumor$condition))) +
      ggnewscale::new_scale("fill") +
      geom_jitter(aes(fill = batch),  colour="black",pch=21, size=2, width = 0.3) +
      scale_fill_manual(values = setNames(unique(spe_tumor$batch_color), unique(spe_tumor$batch))) +
      theme_classic() + xlab("") + ylab(paste0("% of ", celltype1, " in direct contact with ", celltype2)) +
      theme(axis.text.x = element_text(angle = 90))  +
      stat_compare_means(aes(label = after_stat(p.signif)),
                         method = "t.test", ref.group = "MF")
    print(p)
    dev.off()
  }
}

################################################################################
# Map of interactions
################################################################################
general_types = list(
  "Structure" = c(struct_celltype),
  "Myeloid" = APC_types,
  "T_helper1" = "T_helper",
  "Toxic" = c("T_cytotoxic", "NKT"),
  "T_regulatory1" = "T_regulatory"
)

spe_tumor. = spe_tumor[, spe_tumor$cell_type_CellSighter %in% unlist(general_types)]
spe_tumor.$General_Types = gsub(".$","",names(unlist(general_types)))[match(spe_tumor.$cell_type_CellSighter, unlist(general_types))]

# set.seed(47)
# cells = c()
# for(i in unique(spe_tumor.$sample_id)){
#   cells = c(cells, sample(spe_tumor.$cell_id[spe_tumor.$sample_id == i], 1500))
# }
# spe_tumor. = spe_tumor.[,match(cells, spe_tumor.$cell_id)]
# table(spe_tumor.$sample_id)

for(celltype1 in unique(spe_tumor.$cell_type_CellSighter)){
  
  interaction_list = list()
  
  for(samp in unique(spe_tumor.$sample_id)){
    spe_tumor.. = spe_tumor.[,spe_tumor.$sample_id == samp]
    
    cells1 = spe_tumor..$cell_id[which(spe_tumor..$cell_type_CellSighter == celltype1  & spe_tumor..$immune_interaction_community != "infiltrate")]
    
    if(length(cells1) > 1){
      mat = spe_tumor..@assays@data$interaction[cells1,]
      mat[mat>0]=1
      types = as.factor(spe_tumor..$General_Types)
      types_tab = data.frame(
        sample_id = samp,
        condition = unique(spe_tumor..$condition),
        batch = unique(spe_tumor..$batch),
        type = levels(types),
        sum = rowSums(apply(mat, 1, function(i)
          as.data.frame(table(types[i > 0]))[,2]
        )),
        total = as.numeric(table(types))
        )
      interaction_list[[samp]] = types_tab
    }
  }
  
  df = do.call("rbind", interaction_list)
  if(!is.null(df)){

  df = df %>% mutate(weighted_sum = sum / total)
  df = df %>% group_by(sample_id) %>%
    mutate(percent = 100 * weighted_sum / sum(weighted_sum))
  df = df %>% group_by(condition, type) %>%
    summarise(percent = mean(percent))
  
  # for(celltype2 in unique(df$type)){
  #   contingency = data.frame("Type" = c(df$sum[df$type == celltype2 & df$condition == "MF"], df$sum[df$type == celltype2 & df$condition == "SS"]), 
  #                       "NonType" = c( sum(df$sum[df$type != celltype2 & df$condition == "MF"]), sum(df$sum[df$type != celltype2 & df$condition == "SS"])),
  #                       row.names = c("MF", "SS"))
  #   colnames(contingency) = paste0(c("", "Non-"),celltype2)
  #   print(paste0(celltype1," - Number of interactions ",celltype2, " p-value = ", fisher.test(contingency)$p.value))
  #   print(contingency)
  #   write.table(paste0(celltype1," - Number of interactions ",celltype2, " p-value = ", fisher.test(contingency)$p.value), file.path(output, "Percent_in_contact", "contingency_tab_pvalues.csv"), append = TRUE)
  #   write.table(contingency, file.path(output, "Percent_in_contact", "contingency_tab_pvalues.csv"), append = TRUE)
  # }
  # 
  
  for(cond in c("MF", "SS")){
    data = df %>% filter(condition == cond)
    
    # Compute percentages
    data$fraction <- data$percent / sum(data$percent)
    
    # Compute the cumulative percentages (top of each rectangle)
    data$ymax <- cumsum(data$fraction)
    
    # Compute the bottom of each rectangle
    data$ymin <- c(0, head(data$ymax, n=-1))
    
    # Compute label position
    data$labelPosition <- (data$ymax + data$ymin) / 2
    
    # Compute a good label
    data$label <- paste0(data$type, "\n (", round(100*data$fraction,1) ,"%)")
    
    png(file.path(output, "Percent_in_contact", paste0("Interaction_map_", celltype1,"_",cond,"_non_infiltrate.png")), width = 1300, height = 1200, res = 150)
    
    colors = setNames(
      c( "#D72638","#f04410ff", "#0d522aff", "#6B0504","#780C50"),
      c("Toxic", "T_regulatory", "T_helper", "Structure", "Myeloid")
    )
    
    # Make the plot
    print(ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=type)) +
            geom_rect() +
            geom_text( x=2, aes(y=labelPosition, label=label, color=type), size=6) + # x here controls label position (inner / outer)
            scale_fill_manual(values = colors) +
            scale_color_manual(values = colors) +
            annotate("text", x = -1, y = 0, label = celltype1, size = 6) +
            coord_polar(theta="y") +
            xlim(c(-1, 4)) +
            theme_void() +
            theme(legend.position = "none"))
    dev.off()
  }
  
  }
}

################################################################################
# Ratio T-Helper / T-Cytotoxic
################################################################################

library(dplyr)
meta = as.data.frame(colData(spe_tumor))
meta = meta  %>% filter(!condition %in% c("THY", "TON", "LN")) %>%
  filter(cell_type_CellSighter %in% c("T_helper", "T_cytotoxic")) %>%
  dplyr::group_by(sample_id, condition) %>% dplyr::summarise(ratio = length(which(cell_type_CellSighter == "T_helper")) /
                                                               length(which(cell_type_CellSighter == "T_cytotoxic")))

png(file.path(output, paste0("Ratio_T_helper_T_cytotoxic_cells_per_condition.png")), width = 1800, height = 1600, res = 300)
p = meta %>% ggplot(aes(x = condition, y = ratio, fill = condition)) +
  geom_boxplot(outlier.colour =  "white") +
  geom_jitter(aes(color = sample_id),  colour="black",pch=21, size=2, width = 0.3) +
  scale_fill_manual(values = setNames(unique(spe_tumor$condition_color), unique(spe_tumor$condition))) +
  theme_classic() + xlab("") + ylab(paste0("Ratio T helper/T_cytotoxic")) +
  NoLegend()  + ggtitle("Ratio T helper/T_cytotoxic") + theme(axis.text.x = element_text(angle = 90))  +
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "t.test", ref.group = "MF")  +
  geom_hline(yintercept = 1, lwd = 0.35, lty = 2, col = "grey20") 
print(p)
dev.off()


################################################################################
# Molecular profile of T cells interacting with macrophages
################################################################################

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
    df = as.data.frame(colData(spe_tumor[,match(colnames(mat), spe_tumor$cell_id)]))[,c("sample_id", "condition")]
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



