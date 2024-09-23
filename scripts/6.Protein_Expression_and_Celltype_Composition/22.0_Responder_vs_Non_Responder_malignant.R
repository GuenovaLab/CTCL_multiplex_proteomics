# Comparing the diseases at different stage 
# authors: Pacome Prompsy
# contact: pacome.prompsy@chuv.ch
# Guenova Lab
# CHUV (Centre Hospitalier Universitaire Vaudois), Lausanne, Suisse

cat("Comparing the diseases at different stage  (Early, Late ...) \n")

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

output_dir = file.path(args$output, "CTCL_response")
if (!dir.exists(output_dir))
  dir.create(output_dir)

spe = qs::qread("output/SpatialExperiment.qs")
colors_early_late = setNames(c("#EDCC5F", "#941D1B"), c("Early", "Late"))
colors_progression_regression = setNames(c("#7699D4", "#B3A632FA", "#B54A48", "#74C276"),
                                         c("Healthy", "Benign", "Stable/Progression", "Regression"))

######################## ######################## ######################## ########################
# Add the stage information to the SingleCellExperiment ########################
sample_metadata = read.csv(file.path("annotation", "Sample_metadata.csv"), sep = ";", na.strings = "")
spe$stage = sample_metadata$Relative.Severity[match(spe$sample_id, sample_metadata$ROI)]
spe$stage =  gsub("\\/Mild|\\/Sever", "", spe$stage)

spe$CTCL.CAILS = sample_metadata$CTCL.CAILS[match(spe$sample_id, sample_metadata$ROI)]
spe$CTCL.CAILS[which(!spe$CTCL.CAILS %in% c("Stable/Progression", "Regression"))] = NA
spe$CTCL.CAILS[which(spe$condition == "HD")] = "Healthy"
spe$CTCL.CAILS[which(spe$disease == "Benign")] = "Benign"

spe$Inflammaroty.severity = sample_metadata$Inflammarty.severity[match(spe$sample_id, sample_metadata$ROI)]
spe$Inflammaroty.severity = gsub("\\/Mild|Mild\\/|unsecure|Sever\\/| ", "", spe$Inflammaroty.severity)

spe$CTCL.CAILS = factor(spe$CTCL.CAILS, levels = c("Healthy", "Benign", "Stable/Progression", "Regression"))

qs::qsave(spe, "output/SpatialExperiment.qs")
######################## ######################## ######################## ########################

###############################################################################
# Composition changes Late vs Early Stage
###############################################################################

spe. = spe[setdiff(rownames(spe), markers_to_remove), !spe$condition %in% c("THY", "TON", "LN") & !is.na(spe$stage)]
spe. = spe.[, spe.$disease %in% c("Malign")]
spe.= spe.[,!spe.$cell_type_CellSighter %in% c(
  "Keratinocyte",
  "Lymphatic",
  "Endothelial",
  "Unknown-Stroma"
)]

dir.create(file.path(output_dir, "Stage"))
for(i in unique(spe$cell_type_CellSighter)){
  meta = as.data.frame(colData(spe.)) %>% dplyr::group_by(sample_id, condition, stage, disease, batch) %>% 
    dplyr::summarise(percent_CN = 100 * length(
      which(cell_type_CellSighter == i)) / length(cell_type_CellSighter))
  
  png(file.path(output_dir, "Stage", paste0("Boxplot_composition_",i,"_no_temoin_w_stats.png")), width = 1800, height = 1600, res = 300)
  p = grouped_dotplot(meta, y = "percent_CN", categ1 = "stage", categ2 = NULL,
                      ref.group = "Early", add_violin = F,
                      color_categ1 = colors_early_late)
  p = p + xlab("") + ylab(paste0("Percentage of ", gsub("_"," ",i), " (%)")) +
    guides(fill="none", color = "none")
  print(p)
  dev.off()
}


###############################################################################
# Composition changes - CTCL. CAILS evolution
###############################################################################
spe. = spe[setdiff(rownames(spe), markers_to_remove), !spe$condition %in% c("THY", "TON", "LN") & !is.na(spe$CTCL.CAILS)]
spe.= spe.[,!spe.$cell_type_CellSighter %in% c(
  "Keratinocyte",
  "Lymphatic",
  "Endothelial",
  "Unknown-Stroma"
)]

dir.create(file.path(output_dir, "CAILS" ,"Composition"))
for(i in unique(spe$cell_type_CellSighter)){
  meta = as.data.frame(colData(spe.)) %>% dplyr::group_by(sample_id, condition, CTCL.CAILS, disease, batch) %>% 
    dplyr::summarise(percent_CN = 100 * length(
      which(cell_type_CellSighter == i)) / length(cell_type_CellSighter))
  
  png(file.path(output_dir, "CAILS", "Composition", paste0("Boxplot_composition_",i,"_no_temoin_w_stats.png")), width = 1800, height = 1600, res = 300)
  p = grouped_dotplot(meta, y = "percent_CN", categ1 = "CTCL.CAILS", categ2 = NULL,
                      ref.group = "Stable/Progression", add_violin = F,
                      color_categ1 = colors_progression_regression)
  p = p + xlab("") + ylab(paste0("Percentage of ", gsub("_"," ",i), " (% of immune cells)")) +
    guides(fill="none", color = "none")
  print(p)
  dev.off()
}


###############################################################################
# Composition changes - in the IPI
###############################################################################
spe. = spe[setdiff(rownames(spe), markers_to_remove), !spe$condition %in% c("THY", "TON", "LN") & !is.na(spe$CTCL.CAILS)]
spe.= spe.[,!spe.$cell_type_CellSighter %in% c(
  "Keratinocyte",
  "Lymphatic",
  "Endothelial",
  "Unknown-Stroma"
)]
spe. = spe.[,which(!is.na(spe.$immune_interaction_community))]

for(IPI in c("isolated_cells", "infiltrate")){
  dir.create(file.path(output_dir, "CAILS" ,  paste0("Composition_", IPI)))
  spe.. = spe.[,spe.$immune_interaction_community == IPI]
  
  for(i in unique(spe$cell_type_CellSighter)){
    meta = as.data.frame(colData(spe..)) %>% dplyr::group_by(sample_id, condition, CTCL.CAILS, disease, batch) %>% 
      dplyr::summarise(percent_CN = 100 * length(
        which(cell_type_CellSighter == i)) / length(cell_type_CellSighter))
    
    png(file.path(output_dir, "CAILS" , paste0("Composition_", IPI), paste0("Boxplot_composition_",i,"_no_temoin_w_stats.png")), width = 1800, height = 1600, res = 300)
    p = grouped_dotplot(meta, y = "percent_CN", categ1 = "CTCL.CAILS", categ2 = NULL,
                        ref.group = "Stable/Progression", add_violin = F,
                        color_categ1 = colors_progression_regression)
    p = p + xlab("") + ylab(paste0("", gsub("_"," ",i), " (% of immune cells)")) +
      guides(fill="none", color = "none")
    print(p)
    dev.off()
  }
  
}


###############################################################################
# Composition changes - CTCL. CAILS evolution - simplified
###############################################################################
spe. = spe[setdiff(rownames(spe), markers_to_remove), !spe$condition %in% c("THY", "TON", "LN") & !is.na(spe$CTCL.CAILS)]
spe.= spe.[,!spe.$cell_type_CellSighter %in% c(
  "Keratinocyte",
  "Lymphatic",
  "Endothelial",
  "Unknown-Stroma"
)]

dir.create(file.path(output_dir, "CAILS", "Composition"))
for(i in unique(spe$cell_type_simplified)){
  meta = as.data.frame(colData(spe.)) %>% dplyr::group_by(sample_id, condition, CTCL.CAILS, disease, batch) %>% 
    dplyr::summarise(percent_CN = 100 * length(
      which(cell_type_simplified == i)) / length(cell_type_simplified))
  
  png(file.path(output_dir, "CAILS",  "Composition", paste0("Boxplot_composition_",i,"_no_temoin_w_stats.png")), width = 1800, height = 1600, res = 300)
  p = grouped_dotplot(meta, y = "percent_CN", categ1 = "CTCL.CAILS", categ2 = NULL,
                      ref.group = "Stable/Progression", add_violin = F,
                      color_categ1 = colors_progression_regression)
  p = p + xlab("") + ylab(paste0("Percentage of ", gsub("_"," ",i), " (% of immune)")) +
    guides(fill="none", color = "none")
  print(p)
  dev.off()
}

#################################################################################
# Density  - CAILS evolution
#################################################################################
spe. = spe[setdiff(rownames(spe), markers_to_remove), !spe$condition %in% c("THY", "TON", "LN") & !is.na(spe$CTCL.CAILS)]
spe.= spe.[,!spe.$cell_type_CellSighter %in% c(
  "Keratinocyte",
  "Lymphatic",
  "Endothelial",
  "Unknown-Stroma"
)]

# Plot and calculate the density of cell types per region ----------------------
meta = read.csv(file.path(args$output, "pixie", "region_pixel_output_dir", "areas_per_sample.csv"))
cold = as.data.frame(colData(spe.))
df = cold %>% group_by(sample_id, batch, disease, cell_type_CellSighter, CTCL.CAILS) %>%
  dplyr::summarise(n_cells = n())
df = df %>% pivot_wider(names_from = cell_type_CellSighter, values_from = n_cells, values_fill = 0) 

meta = left_join(meta, df, by = "sample_id")
meta = meta[which(!is.na(meta$CTCL.CAILS)),]

dir.create(file.path(output_dir, "CAILS", "Density"))

df = meta %>% filter(!condition %in% c("TON", "THY", "LN"))  %>% select(sample_id, condition, batch, disease)

library(ggbeeswarm)
library(ggpubr)
for(type in unique(spe.$cell_type_CellSighter)){
  meta. = meta %>% dplyr::mutate(density = .data[[type]] / tissue_area) %>%
    dplyr::select(CTCL.CAILS, disease, batch, density, sample_id)
  df[[type]] = 0
  df[[type]] = meta.$density
  
  pdf(file.path(output_dir, "CAILS", "Density", paste0("Density_", type,"_within_tissue.pdf")), width = 4.5, height = 4)
  p = grouped_dotplot(meta., y = "density", categ1 = "CTCL.CAILS", categ2 = NULL,
                      color_categ1 = colors_progression_regression,
                      ref.group = "Stable/Progression", add_violin = F)
  p = p + xlab("") + ylab(paste0("Density of ", gsub("_"," ",type), " in tissue (cells / mm²)")) +
    guides(fill="none", color = "none")
  print(p)
  dev.off()
  
}


#################################################################################
# Interaction Map CAILS Regression vs Non Regression
################################################################################# 
dir.create(file.path(output_dir, "CAILS", "Interaction_Map"))

for(evolution in levels(spe$CTCL.CAILS)){
  
  spe. = spe[,!spe$condition %in% c("LN", "THY", "TON")]
  spe. = spe.[, which(spe.$CTCL.CAILS == evolution)]
  
  interaction_list = list()
  
  for(celltype1 in c("T_helper", "T_cytotoxic", "T_regulatory", "NKT", "Macrophages", "Monocytes", "APC", "Monocytic_Lineage", "B_cell")){
    
    cat(celltype1, ": ")
    for(samp in unique(spe.$sample_id)){
      spe.. = spe.[,spe.$sample_id == samp]
      
      cells1 = spe..$cell_id[which(spe..$cell_type_CellSighter == celltype1)]
      
      if(length(cells1) > 1){
        cat(".")
        mat = spe..@metadata$interaction[cells1,spe..$cell_id]
        mat[mat>0]=1
        types = as.factor(spe..$cell_type_simplified)
        
        types_tab = data.frame(
          sample_id = samp,
          batch = unique(spe..$batch),
          type = levels(types),
          sum = rowSums(apply(mat, 1, function(i)
            as.data.frame(table(types[i > 0]))[,2]
          )),
          total = as.numeric(table(types))
        )
        interaction_list[[samp]] = types_tab
      }
    }
    cat("\n")
    df = do.call("rbind", interaction_list)
    
    if(!is.null(df)){
      
      df = df %>% mutate(weighted_sum = sum / total)
      df = df %>% group_by(sample_id) %>%
        mutate(percent = 100 * weighted_sum / sum(weighted_sum))
      df = df %>% group_by(type) %>%
        summarise(percent = mean(percent))
      
      data = df
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
      
      png(file.path(output_dir, "CAILS", "Interaction_Map", paste0("Interaction_map_non_infiltrate_CAILS_",gsub("/", "_", evolution),"_", celltype1,".png")), width = 1300, height = 1200, res = 150)
      
      colors = cell_type_simplified_color_df
      
      # Make the plot
      print(ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=type)) +
              geom_rect() +
              geom_text( x=2, aes(y=labelPosition, label=label, color=type), size=6) + # x here controls label position (inner / outer)
              scale_fill_manual(values = colors) +
              scale_color_manual(values = colors) +
              annotate("text", x = -1, y = 0, label =paste0(evolution, " - ", celltype1), size = 6) +
              coord_polar(theta="y") +
              xlim(c(-1, 4)) +
              theme_void() +
              theme(legend.position = "none"))
      dev.off()
      
    }
    
  }
}


#################################################################################
# Saturation curve  CAILS Regression vs Non Regression
################################################################################# 
list_distances = qs::qread("output/cell_neighborhood/list_distances.qs")
for(type in struct_celltype){
  
  stage_1 = c("Stable/Progression")
  stage_2 = c("Regression")
  
  spe. = spe[, which(spe$CTCL.CAILS %in% c(stage_1, stage_2))]
  dir.create(file.path(output_dir, "CAILS", "Saturation_Curve", type))
  
  for(type2 in unique(spe.$cell_type_CellSighter)){
    cat(type2, ": ")
    list_type2_to_type = list()
    
    for(samp in unique(spe.$sample_id[which(spe.$CTCL.CAILS %in% c(stage_1, stage_2))])){
      cat(".")
      d = list_distances[[samp]]
      cells2 = spe.$cell_id[which(spe.$sample_id == samp & spe.$cell_type_CellSighter %in% type2)]
      cells = spe.$cell_id[which(spe.$sample_id == samp & spe.$cell_type_CellSighter %in% type)]
      
      if(length(cells) > 1 & length(cells2) > 1){
        d_r = d[cells2, cells] 
        
        list_increasing_distance = list()
        
        for (dist in seq(0,1000,50)) {
          
          within_range = apply(d_r, 2, function(i) which(i < dist), simplify = FALSE)
          
          list_increasing_distance[[as.character(dist)]] = 
            data.frame(n_active = length(which(sapply(within_range, length)>0)),
                       n_tot = length(cells),
                       sample_id = samp, stage = unique(spe.$CTCL.CAILS[spe.$sample_id == samp]),
                       dist = dist)
        }
      }
      list_type2_to_type[[samp]] = do.call("rbind", list_increasing_distance)
      
    }
    cat("\n")
    saturation_df =  do.call("rbind", list_type2_to_type)
    
    saturation_df = saturation_df %>% mutate(percent = 100 * n_active / n_tot)
    saturation_df_sum <- summarySE(saturation_df, measurevar="percent", groupvars=c("dist", "stage"))
    
    png(file.path(output_dir, "CAILS", "Saturation_Curve", type, paste0("Saturation_curve_",type,"within_reach_of_",type2,"_CAILS.png")), width = 2000, height = 1300, res = 300)
    saturation_df_sum = saturation_df_sum  %>% 
      mutate(stage = factor(stage))
    t_test_result <- t.test(percent ~ stage, data = saturation_df_sum, paired = TRUE)
    
    p = saturation_df_sum %>% ggplot(aes(x=dist, y=percent, colour=stage)) + 
      geom_errorbar(aes(ymin=percent-se, ymax=percent+se), width=.1) +
      geom_line() +
      geom_point()  +
      theme_classic() + xlab(paste0("Distance of cells within range")) + ylab(paste0("% of ", type, " within range of ", type2)) +
      ggtitle(paste0("% of ", type, " within range of ", type2)) + theme(axis.text.x = element_text(angle = 90)) + 
      annotate("text", x = mean(saturation_df_sum$dist) + 700, y = max(saturation_df_sum$percent) + 15, 
               label = paste("Wilcoxon Paired Test - p =", signif(t_test_result$p.value, digits = 4)), 
               hjust = 1, vjust = 1, color = "black") + 
      annotate("text", x = min(saturation_df_sum$dist) + 400, y = max(saturation_df_sum$percent) +15, 
               label = paste0(type, " (n) = ", length(which(spe.$cell_type_CellSighter[spe.$CTCL.CAILS %in% c(stage_1, stage_2)] %in% type))), 
               hjust = 1, vjust = 1, color = "black") + 
      annotate("text", x = min(saturation_df_sum$dist) + 400, y = max(saturation_df_sum$percent) +10, 
               label = paste0("Type 2", " (n) = ", length(which(spe.$cell_type_CellSighter[spe.$CTCL.CAILS %in% c(stage_1, stage_2)] %in% type2))), 
               hjust = 1, vjust = 1, color = "black") +
      scale_color_manual(values = colors_progression_regression)
    print(p)
    dev.off()
    
  }
}


#################################################################################
# Saturation curve  CAILS Regression vs Non Regression - Grouping cell types
################################################################################# 

type2_list = list("Mono_Mac" = c("Macrophages", "Monocytes", "Monocytic_Lineage","APC"))

for(type in c("T_helper", "T_cytotoxic", "T_regulatory")){
  
  stage_1 = c("Stable/Progression")
  stage_2 = c("Regression")
  
  spe. = spe[, which(spe$CTCL.CAILS %in% c(stage_1, stage_2))]
  dir.create(file.path(output_dir, "CAILS", "Saturation_Curve", type))
  
  for(type2 in names(type2_list)){
    cat(type2, ": ")
    
    list_type2_to_type = list()
    
    for(samp in unique(spe.$sample_id[which(spe.$CTCL.CAILS %in% c(stage_1, stage_2))])){
      cat(".")
      d = list_distances[[samp]]
      cells2 = spe.$cell_id[which(spe.$sample_id == samp & spe.$cell_type_CellSighter %in% type2_list[[type2]])]
      cells = spe.$cell_id[which(spe.$sample_id == samp & spe.$cell_type_CellSighter %in% type)]
      
      if(length(cells) > 1 & length(cells2) > 1){
        d_r = d[cells2, cells] 
        
        list_increasing_distance = list()
        
        for (dist in seq(0,1000,50)) {
          
          within_range = apply(d_r, 2, function(i) which(i < dist), simplify = FALSE)
          
          list_increasing_distance[[as.character(dist)]] = 
            data.frame(n_active = length(which(sapply(within_range, length)>0)),
                       n_tot = length(cells),
                       sample_id = samp, stage = unique(spe.$CTCL.CAILS[spe.$sample_id == samp]),
                       dist = dist)
        }
      }
      list_type2_to_type[[samp]] = do.call("rbind", list_increasing_distance)
      
    }
    cat("\n")
    saturation_df =  do.call("rbind", list_type2_to_type)
    
    saturation_df = saturation_df %>% mutate(percent = 100 * n_active / n_tot)
    saturation_df_sum <- summarySE(saturation_df, measurevar="percent", groupvars=c("dist", "stage"))
    
    png(file.path(output_dir, "CAILS", "Saturation_Curve", type, paste0("Saturation_curve_",type,"within_reach_of_",type2,"_CAILS.png")), width = 2000, height = 1300, res = 300)
    saturation_df_sum = saturation_df_sum  %>% 
      mutate(stage = factor(stage))
    t_test_result <- t.test(percent ~ stage, data = saturation_df_sum, paired = TRUE)
    
    p = saturation_df_sum %>% ggplot(aes(x=dist, y=percent, colour=stage)) + 
      geom_errorbar(aes(ymin=percent-se, ymax=percent+se), width=.1) +
      geom_line() +
      geom_point()  +
      theme_classic() + xlab(paste0("Distance of cells within range")) + ylab(paste0("% of ", type, " within range of ", type2)) +
      ggtitle(paste0("% of ", type, " within range of ", type2)) + theme(axis.text.x = element_text(angle = 90)) + 
      annotate("text", x = mean(saturation_df_sum$dist) + 700, y = max(saturation_df_sum$percent) + 15, 
               label = paste("Wilcoxon Paired Test - p =", signif(t_test_result$p.value, digits = 4)), 
               hjust = 1, vjust = 1, color = "black") + 
      annotate("text", x = min(saturation_df_sum$dist) + 400, y = max(saturation_df_sum$percent) +15, 
               label = paste0(type, " (n) = ", length(which(spe.$cell_type_CellSighter[spe.$CTCL.CAILS %in% c(stage_1, stage_2)] %in% type))), 
               hjust = 1, vjust = 1, color = "black") + 
      annotate("text", x = min(saturation_df_sum$dist) + 400, y = max(saturation_df_sum$percent) +10, 
               label = paste0("Type 2", " (n) = ", length(which(spe.$cell_type_CellSighter[spe.$CTCL.CAILS %in% c(stage_1, stage_2)] %in% type2))), 
               hjust = 1, vjust = 1, color = "black") +
      scale_color_manual(values = colors_progression_regression)
    print(p)
    dev.off()
    
  }
}



###############################################################################
# Compare markers -  CAILS Regression vs Non Regression
###############################################################################
spe. = spe[, !spe$condition %in% c("THY", "TON", "LN") & !is.na(spe$CTCL.CAILS)]
spe. = spe.[,spe.$CTCL.CAILS %in% c("Stable/Progression", "Regression")]

intensities = as.data.frame(t(spe.@assays@data$CellSighter_marker_mat))
intensities$sample_id = spe.$sample_id
intensities$CTCL.CAILS = spe.$CTCL.CAILS
intensities$batch = spe.$batch
intensities$cell_type_CellSighter = spe.$cell_type_CellSighter
intensities = intensities %>% dplyr::group_by(sample_id, CTCL.CAILS, batch, cell_type_CellSighter) %>% 
  dplyr::summarise(across(Actin:Ki67, mean))

dir.create(file.path(output_dir, "CAILS", "Markers"))
for(celltype in unique(spe.$cell_type_CellSighter)){
  print(celltype)
  
  pdf(file.path(output_dir, "CAILS", "Markers", paste0("Functional_Marker_in_cell_type_CellSighter_",celltype,".pdf")), height = 4, width = 4)
  for(mark in colnames(intensities.)[5:65]){
    intensities. = intensities %>% filter(cell_type_CellSighter %in% celltype)
    intensities.[[mark]] = 100 * intensities.[[mark]]
    
    p = grouped_dotplot(intensities., y =  mark, categ1 = "CTCL.CAILS", categ2 = NULL,
                        ref.group = "Stable/Progression", add_violin = T,
                        color_categ1 = colors_progression_regression)
    p = p + xlab("") + ylab(paste0(celltype, " expressing ", mark, " (%)")) +
      guides(fill="none", color = "none")
    print(p)
  }
  dev.off()
}


##############################################################################
# Signature - CTCL CAILS
##############################################################################

signatures = list(
  migration = c("Galectin3", "CLA", "Galectin9", "CD44", "CD181", "CD65", "CD162"),
  exhaustion = c("CD134", "CD305", "TIM3"),
  galectins = c("Galectin3", "Galectin9", "CLA")
)

# Co-occurence score migration
library(eulerr)
pdf(file.path(output_dir, paste0("Venn_T_cytotoxic_all.pdf")))
mat = spe.@assays@data$CellSighter_marker_mat[signatures[[sig]], which(spe.$cell_type_CellSighter == "T_cytotoxic" &
                                                                         spe.$CTCL.CAILS == "Regression")]
plot(euler(t(mat)))
dev.off()

pdf(file.path(output_dir, paste0("Venn_T_cytotoxic_Inflammatory.pdf")))
mat = spe.@assays@data$CellSighter_marker_mat[signatures[[sig]],which(spe.$cell_type_CellSighter == "T_cytotoxic" & 
                                                                        spe.$CTCL.CAILS == "Inflammatory")]
plot(venn(t(mat)[,1:5]))
dev.off()


pdf(file.path(output_dir, paste0("Venn_T_cytotoxic_Responder.pdf")))
mat = spe.@assays@data$CellSighter_marker_mat[signatures[[sig]],which(spe.$cell_type_CellSighter == "T_cytotoxic" & 
                                                                        spe.$CTCL.CAILS == "Regression")]
plot(venn(t(mat)))
dev.off()

pdf(file.path(output_dir, paste0("Venn_T_cytotoxic_Non_Responder.pdf")))
mat = spe.@assays@data$CellSighter_marker_mat[signatures[[sig]],which(spe.$cell_type_CellSighter == "T_cytotoxic" & 
                                                                        spe.$CTCL.CAILS == "Stable/Progression")]
plot(venn(t(mat)))
dev.off()

# Co-occurence score
mat1 = colSums(spe.@assays@data$CellSighter_marker_mat[signatures[[sig]],which(spe.$cell_type_CellSighter == "T_cytotoxic" & 
                                                                                 spe.$CTCL.CAILS == "Regression")])
mat2 = colSums(spe.@assays@data$CellSighter_marker_mat[signatures[[sig]],which(spe.$cell_type_CellSighter == "T_cytotoxic" & 
                                                                                 spe.$CTCL.CAILS == "Stable/Progression")])
mat3 = colSums(spe.@assays@data$CellSighter_marker_mat[signatures[[sig]],which(spe.$cell_type_CellSighter == "T_cytotoxic" & 
                                                                                 spe.$CTCL.CAILS == "Inflammatory")])

tab = data.frame(coocurrence = c(mat3, mat1, mat2), type = c(rep("Inflammatory", length(mat3) ),
                                                             rep("Regression", length(mat1) ),
                                                             rep("Stable/Progression", length(mat2))) )

tab = tab %>% group_by(coocurrence, type) %>% 
  summarise(occ = n())

pdf(file.path(output_dir, paste0("Coocurence_migration_markers.pdf")))
tab %>% ggplot(aes(y = occ, x = coocurrence, fill = type)) + geom_bar(stat = "identity",
                                                                      position = position_dodge2()) +
  ylab("Co-occurence migration markers") + xlab("N markers") + theme_bw() + 
  scale_fill_manual(values = colors_progression_regression[1:4])
dev.off() 

for(celltype in unique(spe.$cell_type_CellSighter)){
  
  pdf(file.path(output_dir, "Boxplots_Marker_IPI", paste0("Boxplot_",celltype,".pdf")))
  for(sig in names(signatures)){
    meta = as.data.frame(colData(spe.))
    
    mat = spe.@assays@data$CellSighter_marker_mat[signatures[[sig]],]
    meta$marker = (colSums(spe.@assays@data$CellSighter_marker_mat[signatures[[sig]],]))
    
    mat = as.matrix(spe.@assays@data$CellSighter_marker_mat[signatures[[sig]],])
    mat = mat[,which(spe.$CTCL.CAILS == "Stable/Progression" & spe.$cell_type_CellSighter == celltype)]
    plot(euler(t(mat)))
    
    meta =  meta %>% dplyr::filter(cell_type_CellSighter == celltype, !is.na(CTCL.CAILS)) %>%
      dplyr::group_by(sample_id, CTCL.CAILS, batch) %>% dplyr::summarise(marker = mean(marker))
    
    print(meta %>% ggplot(aes(x = CTCL.CAILS, y = marker, fill = CTCL.CAILS)) +
            geom_boxplot(outlier.colour = "white") + 
            scale_fill_manual(values = colors_progression_regression[1:4]) +
            ggnewscale::new_scale_fill()+
            geom_jitter(aes(fill = batch),  colour="black",pch=21, size=2, width = 0.3) + 
            scale_fill_manual(values = setNames(unique(spe.$batch_color), unique(spe.$batch))) +
            theme_classic() + ggtitle(paste0(celltype, " - ", sig)) + 
            ylab(sig) + xlab("") +
            stat_compare_means(aes(label = after_stat(p.signif)),
                               method = "t.test", ref.group = "Stable/Progression")
    )
  }  
  dev.off()
  
}


# Boxplots - celltypes - infiltrate only
spe. = spe[setdiff(rownames(spe), markers_to_remove), !spe$condition %in% c("HD","THY", "TON", "LN")]
spe.= spe.[,!spe.$cell_type_CellSighter %in% c(
  "Keratinocyte",
  "Lymphatic",
  "Endothelial",
  "Unknown-Stroma"
)]
spe. = spe.[,which(spe.$immune_interaction_community == "infiltrate")]

spe.$CTCL.CAILS = as.character(spe.$CTCL.CAILS)
spe.$CTCL.CAILS[which(spe.$condition == "HD")] = "HD"
spe.$CTCL.CAILS[which(spe.$condition %in% c("LP", "AD", "PS", "DAR") )] = "Inflammatory"

for(i in unique(spe.$cell_type_CellSighter)){
  # Immune cells
  spe.. = spe.
  spe..$cell_type_CellSighter[spe..$cell_type_CellSighter != i] = "Aa"
  meta = as.data.frame(colData(spe..))
  
  meta = meta  %>% filter(!is.na(CTCL.CAILS)) %>%
    dplyr::group_by(sample_id, CTCL.CAILS, condition) %>% 
    dplyr::summarise(percent_CN = 100 * length(
      which(cell_type_CellSighter == i)) / length(cell_type_CellSighter))
  
  dir.create(file.path(output_dir, "Boxplots_IPI"))
  png(file.path(output_dir, "Boxplots_IPI", paste0("Composition_changes_CAILS_", i, ".png")), width = 1200, height = 1000, res = 250)
  print(meta %>% ggplot(aes(x = CTCL.CAILS, y = percent_CN, fill = CTCL.CAILS)) +
          geom_boxplot(outlier.colour = "white") + 
          scale_fill_manual(values = colors_progression_regression) +
          ggnewscale::new_scale("fill") +
          geom_jitter(aes(fill = condition),  colour="black",pch=21, size=2, width = 0.3) +
          scale_fill_manual(values = unique(spe$condition_color[match(meta$condition, spe$condition)])) +
          theme_classic() + ggtitle(paste0(" CellType = ", i)) + 
          ylab(paste0("Percent of ",i," (% of immune cells)")) + xlab("") +
          stat_compare_means(aes(label = after_stat(p.signif)),
                             method = "t.test", ref.group = "Stable/Progression") +
          theme(axis.text.x = element_blank())
  )
  dev.off()
}



###############################################################################
# Compare infiltrate - CAILS - Regression vs Stable
###############################################################################

spe. = spe[setdiff(rownames(spe), markers_to_remove), !spe$condition %in% c("THY", "TON", "LN") & !is.na(spe$CTCL.CAILS)]
spe.= spe.[,!spe.$cell_type_CellSighter %in% c(
  "Keratinocyte",
  "Lymphatic",
  "Endothelial",
  "Unknown-Stroma"
)]
library(ggpubr)

# Immune cells
meta = as.data.frame(colData(spe.))

meta = meta  %>%
  dplyr::group_by(sample_id, CTCL.CAILS , batch) %>% 
  dplyr::summarise(percent_CN = 100 * length(
    which(immune_interaction_community == "infiltrate")) / length(immune_interaction_community))

pdf(file.path(output_dir, "Infiltrate_changes_CAILS.pdf"), width = 6, height = 5)
print(meta %>% ggplot(aes(x = CTCL.CAILS, y = percent_CN, fill = CTCL.CAILS)) +
        geom_boxplot(outlier.colour = "white") + 
        scale_fill_manual(values = colors_progression_regression[3:4]) +
        ggnewscale::new_scale("fill") +
        geom_jitter(aes(fill = batch),  colour="black",pch=21, size=2, width = 0.3) +
        scale_fill_manual(values = unique(spe$batch_color[match(meta$batch, spe$batch)])) +
        theme_classic() + ggtitle("Infiltrate changes in Regression vs Stable") + 
        ylab(paste0("Infiltrate (% of immune cells)")) + xlab("") +
        stat_compare_means(aes(label = after_stat(p.signif)),
                           method = "t.test", ref.group = "Stable/Progression")
)
dev.off()



###############################################################################
# Compare markers -  all conditions Late vs Early
###############################################################################

spe. = spe[setdiff(rownames(spe), markers_to_remove), !spe$condition %in% c("THY", "TON", "LN") & !is.na(spe$stage)]

pdf(file.path(output_dir, paste0("DA_Late_vs_Early_per_celltype.pdf")), width = 6, height = 5)
for(celltype in unique(spe.$cell_type_CellSighter)){
  
  spe_early = spe.[,spe.$stage == "Early"  & spe.$cell_type_CellSighter == celltype]
  early_cells = c()
  # subsample 
  for(early in unique(spe_early$sample_id)){
    cells = spe_early$cell_id[spe_early$sample_id == early]
    early_cells = c(early_cells, sample(cells, min(200, length(cells)), replace = FALSE))
    print( min(200, length(cells)))
  }
  
  spe_late = spe.[, spe.$stage == "Late"  & spe.$cell_type_CellSighter == celltype]
  late_cells = c()
  # subsample 
  for(late in unique(spe_late$sample_id)){
    cells = spe_late$cell_id[spe_late$sample_id == late]
    late_cells = c(late_cells, sample(cells, min(200, length(cells)), replace = FALSE))
    print(min(200, length(cells)))
  }
  if(length(early_cells) > 10 & length(late_cells) > 10){
    
    mat1 = spe.@assays@data$CellSighter_marker_mat[,early_cells]
    mat2 = spe.@assays@data$CellSighter_marker_mat[,late_cells]
    
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
    
    df = data.frame(gene = names(values1), value.HD = values1, 
                    value.condition = values2, percent.HD = percent1,
                    percent.cond = percent2, log2FC = log2FC, 
                    p.value.adjusted = p.values, effect = effect)
    
    df = df %>% mutate(categ = ifelse( p.values < 0.01, ifelse(log2FC > 1, "signifplus", ifelse(log2FC < -1, "signifminus", "unsig")), "unsig"))
    
    p = df  %>%
      ggplot(aes(x = log2FC, y = -log10(p.values), color = categ, label = gene)) +
      geom_point() +
      scale_color_manual(values = c("darkgreen","red","black")) +
      theme_classic() + xlab("Percent - log2FC") + ylab("-log10(adjusted p-value)") + 
      ggtitle(paste0(celltype, " - volcano plot - Late vs Early")) + 
      geom_hline(yintercept = 2, lwd = 0.35, lty = 2, col = "grey20") +
      geom_vline(xintercept = -1, lwd = 0.35, lty = 2, col = "grey20") +
      geom_vline(xintercept = 1,  lwd = 0.35, lty = 2, col = "grey20") +
      ggrepel::geom_text_repel()
    print(p)
    
  }    
}
dev.off()


###############################################################################
# Compare markers - celltype by condition - Late vs Early
###############################################################################

spe. = spe[setdiff(rownames(spe), markers_to_remove), !spe$condition %in% c("THY", "TON", "LN") & !is.na(spe$stage)]

pdf(file.path(output_dir, paste0("DA_Late_vs_Early_per_condition_per_celltype.pdf")), width = 6, height = 5)
for(cond in setdiff(spe.$condition, "HD")){
  for(celltype in unique(spe.$cell_type_CellSighter)){
    
    spe_early = spe.[,spe.$stage == "Early" & spe.$condition == cond & spe.$cell_type_CellSighter == celltype]
    early_cells = c()
    # subsample 
    for(early in unique(spe_early$sample_id)){
      cells = spe_early$cell_id[spe_early$sample_id == early]
      early_cells = c(early_cells, sample(cells, min(200, length(cells)), replace = FALSE))
      print( min(200, length(cells)))
    }
    
    spe_late = spe.[, spe.$stage == "Late" &  spe.$condition %in% cond & spe.$cell_type_CellSighter == celltype]
    late_cells = c()
    # subsample 
    for(late in unique(spe_late$sample_id)){
      cells = spe_late$cell_id[spe_late$sample_id == late]
      late_cells = c(late_cells, sample(cells, min(200, length(cells)), replace = FALSE))
      print(min(200, length(cells)))
    }
    if(length(early_cells) > 10 & length(late_cells) > 10){
      
      mat1 = spe.@assays@data$CellSighter_marker_mat[,early_cells]
      mat2 = spe.@assays@data$CellSighter_marker_mat[,late_cells]
      
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
        ggtitle(paste0(cond, " - ", celltype, " - volcano plot - Late vs Early")) + 
        geom_hline(yintercept = 2, lwd = 0.35, lty = 2, col = "grey20") +
        geom_vline(xintercept = -1, lwd = 0.35, lty = 2, col = "grey20") +
        geom_vline(xintercept = 1,  lwd = 0.35, lty = 2, col = "grey20") +
        ggrepel::geom_text_repel()
      print(p)
      
    }    
  }
}
dev.off()



###############################################################################
# Compare infiltrate - all - Late vs Early
###############################################################################

spe. = spe[setdiff(rownames(spe), markers_to_remove), !spe$condition %in% c("THY", "TON", "LN") & !is.na(spe$stage)]
spe.= spe.[,!spe.$cell_type_CellSighter %in% c(
  "Keratinocyte",
  "Lymphatic",
  "Endothelial",
  "Unknown-Stroma"
)]
library(ggpubr)

# Immune cells
meta = as.data.frame(colData(spe.))

meta = meta  %>%
  dplyr::group_by(sample_id, stage, batch) %>% 
  dplyr::summarise(percent_CN = 100 * length(
    which(immune_interaction_community == "infiltrate")) / length(immune_interaction_community))

pdf(file.path(output_dir, "Infiltrate_changes_Late_vs_early.pdf"), width = 6, height = 5)
print(meta %>% ggplot(aes(x = stage, y = percent_CN, fill = stage)) +
        geom_boxplot(outlier.colour = "white") + 
        scale_fill_manual(values = colors_early_late) +
        ggnewscale::new_scale("fill") +
        geom_jitter(aes(fill = batch),  colour="black",pch=21, size=2, width = 0.3) +
        scale_fill_manual(values = unique(spe$batch_color[match(meta$batch, spe$batch)])) +
        theme_classic() + ggtitle("Infiltrate changes in Late vs Early") + 
        ylab(paste0("Infiltrate (% of immune cells)")) + xlab("") +
        stat_compare_means(aes(label = after_stat(p.signif)),
                           method = "t.test", ref.group = "Early")
)
dev.off()


###############################################################################
# Compare infiltrate - Tumor - Late vs Early
###############################################################################

spe. = spe[setdiff(rownames(spe), markers_to_remove), spe$condition %in% c("MF", "SS") & !is.na(spe$stage) ]
spe.= spe.[,!spe.$cell_type_CellSighter %in% c(
  "Keratinocyte",
  "Lymphatic",
  "Endothelial",
  "Unknown-Stroma"
)]
library(ggpubr)

# Immune cells
meta = as.data.frame(colData(spe.))

meta = meta  %>%
  dplyr::group_by(sample_id, stage, condition) %>% 
  dplyr::summarise(percent_CN = 100 * length(
    which(immune_interaction_community == "infiltrate")) / length(immune_interaction_community))

pdf(file.path(output_dir, "Infiltrate_changes_Late_vs_early_tumor.pdf"), width = 6, height = 5)
print(meta %>% ggplot(aes(x = stage, y = percent_CN, fill = stage)) +
        geom_boxplot(outlier.colour = "white") + 
        scale_fill_manual(values = colors_early_late) +
        ggnewscale::new_scale("fill") +
        geom_jitter(aes(fill = condition),  colour="black",pch=21, size=2, width = 0.3) +
        scale_fill_manual(values = setNames(unique(spe$condition_color),unique(spe$condition)) ) +
        theme_classic() + ggtitle("Infiltrate changes in Late vs Early") + 
        ylab(paste0("Infiltrate (% of immune cells)")) + xlab("") +
        stat_compare_means(aes(label = after_stat(p.signif)),
                           method = "t.test", ref.group = "Early")
)
dev.off()


###############################################################################
# Compare infiltrate -Non Tumor - Late vs Early
###############################################################################

spe. = spe[setdiff(rownames(spe), markers_to_remove), spe$condition %in% c("AD", "PS", "DAR", "LP") & !is.na(spe$stage) ]
spe.= spe.[,!spe.$cell_type_CellSighter %in% c(
  "Keratinocyte",
  "Lymphatic",
  "Endothelial",
  "Unknown-Stroma"
)]
library(ggpubr)

# Immune cells
meta = as.data.frame(colData(spe.))

meta = meta  %>%
  dplyr::group_by(sample_id, stage, condition) %>% 
  dplyr::summarise(percent_CN = 100 * length(
    which(immune_interaction_community == "infiltrate")) / length(immune_interaction_community))

pdf(file.path(output_dir, "Infiltrate_changes_Late_vs_early_non_tumor.pdf"), width = 6, height = 5)
print(meta %>% ggplot(aes(x = stage, y = percent_CN, fill = stage)) +
        geom_boxplot(outlier.colour = "white") + 
        scale_fill_manual(values = colors_early_late) +
        ggnewscale::new_scale("fill") +
        geom_jitter(aes(fill = condition),  colour="black",pch=21, size=2, width = 0.3) +
        scale_fill_manual(values = setNames(unique(spe$condition_color),unique(spe$condition)) ) +
        theme_classic() + ggtitle("Infiltrate changes in Late vs Early") + 
        ylab(paste0("Infiltrate (% of immune cells)")) + xlab("") +
        stat_compare_means(aes(label = after_stat(p.signif)),
                           method = "t.test", ref.group = "Early")
)
dev.off()




