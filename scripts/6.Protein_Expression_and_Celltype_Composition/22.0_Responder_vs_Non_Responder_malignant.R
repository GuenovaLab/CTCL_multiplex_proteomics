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




