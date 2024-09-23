# Analyzing distances between cell types and plotting "Saturation Curves"
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
              "ggpubr")
suppressPackageStartupMessages(invisible(lapply(libraries, require, character.only = TRUE)))
setwd("/mnt/RECHERCHE/GUENOVA_LAB/Project_Multiplex_Phenotyping/Miltenyi/Muliplex_Imaging_Pipeline/")

source("scripts/MultiplexImaging_utils.R")
source("scripts/GlobalVars.R")

# Change when running from R:
args = list(output = "output/")

cat("Output = ", args$output, "\n")

output_dir = file.path(args$output, "cell_interactions", "Saturation_Curves")
if (!dir.exists(output_dir))
  dir.create(output_dir)

spe = qs::qread("output/SpatialExperiment.qs")
list_distances = qs::qread("output/cell_neighborhood/list_distances.qs")

################################################################################
# Average distance between T helper and other celltypes
################################################################################

for(type in unique(spe$cell_type_CellSighter)){
  list_T_help_to_type = list()
  for(samp in unique(spe$sample_id)){
    d = list_distances[[samp]]
    if(any(colnames(d) %in% rownames(d))){
      diag(d) = Inf
    }
    T_helpers = spe$cell_id[which(spe$sample_id == samp & spe$cell_type_CellSighter == "T_helper")]
    cells = spe$cell_id[which(spe$sample_id == samp & spe$cell_type_CellSighter == type)]
    if(length(cells) > 1){
      d = d[T_helpers, cells] 
      
      closest = apply(d, 1, which.min)
      
      dist_closest = rep(0, nrow(d))
      for(i in 1:nrow(d)){
        dist_closest[i] = d[i,closest[i]]
      }
      
      list_T_help_to_type[[samp]] = 
        data.frame(dist_closest = dist_closest, sample_id = samp, condition = unique(spe$condition[spe$sample_id == samp]))
    }
  }
  T_help_type_df = do.call("rbind", list_T_help_to_type)
  
  png(file.path(output_dir, paste0("Average_Distance_Closest_T_helper_to_",type,".png")), width = 3500, height = 1600, res = 300)
  p =T_help_type_df %>% ggplot(aes(x = sample_id, y = log2(dist_closest), fill = sample_id)) +
    geom_boxplot(outlier.size = 0.5) + facet_wrap(~condition,  scales = "free_x", nrow = 1) + 
    scale_fill_manual(values = unique(spe$sample_id_color[match(T_help_type_df$sample_id, spe$sample_id)])) +
    theme_classic() + xlab("") + ylab(paste0("Log2(distance from T helper to closest ", type, ")")) +
    NoLegend()  + ggtitle(type) + theme(axis.text.x = element_text(angle = 90))
  print(p)
  dev.off()
  
}

################################################################################
# % of cell type between T helper and others
################################################################################

for(type in unique(spe$cell_type_CellSighter)){
  list_T_help_to_type = list()
  
  for(samp in unique(spe$sample_id)){
    d = list_distances[[samp]]
    T_helpers = spe$cell_id[which(spe$sample_id == samp & spe$cell_type_CellSighter == "T_helper")]
    cells = spe$cell_id[which(spe$sample_id == samp & spe$cell_type_CellSighter == type)]
    if(length(cells) > 1){
      d = d[T_helpers, cells] 
      within_range = apply(d, 2, function(i) which(i < max_dist))
      
      list_T_help_to_type[[samp]] = 
        data.frame(n_active = length(which(sapply(within_range, length)>0)),
                   n_tot = length(cells),
                   sample_id = samp, condition = unique(spe$condition[spe$sample_id == samp]))
    }
  }
  Type_to_T_helper_df = do.call("rbind", list_T_help_to_type)
  
  png(file.path(output_dir, paste0("Percentage_",type,"_within_reach_of_T_helper.png")), width = 3500, height = 1600, res = 300)
  p =Type_to_T_helper_df %>% ggplot(aes(x = condition, y =100 * n_active / n_tot, fill = condition)) +
    geom_boxplot(outlier.size = 0.5) +
    scale_fill_manual(values = unique(spe$condition_color[match(Type_to_T_helper_df$condition, spe$condition)])) +
    theme_classic() + xlab("") + ylab(paste0("% of ", type, " within range of T helper")) +
    NoLegend()  + ggtitle(type) + theme(axis.text.x = element_text(angle = 90))  +
    stat_compare_means(aes(label = after_stat(p.signif)),
                       method = "t.test", ref.group = "HD")
  print(p)
  dev.off()
  
}

################################################################################
# Per condition Saturation curve of type1 to any type2
################################################################################
spe = qs::qread("output/SpatialExperiment.qs")
spe = spe[, spe$condition %in% c("PS", "AD", "MF", "SS")]
dir.create(file.path(output_dir, "AD.PS"))

for(type in c("NKT", "T_cytotoxic", "T_helper")){
  print(type)
  print("")
  
  for(type2 in unique(spe$cell_type_CellSighter)){
    print(type2)
    print("")
    
    list_type2_to_type = list()
    for(samp in unique(spe$sample_id)){
      d = list_distances[[samp]]
      cells2 = spe$cell_id[which(spe$sample_id == samp & spe$cell_type_CellSighter == type2)]
      cells = spe$cell_id[which(spe$sample_id == samp & spe$cell_type_CellSighter == type)]
      if(length(cells) > 1 & length(cells2) > 1){
        d = d[cells2, cells] 
        if(type == type2)
          diag(d) = Inf
        
        list_increasing_distance = list()
        for (dist in seq(0,1000,50)) {
          print(dist)
          
          within_range = apply(d, 2, function(i) which(i < dist), simplify = FALSE)
          
          list_increasing_distance[[as.character(dist)]] = 
            data.frame(n_active = length(which(sapply(within_range, length)>0)),
                       n_tot = length(cells),
                       sample_id = samp, condition = unique(spe$condition[spe$sample_id == samp]),
                       dist = dist)
        }
      }
      list_type2_to_type[[samp]] = do.call("rbind", list_increasing_distance)
    }
    saturation_df =  do.call("rbind", list_type2_to_type)
    
    
    saturation_df = saturation_df %>% mutate(percent = 100 * n_active / n_tot)
    saturation_df_sum <- summarySE(saturation_df, measurevar="percent", groupvars=c("dist", "condition"))
    
    
    png(file.path(output_dir, "AD.PS", paste0("Saturation_curve_percentage_",type,"within_reach_of_",type2,".png")), width = 1500, height = 1200, res = 300)
    p = saturation_df_sum  %>% filter(!condition %in% c("TON", "THY", "LN"))  %>% 
      ggplot(aes(x=dist, y=percent, colour=condition)) + 
      geom_errorbar(aes(ymin=percent-se, ymax=percent+se), width=.1) +
      geom_line() +
      geom_point()  +
      scale_color_manual(values = setNames(unique(spe$condition_color), unique(spe$condition))) +
      theme_classic() + xlab(paste0("Distance of cells within range")) + ylab(paste0("% of ", type, " within range of ", type2)) +
      ggtitle(type) + theme(axis.text.x = element_text(angle = 90))
    print(p)
    dev.off()
    
  }
}

################################################################################
# FOCUS - Per condition Saturation curve of type1 to any type2
################################################################################

condition_1 = c("MF","SS")
condition_2 = c("AD", "PS", "DAR", "LP")
dir.create(file.path(output_dir, "Benign_vs_Malign"))

for(type in c("T_helper", "T_cytotoxic", "T_regulatory", "NKT",
              "Macrophages", "Monocytic_Lineage", "APC", "Monocytes", "B_cell",
              "pDC")){
  for(type2 in unique(spe$cell_type_CellSighter)){
    
    list_type2_to_type = list()
    for(samp in unique(spe$sample_id[which(spe$condition %in% c(condition_1, condition_2))])){
      d = list_distances[[samp]]
      cells2 = spe$cell_id[which(spe$sample_id == samp & spe$cell_type_CellSighter %in% type2)]
      cells = spe$cell_id[which(spe$sample_id == samp & spe$cell_type_CellSighter %in% type)]
      if(length(cells) > 1 & length(cells2) > 1){
        d = d[cells2, cells] 
        
        list_increasing_distance = list()
        for (dist in seq(0,1000,50)) {
          print(dist)
          
          within_range = apply(d, 2, function(i) which(i < dist), simplify = FALSE)
          
          list_increasing_distance[[as.character(dist)]] = 
            data.frame(n_active = length(which(sapply(within_range, length)>0)),
                       n_tot = length(cells),
                       sample_id = samp, condition = unique(spe$condition[spe$sample_id == samp]),
                       dist = dist)
        }
      }
      list_type2_to_type[[samp]] = do.call("rbind", list_increasing_distance)
    }
    
    saturation_df =  do.call("rbind", list_type2_to_type)
    saturation_df = saturation_df %>% mutate(percent = 100 * n_active / n_tot)
    saturation_df$class = "Malign"
    saturation_df$class[saturation_df$condition %in% condition_2] = "Benign"
    saturation_df_sum <- summarySE(saturation_df, measurevar="percent", groupvars=c("dist", "class"))
    
    png(file.path(output_dir, "Benign_vs_Malign", paste0("Saturation_curve_",type, "within_reach_of_", type2, "_Benign_vs_Malign.png")),
        width = 2500, height = 1600, res = 300)
    
    saturation_df_sum = saturation_df_sum  %>% 
      mutate(condition = factor(class))
    t_test_result <- t.test(percent ~ class, data = saturation_df_sum, paired = TRUE)
    
    p = saturation_df_sum %>% ggplot(aes(x=dist, y=percent, colour=condition)) + 
      geom_errorbar(aes(ymin=percent-se, ymax=percent+se), width=.1) +
      geom_line() +
      geom_point()  +
      scale_color_manual(values = setNames(unique(spe$disease_color), unique(spe$disease))) + 
      theme_classic() + xlab(paste0("Distance of cells within range")) + ylab(paste0("% of ", type, " within range of ", type2)) +
      ggtitle(paste0("% of ", type, " within range of ", type2)) + theme(axis.text.x = element_text(angle = 90)) + 
      annotate("text", x = mean(saturation_df_sum$dist) + 700, y = max(saturation_df_sum$percent) + 15, 
               label = paste("Wilcoxon Paired Test - p =", signif(t_test_result$p.value, digits = 4)), 
               hjust = 1, vjust = 1, color = "black") + 
      annotate("text", x = min(saturation_df_sum$dist) + 700, y = max(saturation_df_sum$percent) +15, 
               label = paste0(type, " (n) = ", length(which(spe$cell_type_CellSighter[spe$condition %in% c(condition_1, condition_2)] %in% type))), 
               hjust = 1, vjust = 1, color = "black") + 
      annotate("text", x = min(saturation_df_sum$dist) + 700, y = max(saturation_df_sum$percent) +10, 
               label = paste0("Type 2", " (n) = ", length(which(spe$cell_type_CellSighter[spe$condition %in% c(condition_1, condition_2)] %in% type2))), 
               hjust = 1, vjust = 1, color = "black")
    print(p)
    dev.off()
    
    
  }
}

################################################################################
# FOCUS - All Saturation curve of type1 to any type2 classed by expression
################################################################################

type1 = "Endothelial"

# type2 = unique(spe$cell_type_CellSighter)
# type2 = setdiff(type2, type1)
type2 = "T_helper"

spe. = spe[, !spe$condition %in% c("TON", "LN", "THY")]
spe. = spe.[, spe.$disease %in% c("Malign")]
saturation_df_list = list()
max_d = 1000

for(marker in c("CD5")){ # setdiff(rownames(spe.), markers_to_remove)
  
  list_type2_to_type = list()
  
  for(samp in unique(spe.$sample_id)){
    print(samp)
    d = list_distances[[samp]]
    
    spe_samp = spe.[,spe.$sample_id == samp]
    spe_celltype1 = spe_samp[,spe_samp$cell_type_CellSighter %in% type1]
    spe_celltype2 = spe_samp[,spe_samp$cell_type_CellSighter %in% type2]
    
    cells1 = spe_celltype1$cell_id
    cells2_pos = spe_celltype2$cell_id[which(spe_celltype2@assays@data$CellSighter_marker_mat[marker,] == 1)]
    cells2_neg = spe_celltype2$cell_id[which(spe_celltype2@assays@data$CellSighter_marker_mat[marker,] == 0)]
    
    if(length(cells1) > 1 & length(cells2_pos) > 1 & length(cells2_neg) > 1){
      
      list_increasing_distance = list()
      
      for(categ in c("pos", "neg")){
        
        if(categ == "pos") cells2 = cells2_pos
        if(categ == "neg") cells2 = cells2_neg
        
        d. = d[cells1, cells2] 
        
        for (dist in seq(0,max_d,50)) {
          within_range = apply(d., 2, function(i) which(i < dist), simplify = FALSE)
          
          list_increasing_distance[[paste0(as.character(dist), "_", categ)]] = 
            data.frame(n_active = length(which(sapply(within_range, length)>0)),
                       n_tot = length(cells2),
                       sample_id = samp,
                       condition = unique(spe_samp$condition),
                       dist = dist,
                       categ = categ)
        }
      }
      list_type2_to_type[[samp]] = do.call("rbind", list_increasing_distance)
    }
    
  }
  
  saturation_df =  do.call("rbind", list_type2_to_type)
  saturation_df = saturation_df %>% mutate(percent = 100 * n_active / n_tot)
  saturation_df_sum <- summarySE(saturation_df, measurevar="percent", groupvars=c("dist", "categ"))
  
  if(length(type2) == 1 & length(type1) == 1){
    png(file.path(output_dir,  "Marker_based", 
                  paste0("Saturation_curve_",type2,"within_reach_of_",type1,"_",marker,".png")), width = 2500, height = 1600, res = 300)
    saturation_df_sum = saturation_df_sum  %>% 
      mutate(condition = factor(categ))
    t_test_result <- t.test(saturation_df_sum$percent[saturation_df_sum$categ == "pos"],
                            saturation_df_sum$percent[saturation_df_sum$categ == "neg"],
                            data = saturation_df_sum, paired = TRUE)
    
    p = saturation_df_sum %>% ggplot(aes(x=dist, y=percent, colour=condition)) + 
      geom_errorbar(aes(ymin=percent-se, ymax=percent+se), width=.1) +
      geom_line() +
      geom_point()  +
      theme_classic() + xlab(paste0("Distance of cells within range")) + 
      ylab(paste0("% of ", type2," cells expressing ", marker, " within range of ", type1)) +
      ggtitle(paste0("Distance from ", type1," depending on ", marker, " expression")) + 
      theme(axis.text.x = element_text(angle = 90)) + 
      scale_color_manual(values = c("forestgreen", "red"), name = marker, labels = c("-", "+")) + 
      annotate("text", x = mean(saturation_df_sum$dist) + 700, y = max(saturation_df_sum$percent) + 15, 
               label = paste("Wilcoxon Paired Test - p =", signif(t_test_result$p.value, digits = 4)), 
               hjust = 1, vjust = 1, color = "black") + 
      annotate("text", x = min(saturation_df_sum$dist) + 700, y = max(saturation_df_sum$percent) +15, 
               label = paste0(marker, "- (n) = ", length(which(spe.$cell_type_CellSighter %in% type2 &
                                                                 spe.@assays@data$CellSighter_marker_mat[marker,] == 0))), 
               hjust = 1, vjust = 1, color = "black") + 
      annotate("text", x = min(saturation_df_sum$dist) + 700, y = max(saturation_df_sum$percent) +10, 
               label = paste0(marker, "+ (n) = ", length(which(spe.$cell_type_CellSighter %in% type2 &
                                                                 spe.@assays@data$CellSighter_marker_mat[marker,] == 1))), 
               hjust = 1, vjust = 1, color = "black")
    print(p)
    dev.off()
    
  } else {
    png(file.path(output_dir, "Marker_based", 
                  paste0("Saturation_curve_non",type1,"_within_reach_of_",type1,"_",marker,".png")), width = 2500, height = 1600, res = 300)
    saturation_df_sum = saturation_df_sum  %>% 
      mutate(condition = factor(categ))
    t_test_result <- t.test(percent ~ categ, data = saturation_df_sum, paired = TRUE)
    
    p = saturation_df_sum %>% ggplot(aes(x=dist, y=percent, colour=condition)) + 
      geom_errorbar(aes(ymin=percent-se, ymax=percent+se), width=.1) +
      geom_line() +
      geom_point()  +
      theme_classic() + xlab(paste0("Distance of cells within range")) + 
      ylab(paste0("% of All non-", type1," cells within range of ", type1)) +
      ggtitle(paste0("Distance from ", type1," depending on ", marker, " expression")) + 
      theme(axis.text.x = element_text(angle = 90)) + 
      scale_color_manual(values = c("forestgreen", "red"), name = marker, labels = c("-", "+")) + 
      annotate("text", x = mean(saturation_df_sum$dist) + 700, y = max(saturation_df_sum$percent) + 15, 
               label = paste("Wilcoxon Paired Test - p =", signif(t_test_result$p.value, digits = 4)), 
               hjust = 1, vjust = 1, color = "black") + 
      annotate("text", x = min(saturation_df_sum$dist) + 700, y = max(saturation_df_sum$percent) +15, 
               label = paste0(marker, "- (n) = ", length(which(spe.$cell_type_CellSighter %in% type2 &
                                                                 spe.@assays@data$CellSighter_marker_mat[marker,] == 0))), 
               hjust = 1, vjust = 1, color = "black") + 
      annotate("text", x = min(saturation_df_sum$dist) + 700, y = max(saturation_df_sum$percent) +10, 
               label = paste0(marker, "+ (n) = ", length(which(spe.$cell_type_CellSighter %in% type2 &
                                                                 spe.@assays@data$CellSighter_marker_mat[marker,] == 1))), 
               hjust = 1, vjust = 1, color = "black")
    print(p)
    dev.off()
    
    saturation_df$marker = marker
    saturation_df_list[[paste0(type1, "_", marker)]] = saturation_df
  }
}

dsaturation_df_markers = do.call("rbind", saturation_df_list)
qs::qsave(saturation_df_markers, file.path(output_dir, "Saturation_Curves", "Marker_based", "saturation_df_markers.qs"))


################################################################################
# FOCUS - All Saturation curve of type1 to any type2 classed by expression -
#  Malign vs Benign
################################################################################

type1 = "Lymphocyte"
type2 = "Myeloid"

saturation_df_list = list()
max_d = 1000

for(disease in c("Malign", "Benign")){
  spe. = spe[, !spe$condition %in% c("TON", "LN", "THY", "HD")]
  spe. = spe.[, which(spe.$disease == disease)]
  
  for(marker in c( "CD195", "CD196", "CXCR4", "CXCL12","CD305", "CD65", "CD162", "CLA",  "ZAP70", "HLADR",
                   "HLAABC", "Bcl2", "Slug", "Galectin3", "Galectin9")){ 
    cat(marker, ": ")
    list_type2_to_type = list()
    
    for(samp in unique(spe.$sample_id)){
      cat(".")
      d = list_distances[[samp]]
      
      spe_samp = spe.[,spe.$sample_id == samp]
      if(type1 %in% unique(spe$cell_type_simplified)) 
        spe_celltype1 = spe_samp[,spe_samp$cell_type_simplified %in% type1] else
          spe_celltype1 = spe_samp[,spe_samp$cell_type_CellSighter %in% type1]
      if(type2 %in% unique(spe$cell_type_simplified)) 
        spe_celltype2 = spe_samp[,spe_samp$cell_type_simplified %in% type2] else
          spe_celltype2 = spe_samp[,spe_samp$cell_type_CellSighter %in% type2]
      
      
      cells1 = spe_celltype1$cell_id
      cells2_pos = spe_celltype2$cell_id[which(spe_celltype2@assays@data$CellSighter_marker_mat[marker,] == 1)]
      cells2_neg = spe_celltype2$cell_id[which(spe_celltype2@assays@data$CellSighter_marker_mat[marker,] == 0)]
      
      if(length(cells1) > 1 & length(cells2_pos) > 1 & length(cells2_neg) > 1){
        
        list_increasing_distance = list()
        
        for(categ in c("pos", "neg")){
          
          if(categ == "pos") cells2 = cells2_pos
          if(categ == "neg") cells2 = cells2_neg
          
          d. = d[cells1, cells2] 
          
          for (dist in seq(0,max_d,50)) {
            within_range = apply(d., 2, function(i) which(i < dist), simplify = FALSE)
            
            list_increasing_distance[[paste0(as.character(dist), "_", categ)]] = 
              data.frame(n_active = length(which(sapply(within_range, length)>0)),
                         n_tot = length(cells2),
                         sample_id = samp,
                         condition = unique(spe_samp$condition),
                         dist = dist,
                         categ = categ)
          }
        }
        list_type2_to_type[[samp]] = do.call("rbind", list_increasing_distance)
      }
      
    }
    cat("\n")
    saturation_df =  do.call("rbind", list_type2_to_type)
    saturation_df = saturation_df %>% mutate(percent = 100 * n_active / n_tot)
    saturation_df_sum <- summarySE(saturation_df, measurevar="percent", groupvars=c("dist", "categ"))
    
    dir.create(file.path(output_dir, "Marker_based", disease))
    if(length(type2) == 1 & length(type1) == 1){
      png(file.path(output_dir, "Marker_based", disease, 
                    paste0("Saturation_curve_",type2,"within_reach_of_",type1,"_",marker,".png")), width = 2500, height = 1600, res = 300)
      saturation_df_sum = saturation_df_sum  %>% 
        mutate(condition = factor(categ))
      t_test_result <- t.test(percent ~ categ, data = saturation_df_sum, paired = TRUE)
      
      nmin = length(which( (spe.$cell_type_CellSighter %in% type2 | spe.$cell_type_simplified %in% type2) &
                             spe.@assays@data$CellSighter_marker_mat[marker,] == 0))
      npos =  length(which( (spe.$cell_type_CellSighter %in% type2 | spe.$cell_type_simplified %in% type2) &
                              spe.@assays@data$CellSighter_marker_mat[marker,] == 1)) 
      p = saturation_df_sum %>% ggplot(aes(x=dist, y=percent, colour=condition)) + 
        geom_errorbar(aes(ymin=percent-se, ymax=percent+se), width=.1) +
        geom_line() +
        geom_point()  +
        theme_classic() + xlab(paste0("Distance of cells within range")) + 
        ylab(paste0("% of ", type2," cells expressing ", marker, " within range of ", type1)) +
        ggtitle(paste0("Distance from ", type1," depending on ", marker, " expression - ", disease)) + 
        theme(axis.text.x = element_text(angle = 90)) + 
        scale_color_manual(values = c("forestgreen", "red"), name = marker, labels = c("-", "+")) + 
        annotate("text", x = mean(saturation_df_sum$dist) + 700, y = max(saturation_df_sum$percent) + 15, 
                 label = paste("Wilcoxon Paired Test - p =", signif(t_test_result$p.value, digits = 4)), 
                 hjust = 1, vjust = 1, color = "black") + 
        annotate("text", x = min(saturation_df_sum$dist) + 700, y = max(saturation_df_sum$percent) +15, 
                 label = paste0(marker, "- (n) = ", nmin), 
                 hjust = 1, vjust = 1, color = "black") + 
        annotate("text", x = min(saturation_df_sum$dist) + 700, y = max(saturation_df_sum$percent) +10, 
                 label = paste0(marker, "+ (n) = ", npos), 
                 hjust = 1, vjust = 1, color = "black")
      print(p)
      dev.off()
      
    }
  }
}
saturation_df_markers = do.call("rbind", saturation_df_list)
qs::qsave(saturation_df_markers, file.path(output_dir, "Saturation_Curves", "Marker_based", "saturation_df_markers.qs"))

################################################################################
# FOCUS - Saturation curve - % celltype1 within reach of other cells expressing
#  or not  marker X
################################################################################

type1 = "Endothelial"

type2 = unique(spe$cell_type_CellSighter)
type2 = setdiff(type2, type1)
type2 = T_cells

spe. = spe[, !spe$condition %in% c("TON", "LN", "THY")]
saturation_df_list = list()

for(marker in rownames(spe)){ # setdiff(rownames(spe.), markers_to_remove)
  
  list_type2_to_type = list()
  
  for(samp in unique(spe.$sample_id)){
    print(samp)
    d = list_distances[[samp]]
    
    spe_samp = spe.[,spe.$sample_id == samp]
    spe_celltype1 = spe_samp[,spe_samp$cell_type_CellSighter %in% type1]
    spe_celltype2 = spe_samp[,spe_samp$cell_type_CellSighter %in% type2]
    
    
    n_cells = max(ncol(spe_celltype2), ncol(spe_celltype1))
    
    cells1 = spe_celltype1$cell_id
    cells2_pos = spe_celltype2$cell_id[which(spe_celltype2@assays@data$CellSighter_marker_mat[marker,] == 1)]
    cells2_neg = spe_celltype2$cell_id[which(spe_celltype2@assays@data$CellSighter_marker_mat[marker,] == 0)]
    
    if(length(cells1) > 1 & length(cells2_pos) > 1 & length(cells2_neg) > 1){
      
      list_increasing_distance = list()
      
      for(categ in c("pos", "neg")){
        
        if(categ == "pos") cells2 = cells2_pos
        if(categ == "neg") cells2 = cells2_neg
        
        d. = d[cells2, cells1] 
        
        for (dist in seq(0,2000,50)) {
          within_range = apply(d., 1, function(i) which(i < dist), simplify = FALSE)
          
          list_increasing_distance[[paste0(as.character(dist), "_", categ)]] = 
            data.frame(n_active = length(which(sapply(within_range, length)>0)),
                       n_tot_1 = length(cells1),
                       n_tot_2 = length(cells2),
                       sample_id = samp,
                       condition = unique(spe_samp$condition),
                       dist = dist,
                       categ = categ)
        }
      }
      list_type2_to_type[[samp]] = do.call("rbind", list_increasing_distance)
    }
    
  }
  
  saturation_df =  do.call("rbind", list_type2_to_type)
  saturation_df = saturation_df %>% mutate(percent = 100 * n_active / n_tot_2)
  saturation_df_sum <- summarySE(saturation_df, measurevar="percent", groupvars=c("dist", "categ"))
  
  png(file.path(output_dir, "Marker_based", 
                paste0("Saturation_curve_",type1,"_within_reach_of_non",type1,"_",marker,".png")), width = 2500, height = 1600, res = 300)
  saturation_df_sum = saturation_df_sum  %>% 
    mutate(condition = factor(categ))
  t_test_result <- t.test(percent ~ categ, data = saturation_df_sum, paired = TRUE)
  
  p = saturation_df_sum %>% ggplot(aes(x=dist, y=percent, colour=condition)) + 
    geom_errorbar(aes(ymin=percent-se, ymax=percent+se), width=.1) +
    geom_line() +
    geom_point()  +
    theme_classic() + xlab(paste0("Distance of cells within range")) + 
    ylab(paste0("% of ",marker,"+/- non-",type1," within range of ", type1)) +
    ggtitle(paste0("Distance from ", type1," depending on ", marker, " expression")) + 
    theme(axis.text.x = element_text(angle = 90)) + 
    scale_color_manual(values = c("forestgreen", "red"), name = paste0("Non-", type1, " ", marker), labels = c("-", "+")) + 
    annotate("text", x = mean(saturation_df_sum$dist) + 900, y = max(saturation_df_sum$percent) + 15, 
             label = paste("Wilcoxon Paired Test - p =", signif(t_test_result$p.value, digits = 4)), 
             hjust = 1, vjust = 1, color = "black") + 
    annotate("text", x = min(saturation_df_sum$dist) + 700, y = max(saturation_df_sum$percent) +15, 
             label = paste0(marker, "- (n) = ", length(which(spe.$cell_type_CellSighter %in% type2 &
                                                               spe.@assays@data$CellSighter_marker_mat[marker,] == 0))), 
             hjust = 1, vjust = 1, color = "black") + 
    annotate("text", x = min(saturation_df_sum$dist) + 700, y = max(saturation_df_sum$percent) +10 , 
             label = paste0(marker, "+ (n) = ", length(which(spe.$cell_type_CellSighter %in% type2 &
                                                               spe.@assays@data$CellSighter_marker_mat[marker,] == 1))), 
             hjust = 1, vjust = 1, color = "black")
  print(p)
  dev.off()
  
}

################################################################################
# Marker evolution of type1 across distance to type2, focusing on 2 cell types 
################################################################################

type = "T_cytotoxic"
type2 = "T_helper"
condition_1 = "MF"

dist = c(200)

list_type2_to_type = list()
for(samp in unique(spe$sample_id[spe$condition == condition_1])){
  d = list_distances[[samp]]
  cells2 = spe$cell_id[which(spe$sample_id == samp & spe$cell_type_CellSighter == type2)]
  cells = spe$cell_id[which(spe$sample_id == samp & spe$cell_type_CellSighter == type)]
  if(length(cells) > 1 & length(cells2) > 1){
    d = d[cells2, cells] 
    if(type == type2)
      diag(d) = Inf
    
    list_increasing_distance = list()
    
    
    lower_range = apply(d, 2, function(i) which(i <= 100), simplify = FALSE)
    idx_low = names(which(sapply(lower_range, length)>0))
    higher_range = apply(d, 2, function(i) which(i > 500), simplify = FALSE)
    idx_high =names(which(sapply(higher_range, length)>0))
    
    idx_low = match(idx_low, spe$cell_id)
    pos_mark_low = unlist(lapply(rownames(spe), function(marker){sum(spe@assays@data$CellSighter_marker_mat[marker, idx_low ])}))
    idx_high = match(idx_high, spe$cell_id)
    pos_mark_high = unlist(lapply(rownames(spe), function(marker){sum(spe@assays@data$CellSighter_marker_mat[marker, idx_high ])}))
    names(pos_mark_low) = rownames(spe)
    names(pos_mark_high) = rownames(spe)
    
    list_increasing_distance[[as.character(dist)]] = 
      data.frame(
        marker = rownames(spe),
        n_active_low = pos_mark_low,
        n_active_high = pos_mark_high,
        ntot_low = length(idx_low),
        ntot_high = length(idx_high),
        dist = dist,
        condition = unique(spe$condition[spe$sample_id == samp]),
        sample_id = samp
      )
    
  }
  list_type2_to_type[[samp]] = do.call("rbind", list_increasing_distance)
}

saturation_df =  do.call("rbind", list_type2_to_type)
saturation_df = saturation_df %>% mutate(close = 100 *  n_active_low / (ntot_low),
                                         distant = 100 *  n_active_high / (ntot_high))

saturation_df = saturation_df %>% pivot_longer(cols = c("close","distant"))

pdf(file.path(output_dir, "Saturation_Curves", paste0("Expression_in_",type,"with_increasing_distance_from",type2,"_",condition_1, ".pdf")))
for(mark in unique(saturation_df$marker)){
  df = saturation_df %>% filter(marker == mark) %>%
    filter(!is.nan(value), !is.infinite(value))
  
  saturation_df_sum <- summarySE(df, measurevar="percent", groupvars=c("dist", "condition"))
  
  df = df %>% dplyr::group_by(dist) %>% dplyr::summarise(ntot =  sum(n_active))
  
  cortest = cor.test(saturation_df_sum$percent, y = saturation_df_sum$dist)
  cor = cor(saturation_df_sum$percent, y = saturation_df_sum$dist)
  
  png(file.path(output_dir, "Saturation_Curves", paste0(mark, "_expression_in_",type,"with_increasing_distance_from",type2,"_",condition_1, ".png")), width = 2500, height = 1600, res = 300)
   p = saturation_df_sum %>%
    ggplot(aes(x=dist, y=percent, color = condition)) +
    geom_errorbar(aes(ymin=percent-se, ymax=percent+se), width=.1) +
    geom_line() +
    geom_point()  +
    scale_color_manual(values = unique(spe$condition_color[match(unique(saturation_df_sum$condition), spe$condition)])) +
    theme_classic() + xlab(paste0("Distance from closest ", type2)) + ylab(paste0("% of ", type, " expressing ", mark)) +
    ggtitle(paste0("% of ", type, " expressing ", mark)) + theme(axis.text.x = element_text(angle = 90)) +
     annotate("text", x = mean(saturation_df_sum$dist) + 700, y = max(saturation_df_sum$percent) + 0.4 * max(saturation_df_sum$percent),
              label = paste("Pearson's Correlation Test - p =", signif(cortest$p.value, digits = 3), ", R=", round(cor,3)),
              hjust = 1, vjust = 1, color = "black")
   print(p)
   dev.off()

  p = df %>%
    ggplot(aes(x=name, y=value, color = condition)) +
    geom_boxplot(outlier.colour = "white") +
    geom_jitter(width = 0.05, size = 0.25)  +
    scale_color_manual(values = unique(spe$condition_color[match(unique(df$condition), spe$condition)])) +
    theme_classic() + xlab(paste0("Distance from closest ", type2)) + ylab(paste0("% of ", type, " expressing ", mark)) +
    ggtitle(paste0("% of ", type, " expressing ", mark)) + theme(axis.text.x = element_text(angle = 90)) +
    stat_compare_means(aes(label = after_stat(p.signif)),
                       method = "t.test", ref.group = "close", method.args = list(alternative = "greater"))
  print(p)
  
}
dev.off()
