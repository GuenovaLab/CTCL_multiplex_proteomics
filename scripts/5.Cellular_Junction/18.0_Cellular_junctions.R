# Retrieve and analyse cellular junctions
# authors: Pacome Prompsy
# contact: pacome.prompsy@chuv.ch
# Guenova Lab
# CHUV (Centre Hospitalier Universitaire Vaudois), Lausanne, Suisse

cat("Analyzing inter cellular interactions... \n")

# Loading packages --------------------------------------------------------
libraries = c("argparse",
              "ggplot2",
              "plyr",
              "dplyr",
              "tidyr",
              "Seurat",
              "SpatialExperiment",
              "ggbeeswarm")
suppressPackageStartupMessages(invisible(lapply(libraries, require, character.only = TRUE)))
setwd("/mnt/RECHERCHE/GUENOVA_LAB/Project_Multiplex_Phenotyping/Miltenyi/Muliplex_Imaging_Pipeline/")

source("scripts/MultiplexImaging_utils.R")
source("scripts/GlobalVars.R")
library(reticulate)
np <- import("numpy")


# Change when running from R:
args = list(output = "output/")

cat("Output = ", args$output, "\n")

output_dir = file.path(args$output, "cell_interactions", "enriched_interactions")
if (!dir.exists(output_dir))
  dir.create(output_dir)

spe = qs::qread("output/SpatialExperiment.qs")

marker_metadata = read.csv("annotation/marker_metadata.csv")
markers = marker_metadata$Marker[marker_metadata$PassOverallQuality == "True"]
markers = markers[markers!= "Ki67"]

################################################################################
# Load enriched protein interactions for each sample
################################################################################

matrix_list = list()

interaction_mat = as.sparse(spe@metadata$interaction)

interaction_mat[,] = 0

for(marker in markers){
  matrix_list[marker] = interaction_mat
}

markers_no_cxcr4 = setdiff(markers, "CXCR4")

samples = unique(spe$sample_id)
samples = setdiff(samples, spe$sample_id[spe$condition %in% c("LN", "THY", "TON")])

for(sample in samples){
  print(sample)
  dict_interactions <- np$load(file.path(output_dir, "Dictionnaries_interaction", paste0("dict_",sample,".pkl")), allow_pickle = T)
  
  
  for(cell_1 in names(dict_interactions)){
    print(cell_1)
    list_cells = dict_interactions[[cell_1]]
    for(cell_2 in names(list_cells)){
      df = as.data.frame(list_cells[[cell_2]])
      colnames(df) = c("p_value_1", "p_value_2", "mean_border", "mean_background", "FC")
      
      if(unique(spe$batch[spe$sample_id == sample]) == "TMA1") df$markers = markers
      if(unique(spe$batch[spe$sample_id == sample]) == "TMA2") df$markers = markers_no_cxcr4
      
      df$corrected_pvalue_1 = p.adjust(df$p_value_1, method = "bonferroni")
      df$corrected_pvalue_2 = p.adjust(df$p_value_2, method = "bonferroni")
      df$log2FC = log2(df$FC)
      
      # Protein is enriched at cell border of cell A & B if 
      # both adjusted pvalues are lower than 0.01 & FC is greater than 1.5
      df = df %>% mutate(categ = ifelse(corrected_pvalue_1 < 0.01 & corrected_pvalue_2 < 0.01,
                                        ifelse(log2FC > log2(1.5), "signifplus", 
                                               ifelse(log2FC < -log2(1.5), "signifminus", "unsig")), "unsig"))
      signif_markers  = df$markers[which(df$categ == "signifplus")]
      
      if(length(signif_markers) > 0){
        
        if(paste0(sample,"-",cell_1) %in% rownames(interaction_mat) & 
           paste0(sample,"-",cell_2) %in% colnames(interaction_mat)){
          for(mark in signif_markers){
            matrix_list[[mark]][match(paste0(sample,"-",cell_1), rownames(matrix_list[[mark]])),
                                match(paste0(sample,"-",cell_2), colnames(matrix_list[[mark]]))] = 1 
          }
        }
      }

      p = df  %>%
        ggplot(aes(x = log2FC, y = -log10(corrected_avg_pvalue), color = categ, label = markers)) +
        geom_point() +
        scale_color_manual(values = c("darkgreen","red","black")) +
        theme_classic() + xlab("Percent - log2FC") + ylab("-log10(adjusted p-value)") +
        ggtitle(paste0("Interactions ",cell_1,"-",cell_2,"- volcano plot - Late vs Early")) +
        geom_hline(yintercept = 2, lwd = 0.35, lty = 2, col = "grey20") +
        geom_vline(xintercept = -log2(1.5), lwd = 0.35, lty = 2, col = "grey20") +
        geom_vline(xintercept = log2(1.5),  lwd = 0.35, lty = 2, col = "grey20") +
        ggrepel::geom_text_repel()
      print(p)
    }
  }

}
qs::qsave(matrix_list, file.path(output_dir, "matrix_list.qs"))

################################################################################
# Combine enriched protein interactions and save per celltype
################################################################################

matrix_list = qs::qread(file.path(output_dir, "matrix_list.qs"))

# Select only "active" interactions, e.g. interactions with at least 3 proteins enriched for
all_protein = Reduce("+", matrix_list)
all_protein = all_protein[match( spe$cell_id, rownames(all_protein)), match(spe$cell_id, colnames(all_protein))]
all_protein = all_protein + t(all_protein)
interaction_mat = as.sparse(spe@metadata$interaction)
interaction = interaction_mat[match(spe$cell_id, rownames(interaction_mat)), match(spe$cell_id, colnames(interaction_mat))]

interaction_idx = summary(interaction)
interaction_idx$ij = paste0(interaction_idx$i, "_", interaction_idx$j)
all_protein_idx = summary(all_protein)
all_protein_idx$ij = paste0(all_protein_idx$i, "_", all_protein_idx$j)
all_protein_idx = all_protein_idx %>% filter(x  > 0)
interaction_idx = interaction_idx[match(all_protein_idx$ij, interaction_idx$ij),]
interaction_idx$x = TRUE
interaction_filtered = Matrix::sparseMatrix(i=interaction_idx$i,    #rows to fill in
                                            j=interaction_idx$j,    #cols to fill in
                                            x=interaction_idx$x, #values to use (i.e. the upper values of e)
                                            dims=c(nrow(all_protein),ncol(all_protein)),
                                            dimnames = list(rownames(all_protein),
                                                            colnames(all_protein))) 

# For each cell type, retrieve the number of significantly enriched interaction with T helper for each marker.
for(celltype1 in c( "Keratinocyte", "Endothelial")){
  celltype1_dir = file.path(output_dir, celltype1)
  if(!dir.exists(celltype1_dir)) dir.create(celltype1_dir)
  
  
  list_celltype_dfs = list()
  
  for(celltype2 in unique(spe$cell_type_CellSighter)){
    
    print("############################################")
    print(celltype2)
    print("############################################")
    
    list_df = list()
    
    interaction_mat = interaction_filtered
    
    for(sample in unique(spe$sample_id)){
      print(sample)
      spe. = spe[,spe$sample_id == sample]
      cells1 = spe.$cell_id[spe.$cell_type_CellSighter == celltype1]
      cells2 = spe.$cell_id[spe.$cell_type_CellSighter == celltype2]
      
      values_celltype1 = c()
      values_celltype2 = c()
      percent_interaction = c()
      
      interaction_mat. = interaction_mat[match(cells1, rownames(interaction_mat)),
                                         match(cells2, rownames(interaction_mat))]
      
      if(!is.null(ncol(interaction_mat.)) & !is.null(nrow(interaction_mat.))){
        cells1. = rownames(interaction_mat.)[rowSums(interaction_mat.) > 0]
        cells2. = colnames(interaction_mat.)[colSums(interaction_mat.) > 0]
        
        if(length(cells1.) > 1 & length(cells2.) > 1){
          
          interaction_mat.. = interaction_mat.[which(rownames(interaction_mat.) %in% cells1.),
                                               which(colnames(interaction_mat.) %in% cells2.)]
          n_interactions  = sum(apply(interaction_mat.., 1, function(i) length(which(i > 0))))
          
          for(marker in markers){
            mat = matrix_list[[marker]]
            
            mat_sym = mat[which(rownames(mat) %in% cells1.),
                          which(colnames(mat) %in% cells2.)] +
              t(mat[which(colnames(mat) %in% cells2.),
                    which(rownames(mat) %in% cells1.)])
            idx = which(mat_sym == 1, arr.ind = T)
            
            if(!is.null(nrow(idx))){
              values_celltype1 = c(values_celltype1, 100 * nrow(idx) / (length(cells1.)))
              values_celltype2 = c(values_celltype2, 100 * nrow(idx) / (length(cells2.)))
              percent_interaction = c(percent_interaction, 100 * nrow(idx) / n_interactions)
            } else{
              values_celltype1 = c(values_celltype1, 0)
              values_celltype2 = c(values_celltype2, 0)
              percent_interaction = c(percent_interaction, 0)
            }
            
          }
          
          list_df[[sample]] = data.frame(
            sample_id = sample,
            condition = unique(spe.$condition),
            batch = unique(spe.$batch),
            Th_group = unique(spe.$Th_group),
            marker = markers,
            celltype1 = celltype1,
            celltype2 = celltype2,
            n_celltype1 = length(cells1.),
            n_celltype2 = length(cells2.),
            percent_celltype_1_enriched =  values_celltype1,
            percent_celltype_2_enriched =  values_celltype2,
            mean_percent_celltype = (values_celltype1 + values_celltype2) / 2,
            n_interactions = n_interactions,
            percent_interaction = percent_interaction
          )
        }
        
      }
    }
    
    df = as.data.frame(do.call("rbind", list_df))
    
    list_celltype_dfs[[celltype2]] = df
  }
  
  all_interactions = do.call("rbind", list_celltype_dfs)
  
  all_interactions$disease = "Inflammatory"
  all_interactions$disease[all_interactions$condition == "HD"] = "HD"
  all_interactions$disease[all_interactions$condition %in% c("SS", "MF")] = "Lymphoma"
  
  qs::qsave(all_interactions, file.path(celltype1_dir, paste0(celltype1, "_interactions.qs")))
}



################################################################################
# Load log2FC protein interactions for each sample
################################################################################
dir.create(file.path(output_dir, "intensity"))

markers_no_cxcr4 = setdiff(markers, "CXCR4")
samples = unique(spe$sample_id)
samples = setdiff(samples, spe$sample_id[spe$condition %in% c("LN", "THY", "TON")])
intensity_list_sample = list()

for(celltype1 in c("NKT", "T_cytotoxic", "T_regulatory","pDC", "APC", "Macrophages")){
  
  list_per_celltype_2 = list()
  for(celltype2 in unique(spe$cell_type_CellSighter)){
    
    cat("#################################################################\n")
    cat("#################################################################\n")
    cat("                        ",celltype2, "\n")
    cat("#################################################################\n")
    cat("#################################################################\n")
    list_df = list()
    
    for(sample in samples){
      print(sample)
      spe. = spe[, which(spe$sample_id == sample)]
      dict_interactions <- np$load(file.path(output_dir, "Dictionnaries_interaction", paste0("dict_",sample,".pkl")), allow_pickle = T)
      names(dict_interactions) = paste0(sample,"-", names(dict_interactions))
      
      
      interaction_mat = interaction_filtered
      
      cells1 = spe.$cell_id[spe.$cell_type_CellSighter == celltype1]
      cells2 = spe.$cell_id[spe.$cell_type_CellSighter == celltype2]
      
      values_celltype1 = c()
      values_celltype2 = c()
      percent_interaction = c()
      
      interaction_mat. = interaction_mat[match(cells1, rownames(interaction_mat)),
                                         match(cells2, rownames(interaction_mat))]
      
      
      if(!is.null(ncol(interaction_mat.)) & !is.null(nrow(interaction_mat.))){
        cells1. = rownames(interaction_mat.)[rowSums(interaction_mat.) > 0]
        cells2. = colnames(interaction_mat.)[colSums(interaction_mat.) > 0]
        
        
        
        if(length(cells1.) > 1 & length(cells2.) > 1){
          
          cat("Num cells", celltype1,"=",length(cells1.) ,"\n")
          cat("Num cells", celltype2,"=",length(cells2.) ,"\n")
          gc()
          
          interaction_mat.. = interaction_mat.[which(rownames(interaction_mat.) %in% cells1.),
                                               which(colnames(interaction_mat.) %in% cells2.), drop = F]
          n_interactions  = sum(apply(interaction_mat.., 1, function(i) length(which(i > 0))))
          
          cells = intersect(names(dict_interactions), c(cells1., cells2.))
          
          if(!is.null(cells) & length(cells) > 0){
            
            list_interaction_intensity = list()
            for(cell_x in cells){
              
              list_cells = dict_interactions[[cell_x]]
              names(list_cells) = paste0(sample,"-", names(list_cells))
              cells. = intersect(names(list_cells), c(cells1., cells2.))
              
              if(!is.null(cells.) & length(cells.) > 0){
                
                for(cell_y in cells.){
                  df = as.data.frame(list_cells[[cell_y]])
                  colnames(df) = c("p_value_1", "p_value_2", "mean_border", "mean_background", "FC")
                  
                  if(unique(spe$batch[spe$sample_id == sample]) == "TMA1") df$markers = markers
                  if(unique(spe$batch[spe$sample_id == sample]) == "TMA2") df$markers = markers_no_cxcr4
                  
                  df$log2FC = log2(df$FC)
                  df$effect = abs(df$mean_border - df$mean_background)
                  df$sample_id = sample
                  df$condition = unique(spe.$condition)
                  df$batch = unique(spe.$batch)
                  df$celltype1 = celltype1
                  df$celltype2 = celltype2
                  df$cell1 = cell_x
                  df$cell2 = cell_y
                  
                  list_interaction_intensity[[paste0(cell_x, "_", cell_y)]] = df
                }
              }
            }
            list_df[[sample]] = do.call("rbind", list_interaction_intensity)
          }
          
        }
        
      }
      
    }
    
    list_per_celltype_2[[celltype2]] = as.data.frame(do.call("rbind", list_df))
    
  }  
  
  df = do.call("rbind", list_per_celltype_2)
  qs::qsave(df, file.path(output_dir, "intensity", paste0(celltype1, ".qs")))
  
}


################################################################################
# Boxplots of log2FC intensity of interactions - Malignant vs Benign
################################################################################

library(ggpubr)
list_DA_celltype1 = list()
for(celltype1 in c("T_helper", "NKT", "T_cytotoxic", "T_regulatory","pDC", "APC", "Macrophages")){
  all_interactions = qs::qread(file.path(output_dir, "intensity", paste0(celltype1, ".qs")))
  
  list_DA_celltype = list()
  for(type2 in c("T_helper", "NKT", "T_cytotoxic", "T_regulatory","pDC", "APC", "Macrophages")){
    print(paste0(celltype1, " <--> ", type2))
    
    df = all_interactions %>% dplyr::filter(.data[["celltype2"]] == type2)
    
    df = df[df$condition != "HD", ]
    df$disease = spe$disease[match(df$sample_id, spe$sample_id)]
    df$disease = as.factor(df$disease)
    df = df[which(!is.na(df$batch)),]
    
    dir.create(file.path(output_dir, "intensity", celltype1))
    
    pdf(file.path(output_dir, "intensity", celltype1, paste0(celltype1,"-",type2, "_by_disease.pdf")), width = 3.5, height = 4)
    list_DA = list()
    for(mark in setdiff(markers, "DAPI")){
      
      df. = df %>% dplyr::filter(markers == mark)
      df. = df. %>% dplyr::filter(log2FC > 0)
      
      if(nrow(df.[df.$disease == "Malign",])> 2 & nrow(df.[df.$disease != "Malign",]) > 2){
        
        df. = df. %>% dplyr::filter(!is.infinite(log2FC))
        df.$ID = paste0(df.$cell1, "_", df.$cell2)
        set.seed(47)
        cells = c()
        for(i in unique(df.$sample_id)){
          cells = c(cells, sample(df.$ID[df.$sample_id == i], min(100, nrow(df.[df.$sample_id == i,]))))
        }
        df_equibrilated = df.[match(cells, df.$ID), ]
        
        
        p = df_equibrilated %>% ggplot(aes(x = disease, fill = disease, y = log2FC)) +
          geom_violin(draw_quantiles = 0.5) + geom_jitter(size = 0.05, position = position_quasirandom())+
          scale_fill_manual(values = setNames(unique(spe$disease_color), unique(spe$disease))) +
          stat_compare_means(ref.group = "Benign", method = "t.test") + 
          theme_minimal()
        p = p + xlab("") + ylab(paste0("log2FC(Border/Background) - ", mark, ")")) +
          guides(fill="none", color = "none") +
          ggtitle(paste0(gsub("_"," ",celltype1), " <--> ",gsub("_"," ",type2))) 
        print(p)
        
        df_grouped = df. %>% group_by(sample_id, disease, batch) %>% 
          dplyr::summarise(log2FC = mean(log2FC))
        
        p = grouped_dotplot(df_grouped, y = "log2FC", categ1 = "disease", categ2 = NULL,
                            ref.group = "Benign", add_violin = F,
                            color_categ1 = setNames(unique(spe$disease_color), unique(spe$disease)),
                            size = 3, beeswarm = FALSE)
        p = p + xlab("") + ylab(paste0("log2FC(Border/Background) - ", mark, ")")) +
          guides(fill="none", color = "none") +
          ggtitle(paste0(gsub("_"," ",celltype1), " <--> ",gsub("_"," ",type2))) 
        p = p + stat_summary(fun.data = "mean_cl_boot",  linewidth = 1, size = 2, pch = 18) 
        print(p)  
        
        list_DA[[mark]] = data.frame(
          celltype1 = celltype1,
          celltype2 = type2,
          marker = mark,
          condition = "Malign_vs_Benign",
          diff_log2FC = mean(df.$log2FC[df.$disease == "Malign"]) - mean(df.$log2FC[df.$disease != "Malign"]),
          diff_effect = mean(df.$effect[df.$disease == "Malign"]) - mean(df.$effect[df.$disease != "Malign"]),
          p.value_single_cell = tryCatch({t.test(df.$log2FC[df.$disease == "Malign"], df.$log2FC[df.$disease != "Malign"])$p.value}, error = function(e){return(1)}),
          p.value_by_sample = tryCatch({t.test(df_grouped$log2FC[df_grouped$disease == "Malign"], df_grouped$log2FC[df_grouped$disease != "Malign"])$p.value}, error = function(e){return(1)})
        )
      }
      
    }
    dev.off()
    list_DA_celltype[[type2]] = do.call("rbind", list_DA)
  }
  
  DA = do.call("rbind", list_DA_celltype)
  list_DA_celltype1[[celltype1]] = DA
} 

WriteXLS::WriteXLS(list_DA_celltype1, ExcelFileName = file.path(output_dir, "intensity", "Differential_Enriched_Interactions_Malign_vs_Benign.xlsx"),
                   SheetNames = names(list_DA_celltype1), AdjWidth = T, AutoFilter = T, BoldHeaderRow = T)



################################################################################
# Boxplots of differentially enriched interactions - Malignant vs Benign
################################################################################

library(ggpubr)
list_DA_celltype1 = list()
for(celltype1 in c("T_cytotoxic", "T_regulatory", "Macrophages", "NKT", "Monocytic_Lineage", "B_cell", "pDC", "APC", "Keratinocyte", "Endothelial")){
  
  celltype1_dir = file.path(output_dir, celltype1)
  all_interactions = qs::qread(file.path(celltype1_dir, paste0(celltype1, "_interactions.qs")))
  
  
  interaction_threshold = 5
  list_DA_celltype = list()
  for(type2 in unique(spe$cell_type_CellSighter)){
    print(paste0(celltype1, " <--> ", type2))
    
    df = all_interactions %>% dplyr::filter(.data[["celltype2"]] == type2)
    
    df = df[df$condition != "HD", ]
    df$disease = "Benign"
    df$disease[df$condition %in% c("SS", "MF")] = "Malign"
    df$disease = as.factor(df$disease)
    df$batch = spe$batch[match(df$sample_id, spe$sample_id)]
    df = df[which(!is.na(df$batch)),]
    
    pdf(file.path(celltype1_dir, paste0(celltype1,"-",type2, "_by_disease.pdf")), width = 3.5, height = 4)
    list_DA = list()
    for(mark in markers){
      
      df. = df %>% dplyr::filter(n_interactions >= interaction_threshold) %>%
        dplyr::filter(marker == mark)
      
      if(nrow(df.[df.$disease == "Malign",])> 2 & nrow(df.[df.$disease != "Malign",]) > 2){
        
        p = grouped_dotplot(df., y = "percent_interaction", categ1 = "disease", categ2 = NULL,
                            ref.group = "Benign", add_violin = F,
                            color_categ1 = setNames(unique(spe$disease_color), unique(spe$disease)))
        p = p + xlab("") + ylab(paste0("Interactions enriched in ", mark, " (%)")) +
          guides(fill="none", color = "none") +
          ggtitle(paste0(gsub("_"," ",celltype1), " <--> ",gsub("_"," ",type2)))
        print(p)  
        
        list_DA[[mark]] = data.frame(
          celltype1 = celltype1,
          celltype2 = type2,
          marker = mark,
          condition = "Malign_vs_Benign",
          percent_interaction_Lymphoma = mean(df.$percent_interaction[df.$disease == "Malign"]),
          percent_interaction_Inflammatory = mean(df.$percent_interaction[df.$disease != "Malign"]),
          log2FC = log2(mean(df.$percent_interaction[df.$disease == "Malign"]) / mean(df.$percent_interaction[df.$disease != "Malign"])),
          effect = abs(mean(df.$percent_interaction[df.$disease == "Malign"]) - mean(df.$percent_interaction[df.$disease != "Malign"])),
          p.value = t.test(df.$percent_interaction[df.$disease == "Malign"], df.$percent_interaction[df.$disease != "Malign"])$p.value
        )
      }
      
    }
    dev.off()
    
    list_DA_celltype[[type2]] = do.call("rbind", list_DA)
  }
  
  DA = do.call("rbind", list_DA_celltype)
  DA$log2FC[is.infinite(DA$log2FC)] = 0
  DA$p.value[is.nan(DA$p.value)] = 1
  DA$Significant = FALSE
  DA$Significant[abs(DA$log2FC) > log2(1.5) & DA$p.value < 0.05] = TRUE
  list_DA_celltype1[[celltype1]] = DA
} 

WriteXLS::WriteXLS(list_DA_celltype1, ExcelFileName = file.path(output_dir, "Differential_Enriched_Interactions_Malign_vs_Benign.xlsx"),
                   SheetNames = names(list_DA_celltype1), AdjWidth = T, AutoFilter = T, BoldHeaderRow = T)



################################################################################
# Boxplots of differentially enriched interactions - Malignant vs Benign - simplified
################################################################################

library(ggpubr)
list_DA_celltype1 = list()
for(celltype1 in c("Myeloid", "Lymphocyte")){
  
  celltypes = unique(spe$cell_type_CellSighter[spe$cell_type_CellSighter %in% spe$cell_type_CellSighter[spe$cell_type_simplified == celltype1]])
  
  list_all_interactions = list()
  for(type in celltypes){
    celltype1_dir = file.path(output_dir, type)
    if(file.exists(file.path(celltype1_dir, paste0(type, "_interactions.qs"))))
      list_all_interactions[[type]] = qs::qread(file.path(celltype1_dir, paste0(type, "_interactions.qs")))
  }
  all_interactions = do.call("rbind", list_all_interactions)
  all_interactions$batch[all_interactions$batch == "batch_1"] = "TMA1"
  all_interactions$batch[all_interactions$batch == "batch_2"] = "TMA2"
  all_interactions = all_interactions %>% group_by(celltype2, sample_id, disease, condition, marker, batch) %>%
    summarise(percent_interaction = mean(percent_interaction), n_interactions = sum(n_interactions))
  all_interactions$celltype1 = celltype1
  
  
  interaction_threshold = 5
  list_DA_celltype = list()
  for(type2 in unique(spe$cell_type_CellSighter)){
    print(paste0(celltype1, " <--> ", type2))
    
    df = all_interactions %>% dplyr::filter(.data[["celltype2"]] == type2)
    
    df = df[df$condition != "HD", ]
    df$disease = "Benign"
    df$disease[df$condition %in% c("SS", "MF")] = "Malign"
    df$disease = as.factor(df$disease)
    df$batch = spe$batch[match(df$sample_id, spe$sample_id)]
    df = df[which(!is.na(df$batch)),]
    
    dir.create(file.path(output_dir, celltype1))
    pdf(file.path(output_dir, celltype1, paste0(celltype1,"-",type2, "_by_disease.pdf")), width = 3.5, height = 4)
    list_DA = list()
    for(mark in markers){
      
      df. = df %>% dplyr::filter(n_interactions >= interaction_threshold) %>%
        dplyr::filter(marker == mark)
      
      if(nrow(df.[df.$disease == "Malign",])> 2 & nrow(df.[df.$disease != "Malign",]) > 2){
        
        p = grouped_dotplot(df., y = "percent_interaction", categ1 = "disease", categ2 = NULL,
                            ref.group = "Benign", add_violin = F,
                            color_categ1 = setNames(unique(spe$disease_color), unique(spe$disease)))
        p = p + xlab("") + ylab(paste0("Interactions enriched in ", mark, " (%)")) +
          guides(fill="none", color = "none") +
          ggtitle(paste0(gsub("_"," ",celltype1), " <--> ",gsub("_"," ",type2)))
        print(p)  
        
        list_DA[[mark]] = data.frame(
          celltype1 = celltype1,
          celltype2 = type2,
          marker = mark,
          condition = "Malign_vs_Benign",
          percent_interaction_Lymphoma = mean(df.$percent_interaction[df.$disease == "Malign"]),
          percent_interaction_Inflammatory = mean(df.$percent_interaction[df.$disease != "Malign"]),
          log2FC = log2(mean(df.$percent_interaction[df.$disease == "Malign"]) / mean(df.$percent_interaction[df.$disease != "Malign"])),
          effect = abs(mean(df.$percent_interaction[df.$disease == "Malign"]) - mean(df.$percent_interaction[df.$disease != "Malign"])),
          p.value = t.test(df.$percent_interaction[df.$disease == "Malign"], df.$percent_interaction[df.$disease != "Malign"])$p.value
        )
      }
      
    }
    dev.off()
    
    list_DA_celltype[[type2]] = do.call("rbind", list_DA)
  }
  
  DA = do.call("rbind", list_DA_celltype)
  DA$log2FC[is.infinite(DA$log2FC)] = 0
  DA$p.value[is.nan(DA$p.value)] = 1
  DA$Significant = FALSE
  DA$Significant[abs(DA$log2FC) > log2(1.5) & DA$p.value < 0.05] = TRUE
  list_DA_celltype1[[celltype1]] = DA
} 

WriteXLS::WriteXLS(list_DA_celltype1, ExcelFileName = file.path(output_dir, "Differential_Enriched_Interactions_Malign_vs_Benign_simplified.xlsx"),
                   SheetNames = names(list_DA_celltype1), AdjWidth = T, AutoFilter = T, BoldHeaderRow = T)


################################################################################
# Interaction plots  - celltype synapse specificity stats
################################################################################


# Find protein enriched for specific celltypes
# Plot
interaction_threshold = 5
library(ggpubr)
library(igraph)
dis_top_list = list()


# Calculate celltype specificity stats 
for(dis in unique(spe$disease)[2:3]){
  for(celltype1 in  c("T_helper")){
    print(celltype1)
    celltype1_dir = file.path(output_dir, celltype1)
    all_interactions = qs::qread(file.path(celltype1_dir, paste0(celltype1, "_interactions.qs")))
    all_interactions = all_interactions %>% mutate(n_enriched = percent_interaction * n_interactions / 100)
    all_interactions$disease = spe$disease[match(all_interactions$sample_id, spe$sample_id)]
    all_interactions$score = changeRange(all_interactions$percent_interaction)
    all_interactions = all_interactions %>% filter(disease == dis)
    
    list_stats = list()
    
    for(type2 in unique(spe$cell_type_CellSighter)){
      
      print(type2)
      list_tmp = list()
      dir.create(file.path(output_dir, celltype1))
      
      pdf(file.path(output_dir, celltype1, paste0(celltype1,"_",type2, "_permutation_",dis,".pdf")), width = 5, height = 5.75)
      
      for(mark in setdiff(rownames(spe), markers_to_remove)){
        interactions = all_interactions 
        interactions = interactions %>% mutate(n_enriched = percent_interaction * n_interactions / 100)
        interactions = interactions %>% filter(n_interactions  >= interaction_threshold)
        interactions = interactions %>% filter(marker == mark)
        interactions$permutation = "Real"
        
        # Other celltypes 
        permuted = interactions %>% filter(celltype2 != type2)
        permuted$permutation = "Permuted"
        
        
        df = rbind(permuted, interactions[interactions$celltype2 == type2,])
        df$batch = spe$batch[match(df$sample_id, spe$sample_id)]
        if(nrow(interactions[interactions$celltype2 == type2,]) > 2){
          
          stats = df %>% rstatix::t_test(percent_interaction~permutation, ref.group = "Permuted", alternative = "less")
          stats$celltype1 = celltype1
          stats$celltype2 = type2
          stats$marker = mark
          list_tmp[[mark]] = stats
          if(stats$p < 0.05){
            
            p = grouped_dotplot(df, y = "percent_interaction", categ1 = "permutation", categ2 = NULL,
                                ref.group = "Permuted", add_violin = F, beeswarm = F,
                                color_categ1 = c("Real" = "darkgreen", "Permuted" ="grey"))
            p = p + xlab("") + ylab(paste0("Interactions enriched in ", mark, " (%)")) +
              guides(fill="none", color = "none") +
              ggtitle(paste0(gsub("_"," ",celltype1), " <--> ",gsub("_"," ",type2)))
            p = p + stat_summary(fun.data = "mean_cl_boot", color = "black", linewidth = 1, size = 2, pch = 18) 
            print(p)  
          }
          run_donut = F
          if(run_donut){
            
            interactions = all_interactions %>% filter(celltype1 != celltype2)
            interactions$type = interactions$celltype2
            interactions$disease = spe$disease[match(interactions$sample_id, spe$sample_id)]
            interactions$n_interactions = NULL
            interactions$n_interactions = interactions$n_enriched
            interactions$total_cells = interactions$n_celltype2
            
            .interaction_doughnut_plot(spe, interactions,
                                       celltype = celltype1,
                                       group_by = "cell_type_CellSighter",
                                       plot_by = "disease",
                                       colors = setNames(unique(spe$cell_type_CellSighter_color), unique(spe$cell_type_CellSighter)),
                                       weighted =  "absolute"
            )
            
            
          }
          
        }
      }
      dev.off()
      
      list_stats[[type2]] = do.call("rbind", list_tmp)
    }
    stats = do.call("rbind",list_stats)  
    WriteXLS::WriteXLS(stats , ExcelFileName = file.path(output_dir, celltype1, paste0("stats_permutations_",dis,".xlsx")),
                       AdjWidth = T, AutoFilter = T, BoldHeaderRow = T)
    
  }
}

################################################################################
# Interaction plots - celltype synapse network 
################################################################################

for(dis in unique(spe$disease)[2:3]){
  
  top_list = list()
  for(celltype1 in c("T_helper")){
    print(celltype1)
    celltype1_dir = file.path(output_dir, celltype1)
    all_interactions = qs::qread(file.path(celltype1_dir, paste0(celltype1, "_interactions.qs")))
    all_interactions = all_interactions %>% mutate(n_enriched = percent_interaction * n_interactions / 100)
    all_interactions = all_interactions %>% filter(celltype2 != "Neutrophils") # Not in both diseases
    all_interactions$disease = spe$disease[match(all_interactions$sample_id, spe$sample_id)]
    all_interactions$score = changeRange(all_interactions$percent_interaction)
    other = ifelse(dis == "Benign", "Malign", "Benign" )
    all_interactions_other = all_interactions %>% filter(disease == other)
    all_interactions = all_interactions %>% filter(disease == dis)
    
    meta = as.data.frame(colData(spe[, spe$sample_id %in% unique(all_interactions$sample_id)]))
    meta. = meta %>% group_by(cell_type_CellSighter, sample_id) %>%
      dplyr::summarise(n =  dplyr::n()) %>% group_by(cell_type_CellSighter) %>%
      dplyr::summarise(median_ct =  median(n))
    
    stats = readxl::read_excel(file.path(output_dir, celltype1, paste0("stats_permutations_",dis,".xlsx")))
    stats = stats %>% filter(celltype2 != "Neutrophils") # Not in both diseases
    stats = stats %>% filter(!marker %in% c("CD3", golden_markers)) # Not in both diseases
    stats = stats %>% group_by(celltype1, celltype2) %>% slice_min(p, n = 1, with_ties = F)
    
    top = all_interactions %>% 
      group_by(marker, celltype1, celltype2) %>%
      dplyr::summarise(percent_interaction = mean(percent_interaction),
                       score = mean(score),
                       n_celltype1 = mean(n_celltype1)) 
    top = top[match(paste0(stats$celltype1, stats$celltype2, stats$marker),
                    paste0(top$celltype1, top$celltype2, top$marker)
    ),]
    
    top$specificity_p = stats$p
    
    # Compare to other disease
    other_tab = all_interactions_other %>% 
      group_by(marker, celltype1, celltype2) %>%
      dplyr::summarise(percent_interaction = mean(percent_interaction),
                       score = mean(score),
                       n_celltype1 = mean(n_celltype1)) 
    other_tab = other_tab[match(paste0(stats$celltype1, stats$celltype2, stats$marker),
                                paste0(other_tab$celltype1, other_tab$celltype2, other_tab$marker)
    ),]
    top$effect_with_other = top$percent_interaction - other_tab$percent_interaction
    top$log2FC_with_other = log2(top$percent_interaction / other_tab$percent_interaction)
    top$log2FC_with_other[which(top$log2FC_with_other == Inf)] = max(top$log2FC_with_other[which(top$log2FC_with_other != Inf)]) 
    top$log2FC_with_other[which(top$log2FC_with_other == -Inf)] = min(top$log2FC_with_other[which(top$log2FC_with_other != -Inf)]) 
    
    # Order celltypes
    top = top[match(intersect(celltype_levels, top$celltype2), top$celltype2), ]
    
    top$celltype2[which(top$celltype2 == celltype1)] = paste0(celltype1, " ")
    
    top_list[[celltype1]] = top
    
    
    pdf(file.path(output_dir, celltype1, paste0(celltype1, "_interaction_graph_by_celltype_specific_",dis,".pdf")))
    par(bg = "white", mar = c(1.1, 1.1, 1.1, 1.1), xpd=TRUE)
    
    g = graph_from_data_frame(top[,c(2,3,4,1)], directed = F, vertices = c(celltype1, top$celltype2))
    vertex.color = c(spe$cell_type_CellSighter_color[match(gsub(" ", "", V(g)$name), spe$cell_type_CellSighter)],
                     setNames(unique(spe$cell_type_CellSighter_color[spe$cell_type_CellSighter == celltype1]),
                              paste0(celltype1, " "))
    )
    layout = circle_layout(length(g))
    
    pal = colorRampPalette(c("blue", "grey", "red"))(1000)
    print(plot(g,
               layout = layout,
               vertex.size = sqrt(meta.$median_ct[match(gsub(" ", "", V(g)$name), meta.$cell_type_CellSighter)]) * 1.25, 
               vertex.color = vertex.color, 
               vertex.label.color = "black",
               vertex.label = "",
               edge.width = changeRange(top$score[match(setdiff(V(g)$name,celltype1), top$celltype2)], 1,10) * 2,
               # edge.color = pal[round(changeRange(top$effect_with_other[match(setdiff(V(g)$name,celltype1), top$celltype2)], 1,1000))],
               edge.label = top$marker[match(setdiff(V(g)$name,celltype1), top$celltype2)],
               edge.label.cex = 1.5,
               edge.label.color = "black",
               edge.label.family = "sans",
               edge.curved = c(0.05, 0.1, 0.1, 0.1, 0.1,
                               -0.1, -0.1, -0.1, -0.1,
                               0.05, 0.1, 0.1, 0.1,
                               -0.1, -0.1, -0.1)
    ))
    
    dev.off()
    
    
  }
  dis_top_list[[dis]] = top_list
  
}


################################################################################
# Interaction plots - Direct contact + synapses image overlay 
################################################################################

celltype1 = "T_helper" 
celltype2 =  "T_cytotoxic"
mark =  "CLA"
matrix_list = qs::qread(file.path(output_dir, "matrix_list.qs"))

for(samp in c("ROI-06", "ROI-61", "ROI-22", "ROI-17", "ROI-05", "ROI-63")){
  dir.create(file.path(output_dir, "network"))
  
  spe. = spe[,spe$sample_id == samp]
  spe. = spe.[,which(spe.$centroid.0_nuclear !=0 & spe.$centroid.1_nuclear !=0)]
  mat = spe.@metadata$interaction
  mat = mat[match(spe.$cell_id, rownames(mat)), match(spe.$cell_id, rownames(mat))]
  
  g = igraph::graph_from_adjacency_matrix(mat)
  g = simplify(g)
  g = as.undirected(g)
  
  layout_matrix <- cbind(spe.$centroid.0_nuclear, spe.$centroid.1_nuclear)
  
  mat = matrix_list[[mark]]
  mat = mat[match(spe.$cell_id, rownames(mat)), match(spe.$cell_id, rownames(mat))]
  mat = mat[match(spe.$cell_id[spe.$cell_type_CellSighter == celltype1], rownames(mat)),
            match(spe.$cell_id[spe.$cell_type_CellSighter == celltype2],  colnames(mat))]
  
  # Celltype1 
  cells_1 = rownames(mat)[which(rowSums(mat) > 0)]
  cells_2 = colnames(mat)[which(colSums(mat) > 0)]
  
  colData(spe.)["enriched"] = "white"
  colData(spe.)[match(cells_1, spe.$cell_id), "enriched"] = "darkgreen"
  colData(spe.)[match(cells_2, spe.$cell_id), "enriched"] = "darkred"
  
  png(file.path(output_dir, "network", paste0(samp, "_interaction_graph_",
                                              celltype1,"_",celltype2, "_",mark,".png")),
      width = 2000, height = 2000, res = 300)
  
  par(bg = "transparent", mar = c(1.1, 1.1, 1.1, 1.1), xpd=TRUE)
  print(plot(g,
             layout = layout_matrix,
             vertex.size = ifelse(spe.$enriched == "white", 1, 2), 
             vertex.color = spe.$enriched, 
             vertex.label = "", 
             edge.color = "black"
  ))
  dev.off()
}

################################################################################
# Interaction plots - HEATMAP
################################################################################

library(corrplot)
library(ComplexHeatmap)
for(dis in unique(spe$disease)[2:3]){
  
  top_list = list()
  for(celltype1 in c("T_helper")){
    print(celltype1)
    celltype1_dir = file.path(output_dir, celltype1)
    all_interactions = qs::qread(file.path(celltype1_dir, paste0(celltype1, "_interactions.qs")))
    all_interactions = all_interactions %>% mutate(n_enriched = percent_interaction * n_interactions / 100)
    all_interactions = all_interactions %>% filter(celltype2 != "Neutrophils") # Not in both diseases
    all_interactions$disease = spe$disease[match(all_interactions$sample_id, spe$sample_id)]
    all_interactions$score = changeRange(all_interactions$percent_interaction)
    other = ifelse(dis == "Benign", "Malign", "Benign" )
    all_interactions_other = all_interactions %>% filter(disease == other)
    all_interactions = all_interactions %>% filter(disease == dis)
    
    meta = as.data.frame(colData(spe[, spe$sample_id %in% unique(all_interactions$sample_id)]))
    meta. = meta %>% group_by(cell_type_CellSighter, sample_id) %>%
      dplyr::summarise(n =  dplyr::n()) %>% group_by(cell_type_CellSighter) %>%
      dplyr::summarise(median_ct =  median(n))
    
    stats = readxl::read_excel(file.path(output_dir, celltype1, paste0("stats_permutations_",dis,".xlsx")))
    stats = stats %>% filter(celltype2 != "Neutrophils") # Not in both diseases
    stats = stats %>% filter(!marker %in% c("CD3", golden_markers)) # Not in both diseases
    stats = stats %>% group_by(celltype1, celltype2) %>% slice_min(p, n = 5, with_ties = F)
    
    top = all_interactions %>% 
      group_by(marker, celltype1, celltype2) %>%
      dplyr::summarise(percent_interaction = mean(percent_interaction),
                       score = mean(score),
                       n_celltype1 = mean(n_celltype1)) 
    top = top[match(paste0(stats$celltype1, stats$celltype2, stats$marker),
                    paste0(top$celltype1, top$celltype2, top$marker)
    ),]
    top$specificity_p = stats$p
    
    top = top %>% group_by(celltype1, celltype2) %>%mutate(rank = rank(-specificity_p))
    
    # Compare to other disease
    other_tab = all_interactions_other %>% 
      group_by(marker, celltype1, celltype2) %>%
      dplyr::summarise(percent_interaction = mean(percent_interaction),
                       score = mean(score),
                       n_celltype1 = mean(n_celltype1)) 
    other_tab = other_tab[match(paste0(stats$celltype1, stats$celltype2, stats$marker),
                                paste0(other_tab$celltype1, other_tab$celltype2, other_tab$marker)
    ),]
    top$effect_with_other = top$percent_interaction - other_tab$percent_interaction
    top$log2FC_with_other = log2(top$percent_interaction / other_tab$percent_interaction)
    top$log2FC_with_other[which(top$log2FC_with_other == Inf)] = max(top$log2FC_with_other[which(top$log2FC_with_other != Inf)]) 
    top$log2FC_with_other[which(top$log2FC_with_other == -Inf)] = min(top$log2FC_with_other[which(top$log2FC_with_other != -Inf)])
    
    
    pdf(file.path(output_dir, celltype1, paste0("Heatmap_top_5_synapses_", dis, " _vs_", other, ".pdf")),
        width = 11,
        height = 6)
    p = ggplot(top, aes(x = marker, y = celltype2)) +
      geom_point(aes(size = rank, fill = log2FC_with_other), shape = 21) + 
      scale_size_continuous(name = "Celltype specificity rank", range = c(2,10)) +  # Adjust range for circle size
      scale_fill_gradient2(name = paste0("Log2FC ", dis, " vs ", other), low = "blue",  high = "red") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for better readability
      labs(title = paste0("Synapses with ", gsub("_", " ", celltype1))) +
      xlab("") + ylab("")
    print(p)
    p =  ggplot(top, aes(x = marker, y = celltype2)) +
      geom_point(aes(size = rank, fill = percent_interaction), shape = 21) + 
      scale_size_continuous(name = "Celltype specificity rank", range = c(2,10)) +  # Adjust range for circle size
      scale_fill_gradient(name = "Enriched Synapes (% of interactions)", low = "royalblue4", high = "gold") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for better readability
      labs(title = paste0("Synapses with ", gsub("_", " ", celltype1))) +
      xlab("") + ylab("")
    print(p)
    dev.off()
    
    
  }
}
################################################################################
# Boxplots of differentially enriched interactions - Responder vs Non Responder
################################################################################
colors_progression_regression = setNames(c("#7699D4", "#B3A632FA", "#B54A48", "#74C276"),
                                         c("Healthy", "Benign", "Stable/Progression", "Regression"))
library(ggpubr)
list_DA_celltype1 = list()
output_dir_CTCL = "output/CTCL_response/CAILS/enriched_interactions"
dir.create(output_dir_CTCL)

for(celltype1 in c("T_cytotoxic", "T_regulatory", "Macrophages", "NKT", "Monocytic_Lineage", "B_cell", "pDC", "APC", "Keratinocyte", "Endothelial")){
  dir.create(file.path(output_dir_CTCL, celltype1))
  
  celltype1_dir = file.path(output_dir, celltype1)
  all_interactions = qs::qread(file.path(celltype1_dir, paste0(celltype1, "_interactions.qs")))
  
  interaction_threshold = 5
  list_DA_celltype = list()
  
  for(type2 in unique(spe$cell_type_CellSighter)){
    print(paste0(celltype1, " <--> ", type2))
    
    df = all_interactions %>% dplyr::filter(.data[["celltype2"]] == type2)
    df$CTCL.CAILS = spe$CTCL.CAILS[match(df$sample_id, spe$sample_id)]
    df = df[which(df$CTCL.CAILS %in% c("Stable/Progression", "Regression")),]
    df$batch = spe$batch[match(df$sample_id, spe$sample_id)]
    
    pdf(file.path(output_dir_CTCL, celltype1, paste0(celltype1,"-",type2, "_by_disease.pdf")), width = 5, height = 6)
    list_DA = list()
    for(mark in markers){
      
      df. = df %>% dplyr::filter(n_interactions >= interaction_threshold) %>%
        dplyr::filter(marker == mark)
      
      if(nrow(df.[df.$CTCL.CAILS == "Stable/Progression",])> 1 & nrow(df.[df.$CTCL.CAILS == "Regression",]) > 1){
        
        p = grouped_dotplot(df., y = "percent_interaction", categ1 = "CTCL.CAILS", categ2 = NULL,
                            ref.group = "Stable/Progression", add_violin = F,
                            color_categ1 = colors_progression_regression)
        p = p + xlab("") + ylab(paste0("Interactions enriched in ", mark, " (%)")) +
          guides(fill="none", color = "none") +
          ggtitle(paste0(gsub("_"," ",celltype1), " <--> ",gsub("_"," ",type2)))
        print(p)  
        
        list_DA[[mark]] = data.frame(
          celltype1 = celltype1,
          celltype2 = type2,
          marker = mark,
          condition = "Regression_vs_Progression",
          percent_interaction_Lymphoma = mean(df.$percent_interaction[df.$CTCL.CAILS == "Stable/Progression"]),
          percent_interaction_Inflammatory = mean(df.$percent_interaction[df.$CTCL.CAILS != "Stable/Progression"]),
          log2FC = log2(mean(df.$percent_interaction[df.$CTCL.CAILS == "Stable/Progression"]) / mean(df.$percent_interaction[df.$CTCL.CAILS != "Stable/Progression"])),
          effect = abs(mean(df.$percent_interaction[df.$CTCL.CAILS == "Stable/Progression"]) - mean(df.$percent_interaction[df.$CTCL.CAILS != "Stable/Progression"])),
          p.value = t.test(df.$percent_interaction[df.$CTCL.CAILS == "Stable/Progression"], df.$percent_interaction[df.$CTCL.CAILS != "Stable/Progression"])$p.value
        )
      }
      
    }
    dev.off()
    
    list_DA_celltype[[type2]] = do.call("rbind", list_DA)
  }
  
  DA = do.call("rbind", list_DA_celltype)
  DA$log2FC[is.infinite(DA$log2FC)] = 0
  DA$p.value[is.nan(DA$p.value)] = 1
  DA$Significant = FALSE
  DA$Significant[abs(DA$log2FC) > log2(1.5) & DA$p.value < 0.05] = TRUE
  list_DA_celltype1[[celltype1]] = DA
} 

WriteXLS::WriteXLS(list_DA_celltype1, ExcelFileName = file.path(output_dir_CTCL, "Differential_Enriched_Interactions_Progression_vs_Regression.xlsx"),
                   SheetNames = names(list_DA_celltype1), AdjWidth = T, AutoFilter = T, BoldHeaderRow = T)


################################################################################
# Heatmaps of protein interactions
################################################################################

list_mat = list()
for(dis in c("Lymphoma", "Inflammatory")){
  all_df = all_interactions
  all_df = all_df[all_df$disease != "HD",]
  all_df$stage = spe$stage[match(all_df$sample_id, spe$sample_id)]
  all_df$CTCL.CAILS = spe$CTCL.CAILS[match(all_df$sample_id, spe$sample_id)]
  
  all_df = all_df %>% filter(disease %in% c(dis))
  
  all_df  = all_df %>% dplyr::filter(n_interactions >= 20)
  
  # Mean by sample
  all_df = all_df %>% dplyr::group_by(marker, sample_id, celltype2) %>% 
    dplyr::summarise(percent_interaction = mean(percent_interaction))
  
  # Mean by celltype
  all_df = all_df %>% dplyr::group_by(marker, celltype2) %>% 
    dplyr::summarise(percent_interaction = mean(percent_interaction))
  
  
  all_df = all_df %>% pivot_wider(names_from = celltype2, values_from = percent_interaction, values_fill = NA)
  
  mat = as.matrix(all_df[,3:ncol(all_df)])
  rownames(mat) = all_df$marker
  mat = t(mat)
  
  list_mat[[dis]] = mat
  library(ComplexHeatmap)
  
  pdf(file.path(celltype1_dir, paste0("Protein_interaction_by_celltype_",dis,".pdf")),
      width = 15,
      height = 10)
  
  protein_cor = cor(mat)
  hc_cor = hclust(as.dist(1 - protein_cor), method = "ward.D2")
  
  h = Heatmap( 
    scale(mat),
    cluster_rows = F,
    cluster_columns = hc_cor, 
    name = paste0("Scaled % cells with a given protein \ninteraction by celltype - ", dis),
    column_title = "Proteins",
    row_title = "Celltypes",
    row_dend_side = "right",
    row_dend_width = unit(4, "cm"),
    column_dend_height = unit(5, "cm"),
    show_column_names = T,
    clustering_distance_columns ="pearson",
    clustering_distance_rows = "pearson",
    column_split = 10,
    border = TRUE,
    row_names_gp = gpar(fontsize = 8),
    use_raster = FALSE
  )
  draw(h)
  dev.off()
}

################################################################################
# Interactions per celltypes
################################################################################

for(celltype1 in c("T_helper", "T_cytotoxic", "T_regulatory", "Macrophages", "NKT")){
  celltype1_dir = file.path(output_dir, celltype1)
  all_interactions = qs::qread(file.path(celltype1_dir, paste0(celltype1, "_interactions.qs")))
  
  pdf(file.path(output_dir, paste0(celltype1, "_interactions_per_celltype.pdf")), width = 12, height = 8)
  for(mark in setdiff(rownames(spe), markers_to_remove) ){
    
    interactions = all_interactions %>% filter(marker == mark)
    
    # Plot boxplots:
    interaction_threshold = 5
    
    df = all_interactions
    df = df[which(df$condition %in% c("MF", "SS", "LP", "AD", "PS", "DAR") & df$n_interactions >= interaction_threshold),]
    
    df. = df %>% dplyr::filter(marker == mark)
    p = df. %>% ggplot(aes(x = celltype2, y = percent_interaction, fill = celltype2)) +
      geom_boxplot(outlier.colour = "white") +
      scale_fill_manual(values = setNames(unique(spe$cell_type_CellSighter_color), unique(spe$cell_type_CellSighter))) +
      ggnewscale::new_scale("fill") +
      geom_jitter(aes(fill = sample_id),  colour="black",pch=21, size=2, width = 0.3) +
      scale_fill_manual(values = setNames(unique(spe$sample_id_color), unique(spe$sample_id))) +
      theme_classic() + ggtitle(paste0("Enriched interactions - ", mark)) +
      ylab(paste0("Average % of ",celltype1," -  cells \n wich interaction is enriched for ", mark)) + xlab("") +
      stat_compare_means( aes(label = paste0(after_stat(p.format))), method = "t.test", ref.group = ".all.") + NoLegend() +
      theme(axis.text.x = element_text(angle = 90))
    print(p)
  }
  dev.off()
} 


################################################################################
# Interactions per celltypes, Responder vs Not Responder
################################################################################

list_DA_celltype1 = list()
dir = file.path(args$output, "Early_vs_Late/enriched_interactions")
dir.create(dir)
colors_progression_regression = setNames(c("#7699D4", "#EDCC5F", "#B54A48", "#74C276"),
                                         c("HD", "Inflammatory", "NonResponsive", "Responsive"))
samples = unique(spe$sample_id[which(!is.na(spe$CTCL.CAILS))])

for(celltype1 in c("T_helper", "T_cytotoxic", "T_regulatory", "Macrophages", "NKT")){
  celltype1_dir = file.path(output_dir, celltype1)
  all_interactions = qs::qread(file.path(celltype1_dir, paste0(celltype1, "_interactions.qs")))
  
  # Saving differential enrichments:
  interaction_threshold = 5
  
  list_DA_celltype = list()
  for(type2 in unique(spe$cell_type_CellSighter)){
    print(paste0(celltype1, " <--> ", type2))
    
    df = all_interactions %>% dplyr::filter(.data[["celltype2"]] == type2)
    df = df[which(df$sample_id %in% samples & df$n_interactions >= interaction_threshold),]
    df$CTCL.CAILS = spe$CTCL.CAILS[match(df$sample_id, spe$sample_id)]
    
    list_DA = list()
    for(mark in markers){
      df. = df %>% filter(marker == mark)
      
      if(nrow(df.[df.$CTCL.CAILS == "Responsive",])> 2 & nrow(df.[df.$CTCL.CAILS != "Responsive",]) > 2){
        list_DA[[mark]] = data.frame(
          celltype1 = celltype1,
          celltype2 = type2,
          marker = mark,
          condition = "Responsive_vs_NonResponsive",
          percent_interaction_Responsive = mean(df.$percent_interaction[df.$CTCL.CAILS == "Responsive"]),
          percent_interaction_NonResponsive = mean(df.$percent_interaction[df.$CTCL.CAILS != "Responsive"]),
          log2FC = log2(mean(df.$percent_interaction[df.$CTCL.CAILS == "Responsive"]) / mean(df.$percent_interaction[df.$CTCL.CAILS != "Responsive"])),
          effect = abs(mean(df.$percent_interaction[df.$CTCL.CAILS == "Responsive"]) - mean(df.$percent_interaction[df.$CTCL.CAILS != "Responsive"])),
          p.value = t.test(df.$percent_interaction[df.$CTCL.CAILS == "Responsive"], df.$percent_interaction[df.$CTCL.CAILS != "Responsive"])$p.value
        )
      }
      
    }
    list_DA_celltype[[type2]] = do.call("rbind", list_DA)
  }
  
  DA = do.call("rbind", list_DA_celltype)
  DA$log2FC[is.infinite(DA$log2FC)] = 0
  DA$p.value[is.nan(DA$p.value)] = 1
  DA$Significant = FALSE
  DA$Significant[abs(DA$log2FC) > log2(1.5) & DA$p.value < 0.05] = TRUE
  list_DA_celltype1[[celltype1]] = DA
}

WriteXLS::WriteXLS(list_DA_celltype1, ExcelFileName = file.path(dir, "Differential_Enriched_Interactions_All.xlsx"),
                   SheetNames = names(list_DA_celltype1), AdjWidth = T, AutoFilter = T, BoldHeaderRow = T)

# Plotting Boxplots ---------------------------------------------------------------------
library(ggpubr)
for(celltype1 in c("T_helper", "T_cytotoxic", "T_regulatory", "Macrophages", "NKT")){
  
  celltype1_dir = file.path(output_dir, celltype1)
  all_interactions = qs::qread(file.path(celltype1_dir, paste0(celltype1, "_interactions.qs")))
  
  # Plot boxplots:
  interaction_threshold = 5
  
  pdf(file.path(dir, paste0(celltype1,"-Significant.pdf")), width = 10, height = 8)
  DA = list_DA_celltype1[[celltype1]]
  
  for(i in which(DA$Significant == TRUE)){
    mark = DA$marker[i]
    type2 = DA$celltype2[i]
    
    df = all_interactions %>% dplyr::filter(.data[["celltype2"]] == type2)
    df = df[which((df$sample_id %in% samples | df$condition %in% c("LP", "AD", "PS", "DAR")) & df$n_interactions >= interaction_threshold),]
    df$CTCL.CAILS = as.character(spe$CTCL.CAILS[match(df$sample_id, spe$sample_id)])
    df$CTCL.CAILS[which(df$condition %in% c("LP", "AD", "PS", "DAR") )] = "Inflammatory"
    
    df. = df %>% dplyr::filter(marker == mark)
    p = df. %>% ggplot(aes(x = CTCL.CAILS, y = percent_interaction, fill = CTCL.CAILS)) +
      geom_boxplot(outlier.colour = "white") +
      scale_fill_manual(values = colors_progression_regression) +
      ggnewscale::new_scale("fill") +
      geom_jitter(aes(fill = sample_id),  colour="black",pch=21, size=2, width = 0.3) +
      scale_fill_manual(values = setNames(unique(spe$sample_id_color), unique(spe$sample_id))) +
      theme_classic() + ggtitle(paste0("Enriched interactions - ", mark)) +
      ylab(paste0("Average % of ",celltype1," - ",type2," cells \n wich interaction is enriched for ", mark)) + xlab("") +
      stat_compare_means(method = "t.test", ref.group = "NonResponsive")
    print(p)
  }
  dev.off()
} 

