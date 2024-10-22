# Analyze results from Pixie pixel clustering
# authors: Annotator3 Prompsy
# contact: Annotator3.prompsy@chuv.ch
# Guenova Lab
# CHUV (Centre Hospitalier Universitaire Vaudois), Lausanne, Suisse


cat("Analyzing results from CellSighter cell type detection... \n")

# Loading packages --------------------------------------------------------
libraries = c("argparse",
              "ggplot2",
              "dplyr",
              "tidyr",
              "Seurat",
              "SpatialExperiment",
              "arrow",
              "ChromSCape")
suppressPackageStartupMessages(invisible(lapply(libraries, require, character.only = TRUE)))
setwd("/mnt/RECHERCHE/GUENOVA_LAB/Project_Multiplex_Phenotyping/Miltenyi/Muliplex_Imaging_Pipeline/")

source("scripts/MultiplexImaging_utils.R")
source("scripts/GlobalVars.R")
library(ggplot2)
library(tidyverse)

# Reading in data  -------------------------------------------------------------
args = list(output = "output/")
cat("Output = ", args$output, "\n")

output_dir = file.path(args$output, "CellSighter", "marker")

corresponding_celltype = setNames(c("Negative","Positive"), c(0 , 1))
spe = qs::qread("output/SpatialExperiment.qs")


######################################################################################################
# Predictions 
######################################################################################################
list_deeplearning = list()
list_positive_marker = list()

# Predictions accuracy vs Positive Marker accuracy
predicitions_res = data.frame()
positive_marker_res = data.frame()
list_confusion_predictions = list()

output_dir = file.path("output", "CellSighter", "marker")

for(marker in rownames(spe)[which(!rownames(spe) %in% c("CD270", "CXCR4", "CD209", "Ki67"))]){
  for(sample in unique(spe$sample_id)){
    
      num = as.numeric(gsub(".*-", "", sample))
      best_iter_file = file.path(output_dir,  "marker_classification", marker,  "best_iter.txt")
      
      if(file.exists(best_iter_file)){
        best_iter = read.table(best_iter_file)
      } else {
        best_iter = read.table(file.path(output_dir,  "marker_classification",
                                         marker,  ifelse(num > 28, "batch_2", "batch_1"), "best_iter.txt"))
      }
      pred_file = file.path(output_dir,  "marker_classification", marker,  paste0("val_results_",best_iter[1,1],".csv"))
      pred_file_batch = file.path(output_dir,  "marker_classification", marker, 
                                  ifelse(num > 28, "batch_2", "batch_1"),
                                  paste0("val_results_",best_iter[1,1],".csv"))
      
      if(file.exists(pred_file) | file.exists(pred_file_batch)){
        if(file.exists(pred_file)) labels = read.csv(pred_file)
        if(file.exists(pred_file_batch)) labels = read.csv(pred_file_batch)

        labels$cell_id = paste0(gsub("-Annotator1|-Annotator2|-Annotator3|-Annotator4", "", labels$image_id), "-", labels$cell_id)
        predictions = labels
        
        predictions$marker_pred = corresponding_celltype[match(predictions$pred, as.numeric(names(corresponding_celltype)))]
        predictions$label = labels$label[match(predictions$cell_id, labels$cell_id)]
        predictions$marker_label = corresponding_celltype[match(predictions$label, as.numeric(names(corresponding_celltype)))]
        
        # Deep Learning ------------------------------------------------------------
        # calculate confusion matrix
        confusion_matrix <- table(as.factor(predictions$marker_pred), 
                                  as.factor(predictions$marker_label))
        
        # convert confusion matrix to data frame
        confusion_df <- as.data.frame(confusion_matrix)
        list_confusion_predictions[[marker]]  = confusion_df
        
        # rename columns
        colnames(confusion_df) <- c("True.Class", "Predicted.Class", "Count")
        
        # create heatmap with text labels
        list_deeplearning[[marker]] = (ggplot(data = confusion_df, aes(x = Predicted.Class, y = True.Class, fill = log10(Count+1), label = Count)) + 
                              geom_tile() + 
                              scale_fill_viridis_c(begin = 0.4, direction = -1) +
                              labs(title = marker,
                                   x = "True Class", y = "Predicted Class", fill = "Count") +
                              geom_text(color = "white", size = 5)  + theme(axis.text.x = element_text(angle = 90)) +theme_classic()
        ) + ggtitle(marker)
        
        # Accuracy, Sensitivity, Precision
        TP_p = confusion_df$Count[4]
        FP_p = confusion_df$Count[2]
        TN_p = confusion_df$Count[1]
        FN_p = confusion_df$Count[3]
        
        predicitions_res = rbind(predicitions_res, data.frame(
          marker = marker,
          type = "deeplearning",
          accuracy = (TP_p + TN_p) / (TP_p + FP_p + TN_p + FN_p),
          precision = TP_p / (TP_p + FP_p),
          recall = TP_p / (TP_p + FN_p)
        ))
        
        # Positive Marker  ---------------------------------------------------------
        # calculate confusion matrix
        pos_mat = spe@assays@data$positive_marker
        predictions$positive_marker = pos_mat[marker, match(predictions$cell_id, spe$cell_id)]
        confusion_matrix <- table(as.factor(predictions$positive_marker), 
                                  as.factor(predictions$marker_label))
        
        # convert confusion matrix to data frame
        confusion_df <- as.data.frame(confusion_matrix)
        list_confusion_predictions[[marker]]  = confusion_df
        
        # rename columns
        colnames(confusion_df) <- c("True.Class", "Predicted.Class", "Count")
        
        list_positive_marker[[marker]] = (ggplot(data = confusion_df, aes(x = Predicted.Class, y = True.Class, fill = log10(Count+1), label = Count)) + 
                                            geom_tile() + 
                                            scale_fill_viridis_c(begin = 0.4, direction = -1) +
                                            labs(title = marker,
                                                 x = "True Class", y = "Predicted Class", fill = "Count") +
                                            geom_text(color = "white", size = 5)  + theme(axis.text.x = element_text(angle = 90)) +theme_classic()
        ) + ggtitle(marker)
        
        # Accuracy, Sensitivity, Precision
        TP_m = confusion_df$Count[4]
        FP_m = confusion_df$Count[2]
        TN_m = confusion_df$Count[1]
        FN_m = confusion_df$Count[3]
        
        positive_marker_res = rbind(positive_marker_res, data.frame(
          marker = marker,
          type = "Distribution-based",
          accuracy = (TP_m + TN_m) / (TP_m + FP_m + TN_m + FN_m),
          precision = TP_m / (TP_m + FP_m),
          recall = TP_m / (TP_m + FN_m)
        ))
      }
}
}

results = rbind(predicitions_res, positive_marker_res)
results = results %>% group_by(marker, type) %>%
  summarise(accuracy = mean(accuracy, na.rm = TRUE),
            precision = mean(precision, na.rm = TRUE),
            recall = mean(recall, na.rm = TRUE)
            )
results = results[which(!results$marker %in% results$marker[is.nan(results$accuracy)]),]
results$type = factor(results$type, levels = c("Distribution-based", "deeplearning"))


png(file.path("output/CellSighter/Plots", paste0("positive_marker_vs_pred_accuracy.png")), height = 1200, width = 3000, res = 200)
ggplot(data = results, aes(x = marker, y = accuracy, fill =type)) + 
  geom_bar(stat = "identity", position = "dodge", width = 0.7) + theme_bw() + 
  scale_fill_manual(values = c("#a51d23ff", "#dab600"))+
  labs(title = "Accuracy",
       x = "", y = "Accuracy", fill = "Type") + theme(axis.text.x = element_text(angle = 90))
dev.off()

png(file.path("output/CellSighter/Plots", paste0("positive_marker_vs_pred_precision.png")), height = 1200, width = 3000, res = 200)
ggplot(data = results, aes(x = marker, y = precision, fill =type)) + 
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +  theme_bw() + 
  scale_fill_manual(values = c("#a51d23ff", "#dab600"))+
  labs(title = "Precision",
       x = "", y = "Precision", fill = "Type") + theme(axis.text.x = element_text(angle = 90))
dev.off()

png(file.path("output/CellSighter/Plots", paste0("positive_marker_vs_pred_recall.png")), height = 1200, width = 3000, res = 200)
ggplot(data = results, aes(x = marker, y = recall, fill =type)) + 
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +  theme_bw() + 
  scale_fill_manual(values = c("#a51d23ff", "#dab600"))+
  labs(title = "Recall",
       x = "", y = "Recall", fill = "Type") + theme(axis.text.x = element_text(angle = 90))
dev.off()

library(ggpubr)
png(file.path("output/CellSighter/Plots", paste0("positive_marker_vs_pred_accuracy_overall.png")), height = 1200, width = 1000, res = 200)
ggplot(data = results, aes(x = type, y = accuracy, fill = type)) + 
  geom_boxplot(outlier.colour = "white")  + geom_jitter(size = 0.1, width = 0.25) + theme_bw() + 
  scale_fill_manual(values = c("#a51d23ff", "#dab600"))+
  labs(title = "Accuracy",
       x = "", y = "Accuracy", fill = "Type") + theme(axis.text.x = element_text(angle = 90)) + 
  stat_compare_means(aes(label = after_stat(p.format)),
                     method = "t.test",
                     ref.group = "Distribution-based", paired = TRUE) + ylab("F1 score")
dev.off()

library(ggpubr)
results = results %>%
  mutate(F1 = 2 * (precision * recall) / (precision + recall) )

png(file.path("output/CellSighter/Plots", paste0("positive_marker_vs_pred_precision_overall.png")), height = 1200, width = 750, res = 250)
ggplot(data = results, aes(x = type, y = precision, fill = type)) + 
  geom_boxplot(outlier.colour = "white")  + geom_jitter(size = 0.1, width = 0.25) + theme_classic() + 
  scale_fill_manual(values = c("#a51d23ff", "#dab600"))+
  labs(title = "Precision",
       x = "", y = "Precision", fill = "Type") + theme(axis.text.x = element_text(angle = 90)) + 
  stat_compare_means(aes(label = after_stat(p.format)),
                     method = "t.test", ref.group = "Distribution-based", paired = TRUE) + ylab("Precision") +
  theme(legend.position = "none")
dev.off()

png(file.path("output/CellSighter/Plots", paste0("positive_marker_vs_pred_recall_overall.png")),  height = 1200, width = 750, res = 250)
ggplot(data = results, aes(x = type, y = recall, fill = type)) + 
  geom_boxplot(outlier.colour = "white")  + geom_jitter(size = 0.1, width = 0.25) + theme_classic() + 
  scale_fill_manual(values = c("#a51d23ff", "#dab600"))+
  labs(title = "Recall",
       x = "", y = "Recall", fill = "Type") + theme(axis.text.x = element_text(angle = 90)) + 
  stat_compare_means(aes(label =  after_stat(p.format)),
                     method = "t.test", ref.group = "Distribution-based", paired = TRUE) + ylab("Recall")   +
  theme(legend.position = "none")
dev.off()

png(file.path("output/CellSighter/Plots", paste0("positive_marker_vs_pred_accuracy_overall.png")),  height = 1200, width = 750, res = 250)
ggplot(data = results, aes(x = type, y = accuracy, fill = type)) + 
  geom_boxplot(outlier.colour = "white")  + geom_jitter(size = 0.1, width = 0.25) + theme_classic() + 
  scale_fill_manual(values = c("#a51d23ff", "#dab600"))+
  labs(title = "Accuracy",
       x = "", y = "Accuracy", fill = "Type") + theme(axis.text.x = element_text(angle = 90)) + 
  stat_compare_means(aes(label = after_stat(p.format)),
                     method = "t.test", ref.group = "Distribution-based", paired = TRUE) + ylab("Accuracy") +
  theme(legend.position = "none")
dev.off()

results %>% arrange(desc(accuracy)) %>% head(10)
results %>% arrange((accuracy)) %>% head(50)
results %>% arrange(desc(accuracy)) %>% filter(type != "deeplearning") %>% head(50)


png(file.path("output/CellSighter/Plots", paste0("positive_marker_vs_pred_F1.png")), height = 1200, width = 750, res = 250)
ggplot(data = results, aes(x = type, y = F1, fill = type)) + 
  geom_boxplot(outlier.colour = "white") + geom_jitter(size = 0.1, width = 0.25) + theme_classic() + 
  scale_fill_manual(values = c("#a51d23ff", "#dab600"))+
  labs(title = "F1 score",
       x = "", y = "F1", fill = "Type") + theme(axis.text.x = element_text(angle = 90)) + 
  stat_compare_means(aes(label = after_stat(p.format)),
                     method = "t.test", ref.group = "Distribution-based", paired = TRUE) + ylab("F1 score") +
  theme(legend.position = "none")
dev.off()


mean(results$accuracy[results$type == "deeplearning"])
mean(results$accuracy[results$type == "Distribution-based"], na.rm = T)
mean(results$F1[results$type == "deeplearning"])
mean(results$F1[results$type == "Distribution-based"], na.rm = T)


results %>% arrange(desc(F1)) %>% head(10)
results %>% arrange(desc(F1))  %>% filter(type == "deeplearning") %>% head(10)
results %>% arrange((F1)) %>% head(10)
results %>% arrange(desc(F1))  %>% filter(type != "deeplearning")  %>% head(10)

######################################################################################################
# Comparing Manual Annotations
######################################################################################################

## Compare two annotators ------------------------------------------------------

pdf(file.path(output_dir, "..", "Plots", "ROI-05-Annotator1_vs_Annotator2.pdf"), width = 6, height = 5)
for(marker in rownames(spe)){
  labels = read.csv(file.path(output_dir,  "marker_classification", marker, "val_results_50.csv"))
  labels$cell_id = paste0(gsub("-Annotator1|-Annotator2|-Annotator3|-Annotator4", "", labels$image_id), "-", labels$cell_id)
  
  predictions = labels
  predictions$marker_pred = corresponding_celltype[match(predictions$pred, as.numeric(names(corresponding_celltype)))]
  predictions$marker_label = corresponding_celltype[match(predictions$label, as.numeric(names(corresponding_celltype)))]
  
  lab_chris = predictions[which(predictions$image_id == "ROI-05-Annotator1"),]
  lab_ion = predictions[which(predictions$image_id == "ROI-05-Annotator2"),]
  common = intersect(lab_chris$cell_id, lab_ion$cell_id)
  if(length(common) > 5){
    lab_chris = lab_chris[match(common, lab_chris$cell_id),]
    lab_ion = lab_ion[match(common, lab_ion$cell_id),]
    
    # calculate confusion matrix
    confusion_matrix <- table(as.factor(lab_ion$label), 
                              as.factor(lab_chris$label))
    # convert confusion matrix to data frame
    confusion_df <- as.data.frame(confusion_matrix)
    list_confusion_predictions[[marker]]  = confusion_df
    
    # rename columns
    colnames(confusion_df) <- c("True.Class", "Predicted.Class", "Count")
    
    # create heatmap with text labels
    p = ggplot(data = confusion_df, aes(x = Predicted.Class, y = True.Class, fill = log10(Count+1), label = Count)) + 
      geom_tile() + 
      scale_fill_viridis_c(begin = 0.4, direction = -1) +
      labs(title = marker,
           x = "Annotator2 Class", y = "Annotator1 Class", fill = "Count") +
      geom_text(color = "white", size = 5)  + theme(axis.text.x = element_text(angle = 90)) +theme_classic() +
      ggtitle(marker)
    print(p)
  }
}
dev.off()

######################################################################################################
# Assigning Predicted Positive Marker to Spatial Experiment and save
######################################################################################################

# Assigning predicted positive marker to SpatialExperiment ---------------------------
CellSighter_marker_mat = spe@assays@data$positive_marker
CellSighter_marker_mat[] = 0

for(sample in sort(unique(spe$sample_id))){
  
  for(marker in rownames(CellSighter_marker_mat)){
    pred_file = file.path(output_dir, "Predictions", marker, paste0(sample, "_val_results.csv"))
    samp_num = as.numeric(gsub(".*-0|.*-", "",sample)) 
    pred_file_batch = file.path(output_dir, "Predictions", marker, ifelse(samp_num <= 28 | samp_num >63, "batch_1", "batch_2"), paste0(sample, "_val_results.csv"))
    
    if(file.exists(pred_file) | file.exists(pred_file_batch)){
      
    if(file.exists(pred_file)){
      predictions = read.csv(pred_file)
    } else if(file.exists(pred_file_batch)) {
      predictions = read.csv(pred_file_batch)
    }
    
    predictions$sample_cell_id = paste0(sample, "-", predictions$cell_id)
    predictions = predictions[predictions$image_id == sample,]
    CellSighter_marker_mat[marker,spe$sample_id == sample] = predictions$pred[match(colnames(CellSighter_marker_mat[,spe$sample_id == sample]),
                                                                                    predictions$sample_cell_id)]
    } else{
      cat(marker, " ", sample, "\n")
      CellSighter_marker_mat[marker,spe$sample_id == sample] = -1
    }
  }
}
spe@assays@data$CellSighter_marker_mat = CellSighter_marker_mat
qs::qsave(spe, "output/SpatialExperiment.qs")

######################################################################################################
# Count number of manually annotated markers + and -:
######################################################################################################

marker_dir = list.dirs("output/CellSighter/marker/marker_classification/", recursive = F)
names(marker_dir) = basename(marker_dir)
marker_dir = marker_dir[which(names(marker_dir) %in% rownames(spe))]

annot_images = c("ROI-01-Annotator3",  "ROI-10-Annotator3",  "ROI-20-Annotator3",  "ROI-40-Annotator3", "ROI-50-Annotator3", "ROI-60-Annotator3",
                 "ROI-07-Annotator1", "ROI-15-Annotator3","ROI-05-Annotator1", "ROI-05-Annotator2", "ROI-09-Annotator3",  "ROI-30-Annotator3")

mat_plus = matrix(0, nrow = nrow(spe), ncol = length(annot_images))
colnames(mat_plus) = annot_images
rownames(mat_plus) = rownames(spe)
mat_neg = mat_plus

for(marker in rownames(spe)){
  for(annot_dir in annot_images){
    if(file.exists(file.path("output/CellSighter/manual_annotation_marker/", annot_dir, paste0(marker, "+.csv"))))
      pos = nrow(read.csv(file.path("output/CellSighter/manual_annotation_marker/", annot_dir, paste0(marker, "+.csv"))))
    mat_plus[marker, annot_dir] = pos 
    
    if(file.exists(file.path("output/CellSighter/manual_annotation_marker/", annot_dir, paste0(marker, "-.csv"))))
      neg = nrow(read.csv(file.path("output/CellSighter/manual_annotation_marker/", annot_dir, paste0(marker, "-.csv"))))
    mat_neg[marker, annot_dir] = neg
  }
}

write.csv(mat_plus, "output/CellSighter/manual_annotation_marker/number_of_annotated_cells_markers+.csv")
write.csv(mat_neg, "output/CellSighter/manual_annotation_marker/number_of_annotated_cells_markers-.csv")
write.csv(mat_plus + mat_neg, "output/CellSighter/manual_annotation_marker/number_of_annotated_cells_markers_both.csv")

median(rowSums(mat_plus))
median(rowSums(mat_neg))
median(rowSums(mat_plus))
