# Analyze results from Pixie pixel clustering
# authors: Annotator3 Prompsy
# contact: Annotator3.prompsy@chuv.ch
# Guenova Lab
# CHUV (Centre Hospitalier Universitaire Vaudois), Lausanne, Suisse

cat("Analyzing results from CellSighter cell type detection... \n")
setwd("/mnt/RECHERCHE/GUENOVA_LAB/Project_Multiplex_Phenotyping/Miltenyi/Muliplex_Imaging_Pipeline/")

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

source("scripts/MultiplexImaging_utils.R")
source("scripts/GlobalVars.R")
library(ggplot2)
library(tidyverse)

# Reading in data  -------------------------------------------------------------
args = list(output = "output/")
cat("Output = ", args$output, "\n")

output_dir = file.path(args$output, "CellSighter", "celltype")

config = jsonlite::parse_json(file(file.path(output_dir, "cell_classification", "config.json")))
corresponding_celltype = unlist(config$hierarchy_match)

spe = qs::qread("output/SpatialExperiment.qs")

######################################################################################################
# Validation of training
######################################################################################################

list_confusion_predictions= list()

# Load data
validation_images = c("ROI-07-Annotator1", "ROI-06-Annotator1","ROI-15-Annotator1", "ROI-17-Annotator3")
labels = do.call("rbind", lapply(validation_images, function(i) read.csv(file.path(output_dir,  "Predictions", paste0(i,"_val_results.csv")))))
labels$image_id = gsub("-Annotator1|-Annotator3","",labels$image_id)
labels$cell_id = paste0(labels$image_id,"-", labels$cell_id)
predictions = labels

predictions$marker_pred = corresponding_celltype[match(predictions$pred, as.numeric(names(corresponding_celltype)))]
predictions$marker_label = corresponding_celltype[match(predictions$label, as.numeric(names(corresponding_celltype)))]

# calculate confusion matrix
confusion_matrix <- table(as.factor(predictions$marker_pred),
                          as.factor(predictions$marker_label))

# convert confusion matrix to data frame
confusion_df <- as.data.frame(confusion_matrix)

# rename columns
colnames(confusion_df) <- c("True.Class", "Predicted.Class", "Count")

# confusion_df = scale(confusion_df)
# create heatmap with text labels
# png(file.path("output/CellSighter/Plots", paste0("predictions_vs_labels_marker_",marker,".png")), height = 1200, width = 1200, res = 300)
p = (ggplot(data = confusion_df, aes(x = Predicted.Class, y = True.Class, fill = log10(Count+1), label = Count)) +
       geom_tile() +
       scale_fill_gradient2(midpoint = 0.5, low = "royalblue4", mid = "royalblue4",high = "gold") +
       labs(title = "Celltype confusion matrix",
            x = "True Class", y = "Predicted Class", fill = "log10(Cell Overlap + 1)") +
       geom_text(color = "white", size = 3.5)  + theme(axis.text.x = element_text(angle = 90)) +theme_classic()
) + ggtitle("CellSighter CellType identification ") + theme(axis.text.x = element_text(angle=90))

png(file.path("output/CellSighter/Plots", paste0("positive_marker_vs_pred_accuracy_CellType_all.png")),
    height = 1600, width = 2100, res = 300)
p
dev.off()

# High accuracy only
validation_images = c("ROI-07-Annotator1", "ROI-06-Annotator1","ROI-15-Annotator1", "ROI-17-Annotator3")
labels = do.call("rbind", lapply(validation_images, function(i) read.csv(file.path(output_dir,  "Predictions", paste0(i,"_true_labels.csv")))))
labels$cell_id = paste0(labels$fov,"-", labels$label)
labels$cell_type[labels$cell_type == "Fibroblast"] = "Unknown-Stroma"
labels$cell_type[labels$cell_type == ""] = NA
spe. = spe[,spe$sample_id %in%  c("ROI-07", "ROI-06", "ROI-15", "ROI-17")]
spe.$true_label = labels$cell_type[match(spe.$cell_id, labels$cell_id)]

heatmap(as.matrix(table(spe.$true_label, spe.$cell_type_CellSighter)), Rowv = NA, Colv = NA)

# calculate confusion matrix
confusion_matrix <- table(as.factor(spe..$cell_type_CellSighter),
                          as.factor(spe..$true_label))

# convert confusion matrix to data frame
confusion_df <- as.data.frame(confusion_matrix)

# rename columns
colnames(confusion_df) <- c("True.Class", "Predicted.Class", "Count")

p = (ggplot(data = confusion_df, aes(x = Predicted.Class, y = True.Class, fill = log10(Count+1), label = Count)) +
       geom_tile() +
       scale_fill_gradient2(midpoint = 0.5, low = "royalblue4", mid = "royalblue4",high = "gold") +
       labs(title = "Celltype confusion matrix",
            x = "True Class", y = "Predicted Class", fill = "log10(Cell Overlap + 1)") +
       geom_text(color = "white", size = 3.5)  + theme(axis.text.x = element_text(angle = 90)) +theme_classic()
) + ggtitle("CellSighter CellType identification ") + theme(axis.text.x = element_text(angle=90))

png(file.path("output/CellSighter/Plots", paste0("positive_marker_vs_pred_accuracy_CellType_after_refining.png")),
    height = 1600, width = 2100, res = 300)
p
dev.off()

library(caret)
spe.. = spe.[,!is.na(spe.$true_label) & spe.$true_label != "Unknown-Stroma" &
               spe.$cell_type_CellSighter != "Unknown-Stroma"]
pred_factor = as.factor(spe..$true_label)
pred_label = as.factor(spe..$cell_type_CellSighter)

# Create a confusion matrix
confusion_matrix <- confusionMatrix(pred_factor, pred_label)

# Calculate accuracy
accuracy <- confusion_matrix$overall["Accuracy"]

# Calculate precision, recall (sensitivity), and F1-score for each class
precision <- confusion_matrix$byClass[, "Precision"]
recall <- confusion_matrix$byClass[, "Recall"]
f1_score <- confusion_matrix$byClass[, "F1"]

# Print the results
cat("Confusion Matrix:\n")
print(confusion_matrix$table)

cat("\nAccuracy:", accuracy)
mean(accuracy, na.rm = T)

cat("\n\nPrecision (Per Class):\n")
print(precision)
mean(precision, na.rm = T)

cat("\nRecall (Per Class):\n")
print(recall)
mean(recall, na.rm = T)

cat("\nF1-Score (Per Class):\n")
print(f1_score)
mean(f1_score, na.rm = T)


######################################################################################################
# Training metrics to determine the best training set
######################################################################################################

pdf(file.path("output/CellSighter/", "Plots", paste0("results_marker_classification_0.5_CD195_10_to_150_metrics_used_for_training.pdf")), width = 16, height = 6)

df = data.frame(sample_id = unique(predictions$image_id),
                class = c(rep("CellSighter", 4), rep("Histogram-based",4)),
                accuracy = 0, f1_score = 0)
for(i in unique(predictions$image_id)){

  # Cell Sighter
  predictions. = predictions[predictions$image_id == i,]
  # High accuracy only
  predictions. = predictions.[predictions.$pred_prob > 0.9,]

  pred_factor = as.factor(predictions.$marker_pred)
  pred_label = factor(predictions.$marker_label, levels = levels(pred_factor))
  
  # Create a confusion matrix
  confusion_matrix <- confusionMatrix(pred_factor, pred_label)
  df$accuracy[df$sample_id == i & df$class == "CellSighter"] <- mean(confusion_matrix$overall["Accuracy"], na.rm = T)
  df$f1_score[df$sample_id == i & df$class == "CellSighter"] <- mean(confusion_matrix$byClass[, "F1"], na.rm = T)
  
  # Histogram
  predictions. = predictions[predictions$image_id == i,]
  
  pred_factor = as.factor(spe$cell_type_positive_marker[match(predictions.$cell_id, spe$cell_id)])
  pred_label = factor(predictions.$marker_label, levels = levels(pred_factor))
  
  # Create a confusion matrix
  confusion_matrix <- confusionMatrix(pred_factor, pred_label)
  df$accuracy[df$sample_id == i & df$class == "Histogram-based"] <- mean(confusion_matrix$overall["Accuracy"], na.rm = T)
  df$f1_score[df$sample_id == i& df$class == "Histogram-based"] <- mean(confusion_matrix$byClass[, "F1"], na.rm = T)
}


p = ggplot(data = df, aes(x = class, y = accuracy)) +
  geom_boxplot() +
  labs(title = "CellType Accuracy",
       x = "Training Iteration", y = "TP / (TP + FP)") +
  theme(axis.text.x = element_text(angle = 90)) + theme_classic() + 
  geom_hline(yintercept = 0.85) +
  geom_hline(yintercept = 0.80) +
  ylim(0.5,1) + NoLegend()
print(p)
dev.off()


######################################################################################################
# Assigning predicted celltypes to SpatialExperiment and saving
######################################################################################################

spe = qs::qread("output/SpatialExperiment.qs")
spe$cell_type_CellSighter = "Unknown"
spe$cell_type_CellSighter_prediction_prob = 0

##### Refine Unknown (Unknown-Stroma) using deeplearning cell marker detection
cell_markers = read.csv("annotation/cell_markers.csv")
cell_markers = cell_markers[-which(cell_markers$marker %in% c("CollagenI", "DAPI")),]
cell_markers = cell_markers %>% filter(important == TRUE)
cell_markers = cell_markers[nrow(cell_markers):1,]
cell_type = rep("Unknown-Stroma", ncol(spe))
names(cell_type) = spe$cell_id


pos_mark = spe@assays@data$CellSighter_marker_mat

for(type in unique(cell_markers$cell_type)){
  df = cell_markers[cell_markers$cell_type == type,]
  vec_cell_type = setNames(rep(TRUE, ncol(spe)), colnames(spe))
  for(i in 1:nrow(df)){
    if(df$positive[i]) vec_cell_type = vec_cell_type & (pos_mark[df$marker[i],] == 1)
    else vec_cell_type = vec_cell_type & (pos_mark[df$marker[i],] == 0)
  }
  cell_type[vec_cell_type] = type
  print(table(cell_type))
}

# For each sample assign CellSighter celltype and refine Unknown (Fibroblasts),
# Monocytic Lineage (Monocyte, Macrophages, Neutrophils) and CD45RA+ (T helper,
# T regulatory, T cytotocic)
# Finally, set Cytokeratin+CD16+ to Keratinocytes due to high CD16 background in
# Epidermis.
avg_confident_cells = c()
avg_confident_cells_after_reassignation = c()
for(samp in unique(spe$sample_id)){
  
  # Select sample
  spe. = spe[, spe$sample_id == samp]
  
  # Read the predictions
  predictions = read.csv(file.path(output_dir,  "Predictions", paste0(samp, "_val_results.csv")))
  predictions$sample_cell_id = paste0(predictions$image_id, "-", predictions$cell_id)
  predictions$cell_type_pred = corresponding_celltype[match(predictions$pred, as.numeric(names(corresponding_celltype)))]
  predictions$cell_type_pred[predictions$cell_type_pred=="Fibroblast"] = "Unknown-Stroma"
  
  # Assign the predictions to Spatial Experiment
  spe.$cell_type_CellSighter = predictions$cell_type_pred[match(spe.$cell_id, predictions$sample_cell_id)]
  spe.$CellSighter_prediction_prob = predictions$pred_prob[match(spe.$cell_id, predictions$sample_cell_id)]
  
  # Calculate the number of predictions above the probability threshold
  ntot = nrow(predictions)
  n = nrow(predictions[predictions$pred_prob > 0.9,])
  print(samp)
  print(100 * n/ntot)
  avg_confident_cells = c(avg_confident_cells,  n/ntot)
  
  # Refine the cell types that have probability below threshold based on 
  # CellSighter positive marker identification:
  # e.g. If a T_helper cell has prob < threshold but is assigned CD4+, we can 
  # assign it to T_helper with confidence
  below_threshold_cells = predictions$sample_cell_id[predictions$pred_prob < 0.9]
  spe_low_prob = spe.[,spe.$cell_id %in% below_threshold_cells]
  
  for(type in unique(cell_markers$cell_type)){
    spe_low_prob. = spe_low_prob[,spe_low_prob$cell_type_CellSighter == type]
    pos_mark = spe_low_prob.@assays@data$CellSighter_marker_mat
    
    if(nrow(spe_low_prob.) > 0){
      df = cell_markers[cell_markers$cell_type == type,]
      vec_cell_type = setNames(rep(TRUE, ncol(spe_low_prob.)), colnames(spe_low_prob.))
      for(i in 1:nrow(df)){
        if(df$positive[i]) vec_cell_type = vec_cell_type & (pos_mark[df$marker[i],] == 1)
        else vec_cell_type = vec_cell_type & (pos_mark[df$marker[i],] == 0)
      }
      
      # If any cells with low confidence has not the marker, assign it to "Unknown" 
      if(length(which(vec_cell_type == FALSE)) > 0){
        vec_cell_type = vec_cell_type[which(vec_cell_type == FALSE)]
        spe.$cell_type_CellSighter[match(names(vec_cell_type), spe.$cell_id)] = "Unknown"
      }
    }
  }
  
  # Refine the Unknown and Stroma cells
  unknown_or_stroma = spe.$cell_id[which(spe.$cell_type_CellSighter == "Unknown-Stroma" | spe.$cell_type_CellSighter == "Unknown")]
  
  spe.$cell_type_CellSighter[match(unknown_or_stroma,spe.$cell_id)] = cell_type[unknown_or_stroma]
 
  # Calculate the rate of overall assigned cells after re-assignation of below threshold cells
  ratio = length(which(spe.$cell_type_CellSighter != "Unknown-Stroma"))  / ncol(spe.)
  avg_confident_cells_after_reassignation = c(avg_confident_cells_after_reassignation, ratio)
  
  # Refine Monocytic Lineage
  spe_mono = spe.[,which(spe.$cell_type_CellSighter == "Monocytic_Lineage")]
  cell_type_mono = spe_mono$cell_type_CellSighter
  names(cell_type_mono) = spe_mono$cell_id
  pos_mark = spe_mono@assays@data$CellSighter_marker_mat
  
  for(type in c("Macrophages", "Monocytes", "Neutrophils") ){
    df = cell_markers[cell_markers$cell_type == type,]
    vec_cell_type = setNames(rep(TRUE, ncol(spe_mono)), colnames(spe_mono))
    for(i in 1:nrow(df)){
      if(df$positive[i]) vec_cell_type = vec_cell_type & (pos_mark[df$marker[i],] == 1)
      else vec_cell_type = vec_cell_type & (pos_mark[df$marker[i],] == 0)
    }
    cell_type_mono[vec_cell_type] = type
  }
  
  mono = spe.$cell_id[which(spe.$cell_type_CellSighter == "Monocytic_Lineage")]
  spe.$cell_type_CellSighter[match(mono,spe.$cell_id)] = cell_type_mono[mono]
  
  # Refine Leukocyte
  spe_leukocyte = spe.[,which(spe.$cell_type_CellSighter == "Leukocyte")]
  cell_type_leukocyte = spe_leukocyte$cell_type_CellSighter
  names(cell_type_leukocyte) = spe_leukocyte$cell_id
  pos_mark = spe_leukocyte@assays@data$CellSighter_marker_mat
  
  for(type in c("T_helper", "T_regulatory", "T_cytotoxic", 
                "Neutrophils",  "Monocytes", "Macrophages",
                "Basophil", "pDC", "NKT", "B_cell")){
    df = cell_markers[cell_markers$cell_type == type,]
    vec_cell_type = setNames(rep(TRUE, ncol(spe_leukocyte)), colnames(spe_leukocyte))
    for(i in 1:nrow(df)){
      if(df$positive[i]) vec_cell_type = vec_cell_type & (pos_mark[df$marker[i],] == 1)
      else vec_cell_type = vec_cell_type & (pos_mark[df$marker[i],] == 0)
    }
    cell_type_leukocyte[vec_cell_type] = type
  }
  
  leukocyte = spe.$cell_id[which(spe.$cell_type_CellSighter == "Leukocyte")]
  spe.$cell_type_CellSighter[match(leukocyte,spe.$cell_id)] = cell_type_leukocyte[leukocyte]
  
  # Assign CD16+Cytokeratin+ neutrophils cells to keratinocytes
  neutro = spe.$cell_id[which(spe.$cell_type_CellSighter == "Neutrophils")]
  spe_neutro = spe.[,neutro]
  pos_mark = spe_neutro@assays@data$CellSighter_marker_mat
  spe_neutro$cell_type_CellSighter[pos_mark["Cytokeratin",] == 1] = "Keratinocyte"
  spe.$cell_type_CellSighter[match(neutro,spe.$cell_id)] = spe_neutro$cell_type_CellSighter
  
  # Print total reassinged
  ntot = ncol(spe.)
  n = length(which(!is.na(spe.$cell_type_CellSighter) & (spe.$cell_type_CellSighter != "Unknown-Stroma")))
  print(n/ntot)
  
  # Assign to all spe
  spe$cell_type_CellSighter[match(spe.$cell_id, spe$cell_id)] = spe.$cell_type_CellSighter
  spe$cell_type_CellSighter_prediction_prob[match(spe.$cell_id, spe$cell_id)] = spe.$CellSighter_prediction_prob
  
  # Plotting
  cell_overlay_mat = getImageAsMatrix(file.path("output", "segmentation", paste0(samp,"_whole_cell.tiff")))
  celltype_img_mat = get_metadata_image_overlay(spe., cell_overlay_mat, sample = samp, metadata = "cell_type_CellSighter", levels =  cell_type_color_df$cell_type_CellSighter)
  
  pdf(file.path(output_dir, "Predictions", "Masks", paste0(samp,"_cell_type.pdf")), height = 5, width = 6) # output_dir, "Predictions", "Masks",
  rgb_color_cell_raster(celltype_img_mat, color_df = cell_type_color_df, samp) 
  dev.off()
  
  tiff::writeTIFF(as.matrix(celltype_img_mat),
                  file.path(output_dir, "Predictions", "Masks", paste0(samp,"_cell_type.tiff")),
                  bits.per.sample = 8)

}

confidence_df = data.frame("avg_confident_cells" = avg_confident_cells,
                           "avg_confident_cells_after_reassignation" = avg_confident_cells_after_reassignation,
                           sample_id = unique(spe$sample_id),
                           condition = spe$condition[match(unique(spe$sample_id), spe$sample_id)]
                           )

png(file.path(output_dir, "Confidence_prediction.png"),  width = 3000, height = 1300, res = 300)
confidence_df %>% pivot_longer(cols = c(avg_confident_cells, avg_confident_cells_after_reassignation), values_to = "Ratio", 
                               names_to = "Assignation") %>%
  ggplot(aes(x = sample_id, y = 100 * Ratio, fill =  Assignation)) + theme_classic() + 
  geom_bar(stat = "identity", position = position_dodge2()) + theme(axis.text.x = element_text(angle = 90)) +
  labs(y = "% Cells with prediction probability > 0.9", x = "") +
  geom_hline(yintercept =  100* median(confidence_df$avg_confident_cells), col = "#D66462", lty = 2) +
  geom_hline(yintercept =  100* median(confidence_df$avg_confident_cells_after_reassignation), col = "#36ADBA", lty = 2)
dev.off()

mean(avg_confident_cells) # 0.75
mean(avg_confident_cells_after_reassignation) # 0.87
spe = colors_scExp(spe, annotCol = "cell_type_CellSighter", color_by = "cell_type_CellSighter", color_df = cell_type_color_df)

par(mai = c(1,2,1,1))
tab = table(spe$cell_type_CellSighter)
tab = tab[c("Unknown-Stroma","Leukocyte", "B_cell", "T_helper", "T_cytotoxic", "T_regulatory", 
            "pDC", "Basophil",  "NKT", "APC", "Monocytes", "Monocytic_Lineage",
            "Neutrophils", "Macrophages", "Keratinocyte", "Endothelial", "Lymphatic")]
barplot(tab,  col = cell_type_color_df$cell_type_CellSighter_color[match(names(tab), cell_type_color_df$cell_type_CellSighter)], horiz = T, las = 1)

qs::qsave(spe, "output/SpatialExperiment.qs")

# calculate confusion matrix
confusion_matrix <- table(as.factor(spe$cell_type_CellSighter), 
                          as.factor(spe$cell_type_positive_marker))
confusion_matrix = scale(t(confusion_matrix), center = F)

# convert confusion matrix to data frame
confusion_df <- as.data.frame(confusion_matrix)

# rename columns
colnames(confusion_df) <- c("cell_type_CellSighter", "cell_type_positive_marker", "Count")
  
# create heatmap with text labels
png(file.path(output_dir, "..", "Plots", "Positive_Marker_vs_CellSighter_All.png"), res = 150, width = 1200, height = 1100)
ggplot(data = confusion_df, aes(x = cell_type_CellSighter,
                                y = cell_type_positive_marker,
                                fill =Count, label = round(Count,2))) + 
  geom_tile() + 
  scale_fill_viridis_c(begin = 0, direction = -1) +
  labs(title = "Confusion Matrix Heatmap",
       x = "Positive Marker", y = "CellSighter") + 
    geom_text(color = "white", size = 3) + 
  theme(axis.text.x = element_text(angle = 90)) + guides(fill=guide_legend(title="Scaled Overlap"))
dev.off()

# Inter annotator accuracy---------------------------------------------------------
library(reticulate)

# 1 # Before guidelines - ROI-17
labels_Annotator1 = read.csv("~/Documents/manual_annotation_safe/old/ROI-17-TL-Annotator1/cells_df.csv")
labels_Annotator3 = read.csv("~/Documents/manual_annotation_safe/old/ROI-17-TL-Annotator3/cells_df.csv")

labels_Annotator1 = labels_Annotator1[which(!labels_Annotator1$cell_type %in% c("Artifact","")),]
labels_Annotator1$cell_id = paste0(labels_Annotator1$fov, "-", labels_Annotator1$label)
labels_Annotator3 = labels_Annotator3[which(!labels_Annotator3$cell_type %in% c("Artifact","")),]
labels_Annotator3$cell_id = paste0(labels_Annotator3$fov, "-", labels_Annotator3$label)

cells = intersect(labels_Annotator1$cell_id, labels_Annotator3$cell_id)

labels_Annotator1 = labels_Annotator1[match(cells, labels_Annotator1$cell_id),]
labels_Annotator3 = labels_Annotator3[match(cells, labels_Annotator3$cell_id),]


confusion_df <- as.data.frame(table(labels_Annotator1$cell_type, labels_Annotator3$cell_type))

# rename columns
colnames(confusion_df) <- c("Annotator_1", "Annotator_2", "Count")

# create heatmap with text labels
png(file.path(output_dir, "..", "Plots", "InterAnnotator_agreement_heatmap_without_guidelines.png"), res = 150, width = 1200, height = 1100)
ggplot(data = confusion_df, aes(x = Annotator_1, y = Annotator_2, fill = log10(Count+1), label = Count)) + 
  geom_tile() + 
  scale_fill_gradient2(midpoint = 0.5, low = "royalblue4", mid = "royalblue4",high = "gold") +
  labs(title = "Inter-Annotator Agreement - with guidelines",
       x = "Annotator #1", y = "Annotator #2", fill = "Count") +
  geom_text(color = "white", size = 5)  + theme(axis.text.x = element_text(angle = 90))
dev.off()


labels_Annotator3_cell_type = factor(labels_Annotator3$cell_type)
labels_Annotator1_cell_type = factor(labels_Annotator1$cell_type, levels =  
                                      levels(labels_Annotator3_cell_type))
confusion_matrix <- confusionMatrix(labels_Annotator1_cell_type, labels_Annotator3_cell_type)

# Calculate accuracy
accuracy <- confusion_matrix$overall["Accuracy"]

# Calculate precision, recall (sensitivity), and F1-score for each class
f1_score <- mean(confusion_matrix$byClass[, "F1"], na.rm =T)

# 2 # After guidelines - ROI-17
np <- import("numpy")
labels_Annotator1 <- np$load("output/CellSighter/celltype/cell_classification/CellTypes/cells2labels/ROI-17-Annotator2.npz", allow_pickle = T)["data"]
cells_Annotator1 <- np$load("output/CellSighter/celltype/cell_classification/CellTypes/cells/ROI-17-Annotator2.npz", allow_pickle = T)["data"]
labels_Annotator3 <- np$load("output/CellSighter/celltype/cell_classification/CellTypes/cells2labels/ROI-17-Annotator3.npz", allow_pickle = T)["data"]

labels_Annotator1 <- np$load("output/CellSighter/celltype/cell_classification/CellTypes/cells2labels/ROI-17-Annotator2.npz", allow_pickle = T)["data"]
labels_Annotator3 <- np$load("output/CellSighter/celltype/cell_classification/CellTypes/cells2labels/ROI-17-Annotator3.npz", allow_pickle = T)["data"]

table(labels_Annotator1, labels_Annotator3)

labels_Annotator1 = corresponding_celltype[match(labels_Annotator1, as.numeric(names(corresponding_celltype)))]
labels_Annotator3 = corresponding_celltype[match(labels_Annotator3, as.numeric(names(corresponding_celltype)))]

confusion_df <- as.data.frame(table(labels_Annotator1, labels_Annotator3))

# rename columns
colnames(confusion_df) <- c("Annotator_1", "Annotator_2", "Count")

# create heatmap with text labels
png(file.path(output_dir, "..", "Plots", "InterAnnotator_agreement_heatmap_with_guidelines.png"), res = 150, width = 1200, height = 1100)
ggplot(data = confusion_df, aes(x = Annotator_1, y = Annotator_2, fill = log10(Count+1), label = Count)) + 
  geom_tile() + 
  scale_fill_gradient2(midpoint = 0.5, low = "royalblue4", mid = "royalblue4",high = "gold") +
  labs(title = "Inter-Annotator Agreement - with guidelines",
       x = "Annotator #1", y = "Annotator #2", fill = "Count") +
  geom_text(color = "white", size = 5)  + theme(axis.text.x = element_text(angle = 90))
dev.off()

labels_Annotator3_cell_type = factor(labels_Annotator3)
labels_Annotator1_cell_type = factor(labels_Annotator1)
confusion_matrix <- confusionMatrix(labels_Annotator1_cell_type, labels_Annotator3_cell_type)

# Calculate accuracy
accuracy2 <- confusion_matrix$overall["Accuracy"]

# Calculate precision, recall (sensitivity), and F1-score for each class
f1_score2 <- mean(confusion_matrix$byClass[, "F1"], na.rm =T)



df = data.frame("With_Guidelines" = c(FALSE, TRUE),
                "Accuracy" = c(accuracy, accuracy2),
                "F1" = c(f1_score, f1_score2))
df = df %>% pivot_longer(-With_Guidelines)
write.csv(df, file.path(output_dir, "..", "Plots", "InterAnnotator_agreement_improvement_bargraph.csv"))
png(file.path(output_dir, "..", "Plots", "InterAnnotator_agreement_improvement_bargraph.png"), res = 300, width = 1200, height = 1100)
ggplot(data = df, aes(x = With_Guidelines, y = value, fill = name)) + 
  geom_bar(stat = "identity", position = position_dodge2()) + 
  scale_fill_manual(values = c("#51426E", "#71AD7C")) +
  scale_y_continuous(name = "Accuracy",sec.axis = sec_axis( trans=~.*1, name="F1")) + 
  labs(title = "ROI-17 Inter-Annotator Agreement") +
  theme(axis.text.x = element_text(angle = 90)) + theme_classic()
dev.off()

################################################################################
# Count number of manually annotated markers + and -:
################################################################################

annot_images = c("ROI-14-Annotator4","ROI-10-Annotator2", "ROI-20-Annotator2", "ROI-22-Annotator1", "ROI-12-Annotator3","ROI-08-Annotator3", "ROI-13-Annotator3",
                 "ROI-07-Annotator1", "ROI-06-Annotator1","ROI-15-Annotator1", "ROI-15-Annotator3",  "ROI-17-Annotator3")

celltypes = unique(spe$cell_type_CellSighter)
celltypes[celltypes == "Unknown-Stroma"] = "Fibroblast"
celltypes_annotated = setNames(rep(0, length(celltypes)), celltypes)
for(celltype in celltypes){
  print(celltype)
  for(annot_dir in annot_images){
    print(annot_dir)
    file = list.files(file.path("output/CellSighter/manual_annotation_celltype/", annot_dir), paste0(celltype, ".*.csv"), full.names = T)

    if(length(file) > 0){
      pos = nrow(read.csv(file))
      celltypes_annotated[celltype] = celltypes_annotated[celltype] + pos
      print(pos)
    }
  }
}

write.csv(mat_plus, "output/CellSighter/manual_annotation_celltype/number_of_annotated_cells_celltypes+.csv")
write.csv(mat_neg, "output/CellSighter/manual_annotation_celltype/number_of_annotated_cells_celltypes-.csv")
write.csv(mat_plus + mat_neg, "output/CellSighter/manual_annotation_celltype/number_of_annotated_cells_celltypes_both.csv")

median(rowSums(mat_plus))
median(rowSums(mat_neg))
median(rowSums(mat_plus))