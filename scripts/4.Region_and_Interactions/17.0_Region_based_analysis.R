# Analyzis of regions (results from Pixie pixel clustering)
# authors: Pacome Prompsy
# contact: pacome.prompsy@chuv.ch
# Guenova Lab
# CHUV (Centre Hospitalier Universitaire Vaudois), Lausanne, Suisse

cat("Analyzing differences in the regions of different diseases... \n")
setwd("/mnt/RECHERCHE/GUENOVA_LAB/Project_Multiplex_Phenotyping/Miltenyi/Muliplex_Imaging_Pipeline/")

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
              "ggbeeswarm",
              "rstatix")
suppressPackageStartupMessages(invisible(lapply(libraries, require, character.only = TRUE)))

source("scripts/MultiplexImaging_utils.R")
source("scripts/GlobalVars.R")

#################################################################################
# Reading in data
#################################################################################
args = list(output = "output/")
cat("Output = ", args$output, "\n")

cell_output_dir = file.path(args$output, "pixie", "region_cell_output_dir")
pixel_output_dir = file.path(args$output, "pixie", "region_pixel_output_dir")
spe = qs::qread("output/SpatialExperiment.qs")

# Mapping 
mapping = setNames(c(6,7,10,21,23,24,25,26), c(
  "Epidermis Outer",
  "Blood Vessels Outer",
  "Lymphatic Vessels",
  "Blood Vessels Inner",
  "Dermis",
  "Epidermis Inner",
  "Background",
  "Immune"
))

# Mapping 
mapping_rename = setNames(c(1,2,3,4,5,6,7,8), c(
  "Epidermis Outer",
  "Epidermis Inner",
  "Blood Vessels Outer",
  "Lymphatic Vessels",
  "Blood Vessels Inner",
  "Dermis",
  "Immune",
  "Background"
))

colors_region = data.frame(
  region = 
    c(  "Epidermis Outer",
        "Epidermis Inner",
        "Blood Vessels Outer",
        "Lymphatic Vessels",
        "Blood Vessels Inner",
        "Dermis",
        "Immune",
        "Background"),
  region_color =
    c(
      "#930805ff",
      "#4a0100ff",
      "#C23EB5",
      "#fec911ff",
      "#78108F",
      "grey85",
      "#0d522aff",
      "black")
)

#################################################################################
# Image region overlay plots
#################################################################################

# Masks
spe$region = "Background"
df = arrow::read_feather(file = file.path(cell_output_dir, "cluster_counts_size_norm.feather"))
df$cell_id = paste0(df$fov, "-", df$segmentation_label)
spe$region = df$cell_meta_cluster[match(spe$cell_id, df$cell_id)]
spe$region = names(mapping_rename)[match(spe$region, mapping_rename)]
table(spe$region)
spe = colors_scExp(spe, color_by = "region", annotCol = "region", color_df = colors_region)
qs::qsave(spe, "output/SpatialExperiment.qs")

# Plot mask
for(samp in unique(spe$sample_id)){
  pixel_overlay_mat = getImageAsMatrix(file.path(pixel_output_dir, "pixel_masks", paste0(samp,"_pixel_mask.tiff")))
  
  colors_to_pixel_cluster = setNames(colors_region$region_color[match(names(mapping), colors_region$region)], mapping)
  
  pdf(file.path(pixel_output_dir,paste0(samp,"_pixel_mask.pdf")))
  img_mat = matrix(colors_to_pixel_cluster[as.character(pixel_overlay_mat)], nrow = nrow(pixel_overlay_mat), ncol = ncol(pixel_overlay_mat))
  r = as.raster(img_mat)
  par_save = par()
  par(bg = "black", mar = c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  plot( c(0,ncol(r)), c(0,nrow(r)), axes = FALSE, type = "n", xlab = "", ylab = "")
  title(main = samp, adj = 0, line = 0,  cex.main = 0.6, col.main = "white")
  graphics::rasterImage(r, xleft = 0, ytop = nrow(r), xright = ncol(r), ybottom = 0, interpolate = FALSE)
  legend(x = "topright", inset = c(-0.35, 0.1),  cex = 0.6, title = "Legend", title.col = "white", text.col = "white", ncol = 1,
         legend = names(mapping),
         fill = as.character(colors_to_pixel_cluster), col = "white")
  par(par_save)
  dev.off()
}

#################################################################################
# Epidermis Invasion calculation
#################################################################################

meta = as.data.frame(colData(spe))
meta = meta %>% group_by(sample_id, condition) %>% 
  mutate(total_immune = length(which(!cell_type_CellSighter %in% struct_celltype))) %>% 
  dplyr::summarise(T_cell_invasion = length(
    which(
      (cell_type_CellSighter == "T_cytotoxic") &
        (region == "Epidermis Inner" | region == "Epidermis Outer") )
  ) / mean(total_immune))

png(file.path(pixel_output_dir, "Lymphocyte T invasion in Epidermis (count).png"), width = 1800, height = 1600, res = 300)
meta %>% ggplot(aes(x = condition, y = T_cell_invasion, fill = condition)) +
  geom_boxplot(outlier.colour = "white") + geom_jitter() +  scale_fill_manual(values = unique(spe$condition_color[order(spe$condition)])) +
  theme_classic() + ggtitle("Lymphocyte T invasion in Epidermis") + 
  ylab("Number of T helper/cytotoxic in epidermis") + xlab("") + 
  geom_text(aes(label = sample_id), position = position_jitterdodge())

p = meta %>% ggplot(aes(x = condition, y = type, fill = condition)) +
  geom_boxplot(outlier.colour =  "white") +
  geom_jitter(aes(color = sample_id),  colour="black",pch=21, size=2, width = 0.3) +
  scale_fill_manual(values = unique(spe$condition_color[match(levels(meta$condition), spe$condition)])) +
  theme_classic() + xlab("") + ylab(paste0(type," of ", celltype)) +
  NoLegend()  + ggtitle(paste0(type," of ", celltype)) +
  theme(axis.text.x = element_text(angle = 90))  +
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "t.test", ref.group = "HD") +
  geom_hline(yintercept = median(meta.$type), lwd = 0.35, lty = 2, col = "grey20")
print(p)
dev.off()

library(dittoSeq)
spe. = spe[,spe$region %in% c("Epidermis_Inner", "Epidermis_Outer")]
d = dittoBarPlot(spe., "cell_type_CellSighter", group.by = "sample_id",
                 color.panel = unique(spe$cell_type_CellSighter_color)[order(spe$cell_type_CellSighter)],  main = "Sample x Clusters",
                 split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)

png(file.path(output, "Plots", "Sample x CellType x Condition.png"), width = 3500, height = 1600, res = 300)
colors =  setNames(unique(spe$cell_type_CellSighter_color), unique(spe$cell_type_CellSighter))
dittoBarPlot(spe, "cell_type_CellSighter", group.by = "sample_id", color.panel = colors[levels(d$data$label)],  main = "Sample x Clusters",
             split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)
dev.off()

# Plot %
for(celltype in unique(spe$cell_type_CellSighter)){
  meta = as.data.frame(colData(spe))
  meta = meta %>% filter(!condition %in% c("TON", "THY", "LN")) %>% 
    group_by(sample_id, condition) %>% 
    summarise(invasion = 100 * length(
      which(
        (cell_type_CellSighter == celltype) &
          (region == "Epidermis_Inner" | region == "Epidermis_Outer") )
    ) / length(which((cell_type_CellSighter == "Keratinocyte")))  # / length(which((region == "Epidermis_Inner" | region == "Epidermis_Outer"))) 
    ) %>% mutate(invasion = ifelse(is.nan(invasion),0, invasion))
  
  png(file.path(pixel_output_dir, paste0(celltype, " invasion in Epidermis (percent).png")), width = 1800, height = 1600, res = 300)
  p = meta %>% ggplot(aes(x = condition, y = invasion, fill = condition)) +
    geom_boxplot(outlier.colour =  "white") +
    geom_jitter(aes(color = sample_id),  colour="black",pch=21, size=2, width = 0.3) +
    scale_fill_manual(values = unique(spe$condition_color[match(levels(meta$condition), spe$condition)])) +
    theme_classic() + ylab(paste0("Percentage of epidermis invaded by ", celltype, " (% of Keratinocytes)")) + xlab("")  +
    NoLegend()  + ggtitle(paste0("Percentage of epidermis invaded by ", celltype, " (% of Keratinocytes)")) +
    theme(axis.text.x = element_text(angle = 90))  +
    stat_compare_means(aes(label = after_stat(p.signif)),
                       method = "t.test", ref.group = "HD") +
    geom_hline(yintercept = median(meta$invasion), lwd = 0.35, lty = 2, col = "grey20") +
    ylim(c(0, quantile(meta$invasion, 0.98)))
  print(p)
  dev.off()
  
}

#################################################################################
# Area calculations
#################################################################################

meta = data.frame(sample_id = unique(spe$sample_id))
meta$total_area = 0
meta$area_epidermis = 0
meta$area_lymphatic = 0
meta$area_vasculature = 0
meta$area_dermis = 0
meta$tissue_area = 0

for(sample in unique(spe$sample_id)){
  print(sample)
  mask = tiff::readTIFF(file.path(pixel_output_dir, "pixel_masks", paste0(sample, "_pixel_mask.tiff")), as.is=TRUE)
  meta$total_area[meta$sample_id == sample] = ncol(mask) * nrow(mask) * 0.17^2
  meta$area_epidermis[meta$sample_id == sample] = length(which(mask == mapping["Epidermis Outer"] |
                                                                 mask == mapping["Epidermis Inner"])) * 0.17^2
  meta$area_lymphatic[meta$sample_id == sample] = length(which(mask == mapping["Lymphatic Vessels"])) * 0.17^2
  meta$area_immune[meta$sample_id == sample] = length(which(mask == mapping["Immune"])) * 0.17^2
  meta$area_vasculature[meta$sample_id == sample] = length(which(mask == mapping["Blood Vessels Outer"] |
                                                                   mask == mapping["Blood Vessels Inner"])) * 0.17^2
  meta$area_dermis[meta$sample_id == sample] = length(which(mask == mapping["Dermis"])) * 0.17^2
}
meta$tissue_area = meta$area_epidermis + meta$area_lymphatic + meta$area_vasculature + meta$area_dermis + meta$area_immune

meta$condition = spe$condition[match(meta$sample_id, spe$sample_id)]
write.csv(meta, file.path(pixel_output_dir, "areas_per_sample.csv"))

meta = read.csv(file.path(pixel_output_dir, "areas_per_sample.csv"))

png(file.path(pixel_output_dir, "Total Area per Condition.png"), width = 1800, height = 1600, res = 300)
meta  %>% filter(!condition %in% c("TON", "THY", "LN")) %>%
  mutate(total_area = total_area / 1000^2)  %>% 
  ggplot(aes(x = condition, y = total_area, fill = condition)) +
  geom_boxplot() + scale_fill_manual(values = unique(spe$condition_color[order(spe$condition)])) +
  theme_classic() + ggtitle("Total Area")  +
  ylab("Area (mm²)") + xlab("")
dev.off()

meta =  meta %>% filter(!condition %in% c("TON", "THY", "LN"))
meta$condition = factor(meta$condition, levels = levels(spe$condition))
meta$batch = spe$batch[match(meta$sample_id, spe$sample_id)]
meta$disease = spe$disease[match(meta$sample_id, spe$sample_id)]

list_stats = list()
for(area_name in colnames(meta)[3:9]){
  png(file.path(pixel_output_dir, paste0(area_name, " per Condition.png")), width = 1900, height = 1300, res = 300)
  meta. = meta %>% mutate(area = .data[[area_name]] / 1000^2) 
  
  p = grouped_dotplot(meta., y = "area", categ1 = "condition", categ2 = "disease",
                      color_categ1 = setNames(colCondition$condition_color, colCondition$condition))
  p = p + xlab("") + ylab(paste0("Area of ", area_name, " (mm²)")) +
    guides(fill="none", color = "none")
  print(p)
  dev.off()
  
  stats =  meta. %>%
    pairwise_t_test(area ~ disease, paired = FALSE, p.adjust.method = "none")
  stats$area = area_name
  list_stats[[area_name]] = stats
}
stats_all = do.call("rbind", list_stats)
write.csv(stats_all, file.path(pixel_output_dir, paste0("stats_area.csv")))


#################################################################################
# Density calculations
#################################################################################

### Complex cell types

# Plot and calculate the density of cell types per region ----------------------
meta = read.csv(file.path(pixel_output_dir, "areas_per_sample.csv"))
cold = as.data.frame(colData(spe))
df = cold %>% group_by(sample_id, batch, disease, cell_type_CellSighter) %>%
  dplyr::summarise(n_cells = n())
df = df %>% pivot_wider(names_from = cell_type_CellSighter, values_from = n_cells, values_fill = 0) 

meta = left_join(meta, df, by = "sample_id")
meta$condition = factor(meta$condition, levels = levels(spe$condition))
meta = meta[which(!is.na(meta$APC)),]

dir.create(file.path(cell_output_dir, paste0("Density")))
list_stats = list()

df = meta %>% filter(!condition %in% c("TON", "THY", "LN"))  %>% select(sample_id, condition, batch, disease)

for(type in unique(spe$cell_type_CellSighter)){
  meta. = meta %>% dplyr::mutate(density = .data[[type]] / tissue_area) %>%
    dplyr::filter(!condition %in% c("TON", "THY", "LN")) %>%
    dplyr::select(condition, disease, batch, density, sample_id)
  df[[type]] = 0
  df[[type]] = meta.$density
  
  pdf(file.path(cell_output_dir, paste0("Density"), paste0("Density_", type,"_within_tissue.pdf")), width = 6, height = 4,
      pointsize = 22)
  p = grouped_dotplot(meta., y = "density", categ1 = "condition", categ2 = "disease", stats = FALSE,
                      color_categ1 = setNames(colCondition$condition_color, colCondition$condition))
  p = p + xlab("") + ylab(paste0("Density of ", gsub("_"," ",type), " in tissue (cells / mm²)")) +
    guides(fill="none", color = "none")
  print(p)
  dev.off()
  
  stats = ttest_results <- meta. %>%
    pairwise_t_test(density ~ disease, paired = FALSE, p.adjust.method = "none")
  stats$celltype = type
  list_stats[[type]] = stats
}

stats_all = do.call("rbind", list_stats)
write.csv(stats_all, file.path(cell_output_dir, paste0("Density"), paste0("stats.csv")))


# Cor Plot densities
dir.create(file.path(cell_output_dir, paste0("Density_Cor")))

pdf(file.path(cell_output_dir, paste0("Density_Cor"), "Density_Correlations.pdf"), height = 4, width = 4)
for(celltype1 in c("APC", "T_helper")){
  print(celltype1)
  for(celltype2 in c("T_cytotoxic", "NKT")){
    if(celltype1 != celltype2){
      
      COR = cor.test(df[[celltype1]], df[[celltype2]])
      R = summary(lm(df[[celltype1]] ~ df[[celltype2]]))$r.squared
      p = df %>%
        ggplot(aes(x = .data[[celltype1]], y = .data[[celltype2]])) +
        geom_point(aes(fill = condition, shape = batch), size = 3) + 
        scale_fill_manual(values = setNames(colCondition$condition_color, colCondition$condition)) +
        scale_shape_manual(values = setNames(c(21, 24), c("TMA1", "TMA2"))) +
        guides(shape = guide_legend("TMA")) +
        geom_smooth(method = "lm", color = "black") + 
        theme_classic() +
        theme(
          axis.text.x = element_text(angle = 90), strip.background = element_blank(),
          panel.background = element_blank(),
          panel.ontop = F,
          panel.spacing = unit(0.75, "lines"), strip.text.x = element_text(size = 12.5))
      p = p +  annotate("text", x = quantile(df[[celltype1]],0.85), y = max(df[[celltype2]]), size = 3.5, 
                    label = paste0("Pearson Cor Test, p = ", round(COR$p.value,5))) +
        annotate("text", x = quantile(df[[celltype1]],0.65), y = max(df[[celltype2]]) /2 , size = 3.5, 
                    label = paste0("Pearson Cor = ", round(COR$estimate,5)))
      
      print(p + NoLegend())
    }
  }
}
dev.off()

### Simplified cell types

# Plot and calculate the density of cell types per region ----------------------
meta = read.csv(file.path(pixel_output_dir, "areas_per_sample.csv"))
cold = as.data.frame(colData(spe))
df = cold %>% group_by(sample_id, batch, disease, cell_type_simplified) %>%
  summarise(n_cells = n())
df = df %>% pivot_wider(names_from = cell_type_simplified, values_from = n_cells, values_fill = 0) 

meta = left_join(meta, df, by = "sample_id")
meta$condition = factor(meta$condition, levels = levels(spe$condition))
meta = meta[which(!is.na(meta$Dermis)),]

dir.create(file.path(cell_output_dir, paste0("Density_simplified")))
list_stats = list()
for(type in unique(spe$cell_type_simplified)){
  meta. = meta %>% mutate(density = .data[[type]] / tissue_area) %>%
    filter(!condition %in% c("TON", "THY", "LN")) %>%
    select(condition, disease, batch, density)
  
  png(file.path(cell_output_dir, paste0("Density_simplified"), paste0("Density_", type,"_within_tissue.png")), width = 1900, height = 1300, res = 300)
  p = grouped_dotplot(meta., y = "density", categ1 = "condition", categ2 = "disease",
                      color_categ1 = setNames(colCondition$condition_color, colCondition$condition))
  p = p + xlab("") + ylab(paste0("Density of ", gsub("_"," ",type), " in tissue (cells / mm²)")) +
    guides(fill="none", color = "none")
  print(p)
  dev.off()
  
  stats = ttest_results <- meta. %>%
    pairwise_t_test(density ~ disease, paired = FALSE, p.adjust.method = "none")
  stats$celltype = type
  list_stats[[type]] = stats
}
stats_all = do.call("rbind", list_stats)
write.csv(stats_all, file.path(cell_output_dir, paste0("Density_simplified"), paste0("stats.csv")))

# Total cell density
meta. = meta
meta.$all = meta.$Dermis + meta.$Epidermis + meta.$Lymphocyte + meta.$Myeloid + meta.$Vessels
meta. = meta. %>% mutate(density = .data[["all"]] / tissue_area)

png(file.path(cell_output_dir, paste0("Density"), paste0("Density_", "all","_within_tissue.png")), width = 1900, height = 1300, res = 300)
meta. = meta. %>% mutate(density = .data[["all"]] / tissue_area) %>%
  filter(!condition %in% c("TON", "THY", "LN")) %>%
  select(condition, disease, batch, density)

png(file.path(cell_output_dir, paste0("Density"), paste0("Density_all_within_tissue.png")), width = 1900, height = 1300, res = 300)
p = grouped_dotplot(meta., y = "density", categ1 = "condition", categ2 = "disease",
                    color_categ1 = setNames(colCondition$condition_color, colCondition$condition))
p = p + xlab("") + ylab(paste0("Density of cells in tissue (cells / mm²)")) +
  guides(fill="none", color = "none")
print(p)
dev.off()

#################################################################################
# Basal cell determination
#################################################################################

# Read the contour to find basal layer -----------------------------------------
spe$basal = FALSE
for(sample in unique(spe$sample_id)[1:length(unique(spe$sample_id))]){
  spe. = spe[,spe$sample_id == sample]
  print(sample)
  
  contour = tiff::readTIFF(paste0("output/Regions/", sample, "_contour.tiff"), as.is = TRUE)
  
  whole_cell = tiff::readTIFF(paste0("output/segmentation//", sample, "_whole_cell.tiff"), as.is = TRUE)
  whole_cell_contour_exclude  =  whole_cell
  whole_cell_contour_exclude[which(contour == 0, arr.ind = TRUE)] = 0
  cells = (unique(Matrix::Matrix(whole_cell_contour_exclude, sparse = TRUE)@x))
  
  whole_cell. = Matrix::Matrix(whole_cell, sparse = TRUE)
  keratinocytes = intersect(cells, as.numeric(gsub(".*-","", spe.$cell_id[spe.$cell_type_CellSighter == "Keratinocyte"])))
  whole_cell.@x[!whole_cell.@x %in% keratinocytes] = 0
  
  spe.$basal = FALSE
  
  print(length(keratinocytes))
  n = 0
  ratios = rep(0, length(keratinocytes))
  for(i in keratinocytes){
    if (n %% 50 == 0) print(n)
    n = n+1
    coord = which(whole_cell. == i, arr.ind = TRUE)
    ratio = length(which(contour[coord] != 0)) / length(which(contour[coord] > -1))
    ratios[n] = ratio
    if(ratio > 0.25)
      spe.$basal[spe.$cell_id == paste0(sample, "-", i)] = TRUE
  }
  
  # Plotting
  cell_overlay_mat = getImageAsMatrix(file.path("output", "segmentation", paste0(sample,"_whole_cell.tiff")))
  celltype_img_mat = get_metadata_image_overlay(spe., cell_overlay_mat, sample = sample, metadata = "basal",
                                                levels = c(FALSE, TRUE))
  
  tiff::writeTIFF(as.matrix(celltype_img_mat),
                  file.path("output/Regions",  paste0(sample,"_basal.tiff")),
                  bits.per.sample = 8)
  
  spe$basal[match(spe.$cell_id, spe$cell_id)] = spe.$basal
}

dim(spe)
qs::qsave(spe, "output/SpatialExperiment.qs")


## Basal cell analysis
spe = qs::qread("output/SpatialExperiment.qs")


#################################################################################
# Comparing basal cells disease vs  healthy, by disease 
################################################################################# 
spe. = spe[setdiff(rownames(spe), markers_to_remove), !spe$condition %in% c("THY", "TON", "LN")]

pdf(file.path(cell_output_dir, paste0("DA_basal_disease_vs_HD_per_condition.pdf")), width = 6, height = 5)
for(cond in setdiff(spe.$condition, "HD")){
  
  spe_healthy = spe.[,spe.$condition == "HD" & spe.$basal == TRUE]
  healthy_cells = c()
  # subsample 
  for(healthy in unique(spe_healthy$sample_id)){
    cells = spe_healthy$cell_id[spe_healthy$sample_id == healthy]
    healthy_cells = c(healthy_cells, sample(cells, min(50, length(cells)), replace = FALSE))
    print( min(50, length(cells)))
  }
  
  spe_disease = spe.[,spe.$condition %in% cond  & spe.$basal == TRUE]
  disease_cells = c()
  # subsample 
  for(disease in unique(spe_disease$sample_id)){
    cells = spe_disease$cell_id[spe_disease$sample_id == disease]
    disease_cells = c(disease_cells, sample(cells, min(50, length(cells)), replace = FALSE))
    print(min(50, length(cells)))
  }
  
  
  if(length(healthy_cells) > 10 & length(disease_cells) > 10){
    
    mat1 = spe.@assays@data$CellSighter_marker_mat[,healthy_cells]
    mat2 = spe.@assays@data$CellSighter_marker_mat[,disease_cells]
    
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
      ggtitle(paste0(cond, "  vs HD - Basal - volcano plot")) + 
      geom_hline(yintercept = 2, lwd = 0.35, lty = 2, col = "grey20") +
      geom_vline(xintercept = -1, lwd = 0.35, lty = 2, col = "grey20") +
      geom_vline(xintercept = 1,  lwd = 0.35, lty = 2, col = "grey20") +
      ggrepel::geom_text_repel()
    print(p)
    
  }    
  
}
dev.off()

#################################################################################
# Comparing basal cells disease vs  healthy, by disease 
################################################################################# 
spe. = spe[setdiff(rownames(spe), markers_to_remove), !spe$condition %in% c("THY", "TON", "LN")]

pdf(file.path(cell_output_dir, paste0("DA_basal_disease_vs_HD_per_condition.pdf")), width = 6, height = 5)
for(cond in setdiff(spe.$condition, "HD")){
  
  spe_healthy = spe.[,spe.$condition == "HD" & spe.$basal == TRUE]
  healthy_cells = c()
  # subsample 
  for(healthy in unique(spe_healthy$sample_id)){
    cells = spe_healthy$cell_id[spe_healthy$sample_id == healthy]
    healthy_cells = c(healthy_cells, sample(cells, min(50, length(cells)), replace = FALSE))
    print( min(50, length(cells)))
  }
  
  spe_disease = spe.[,spe.$condition %in% cond  & spe.$basal == TRUE]
  disease_cells = c()
  # subsample 
  for(disease in unique(spe_disease$sample_id)){
    cells = spe_disease$cell_id[spe_disease$sample_id == disease]
    disease_cells = c(disease_cells, sample(cells, min(50, length(cells)), replace = FALSE))
    print(min(50, length(cells)))
  }
  
  
  if(length(healthy_cells) > 10 & length(disease_cells) > 10){
    
    mat1 = spe.@assays@data$CellSighter_marker_mat[,healthy_cells]
    mat2 = spe.@assays@data$CellSighter_marker_mat[,disease_cells]
    
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
      ggtitle(paste0(cond, "  vs HD - Basal - volcano plot")) + 
      geom_hline(yintercept = 2, lwd = 0.35, lty = 2, col = "grey20") +
      geom_vline(xintercept = -1, lwd = 0.35, lty = 2, col = "grey20") +
      geom_vline(xintercept = 1,  lwd = 0.35, lty = 2, col = "grey20") +
      ggrepel::geom_text_repel()
    print(p)
    
  }    
  
}
dev.off()


#################################################################################
# Interaction Maps of basal layer, condition by condition
################################################################################# 

general_types = list(
  "Myeloid" = APC_types,
  "T_helper1" = "T_helper",
  "Toxic" = c("T_cytotoxic", "NKT"),
  "T_regulatory1" = "T_regulatory"
)

for(cond in setdiff(spe$condition, c("LN","THY","TON")) ){
  
  spe. = spe[, (spe$cell_type_CellSighter %in% unlist(general_types) | spe$basal == TRUE) & spe$condition == cond]
  spe.$General_Types = gsub(".$","",names(unlist(general_types)))[match(spe.$cell_type_CellSighter, unlist(general_types))]
  interaction_list = list()
  
  for(samp in unique(spe.$sample_id)){
    spe.. = spe.[,spe.$sample_id == samp]
    
    cells1 = spe..$cell_id[which(spe..$basal == TRUE)]
    
    if(length(cells1) > 1){
      mat = spe..@assays@data$interaction[cells1,]
      mat[mat>0]=1
      types = as.factor(spe..$General_Types)
      
      types_tab = data.frame(
        sample_id = samp,
        disease = unique(spe..$condition),
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
  
  df = do.call("rbind", interaction_list)
  if(!is.null(df)){
    
    # for(celltype2 in unique(df$type)){
    #   contingency = data.frame("Type" = c(df$sum[df$type == celltype2 & df$disease == "MF"], df$sum[df$type == celltype2 & df$disease == "SS"]),
    #                       "NonType" = c( sum(df$sum[df$type != celltype2 & df$disease == "MF"]), sum(df$sum[df$type != celltype2 & df$disease == "SS"])),
    #                       row.names = c("MF", "SS"))
    #   colnames(contingency) = paste0(c("", "Non-"),celltype2)
    #   print(paste0(celltype1," - Number of interactions ",celltype2, " p-value = ", fisher.test(contingency)$p.value))
    #   print(contingency)
    #   write.table(paste0(celltype1," - Number of interactions ",celltype2, " p-value = ", fisher.test(contingency)$p.value), file.path(output, "Percent_in_contact", "contingency_tab_pvalues.csv"), append = TRUE)
    #   write.table(contingency, file.path(output, "Percent_in_contact", "contingency_tab_pvalues.csv"), append = TRUE)
    # }
    
    df = df %>% mutate(weighted_sum = sum / total)
    df = df %>% group_by(sample_id) %>%
      mutate(percent = 100 * weighted_sum / sum(weighted_sum))
    df = df %>% group_by(disease, type) %>%
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
    
    png(file.path(cell_output_dir, paste0("Interaction_map_basal_",cond,".png")), width = 1300, height = 1200, res = 150)
    
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
            annotate("text", x = -1, y = 0, label = "Basal", size = 6) +
            coord_polar(theta="y") +
            xlim(c(-1, 4)) +
            theme_void() +
            theme(legend.position = "none"))
    dev.off()
    
    
  }
}

#################################################################################
# Proliferation in epidermis and basal layer
#################################################################################

meta = as.data.frame(colData(spe[,spe$batch == "batch_2"]))
meta$Ki67 = spe@assays@data$CellSighter_marker_mat["Ki67",match(meta$cell_id, colnames(spe))]

pdf(file.path(cell_output_dir, paste0("Proliferation_per_celltype.pdf")), width = 6, height = 5)
for(celltype in c(unique(spe$cell_type_CellSighter, "Basal"))){
  if(celltype == "Basal"){
    meta. = meta %>% filter(basal == TRUE) %>% 
      group_by(sample_id, condition) %>% summarise(percent_prolif = 100 * mean(Ki67))
  } else{
    meta. = meta %>% filter(cell_type_CellSighter == celltype) %>% 
      group_by(sample_id, condition) %>% summarise(percent_prolif = 100 * mean(Ki67))
  }
  
  p = meta. %>%
    ggplot(aes(x = condition, y = percent_prolif, fill = condition)) +
    geom_boxplot(outlier.colour = "white") +
    geom_jitter(aes(color = sample_id),  colour="black",pch=21, size=2, width = 0.3) +
    scale_fill_manual(values =setNames(unique(spe$condition_color), unique(spe$condition))) +
    theme_classic() + xlab("") + ylab("Ki67+ (%)") + 
    ggtitle(paste0(celltype, "  Proliferation")) + 
    ggpubr::stat_compare_means(aes(label = after_stat(p.signif)),
                               method = "t.test", ref.group = "Inflammatory")
  print(p)
}
dev.off()


meta. = meta %>% filter(cell_type_CellSighter == "Keratinocyte" & basal == F) %>% 
  group_by(sample_id, condition) %>% dplyr::summarise(percent_prolif = 100 * mean(Ki67))

meta.$disease = as.character(meta.$condition)
meta.$disease[meta.$condition %in% c("MF", "SS")] = "Lymphoma"
meta.$disease[meta.$condition %in% c("LP", "DAR", "PS","AD")] = "Inflammatory"
meta.$disease = factor(meta.$disease, levels = c("HD", "Inflammatory", "Lymphoma", "THY", "TON", "LN"))

png(file.path("output/CellSighter/Plots/Tumor_vs_Rest/", paste0("Proliferation_in_Keratinocytes.png")), width = 1200, height = 1000, res = 200)
p = meta. %>%
  ggplot(aes(x = disease, y = percent_prolif, fill = disease)) +
  geom_boxplot(outlier.colour = "white") +
  geom_jitter(aes(color = sample_id),  colour="black",pch=21, size=2, width = 0.3) +
  scale_fill_manual(values = c(setNames(c("#9E2069", "#B3A632FA"), c("Lymphoma", "Inflammatory")),
                               setNames(unique(spe$condition_color), unique(spe$condition)))) +
  theme_classic() + xlab("") + ylab("Ki67+ (%)") + 
  ggtitle(paste0(celltype, "  Proliferation")) +  theme(axis.text.x = element_text(angle = 90)) +
  ggpubr::stat_compare_means(aes(label = after_stat(p.signif)),
                             method = "t.test", ref.group = "HD")
print(p)
dev.off()



meta. = meta %>% filter(basal == TRUE) %>% 
  group_by(sample_id, condition) %>% dplyr::summarise(percent_prolif = 100 * mean(Ki67))

meta.$disease = as.character(meta.$condition)
meta.$disease[meta.$condition %in% c("MF", "SS")] = "Lymphoma"
meta.$disease[meta.$condition %in% c("LP", "DAR", "PS","AD")] = "Inflammatory"
meta.$disease = factor(meta.$disease, levels = c("HD", "Inflammatory", "Lymphoma", "THY", "TON", "LN"))

png(file.path("output/CellSighter/Plots/Tumor_vs_Rest/", paste0("Proliferation_in_Basal.png")), width = 1200, height = 1000, res = 200)
p = meta. %>%
  ggplot(aes(x = disease, y = percent_prolif, fill = disease)) +
  geom_boxplot(outlier.colour = "white") +
  geom_jitter(aes(color = sample_id),  colour="black",pch=21, size=2, width = 0.3) +
  scale_fill_manual(values = c(setNames(c("#9E2069", "#B3A632FA"), c("Lymphoma", "Inflammatory")),
                               setNames(unique(spe$condition_color), unique(spe$condition)))) +
  theme_classic() + xlab("") + ylab("Ki67+ (%)") + 
  ggtitle(paste0("Basal", "  Proliferation")) +  theme(axis.text.x = element_text(angle = 90)) +
  ggpubr::stat_compare_means(aes(label = after_stat(p.signif)),
                             method = "t.test", ref.group = "HD")
print(p)
dev.off()



