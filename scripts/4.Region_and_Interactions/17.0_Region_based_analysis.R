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


