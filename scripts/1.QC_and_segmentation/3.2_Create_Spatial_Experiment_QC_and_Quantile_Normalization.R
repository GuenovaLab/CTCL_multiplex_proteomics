# Creation of Spatial Experiment based on segmentation
# Sample and cell filtering
# Quantile Normalization 
# authors: Pacome Prompsy
# contact: pacome.prompsy@chuv.ch
# Guenova Lab
# CHUV (Centre Hospitalier Universitaire Vaudois), Lausanne, Suisse

cat("Creating Spatial Experiment from segmentation matrices and filtering out
low quality samples and low intensity cells... \n")

# Loading packages -------------------------------------------------------------
libraries = c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "ggspavis",
  "SpatialExperiment",
  "ComplexHeatmap",
  "seriation",
  "EBImage"
)
suppressPackageStartupMessages(invisible(lapply(libraries, require, character.only = TRUE)))

setwd("/mnt/RECHERCHE/GUENOVA_LAB/Project_Multiplex_Phenotyping/Miltenyi/Muliplex_Imaging_Pipeline/")

source("scripts/MultiplexImaging_utils.R")
source("scripts/GlobalVars.R")
marker_metadata = read.csv("annotation/marker_metadata.csv")

# Directorires --------------------------------------------
output = "./output/" 

output_dir = file.path(output, "markers")
if (!dir.exists(output_dir))
  dir.create(output_dir)


# Loading the dataset ----------------------------------------------------------
whole_cell_matrix = read.csv(file.path(output, "cell_table/cell_table_size_normalized.csv.gz"))
whole_cell_normalized_matrix = read.csv(file.path(output,"cell_table/cell_table_arcsinh_transformed.csv.gz"))

samples = unique(whole_cell_matrix$fov)

names =  c(
  "cell_size",
  "area",
  "eccentricity",
  "major_axis_length",
  "minor_axis_length",
  "perimeter",
  "convex_area",
  "equivalent_diameter",
  "centroid.0",
  "centroid.1",
  "major_minor_axis_ratio",
  "perim_square_over_area",
  "major_axis_equiv_diam_ratio",
  "convex_hull_resid",
  "centroid_dif",
  "num_concavities",
  "nc_ratio",
  "label"
)
metadata = DataFrame(whole_cell_matrix[,c(names, paste0(names, "_nuclear"), "fov")])

metadata$cell_id = paste0(metadata$fov, "-", metadata$label)
whole_cell_matrix = as.matrix(whole_cell_matrix[, setdiff(colnames(whole_cell_matrix), colnames(metadata))])
whole_cell_normalized_matrix = as.matrix(whole_cell_normalized_matrix[, setdiff(colnames(whole_cell_normalized_matrix), colnames(metadata))])

whole_cell_matrix[is.na(whole_cell_matrix)] = 0
whole_cell_normalized_matrix[is.na(whole_cell_normalized_matrix)] = 0

rownames(metadata) = metadata$cell_id
rownames(whole_cell_matrix) = metadata$cell_id
rownames(whole_cell_normalized_matrix) = metadata$cell_id
whole_cell_matrix = whole_cell_matrix[,grep("original|spots|mask|smooth|fill", colnames(whole_cell_matrix), invert = TRUE)]
whole_cell_normalized_matrix = whole_cell_normalized_matrix[,grep("original|spots|mask|smooth|fill", colnames(whole_cell_normalized_matrix), invert = TRUE)]

nuclear_matrix = whole_cell_matrix[,grep("_nuclear", colnames(whole_cell_matrix))]
colnames(nuclear_matrix) = gsub("_nuclear", "",colnames(nuclear_matrix))
nuclear_matrix[nuclear_matrix < 1e-10] = 0
nuclear_normalized_matrix = whole_cell_normalized_matrix[,grep("_nuclear", colnames(whole_cell_normalized_matrix))]
colnames(nuclear_normalized_matrix) = gsub("_nuclear", "",colnames(nuclear_normalized_matrix))
nuclear_normalized_matrix[nuclear_normalized_matrix < 1e-10] = 0
nuclear_normalized_matrix[is.infinite(nuclear_normalized_matrix)] = quantile(nuclear_normalized_matrix, 0.95)

whole_cell_matrix = whole_cell_matrix[,grep("_nuclear", colnames(whole_cell_matrix), invert = TRUE)]
whole_cell_matrix = whole_cell_matrix[,-1]
whole_cell_normalized_matrix = whole_cell_normalized_matrix[,grep("_nuclear", colnames(whole_cell_normalized_matrix), invert = TRUE)]
whole_cell_normalized_matrix = whole_cell_normalized_matrix[,-1]

cytoplasm_matrix = whole_cell_matrix - nuclear_matrix
cytoplasm_matrix[cytoplasm_matrix<1e-10] = 0
cytoplasm_normalized_matrix = whole_cell_normalized_matrix - nuclear_normalized_matrix
cytoplasm_normalized_matrix[cytoplasm_normalized_matrix<1e-10] = 0


write.csv(whole_cell_matrix, file.path(output, "cell_table", "whole_cell_matrix.csv"))
write.csv(whole_cell_normalized_matrix, file.path(output, "cell_table", "whole_cell_normalized_matrix.csv"))

write.csv(nuclear_matrix, file.path(output, "cell_table", "nuclear_matrix.csv"))
write.csv(nuclear_normalized_matrix, file.path(output, "cell_table", "nuclear_normalized_matrix.csv"))

write.csv(cytoplasm_matrix, file.path(output, "cell_table", "cytoplasm_matrix.csv"))
write.csv(cytoplasm_normalized_matrix, file.path(output, "cell_table", "cytoplasm_normalized_matrix"))

# Sample Metadata
sample_metadata = read.csv(file = file.path("annotation", "Sample_metadata.csv"), sep = ";")


# Creating the SpatialExperiment object ----------------------------------------
cat("Creating SpatialExperiment...\n")

imageSources = unlist(lapply(list.dirs(args$input), function(dir) {
  list.files(dir, pattern = ".tiff", full.names = TRUE)
}))
names(imageSources) = paste0(gsub("/.*", "", gsub(".*ROI", "ROI", imageSources)), "-", gsub(".tiff", "", basename(imageSources)))

spe <- SpatialExperiment(
  assay = t(whole_cell_matrix),
  colData = metadata,
  spatialCoordsNames = c("centroid.0", "centroid.1"),
  sample_id = metadata$fov,
)
spe@int_metadata$imgSources = imageSources

# Remove samples ---------------------------------------------------------------
# Low quality samples or misdiagnosed:
#  - are too small (#46 an #51 containing only 2 frames and only captured epidermis, not enough cells in dermis)
#  - were misclassified (#30 not Darier & #49 was misclassified as PS but is AD that had change of treatment and #58 not MF but myleoid lymphoma)
#  - were the same patient (#52 is actually the same patient as #19. #52 is removed)
# From 63 patients with 343,774 cells to 57 patients with 319,743 cells.
spe = spe[,which(spe$sample_id %in% sample_metadata$ROI[which(sample_metadata$KeepForAnalysis == TRUE)])]

# Give color codes to samples
cols = setNames(discrette_colors_50, unique(spe$sample_id))
spe$sample_id_color = cols[match(spe$sample_id, names(cols))]

spe@assays@data[["counts"]] = spe@assays@data[[1]]
spe@assays@data[["normcounts"]] = t(whole_cell_normalized_matrix)

spe$log10_area_normalized_intensity = log10(colSums(counts(spe)))
spe$batch = "TMA_1"
spe$batch[spe$sample_id %in% paste0("ROI-",29:63)] = "TMA_2"



# Filter out cells with very low total counts
area_normalized_avg_intensity_th = 15000
png(
  file.path(output_dir, "Histogram_total_counts_before.png"),
  width = 1200,
  height = 800,
  res = 200
)
hist(
  spe$log10_area_normalized_intensity,
  xlab = "Area-Normalized Avg. Intensity",
  main = "Distribution of area-normalized intensities",
  breaks = 150
)
abline(
  v = log10(area_normalized_avg_intensity_th),
  lty = 2,
  col = "darkred"
)
dev.off()

spe = spe[, spe$log10_area_normalized_intensity > log10(area_normalized_avg_intensity_th)]

png(
  file.path(output_dir, "Histogram_total_counts_after.png"),
  width = 1200,
  height = 800,
  res = 200
)
hist(
  spe$log10_area_normalized_intensity,
  xlab = "Area-Normalized Avg. Intensity",
  main = "Distribution of area-normalized intensities",
  breaks = 150
)
abline(
  v = log10(area_normalized_avg_intensity_th),
  lty = 2,
  col = "darkred"
)
dev.off()

# Filter out large cells with low DAPI
spe$DAPI = (counts(spe)["DAPI", ])
spe = spe[-which(rownames(spe) == "DAPI"), ]

Area_th = 10000
DAPI_th = 10000
png(
  file.path(output_dir, "Filter_Large_Blobs_with_Low_DAPI_before.png"),
  width = 1200,
  height = 1200,
  res = 150
)
plot(
  spe$area,
  spe$DAPI,
  xlab = "Area",
  main = "DAPI vs Area",
  pch = 20,
  col = "lightblue"
)
abline(v = Area_th,
       lty = 2,
       col = "darkred")
abline(h = DAPI_th,
       lty = 2,
       col = "darkred")
dev.off()

spe = spe[, spe$DAPI > 10000 | spe$area < 10000]

png(
  file.path(output_dir, "Filter_Large_Blobs_with_Low_DAPI_after.png"),
  width = 1200,
  height = 1200,
  res = 150
)
plot(
  spe$area,
  spe$DAPI,
  xlab = "Area",
  main = "DAPI vs Area",
  pch = 20,
  col = "lightblue"
)
abline(v = Area_th,
       lty = 2,
       col = "darkred")
abline(h = DAPI_th,
       lty = 2,
       col = "darkred")
dev.off()


# Normalize using normalizeQuantile function from limma ------------------------
normcounts(spe) = limma::normalizeQuantiles(counts(spe))

# Saving  ----------------------------------------------------------------------
qs::qsave(spe, file.path(output, "SpatialExperiment.qs"))

# Plotting distribution --------------------------------------------------------
# Per sample distribution of raw intensities
intensities = as.data.frame(t(counts(spe)))
intensities$sample_id = spe$sample_id
intensities = intensities %>% tidyr::pivot_longer(-sample_id, names_to = c("marker"))
intensities$value = log2(intensities$value+1)

pdf(
  file.path(
    output_dir,
    "Marker_log2_Raw_Intensities_Distribution_All.pdf"
  ),
  width = 12,
  height = 5
)
quant0.1 = intensities %>% summarise(quant = quantile(value, 0.1))
quant0.9 = intensities %>% summarise(quant = quantile(value, 0.9))
intensities %>%
  ggplot(aes(y = value, x = sample_id, fill = sample_id)) +
  geom_violin(draw_quantiles = c(0.5)) +
  theme_classic() + theme(legend.position = 'none') +
  ylab("Normalized Average Intensities per Cell") + xlab("") +
  geom_hline(
    data = quant0.1,
    aes(yintercept = quant),
    lty = 2,
    col = "grey15"
  ) +
  geom_hline(
    data = quant0.9,
    aes(yintercept = quant),
    lty = 2,
    col = "grey15"
  ) + ggtitle("All markers")
dev.off()

pdf(
  file.path(output_dir, "Marker_log2_Raw_Intensities_Distribution.pdf"),
  width = 12,
  height = 150
)

quant0.1 = intensities %>% group_by(marker) %>% summarise(quant = quantile(value, 0.1))
quant0.9 = intensities %>% group_by(marker) %>% summarise(quant = quantile(value, 0.9))
intensities %>%
  ggplot(aes(y = value, x = sample_id, fill = sample_id)) +
  geom_violin() +
  
  theme_classic() + theme(legend.position = 'none') +
  facet_wrap( ~ marker, ncol = 1)  +
  ylab("Normalized Average Intensities per Cell") + xlab("") +
  geom_hline(
    data = quant0.1,
    aes(yintercept = quant),
    lty = 2,
    col = "grey15"
  ) +
  geom_hline(
    data = quant0.9,
    aes(yintercept = quant),
    lty = 2,
    col = "grey15"
  )
dev.off()

pdf(
  file.path(output_dir, "Number_of_cells_per_sample.pdf"),
  width = 12,
  height = 5
)
meta = as.data.frame(colData(spe))
meta %>%
  ggplot(aes(x = sample_id, fill = sample_id)) +
  geom_bar() +
  theme_classic() + theme(legend.position = 'none',
                          axis.text.x = element_text(angle = 90))
dev.off()

meta = meta %>% group_by(sample_id) %>% summarise(ncell = n())
mean_int_df = intensities %>% group_by(sample_id) %>% summarise(mean_intensities = mean(value))
meta$mean_log2_raw_intensity = mean_int_df$mean_intensities[match(meta$sample_id, mean_int_df$sample_id)]

png(
  file.path(output_dir, "Correlation_cells_per_sample.png"),
  width = 800,
  height = 800,
  res = 200
)
meta %>%
  ggplot(aes(x = ncell, y = mean_log2_raw_intensity)) +
  geom_point() +
  theme_classic() + theme(legend.position = 'none',
                          axis.text.x = element_text(angle = 90))
dev.off()

# Per sample distribution normalized intensities
intensities = as.data.frame(t(spe@assays@data$normcounts))
intensities$sample_id = spe$sample_id
intensities = intensities %>% tidyr::pivot_longer(-sample_id, names_to = c("marker"))

pdf(
  file.path(
    output_dir,
    "Marker_Normalized_Intensities_Distribution_All.pdf"
  ),
  width = 12,
  height = 5
)
quant0.1 = intensities %>% summarise(quant = quantile(value, 0.1))
quant0.9 = intensities %>% summarise(quant = quantile(value, 0.9))
intensities %>%
  ggplot(aes(y = value, x = sample_id, fill = sample_id)) +
  geom_violin(draw_quantiles = c(0.5)) +
  theme_classic() + theme(legend.position = 'none') +
  ylab("Normalized Average Intensities per Cell") + xlab("") +
  geom_hline(
    data = quant0.1,
    aes(yintercept = quant),
    lty = 2,
    col = "grey15"
  ) +
  geom_hline(
    data = quant0.9,
    aes(yintercept = quant),
    lty = 2,
    col = "grey15"
  ) + ggtitle("All markers")
dev.off()

pdf(
  file.path(
    output_dir,
    "Marker_Normalized_Intensities_Distribution.pdf"
  ),
  width = 12,
  height = 150
)
quant0.1 = intensities %>% group_by(marker) %>% summarise(quant = quantile(value, 0.1))
quant0.9 = intensities %>% group_by(marker) %>% summarise(quant = quantile(value, 0.9))
intensities %>%
  ggplot(aes(y = value, x = sample_id, fill = sample_id)) +
  geom_violin() +
  
  theme_classic() + theme(legend.position = 'none') +
  facet_wrap( ~ marker, ncol = 1)  +
  ylab("Normalized Average Intensities per Cell") + xlab("") +
  geom_hline(
    data = quant0.1,
    aes(yintercept = quant),
    lty = 2,
    col = "grey15"
  ) +
  geom_hline(
    data = quant0.9,
    aes(yintercept = quant),
    lty = 2,
    col = "grey15"
  )
dev.off()


