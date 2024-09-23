qn <- function(.data){
  data_sort <- apply(.data, 2, sort)
  row_means <- rowMeans(data_sort)
  data_sort <- matrix(row_means, 
                      nrow = nrow(data_sort), 
                      ncol = ncol(data_sort), 
                      byrow = TRUE
  )
  index_rank <- apply(.data, 2, order)
  normalized_data <- matrix(nrow = nrow(.data), ncol = ncol(.data))
  for(i in 1:ncol(.data)){
    normalized_data[,i] <- data_sort[index_rank[,i], i]
  }
  return(normalized_data)
}


circle_layout <- function(num_nodes = 10, diameter = 2){
  
  # Calculate the angle increment for equally spacing nodes along the circle
  angle_increment <- 2 * pi / (num_nodes - 1)
  
  # Initialize an empty matrix to store the coordinates
  layout_matrix <- matrix(nrow = num_nodes, ncol = 2)
  
  # Central node (0, 0)
  layout_matrix[1, ] <- c(0, 0)
  
  for (i in 1:(num_nodes - 1)) {
    angle <- i * angle_increment
    x <- diameter / 2 * cos(angle)
    y <- diameter / 2 * sin(angle)
    layout_matrix[i + 1, ] <- c(x, y)
  }
  
  return(layout_matrix)
}

changeRange = function (v, newmin = 0, newmax = 1) 
   {
         oldmin <- min(v, na.rm = TRUE)
         oldmax <- max(v, na.rm = TRUE)
         newmin + ((newmax - newmin) * (v - oldmin)/(oldmax - oldmin))
    }

interaction_doughnut_plot <- function(spe, celltype = "T_helper", group_by = "cell_type_CellSighter",
                                      enriched_protein = NULL, plot_by = "disease",
                                      colors = NULL,  weighted = c("by_celltype_weighted", "by_celltype", "absolute")[1]){
  prepare_df = .interaction_prepare(spe, celltype = celltype, group_by = group_by, enriched_protein = enriched_protein)
  .interaction_doughnut_plot(spe, prepare_df, celltype = celltype, group_by = group_by, plot_by = plot_by, colors = colors, weighted = weighted)
}

.interaction_prepare <- function(spe, celltype = "T_helper", group_by = "cell_type_CellSighter",
                                 enriched_protein = NULL){
  interaction_list = list()
  
  if(!is.null(enriched_protein)){
    all_interactions = qs::qread(file.path("output/cell_interactions/enriched_interactions/", celltype, paste0(celltype, "_interactions.qs")))
    interactions = all_interactions %>% filter(marker == enriched_protein)
    interactions = interactions %>% mutate(n_interactions_enriched = n_interactions * percent_interaction / 100)
  }
  
  for(samp in unique(spe$sample_id)){
    spe. = spe[,spe$sample_id == samp]
    cells1 = spe.$cell_id[which(spe.$cell_type_CellSighter == celltype)]
    
    if(length(cells1) > 1){
      print(samp)
      cells2 = spe.$cell_id[which(spe.$sample_id == samp)]
      
      mat = spe.@metadata$interaction[cells2,cells1]
      mat[mat>0]=1
      mat = as.data.frame(mat)
      types = as.factor(spe.[[group_by]])  
      mat$type = types
      mat = mat %>% group_by(type) %>% summarise_all(sum)
      n_interactions = rowSums(mat[,-1])
      

      types_tab = data.frame(
        sample_id = samp,
        disease = unique(spe.$disease),
        batch = unique(spe.$batch),
        type = levels(types),
        n_interactions = n_interactions,
        total_cells = as.numeric(table(types))
      )
      
      if(!is.null(enriched_protein)){
        mat = interactions %>% filter(sample_id == samp)
        n_interactions_enriched = mat$n_interactions_enriched
        types_tab = types_tab[match(mat$celltype2, types_tab$type),]
        types_tab$n_interactions_enriched = n_interactions_enriched
        
      }
      
      interaction_list[[samp]] = types_tab
    }
  }
  
  df = do.call("rbind", interaction_list)
  return(df)
}


.interaction_doughnut_plot <- function(spe, prepare_df, celltype =  "T_helper", group_by = "cell_type_CellSighter", plot_by = "disease",
                                       colors = setNames(unique(spe$cell_type_CellSighter_color), unique(spe$cell_type_CellSighter)), 
                                       weighted = c("pct_weighted_celltype", "pct", "absolute")[1]
){
  if(!is.null(prepare_df)){
    
    
    if(weighted == "pct_weighted_celltype"){
      df = prepare_df %>% mutate(weighted_sum = n_interactions / total_cells)
      df = df %>% group_by(sample_id) %>%
        mutate(percent = 100 * weighted_sum / sum(weighted_sum))
    } else if(weighted == "pct"){
      df = prepare_df %>% mutate(sum = n_interactions)
      df = df %>% group_by(sample_id) %>%
        mutate(percent = 100 * sum / sum(sum))
    } else if(weighted == "absolute"){
      df = prepare_df %>% group_by(sample_id) %>%
        mutate(percent = n_interactions)
    }
    df = df %>% group_by(disease, type) %>%
      summarise(percent = mean(percent))
print(head(df))
    
for(cond in unique(df[[plot_by]])){
      data = df[which(df[[plot_by]] == cond),]
      
      # Compute percentages
      if(weighted == "absolute"){
        data$fraction <- data$percent
      } else{
        data$fraction <- data$percent / sum(data$percent)
      }
      
      # Compute the cumulative percentages (top of each rectangle)
      data$ymax <- cumsum(data$fraction)
      
      # Compute the bottom of each rectangle
      data$ymin <- c(0, head(data$ymax, n=-1))
      
      # Compute label position
      data$labelPosition <- (data$ymax + data$ymin) / 2
      
      # Compute a good label
      # Compute percentages
      if(weighted == "absolute"){
        data$label <- paste0(data$type, "\n (",round(data$fraction,1) ,")")
      } else{
        data$label <- paste0(data$type, "\n (", round(100*data$fraction,1) ,"%)")
        
      }
      
      if(is.null(colors)){
        set.seed(47)
        colors = setNames(sample(discrette_colors_50, length(unique(spe[[group_by]]))), unique(spe[[group_by]]))
      }
      
      # Make the plot
      p = (ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=type)) +
             geom_rect() +
             geom_text( x=2, aes(y=labelPosition, label=label, color=type), size=6) + # x here controls label position (inner / outer)
             scale_fill_manual(values = colors) +
             scale_color_manual(values = colors) +
             annotate("text", x = -1, y = 0, label = celltype, size = 6) +
             coord_polar(theta="y") +
             xlim(c(-1, 4)) +
             theme_void() +
             theme(legend.position = "none") +
             ggtitle(cond)) 
      print(p)
}

  }
}


grouped_dotplot <- function(metadata, y, categ1, categ2 = NULL,
                            color_categ1, ref.group = "HD",
                            shape_by = "batch", add_violin = FALSE, 
                            stats = TRUE, paired = FALSE, beeswarm = FALSE,
                            size = 3, add_points = TRUE){
  
  stopifnot(!is.null(metadata), !is.null(categ1),  !is.null(y),
            categ1 %in% colnames(metadata), y %in% colnames(metadata),
            ((stats == TRUE && ref.group %in% metadata[[categ1]]) | (stats == F))
  )
  metadata = metadata %>% ungroup
  
  p = metadata   %>%
    ggplot(aes(x = .data[[categ1]], y = .data[[y]], fill = .data[[categ1]], color = .data[[categ1]])) +
    scale_fill_manual(values = color_categ1) +
    scale_color_manual(values = color_categ1) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 90), strip.background = element_blank(),
      panel.background = element_blank(),
      panel.ontop = F,
      panel.spacing = unit(0.75, "lines"), strip.text.x = element_text(size = 12.5))
  
  if(add_violin == TRUE){
    p = p + geom_violin(color = "black", draw_quantiles = 0.5)
  } else {
    p = p + stat_summary(fun.data = "mean_cl_boot",  linewidth = 1, size = 2, pch = 18) 
  }
  
  if(add_points == TRUE){
    
    if(beeswarm == TRUE){
      if( !is.null(shape_by)){
        p = p + geom_beeswarm(aes(shape = .data[[shape_by]]), color = "black", size = size)
        if(shape_by == "batch") p = p + scale_shape_manual(values = setNames(c(21, 24), c("TMA1", "TMA2"))) +
            guides(shape = guide_legend("TMA"))
      }  else{
        p = p + geom_beeswarm( color = "black", size = size, pch = 21)
      }
    } else{
      if(!is.null(shape_by)){
        p = p + geom_quasirandom(aes(shape = .data[[shape_by]]), color = "black", size = size)
        
        if(shape_by == "batch") p = p + scale_shape_manual(values = setNames(c(21, 24), c("TMA1", "TMA2"))) +
            guides(shape = guide_legend("TMA"))  
        
      }  else{
        p = p + geom_quasirandom(color = "black", size = size, pch = 21)
      }
    }
  
    
  }
  
  if(!is.null(categ2)){
    p = p + facet_grid(formula(paste0("~",categ2)), scales = "free_x",
                       space = "free_x", switch="both") 
    
    # facet_grid(formula(paste0("~",categ2)), scales = "free_x", axes = "all_x", 
    #            axis.labels = "margins", space = "free_x", switch="both") 
  }
  
  if(stats == TRUE){
    if(is.null(categ2)){
      p = p + stat_compare_means(aes(label = after_stat(p.format)),
                                 method = "t.test", ref.group = ref.group, paired = paired) 
    } else{
      ttest_results <- metadata %>% ungroup %>%
        pairwise_t_test(formula(paste0(y, " ~ ", categ1)), 
                        paired = paired, p.adjust.method = "none", ref.group = ref.group)
      
      colnames(ttest_results)[3] = categ1
      ttest_results = ttest_results %>% select(.data[[categ1]], p)
      metadata$p  = ttest_results$p[match(metadata[[categ1]], ttest_results[[categ1]])]
      metadata$p[is.na(metadata$p)] = ""
      
      p = p + geom_text(y = max(metadata[[y]]), size = 3.5, check_overlap = T, color = "black", 
                        label = metadata$p, hjust = 0.5, vjust = -0.5) 
    }
    
  }
  
  return(p)
}



plot_summary_Seurat <- function(object, reduction, cluster_name){
  print(DimPlot(object, reduction = reduction, group.by = "study"))
  print(DimPlot(object, reduction = reduction,  group.by = "condition"))
  print(DimPlot(object, reduction = reduction, group.by = "sample_id"))
  print(DimPlot(object, reduction = reduction, group.by = "patient"))
  print(DimPlot(object, reduction = reduction, group.by = cluster_name))
  
  tab = table(object$seurat_clusters)
  
  print(dittoSeq::dittoBarPlot(object[,object$seurat_clusters %in% names(tab[tab>300])], "study", group.by = cluster_name,
                               main = "Cluster x Study"))
  
  for(study in unique(object$study)){
    object$is_study = FALSE
    object$is_study[object$study == study] = TRUE
    
    print(DimPlot(object, reduction = reduction, group.by = "is_study", order = TRUE, cols = c("grey", "red")) + ggtitle(study))
    
  }
  
  markers = c(
    "FABP4", "APOLD1", "EMCN", "PDPN", # Endothelial
    "COL3A1", "DCN", "DPT", "LUM", "PDGFRA", # Fibroblast
    "KRT16", "KRT14", "KRT6A", # Epithelial
    "PTPRC", # Immune
    "CD3E", "CD4", "IL7R", # T helper cells
    "FOXP3", "IL2RA", # T regs
    "CD8A", "CD28", # T cytotoxic
    "MS4A1", "CD79A",# B cells
    "NCAM1",  "B3GAT1", # NK
    "HLA-DRA", "HLA-DRB1", # Myeloid
    "IL3RA", # pDC
    "CD68",  "ITGAM",  # Macrophages
    "FCGR3A" ,# Neutrophils
    "CD14" # Monocytes
  )
  
  for(mark in markers){
    print(FeaturePlot(object, reduction = reduction, features = mark, order = TRUE, cols = c("grey", "red")))
  }
  
}


scaleMat <- function(x, scale = c("row", "column", "both")){
  if (scale == "row") {
    x <- sweep(x, 1L, rowMeans(x, na.rm = T), check.margin = FALSE)
    sx <- apply(x, 1L, sd, na.rm = T)
    x <- sweep(x, 1L, sx, `/`, check.margin = FALSE)
  }
  else if (scale == "column") {
    x <- sweep(x, 2L, colMeans(x, na.rm = T), check.margin = FALSE)
    sx <- apply(x, 2L, sd, na.rm = T)
    x <- sweep(x, 2L, sx, `/`, check.margin = FALSE)
  } else if (scale == "both") {
    colM = colMeans(x, na.rm = T)
    x <- sweep(x, 1L, rowMeans(x, na.rm = T), check.margin = FALSE)
    x <- sweep(x, 2L,  colMeans(x, na.rm = T), check.margin = FALSE)
    
    sx <- apply(x, 2L, sd, na.rm = T)
    x <- sweep(x, 2L, sx, `/`, check.margin = FALSE)
  }
  
}


get_ellipse <- function(spe, cells){
  
  points_matrix = as.matrix(spe@int_colData$spatialCoords[cells,])
  best_ellipse <- ellipsoidhull(points_matrix )
  
  mean_values <- colMeans(points_matrix)
  eig <- eigen(best_ellipse$cov)
  angle <- atan2(eig$vectors[2, 1], eig$vectors[1, 1])
  major_axis <- sqrt(qchisq(0.95, 2)) * sqrt(eig$values[1])
  minor_axis <- sqrt(qchisq(0.95, 2)) * sqrt(eig$values[2])
  
  best_ellipse = predict(best_ellipse)
  
  ellipse_area <- pi * major_axis * minor_axis
  ellipse_perimeter <- pi * (3 * (major_axis + minor_axis) - sqrt((3 * major_axis + minor_axis) * (major_axis + 3 * minor_axis)))
  normalized_diameter= ellipse_area / ellipse_perimeter
  
  out = list(ellipse_points = best_ellipse, major_axis = major_axis, minor_axis = minor_axis, 
             ellipse_area = ellipse_area, ellipse_perimeter = ellipse_perimeter,
             normalized_diameter = normalized_diameter)
  return(out)
}



get_convex_hull <- function(spe, cells){
  
  points_matrix = as.matrix(spe@int_colData$spatialCoords[cells,])
  
  # Calculate the convex hull of your points
  convex_hull_indices <- chull(points_matrix)
  
  # Extract the points that form the convex hull
  convex_hull_points <- points_matrix[convex_hull_indices, ]
  convex_hull_points_sp <- SpatialPoints(convex_hull_points)
  
  # Create a buffered polygon around the convex hull
  buffered_polygon <- gBuffer(convex_hull_points_sp, byid = TRUE)
  buffered_polygon =  SpatialPoints(buffered_polygon)
  convex_hull_indices <- chull(buffered_polygon@coords)
  convex_hull_points <- buffered_polygon@coords[convex_hull_indices, ]
  
  # centroid <- colMeans(convex_hull_points)convex_hull_points
  # transformation_matrix <- matrix(c(scaling_factor, 0, 0, scaling_factor), ncol = 2)
  # translation_vector <- centroid - centroid * scaling_factor
  
  # Apply the affine transformation to each vertex to create the larger polygon
  # larger_polygon <- convex_hull_points %*% t(transformation_matrix) + matrix(1,  nrow  = nrow(convex_hull_points)) %*% translation_vector 
  
  xlim_min <- min(convex_hull_points[,1])
  xlim_max <- max(convex_hull_points[,1])
  ylim_min <- min(convex_hull_points[,2])
  ylim_max <- max(convex_hull_points[,2])
  
  # Plot the original points, convex hull, and enlarged polygon
  plot(points_matrix[,1], points_matrix[,2], pch = 19, col = "blue", 
       xlab = "X", ylab = "Y", main = "Enlarged Convex Hull",
       xlim = c(xlim_min-20, xlim_max+20), ylim = c(ylim_min-20, ylim_max+20))
  polygon(convex_hull_points, col = "transparent", border = "red")
  # polygon(larger_polygon, col = "transparent", border = "red")
  
  return(convex_hull_points)
}


summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- plyr::ddply(data, groupvars, .drop=.drop,
                       .fun = function(xx, col) {
                         c(N    = length2(xx[[col]], na.rm=na.rm),
                           mean = mean   (xx[[col]], na.rm=na.rm),
                           sd   = sd     (xx[[col]], na.rm=na.rm)
                         )
                       },
                       measurevar
  )
  
  # Rename the "mean" column    
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

QuPath_class_to_dataframe <- function(class_vector, sep = ": "){
  separated = sapply(class_vector, function(i) strsplit(i, sep))
  mat = matrix(0, ncol = length(unique(unlist(separated))), nrow =length(class_vector))
  colnames(mat) = unique(unlist(separated))
  mat = mat[,-which(is.na(colnames(mat)))]
  
  for(i in 1:length(separated)){
    vec = separated[[i]]
    for(j in 1:length(vec)){
      if(!is.na(vec[j])) mat[i,vec[j]] = 1
    }
  }
  return(as.data.frame(mat))
}


QuPath_filter_measurements <- function(mat){
  separated = sapply(class_vector, function(i) strsplit(i, sep))
  mat = matrix(0, ncol = length(unique(unlist(separated))), nrow =length(class_vector))
  colnames(mat) = unique(unlist(separated))
  mat = mat[,-which(is.na(colnames(mat)))]
  
  for(i in 1:length(separated)){
    vec = separated[[i]]
    for(j in 1:length(vec)){
      if(!is.na(vec[j])) mat[i,vec[j]] = 1
    }
  }
  return(as.data.frame(mat))
}

normalized_marker_sample <- function(mat, spe, ignore_markers = c("Cytokeratin") ) {
  mat. = mat
  intensities = as.data.frame(t(mat))
  intensities$sample_id = spe$sample_id
  intensities$cell_id = rownames(intensities)
  intensities = intensities %>% tidyr::pivot_longer(-c(sample_id,cell_id), names_to = c("marker"))
  
  median_markers = intensities %>% group_by(sample_id, marker) %>% 
    summarise(median_marker = median(value))
  median_markers$median_marker[median_markers$median_marker == 0] = quantile(median_markers$median_marker, 0.25)
  intensities = left_join(intensities, median_markers, by = c("sample_id", "marker"))
  intensities = intensities %>% mutate(value_normalized_by_median_marker_sample = value / median_marker)
  
  mat = intensities %>% select(cell_id, marker, value_normalized_by_median_marker_sample) %>% 
    pivot_wider(names_from = c(cell_id), values_from = value_normalized_by_median_marker_sample) %>% 
    as.data.frame
  rownames(mat) = mat$marker
  mat = mat[,-1]
  mat = as.matrix(mat)
  
  if(!is.null(ignore_markers)){
    mat[ignore_markers,] = mat.[ignore_markers,] / median(mat.[ignore_markers,])
  }
  return(mat)
}

marker_identification <- function(values, marker){
  # values. = MetaCyto:::baselineCut(values)
  
  q0.01 = quantile(values, 0.01) 
  q0.99 = quantile(values, 0.99) 
  values. = values[values  > q0.01 & values < q0.99]
  
  
  # perform Silverman test with Hall-York calibration:
  # p  = silverman.test(values., k=1, R=1000, adjust=T)
  
  # if(p@p_value > 0.01){
  #   cat("Mean + Sd\n")
  #   if(length(which(values < med)) > length(values)/2)
  #     return(mean(values) + sd(values)) else
  #       return(mean(values) - sd(values))
  # } else{
  cat("findCutoff\n")
  set.seed(48)
  cutoff = findCutoffGaussian(values., marker = marker)
  return(cutoff)
  # }
}

library(mclust)
findCutoffGaussian <- function(intensity, marker){
  
  gmm <- Mclust(intensity, G = 2)
  means <- gmm$parameters$mean
  variances <- gmm$parameters$variance$sigmasq
  variance2 = ifelse(length(variances) > 1, variances[2], variances[1])
  threshold <- max(c(means[1] + sqrt(variances[1]), 
                     means[2] - sqrt(variance2)))
  
  df <- data.frame(intensity)
  
  # Plot histogram
  histogram <- ggplot(df, aes(x = intensity)) +
    geom_histogram(aes(y = ..density..), alpha  = 0.3, bins = 100, fill = "lightblue", color = "grey30") +
    theme_minimal() +
    labs(x = "Intensity", y = "Density") + theme(legend.position = "none")
  
  
  # Add Gaussian curves
  x <- seq(min(intensity), max(intensity), length.out = 100)
  density1 <- dnorm(x, mean = means[1], sd = sqrt(variances[1]))
  density2 <- dnorm(x, mean = means[2], sd = sqrt(variance2))
  
  df_gaussian <- data.frame(x = rep(x, 2), density = c(density1, density2),
                            component = rep(c("Component 1", "Component 2"), each = length(x)))
  
  histogram <- histogram +
    geom_line(data = df_gaussian, aes(x = x, y = density, color = component), linetype = "dashed",  size = 1) + 
    scale_color_manual(values = c("darkred","darkgreen"))
  
  histogram <- histogram +
    geom_vline(xintercept = threshold, color = "#DEBA28", linetype = "dashed", size = 1) + ggtitle(marker)
  
  # Display the histogram
  print(histogram)
  return(threshold)
}


findCutoff = function (x, returnSil = FALSE,  minX = 0) {
  if (length(x) > 10000) {
    x = sample(x, 10000)
  }
  valley = seq(min(x), max(x), length.out = 100)
  if (!is.null(minX)) {
    valley = valley[valley > minX]
  }
  D = distances(x)
  sil = sapply(valley, function(v) {
    cluster = 1 * (x > v) + 1
    if (length(unique(cluster)) < 2) {
      return(-2)
    }
    ss <- cluster::silhouette(cluster, D)
    return(mean(ss[, 3]))
  })
  valley = valley[which.max(sil)]
  if (returnSil == FALSE) {
    return(valley)
  }  else {
    return(c(cuoff = valley, sil = max(sil)))
  }
}

cell_marker_identification <- function(spe, ignore = NULL){
  
  if(!is.null(ignore)){
    df = data.frame(t(normcounts(spe)))
    df$sample_id = spe$sample_id
    df = df %>% tidyr::gather(key = "marker", value = "value", -sample_id)
    for(mark in names(ignore)){
      df = df[which(df$marker != mark | !(df$sample_id %in% ignore[[mark]])),]
    }  
  } else {
    df = data.frame(t(normcounts(spe)))
    df$sample_id = spe$sample_id
    df = df %>% tidyr::gather(key = "marker", value = "value", -sample_id)
  }
  # MetaCyto cell average cutoff
  cyto_cutoffs = sapply(unique(df$marker), function(marker){
    print(marker)
    values = df$value[df$marker == marker]
    
    return(marker_identification(values, marker))
  }, USE.NAMES = TRUE)
  return(cyto_cutoffs)
}

cell_marker_identification_per_sample <- function(spe, ignore){
  df = data.frame(t(spe@assays@data$normcounts))
  df$sample_id = spe$sample_id
  df = df %>% tidyr::gather(key = "marker", value = "value", -sample_id)
  # for(mark in names(ignore)){
  #   df = df[-which(df$marker == mark & df$sample_id %in% ignore[[mark]]),]
  # }  
  
  list_cutoffs = list()
  for(samp in unique(spe$sample_id)){
    print(samp)
    samp_df = df %>% dplyr::filter(sample_id == samp)
    
    # Otsu threshold
    # otsu_cutoffs = sapply(unique(samp_df$marker), function(marker){
    #   img = readImage(spe@int_metadata$imgSources[paste0(samp,"-", marker)])
    #   img@.Data = round(img@.Data * 2^16)
    #   threshold = otsu(img, range = c(0, 256^2), levels = 2^16)
    #   return(threshold)
    # }, USE.NAMES = TRUE)
    
    # MetaCyto cell average cutoff
    med = median(samp_df$value)
    cutoffs = sapply(unique(samp_df$marker), function(marker){
      values = samp_df$value[samp_df$marker == marker]
      q95 = quantile(values, 0.95)
      q5 = quantile(values, 0.05)
      values. = values[values > q5 & values < q95]
      
      # perform Silverman test with Hall-York calibration:
      p  = silverman.test(values., k=1, R=1000, adjust=T)
      if(p@p_value > 0.01){
        
        if(length(which(values < med)) > length(values)/2)
          cutoff = mean(values) + sd(values) else
            cutoff =mean(values) - sd(values)
          
          print(hist(values, breaks = 250, cex.main = 0.9, main = paste0(samp, " - ", marker, " - Mean +- Sd", " - pval =", as.numeric(p@p_value))))
          print(abline(v = cutoff, col = "red", lty = 3))
          return(cutoff)
      } else{
        cutoff = MetaCyto::findCutoff(values, returnSil = FALSE, minX = 0)
        print(hist(values, breaks = 250, cex.main = 0.9, main = paste0(samp, " - ", marker, " - MetaCyto",  " - pval =", as.numeric(p@p_value))))
        print(abline(v = cutoff, col = "red", lty = 3))
        return(cutoff)
      }
    }, USE.NAMES = TRUE)
    
    list_cutoffs[[samp]] = cutoffs
  }
  
  return(list_cutoffs)
}

plotSPE <- function(spe, assay = c("normcounts","counts","positive_marker","metadata")[1], feature = NULL, x_coord = NULL, y_coord = NULL, 
                    sample_id = "sample_id", palette = c("gray90", "navy"), size = 0.3, log = F) 
{
  if (!is.null(feature)) 
    stopifnot(is.character(feature))
  if (is.null(x_coord)) 
    x_coord <- colnames(spatialCoords(spe))[1]
  if (is.null(y_coord)) 
    y_coord <- colnames(spatialCoords(spe))[2]
  
  n_samples <- length(table(colData(spe)[, sample_id]))
  
  if(assay == "normcounts") counts <- as.numeric(normcounts(spe)[feature, ])
  if(assay == "counts") counts <- as.numeric(counts(spe)[feature, ])
  if(assay == "positive_marker") counts <- as.character(spe@assays@data$positive_marker[feature, ])
  if(assay == "metadata") counts <- colData(spe)[,feature]
  
  if(log && is.numeric(counts)) counts = log2(counts + 1)
  stopifnot(length(counts) == ncol(spe))
  #palette <- .get_pal(palette)
  df <- cbind.data.frame(spatialCoords(spe), sum = counts, sample_id = spe$sample_id)
  if(length(unique(df$counts)) == 2){
    df$counts = as.factor(df$counts)
    levels(df$counts) = c(1,0)
  }
  p <- ggplot(df, aes_string(x = x_coord, y = y_coord, color = "sum")) + 
    geom_point(size = size) +  coord_fixed() + 
    ggtitle(feature) + theme_void()
  if(is.numeric(counts)) p = p + scale_color_gradientn(colours = c("gray90", rev(viridis::inferno(100))))
  if(assay == "positive_marker") p = p + scale_color_manual(values = c("gray90", "navy"), name = "Positivity")
  if (n_samples > 1) {
    p <- p + facet_wrap(~sample_id)
  }
  p
}

get_metadata_image_overlay <- function(spe, img_mat, sample, metadata, levels = NULL) {
  if(!is.null(levels)){
    cell_type = factor(colData(spe)[,metadata], level = levels)
  } else {
    cell_type = as.factor(colData(spe)[,metadata])
  }
  
  cell_ids = spe$sample_id == sample
  mapping = as.numeric(cell_type[cell_ids])
  names(mapping) = gsub(".*-","",spe$cell_id[cell_ids])
  # img_mat = as.sparse(img_mat)
  img_mat = as(img_mat, "TsparseMatrix")
  img_mat@x = mapping[match(img_mat@x, as.numeric(names(mapping)))]
  img_mat@x[which(is.na(img_mat@x))] = 0
  return(img_mat)
}

get_numeric_image_overlay <- function(spe, img_mat, sample, metadata, levels = NULL) {
  
  range = c(0, 255)
  
  cell_ids = spe$sample_id == sample
  mapping = as.numeric(spe[[metadata]])
  names(mapping) = gsub(".*-","",spe$cell_id[cell_ids])

  img_mat = as(img_mat, "TsparseMatrix")
  img_mat@x = mapping[match(img_mat@x, as.numeric(names(mapping)))]
  img_mat@x[which(is.na(img_mat@x))] = 0
  return(img_mat)
}



getImageAsMatrix <- function(file) {
  img = tiff::readTIFF(file.path(file),  as.is = TRUE)
  img = as.matrix(img)
  return(img)
}

rainbow_cell_raster <- function(img_mat, main, cex.main = 20, line = -10) {
  
  colors = sample(rainbow(max(img_mat)))
  colors = c("#000000", colors)
  img_mat = matrix(colors[as.vector(img_mat)+1], nrow = nrow(img_mat), ncol = ncol(img_mat))
  r = as.raster(img_mat)
  par(bg = "black")
  plot( c(0,ncol(r)), c(0,nrow(r)), axes = FALSE, type = "n", xlab = "", ylab = "")
  title(main = feature, adj = 0, line = line,  cex.main = cex.main, col.main = "white")
  graphics::rasterImage(r, xleft = 0, ytop = nrow(r), xright = ncol(r), ybottom = 0, interpolate = FALSE)
  par(bg = "white")
  
}

two_color_cell_raster <- function(img_mat, main , cex.main = 20, line = -10) {
  colors = c("grey", "darkblue")
  colors = c("#000000", colors)
  img_mat = matrix(colors[as.vector(img_mat)+1], nrow = nrow(img_mat), ncol = ncol(img_mat))
  r = as.raster(img_mat)
  par(bg = "black")
  plot( c(0,ncol(r)), c(0,nrow(r)), axes = FALSE, type = "n", xlab = "", ylab = "")
  title(main = main, adj = 0, line = line,  cex.main = cex.main, col.main = "white")
  graphics::rasterImage(r, xleft = 0, ytop = nrow(r), xright = ncol(r), ybottom = 0, interpolate = FALSE)
  par(bg = "white")
}


rgb_color_cell_raster <- function(img_mat, color_df, main , cex.main = 0.6, line = 0) {
  colors = unique(color_df[,2])
  colors = c("#FFFFFF", colors)
  img_mat = matrix(colors[as.vector(img_mat) + 1], nrow = nrow(img_mat), ncol = ncol(img_mat))
  r = as.raster(img_mat)
  par_save = par()
  par(bg = "black", mar = c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  plot( c(0,ncol(r)), c(0,nrow(r)), axes = FALSE, type = "n", xlab = "", ylab = "")
  title(main = main, adj = 0, line = line,  cex.main = cex.main, col.main = "white")
  graphics::rasterImage(r, xleft = 0, ytop = nrow(r), xright = ncol(r), ybottom = 0, interpolate = FALSE)
  legend(x = "topright", inset = c(-0.35, 0.1),  cex = 0.6, title = "Legend", title.col = "white", text.col = "white", ncol = 1,
         legend = c("", color_df[,1]), col = color_df[,2],
         fill = colors)
  par(par_save)
}
