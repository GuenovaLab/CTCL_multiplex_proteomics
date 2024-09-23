# Calculate distances between each and every cell
# authors: Pacome Prompsy
# contact: pacome.prompsy@chuv.ch
# Guenova Lab
# CHUV (Centre Hospitalier Universitaire Vaudois), Lausanne, Suisse

matrix_list = qs::qread(file.path(output_dir, "matrix_list.qs"))

samps = c( "ROI-15", "ROI-04", "ROI-08", "ROI-24")

samps = c("ROI-06", "ROI-23", "ROI-12", "ROI-27", "ROI-05", "ROI-11", "ROI-17", "ROI-22")
samps = c("ROI-13")
markers = c("CD2", "CD5", "CD1c", "Bcl2", "CD195", "CD196", "CD11b", "HLAABC", "HLADR", "Galectin9", "Galectin3")
interaction_filtered = spe@metadata$interaction
  
for(celltype1 in c( "T_helper")){
  celltype1_dir = file.path(output_dir, celltype1)
  if(!dir.exists(celltype1_dir)) dir.create(celltype1_dir)
  
  
  list_celltype_dfs = list()
  
  for(celltype2 in  "T_cytotoxic"){
    
    print("############################################")
    print(celltype2)
    print("############################################")
    
    list_df = list()
    
    interaction_mat = interaction_filtered
    
    for(sample in samps){
      print(sample)
      spe. = spe[,spe$sample_id == sample]
      cells1 = spe.$cell_id[spe.$cell_type_CellSighter == celltype1]
      
      if(celltype2 == "immune"){
        cells2 = spe.$cell_id[!spe.$cell_type_CellSighter %in% struct_celltype]
      } else{
        cells2 = spe.$cell_id[spe.$cell_type_CellSighter == celltype2]
      }
      
      
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
          
          
          cell_overlay_mat = getImageAsMatrix(file.path("output", "segmentation", paste0(sample,"_whole_cell.tiff")))
          
          
          for(marker in markers){
            mat = matrix_list[[marker]]
            mat_sym = mat[which(rownames(mat) %in% cells1.),
                          which(colnames(mat) %in% cells2.)] +
              t(mat[which(colnames(mat) %in% cells2.),
                    which(rownames(mat) %in% cells1.)])
            
            cells_interacting = c(names(which(rowSums(mat_sym) != 0)),
              names(which(colSums(mat_sym) != 0)))
            
            
            spe.$interacting = 0
            spe.$interacting[match(cells_interacting, spe.$cell_id)] = 1
            
            celltype_img_mat = get_metadata_image_overlay(spe = spe., img_mat = cell_overlay_mat,
                                                          sample = sample, metadata = "interacting", levels = c(0, 1, 2))

            
            tiff::writeTIFF(as.matrix(celltype_img_mat),
                              file.path(output_dir, "Snapshots", paste0(sample,"_",marker,"_mask.tiff")),
                            bits.per.sample = 8)
          }
        }
      }
    }
  }
}