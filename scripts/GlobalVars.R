discrette_colors_50 = toupper(c(
  "#765aa4ff",
  "#ad2070ff",
  "#fdd62eff",
  "#96c9e6ff",
  "#f48121ff",
  "#68c3a5ff",
  "#ef3e2cff",
  "#0d522aff",
  "#42b649ff",
  "#660a17ff",
  "#3f7fc1ff",
  "#189cd5ff",
  "#0e8342ff",
  "#f9ae1dff",
  "#552e8aff",
  "#8b87c0ff",
  "#984d9dff",
  "#fec964ff",
  "#126badff",
  "#a51d23ff",
  "#e5569fff",
  "#eae71cff",
  "#B7C3F3",
  "#6B0504",
  "#CD8B76",
  "#FF785A",
  "#5FAD56",
  "#47A8BD", 
  "#A27E6F",
  "#FF6978",
  "#D72638",
  "#CE8964",
  "#FA8334",
  "#5E0B15",
  "#D1B1C8",
  "#F8BDC4",
  "#6BAA75",
  "#2D2327",
  "#EB9486",
  "#00E5E8",
  "#A5C4D4",
  "#D0CFEC",
  "#DDC4DD",
  "#9AA899",
  "#3BA99C",
  "#E58F65",
  "#92B4F4",
  "#030027",
  "#080066",
  "#00ABE7",
  "#91A6FF",
  "#FF88DC",
  "#FAFF7F",
  "#FFFFFF",
  "#1F2421",
  "#216869",
  "#49A078",
  "#9CC5A1",
  "#DCE1DE",
  "#FF9B42",
  "#F7A072",
  "#EDDEA4",
  "#D9E5D6",
  "#0FA3B1",
  "#861657",
  "#AA4465",
  "#FFA69E",
  "#93E1D8",
  "#DDFFF7",
  "#9883E5",
  "#FCD3DE",
  "#72A1E5",
  "#50C9CE",
  "#251101",
  "#470024",
  "#5B1865",
  "#2C5784",
  "#5688C7"
))

cell_cluster_color_df = data.frame(
  "cell_type" = sort(c(
    "Fibroblast",
    "Leukocyte",
    "B_cell",
    "T_helper",
    "T_cytotoxic",
    "T_regulatory",
    "pDC",
    "Basophil",
    "NKT",
    "Monocytes",
    "Neutrophils",
    "Macrophages",
    "Keratinocyte",
    "Endothelial",
    "Lymphatic"
  )),
  "cell_type_color" = c(
    "#401F7A",
             "#FF6978",
             "#e5569fff",
             "#DDC4DD",
             "#6B0504",
             "#189cd5ff",
             "#fec911",
             "#42b649ff",
             "#eae71cff",
             "#3BA99C",
             "#8A6C0B",
             "#00E5E8",
             "#D72638",
             "#0d522aff",
             "#f9ae1dff"
  ))


# cell_cluster_color_pixel_df = data.frame(
#   "cell_type" = sort(c(
#     "Fibroblast",
#     "Leukocyte",
#     "B_cell",
#     "T_helper",
#     "T_cytotoxic",
#     "T_regulatory",
#     "pDC",
#     "Basophil",
#     "NKT",
#     "Monocytes",
#     "Neutrophils",
#     "Macrophages",
#     "Keratinocyte-Spinosum",
#     "Keratinocyte-Spinosum",
#     "Endothelial",
#     "Lymphatic"
#   )),
#   "cell_type_color" = c(
#     "#401F7A",
#     "#FF6978",
#     "#e5569fff",
#     "#DDC4DD",
#     "#6B0504",
#     "#189cd5ff",
#     "#fec911",
#     "#42b649ff",
#     "#eae71cff",
#     "#3BA99C",
#     "#8A6C0B",
#     "#00E5E8",
#     "#D72638",
#     "#0d522aff",
#     "#f9ae1dff"
#   ))

cell_type_color_df = data.frame(
  cell_type_CellSighter = c(
    "Macrophages",
    "APC",
    "B_cell",
    "T_regulatory",
    "Monocytes",
    "Neutrophils",
    "NKT",
    "Keratinocyte",
    "Leukocyte",
    "Endothelial",
    "T_cytotoxic",
    "pDC",
    "T_helper",
    "Unknown-Stroma",
    "Lymphatic",
    "Basophil",
    "Monocytic_Lineage"
  ),
  cell_type_CellSighter_color = c(
    "#42b649ff",
    "#780C50",
    "#401F7A",
    "#f04410ff",
    "#eae71cff",
    "#3BA99C",
    "#8A6C0B",
    "#6B0504",
    "#189cd5ff",
    "#e5569fff",
    "#D72638",
    "#00E5E8",
    "#0d522aff",
    "#E5E6E4",
    "#fec911",
    "#FF6978",
    "#25303B"
  )
)


cell_type_simplified_color_df = setNames(
  c("#BD86D1",
    "#3e9488",
    "#E5E6E4",
    "#f29158",
    "#6B0504"
  ),
  c("Vessels",
    "Myeloid",
    "Dermis",
    "Lymphocyte",
    "Epidermis")
)

colCondition = data.frame(condition = c("TON","THY", "LN", "HD","LP","AD","PS","DAR", "MF","SS"),
                          condition_color = c( "#E5E6E4", "#A6A2A2", "#5A5A66", "#7699D4","#fdd62eff", "#FFB8F2", "#765aa4ff", "#39001E", "#a12070ff", "#D11E3Bff"))

cell_type_markers = read.csv("annotation/cell_markers.csv")
golden_markers = unique(cell_type_markers$marker)
golden_markers = golden_markers[c(1,17,3,6,7,8,9,11,12,13,14,15,16,10,5,18)]


struct_celltype = c("Keratinocyte", "Endothelial",  "Lymphatic", "Unknown-Stroma")
celltype_levels = c("Keratinocyte","Endothelial","Lymphatic","Leukocyte","APC","Monocytic_Lineage","Monocytes",
                    "Macrophages","Neutrophils","Basophil","pDC","B_cell","NKT","T_regulatory",
                    "T_cytotoxic","T_helper", "Unknown-Stroma")
APC_types = c("APC","Monocytic_Lineage","Monocytes",
              "Macrophages","Neutrophils","Basophil","pDC")
T_cells = c("NKT","T_regulatory",
            "T_cytotoxic","T_helper")

markers_to_remove = c("Ki67","CD209", "CD270", "CXCR4")




