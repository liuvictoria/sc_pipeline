################# LOAD UTILS + SEURATOBJ  ##############
source("~/Documents/victoria_liu/matching_patients/R_Code/utils.R")

# capture session info, versions, etc.
write_experimental_configs(code_file = "standard_viz")

# this object is fully pre-processed for GEX and/or ADT
if (ObjName != "WNN") {
  tmpObjName <- ObjName
} else {
  # ADT and WNN are preprocessed in the same object, called ADT
  tmpObjName <- "ADT"
}

RDS_filename <- paste0(
  RobjDirectory, tmpObjName, Subset, 
  "_resAll.rds"
)
if (! file.exists(RDS_filename)) {
  print ("no Seurat object to analyze")
}
print(paste0("reading SeuratObj from file ", RDS_filename))
SeuratObj <- readRDS(RDS_filename)

# if this is a denovo object, it may contain type info
# inherited from its superset parent
# remove, so that clustifyr has space to add new info, when called
if (analyses$denovo) {
  SeuratObj$r <- NULL
  SeuratObj$type <- NULL
}

# loads Clusterspecificgenes
Clusterspecificgenes <- master[[
  paste0(analyses[["denovo_lineage"]], "_lineage_markers")
]]

############### ADD CLUSTER METADATA  #################
SeuratObj <- add_cluster_metadata(SeuratObj)

# remove unnecessary metadata; focusing on one resolution at a time
for (metadata_col in colnames(SeuratObj@meta.data)) {
  if (
    grepl("snn_res.", metadata_col, fixed = TRUE) & 
    ! grepl(RESOLUTION, metadata_col, fixed = TRUE)
  ) {
    print (paste0("removing ", metadata_col, " from metadata"))
    SeuratObj[[metadata_col]] <- NULL
  }
}


################# CLUSTER / SAMPLE DIMPLOTS ##################
U2 <- plot_umap(
  seurat_object = SeuratObj,
  group_by = "Sample",
  color_scheme = analyses$hue_palette,
  subtype_by = analyses$chromatose_subtype_by,
  reduction = paste0("umap", analyses$viz_clustering),
  title = "Samples", xlab = "UMAP1", ylab = "UMAP2",
  legend_position = "bottom",
  ncol_guide = 3, shuffle = TRUE
)

pdf(paste0(
  UMAPDirectory, ObjName, Subset,
  "_res", RESOLUTION,
  "_Samples_UMAP.pdf"
), width = 20, height = 25, family = FONT_FAMILY
)
print(U2)
dev.off()





################# ASSIGNMENT BARPLOTS ###################
P5 <- plot_bargraph (
  seurat_object = SeuratObj, aesX = "Sample", fill = "Assignment",
  y_label = "Composition (Number of cells)", x_label = NULL,
  y_lower_limit = 0, y_break = 1000, 
  title = paste0(analyses$which_assignment, " Assignment")
)

pdf(paste0(
  densityplotDirectory, ObjName, Subset, 
  "res", RESOLUTION,
  "_number of cells per sample and Assignment (", 
  analyses$which_assignment,
  ") barplot.pdf"
), width = 35, height = 20, family = FONT_FAMILY
)
print(P5)
dev.off()



P6 <- plot_bargraph (
  seurat_object = SeuratObj, aesX = "Assignment", fill = "Sample",
  y_label = "Composition (percentage of cells)", x_label = NULL, 
  y_lower_limit = 0, y_break = 0.2,
  position = "fill",
  plot_margin = unit(c(0.2, 0.5, 0.2, 0.5), "cm"),
  title = paste0(analyses$which_assignment, " Assignment")
)
pdf(paste0(
  densityplotDirectory, ObjName, Subset, 
  "res", RESOLUTION, 
  "percent of cells per sample and Assignment (", 
  analyses$which_assignment,
  ") barplot.pdf"
), width = 35, height = 20, family = FONT_FAMILY
)
print(P6)
dev.off()


P7 <- plot_bargraph (
  seurat_object = SeuratObj, aesX = "Assignment", fill = "Sample",
  y_label = "Composition (Number of cells)", x_label = NULL, 
  y_lower_limit = 0, y_break = 1000,
  title = paste0(analyses$which_assignment, " Assignment")
)

pdf(paste0(
  densityplotDirectory, ObjName, Subset,
  "res", RESOLUTION, 
  "_percentage of cells per major population and Sample (", 
  analyses$which_assignment,
  ") barplot.pdf"
), width = 35, height = 20, family = FONT_FAMILY
)
print(P7)
dev.off()

P8 <- plot_bargraph (
  seurat_object = SeuratObj, aesX = "Sample", fill = "Assignment",
  y_label = "Composition (percentage of cells)", x_label = NULL,
  y_lower_limit = 0, y_break = 1000, position = "fill",
  title = paste0(analyses$which_assignment, " Assignment")
)

pdf(paste0(
  densityplotDirectory, ObjName, Subset, 
  "res", RESOLUTION,
  "_percentage of cells per sample and Assignment (", 
  analyses$which_assignment,
  ") barplot.pdf"
), width = 35, height = 20, family = FONT_FAMILY
)
print(P8)
dev.off()



################# ASSIGNMENT DIMPLOTS ###################
U4 <- plot_umap(
  seurat_object = SeuratObj, group_by = "Assignment",
  reduction = paste0("umap", analyses$viz_clustering),
  title = paste0(analyses$which_assignment, " Assignment"), 
  xlab = "UMAP1", ylab = "UMAP2",
  legend_position = "bottom",
  ncol_guide = 4,
  label_clusters = TRUE, 
  label_size = 5,
  color_reverse = FALSE,
  repel_labels = TRUE
)

pdf(paste0(
  UMAPDirectory, ObjName, Subset, 
  "_res", RESOLUTION, 
  "_", analyses$viz_clustering,
  "Clusters_Assignment (", 
  analyses$which_assignment,
  ") UMAP.pdf"
), width = 20, height = 25, family = FONT_FAMILY
)
print(U4)
dev.off()


################# CELL CYCLE SORTING ###################

U7 <- plot_umap (
  seurat_object = SeuratObj, group_by = "Phase",
  reduction = paste0("umap", analyses$viz_clustering),
  title = "Seurat Cell Cycle scoring", xlab = "UMAP1", ylab = "UMAP2",
  legend_position = "bottom",
  ncol_guide = 3,
  title_font_size = 16, x_font_size = 16, y_font_size = 16, 
  pt_size = 0.2,
  label_clusters = TRUE, repel_labels = TRUE, label_size = 3,
  shuffle = TRUE
)


pdf(paste0(
  UMAPDirectory, ObjName, Subset, 
  "CellCycle UMAP.pdf"
), width = 20, height = 25, family = FONT_FAMILY
)
print(U7)
dev.off()


############## SAVE SEURAT ################
saveRDS(
  SeuratObj,
  file = paste0(
    RobjDirectory, ObjName, Subset,
    "_res", RESOLUTION, ".rds"
  )
)