################# LOAD UTILS + SEURATOBJ ##############
source("/projects/compsci/USERS/alizae/GBM/matching_patients/R_Code/sc_pipeline/utils.R")

# capture session info, versions, etc.
write_experimental_configs(
  suffix = "_ADT", code_file = "ADT_preprocessing"
  )

# this object is fully pre-processed for GEX
# it must include ADT data, but the ADT data is not yet preprocessed
# QC and doublet removal have already been done, courtesy of GEX
# it includes sample, group, patient info
# important dimreducs: harmonyRNA, umapRNA
# important metadata: ClusterRNA
# important assays: RNA
# for full list, refer to Notion guide
RDS_filename <- paste0(
  RobjDirectory, "GEX", Subset, 
  "_resAll.rds"
)
if (! file.exists(RDS_filename)) {
  print(paste0(RDS_filename, " Seurat obj does not exist"))
} else {
  print(paste0("reading Seurat obj: ", RDS_filename))
}

SeuratObj <- readRDS(
  RDS_filename
)


temp=readRDS(paste0(
   RobjDirectory, "GEX", Subset,
   "_res0.4.rds"))

# need to manually load temp object
# temp object is basically RDS_filename except res that has been visualized

SeuratObj$GEX_Assignment <- temp$Assignment#????

######## QC PLOTS ########

x5 <- plot_densitygraph (
  seurat_object = SeuratObj, aesX = "nCount_ADT", fill = "Sample",
  color_by = "Sample", scale_x_log10 = TRUE
)

x6 <- plot_densitygraph (
  seurat_object = SeuratObj, aesX = "nFeature_ADT", fill = "Sample",
  color_by = "Sample"
)

XX <- grid_arrange_shared_legend(
  x5 + theme_min(), 
  x6 + theme_min(), 
  position = "right", 
  ncol = 2, nrow = 1
)

#violin plots
X <- VlnPlot(
  SeuratObj, 
  features = c("nFeature_ADT", "nCount_ADT", "mito_ratio"), 
  ncol = 3, pt.size = 0, 
  cols = get_colors(seurat_object = SeuratObj, color_by = "Sample"), 
  group.by = "Sample"
) +
  RotatedAxis()
x5 = ggarrange(
  X[[1]] + theme_min() + NoLegend() + RotatedAxis(), 
  X[[2]] + theme_min() + NoLegend() + RotatedAxis(), 
  X[[3]] + theme_min() + NoLegend() + RotatedAxis(), 
  nrow = 1
)

pdf(paste0(
  QCDirectory, ObjName, Subset, 
  " ADT feature and count.pdf"
), width = 12, height = 10, family = FONT_FAMILY
)
grid.arrange(XX, x5, heights = c(0.8, 0.5))
dev.off()



########### BATCH INTEGRATION (ADT) ###########
# this function is also used in de novo clustering
SeuratObj <- RenameAssays(SeuratObj, ADT = "ADTpreInt")
SeuratObj <- ADT_integrate(SeuratObj)



######## ADT RUN UMAP ########
# ADT
SeuratObj <- RunUMAP(
  SeuratObj, 
  reduction = "pcaADT", 
  dims = c(1:analyses$pcaADT_dims),
  reduction.name = "umapADT",
  reduction.key = "umapADT_"
)
######## LOUVAIN CLUSTERING ########
# ADT
for (resolution in config$RESOLUTIONS) {
  SeuratObj <- ADT_louvain(
    SeuratObj, resolution = resolution
  )
}
# for WNN, we need to run louvain before umap
# WNN
for (resolution in config$RESOLUTIONS) {
  SeuratObj <- WNN_louvain(SeuratObj, resolution = resolution)
}

######## WNN RUN UMAP ########

# WNN
SeuratObj <- RunUMAP(
  SeuratObj, 
  nn.name = "weighted.nn",
  reduction.name = "umapWNN",
  reduction.key = "UMAPWNN_"
)

# make sure there are no factors
SeuratObj@meta.data[, ] <- lapply(
  SeuratObj@meta.data, 
  function(x) type.convert(as.character(x), as.is = TRUE)
)

# no default resolution!
SeuratObj$seurat_clusters <- NULL

################ SAMPLE VISUALIZATION #############
U0 <- plot_umap(
  seurat_object = SeuratObj,
  group_by = paste0("Assignment"),
  reduction = paste0("umapADT"),
  title = "Clusters", xlab = "UMAP1", ylab = "UMAP2",
  legend_position = "bottom",
  label_size = 5,
  ncol_guide = 5,
  color_reverse = TRUE, label_clusters = TRUE,
  shuffle = TRUE
)
U0

############## SAVE SEURAT AND SESSION INFO, LOOSE ENDS ################
saveRDS(
  SeuratObj, 
  file = paste0(
    RobjDirectory, ObjName, Subset, 
    "_resAll.rds"
  )
)
