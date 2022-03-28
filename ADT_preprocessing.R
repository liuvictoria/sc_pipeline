################# LOAD UTILS  ##############
source("~/Box/Yun lab projects/victoria_liu/matching_patients/R_Code/utils.R")

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
SeuratObj <- RenameAssays(SeuratObj, ADT = "ADTpreInt")
DefaultAssay(SeuratObj) <- "ADTpreInt"
# normalize per sample
# using CLR normalization based on: 
# https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html
samples_list <- SplitObject(SeuratObj, split.by = config$batch_norm_by) %>%
  lapply(
    FUN = function(x) {
      x <- NormalizeData(
        x, normalization.method = "CLR", margin = 2, assay = "ADTpreInt"
      )
    }
  )
# select features that are repeatedly variable across data sets
features <- SelectIntegrationFeatures(
  object.list = samples_list
)

samples_list <- lapply(
  X = samples_list, 
  FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
  })

# find anchor cells (this step stakes a while)
# use rpca; reference:
# https://satijalab.org/seurat/archive/v3.2/integration.html
cell_anchors <- FindIntegrationAnchors(
  object.list = samples_list, anchor.features = features, reduction = "rpca"
)
# combine into integrated Seurat object
SeuratObj_ADT <- IntegrateData(
  anchorset = cell_anchors, 
  normalization.method = "LogNormalize",
  new.assay.name = "ADT"
)

# scale data & dimreduc PCA, in prep for WNN
DefaultAssay(SeuratObj_ADT) <- "ADT"
SeuratObj_ADT <- SeuratObj_ADT %>%
  ScaleData() %>%
  RunPCA(reduction.name = "pcaADT")

# add ADT data back into original SeuratObj
# necessary bc IntegrateData destroys scale.data for all assays
# ignore any warnings about offending keys
# ref: https://github.com/satijalab/seurat/issues/3843
SeuratObj[["ADT"]] <- SeuratObj_ADT[["ADT"]]
SeuratObj[["pcaADT"]] <- SeuratObj_ADT[["pcaADT"]]

E2 <- ElbowPlot(SeuratObj, ndims = 15, reduction = "pcaADT")
pdf(paste0(
  ElbowDirectory, ObjName, Subset, 
  " elbow plot after CC scaling integration and pca (ADT).pdf"
), width = 12, height = 4.5, family = FONT_FAMILY
)
E2
dev.off()

######## LOUVAIN CLUSTERING ########
analyses <- fromJSON(file = here("analysis.json"))

for (resolution in config$RESOLUTIONS) {
  SeuratObj <- ADT_louvain(SeuratObj, resolution)
}

# no default resolution!
SeuratObj$seurat_clusters <- NULL

######## RUN UMAP ########
SeuratObj <- RunUMAP(
  SeuratObj, 
  nn.name = "weighted.nn",
  reduction.name = "umapWNN",
  reduction.key = "UMAPWNN_"
)



############## SAVE SEURAT AND SESSION INFO, LOOSE ENDS ################
saveRDS(
  SeuratObj, 
  file = paste0(
    RobjDirectory, ObjName, Subset, 
    "_resAll.rds"
  )
)

# capture session info, versions, etc.
writeLines(
  capture.output(sessionInfo()), 
  paste0(ConfigDirectory, ObjName, "_", Subset, "_sessionInfo_ADT.txt")
)
file.copy(
  from = here("config.json"), 
  to = paste0(ConfigDirectory, ObjName, "_", Subset, "_config_params_ADT.json")
)
file.copy(
  from = here("analysis.json"), 
  to = paste0(ConfigDirectory, ObjName, "_", Subset, "_analysis_params_ADT.json")
)
