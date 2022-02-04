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
SeuratObj <- readRDS(
  paste0(RobjDirectory, ObjName, Subset, "_GEX", ".rds")
)

######## QC PLOTS ########

x5 <- plot_densitygraph (
  seurat_object = SeuratObj, aesX = "nCount_ADT", fill = "Sample",
  color_by = "Sample", xintercept = 0.8, scale_x_log10 = TRUE
)

x6 <- plot_densitygraph (
  seurat_object = SeuratObj, aesX = "nFeature_ADT", fill = "Sample",
  color_by = "Sample", xintercept = 0.8, scale_x_log10 = TRUE
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



######## NORMALIZATION (ADT) ########
# Normalize for ADT 
# not a lot of features, so use all antibody rows
# scaling happens after CC regression and integration
DefaultAssay(SeuratObj) <- "ADT"
SeuratObj <- NormalizeData(
  SeuratObj, normalization.method = "CLR", margin = 2
)
VariableFeatures(SeuratObj) <- rownames(SeuratObj[["ADT"]])

# plot variable features
top10 <- head(VariableFeatures(SeuratObj), 10)
plot1 <- VariableFeaturePlot(SeuratObj)
plot2 <- LabelPoints(plot = plot1, points = top10, size = 3)
pdf(paste0(
  QCDirectory, ObjName, Subset, 
  " ADT variable freatures.pdf"
), width = 8, height = 5.5, family = FONT_FAMILY
)
plot2 + theme_min()
dev.off()


########### CELL CYCLE REGRESSION (ADT) #############
# cell cycle scoring already previously calculated
# ADT
DefaultAssay(SeuratObj) <- "ADT"
SeuratObj <- ScaleData(
  SeuratObj,
  vars.to.regress = c("S.Score", "G2M.Score"),
  features = rownames(SeuratObj)
)

# regressing is a really long process, so save for future
saveRDS(
  SeuratObj, 
  file = paste0(RobjDirectory, ObjName, Subset, "_ADT", ".rds")
)

########### BATCH INTEGRATION (ADT) ###########
DefaultAssay(SeuratObj) <- "ADT"
# normalize per sample
# using CLR normalization based on: 
# https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html
samples_list <- SplitObject(SeuratObj, split.by = config$batch_norm_by) %>%
  lapply(
    FUN = function(x) {
      x <- NormalizeData(
        x, normalization.method = "CLR", margin = 2, assay = "ADT"
      )
    }
  )
# select features that are repeatedly variable across data sets
features <- SelectIntegrationFeatures(
  object.list = samples_list
)

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
  new.assay.name = "integratedADT"
)

# scale data & dimreduc PCA, in prep for WNN
DefaultAssay(SeuratObj_ADT) <- "integratedADT"
SeuratObj_ADT <- SeuratObj_ADT %>%
  ScaleData() %>%
  RunPCA(reduction.name = "pcaADT")

# add ADT data back into original SeuratObj
# necessary bc IntegrateData destroys scale.data for all assays
# ignore any warnings about offending keys
# ref: https://github.com/satijalab/seurat/issues/3843
SeuratObj[["integratedADT"]] <- SeuratObj_ADT[["integratedADT"]]
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

SeuratObj <- SeuratObj %>%
  FindMultiModalNeighbors(
    reduction.list = list("harmonyRNA", "pcaADT"), 
    dims.list = list(
      1 : analyses$harmonyRNA_dims, 1 : analyses$pcaADT_dims
    ), 
    modality.weight.name = "RNA.weight"
)

CellInfo <- SeuratObj@meta.data
# Rename Clusters
cluster_count <- length(
  levels(as.factor(SeuratObj@meta.data$seurat_clusters))
)
for(j in 1 : cluster_count){
  if (j < 10){
    CellInfo$ClusterWNN[CellInfo$seurat_clusters == j - 1] <- paste0("WC0", j)
  }
  else {
    CellInfo$ClusterWNN[CellInfo$seurat_clusters == j - 1] <- paste0("WC", j)
  }
}
SeuratObj@meta.data <- CellInfo
Idents(SeuratObj) <- CellInfo$ClusterWNN

# Get number of cells per cluster and per Sample
write.csv(
  as.matrix(table(SeuratObj@meta.data$ClusterWNN, SeuratObj@meta.data$Sample)),
  file = paste0(
    OutputDirectory, ObjName, Subset, 
    " number of cells per cluster and sample.csv"
  )
)

######## RUN UMAP ########
SeuratObj <- RunUMAP(
  SeuratObj, 
  nn.name = "weighted.nn",
  reduction.name = "umapWNN",
  reduction.key = "UMAPWNN_"
)

################ FIND CLUSTER BIOMARKERS ################

Idents(SeuratObj) <- CellInfo$ClusterWNN
markers1 <- FindAllMarkers(
  SeuratObj,
  assay = "RNA",
  logfc.threshold = 0.5, 
  test.use = "wilcox", 
  only.pos = TRUE
)

markers2 <- FindAllMarkers(
  SeuratObj,
  assay = "ADT",
  logfc.threshold = 0.5, 
  test.use = "wilcox", 
  only.pos = TRUE
)

# save, because it takes a little time to calculate
write.csv(
  markers1, 
  paste0(
    OutputDirectory, ObjName, Subset, 
    " RNA cluster markers (by WNN)", "res", RESOLUTION, ".csv"
  )
)

write.csv(
  markers2, 
  paste0(
    OutputDirectory, ObjName, Subset, 
    " ADT cluster markers (by WNN)", "res", RESOLUTION, ".csv"
  )
)

############## SAVE SEURAT AND SESSION INFO, LOOSE ENDS ################
saveRDS(
  SeuratObj, 
  file = paste0(RobjDirectory, ObjName, Subset, "_ADT", ".rds")
)


# capture session info, versions, etc.
writeLines(
  capture.output(sessionInfo()), 
  paste0(ConfigDirectory, "sessionInfo_ADT.txt")
)
file.copy(
  from = here("config.json"), 
  to = paste0(ConfigDirectory, "config_params_ADT.json")
)
file.copy(
  from = here("analysis.json"), 
  to = paste0(ConfigDirectory, "analysis_params_ADT.json")
)