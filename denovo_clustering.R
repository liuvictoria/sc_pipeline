################# LOAD UTILS  ##############
source("~/Documents/victoria_liu/matching_patients/R_Code/utils.R")



############# SUBSET CELLS ############
superset_assay <- ifelse (
  USE_ADT,
  "ADT",
  "GEX"
)
# this object is fully pre-processed for GEX
SeuratObj_filename <- paste0(
  RobjDir, analyses$denovo_superset, "/",
  superset_assay, analyses$denovo_superset, 
  "_res", analyses$denovo_superset_resolution, ".rds"
)

print(paste0("loading SeuratObj from file: ", SeuratObj_filename))
# change the loading info in the config file
SeuratObj <- readRDS(SeuratObj_filename)

# subset cells / remove unwanted clusters
# unwanted clusters can be visually identified in parent UMAPs
Idents(SeuratObj) <- analyses$denovo_subset_col
SeuratObj <- subset(
  SeuratObj, idents = analyses$denovo_subset
)
for (cluster_remove in analyses$denovo_clusters_remove) {
  SeuratObj <- subset(
    SeuratObj, subset = ClusterRNA != cluster_remove
    )
  SeuratObj <- subset(
    SeuratObj, subset = ClusterADT != cluster_remove
  )
}

######## NORMALIZATION ########
SeuratObj <- GEX_normalization(SeuratObj)

# includes (r)PCA for ADT
if (USE_ADT & ObjName == "ADT") {
  SeuratObj <- ADT_integrate(SeuratObj)
}

# #################  CELL POPULATION: SINGLER ###################
# ref_singler <- celldex::BlueprintEncodeData()
# 
# SeuratObj <- SingleR_wrapper(SeuratObj, ref_singler)
# 
# ################# (T CELLS ONLY) ASSIGNMENT ####################
{
#   ################# CELL POPULATION: PROJECTILS NO GATING ###################
#   # load murine TIL reference
#   ref_projectils_filename <- paste0(
#     dataDirectory, "projecTILs/", analyses$projecTILs_ref,
#     "_processed_mouse2human_reference.rds"
#   )
#   print(paste0(
#     "trying to read projecTILs reference: ",
#     ref_projectils_filename
#   ))
#   # ideally, it has already been pre-processed from murine to human
#   if (file.exists(ref_projectils_filename)) {
#     ref_projectils <- readRDS(ref_projectils_filename)
#   } else {
#     print("projecTILs reference has not been processed, processing now")
#     # we will preprocess from murine to human
#     # also includes recomp of umap / pca according to projecTILs spec
#     # first read in the original reference
#     ref_projectils <- readRDS(paste0(
#       dataDirectory, "projecTILs/", analyses$projecTILs_ref,
#       "_mouse_atlas.rds"
#     ))
#     ref_projectils <- projectils_ref_full_translation (
#       ref_projectils, SeuratObj
#     )
#   }
#   
#   SeuratObj <- projecTILs_wrapper(SeuratObj, ref_projectils, gating = TRUE)
#   
#   
#   ##################### REMOVE NK CELLS ####################
#   SeuratObj <- subset_Tcells(SeuratObj)
#   
#   ################# PROJECTILS NO GATING #################
#   SeuratObj <- projecTILs_wrapper(SeuratObj, ref_projectils, gating = FALSE)
}
######## PCA (GEX) ########
SeuratObj <- GEX_pca(
  SeuratObj, 
  paste0("Denovo ", config$Subset), 
  specific_PCA_features = T
  )
E1 <- ElbowPlot(SeuratObj)
pdf(paste0(
  ElbowDirectory, ObjName, Subset, 
  " pca.pdf"
), width = 12, height = 4.5, family = FONT_FAMILY
)
print(E1)
dev.off()


if (analyses$denovo_run_harmony) {
  SeuratObj@reductions$harmonyRNA <- NULL
  SeuratObj <- SeuratObj %>% 
    RunHarmony(
      group.by.vars = config$batch_norm_by,
      reduction.save = "harmonyRNA"
    )
  E2 <- ElbowPlot(SeuratObj, ndims = 50, reduction = "harmonyRNA")

  pdf(paste0(
    ElbowDirectory, ObjName, Subset,
    " elbow plot after harmony.pdf"
  ), width = 12, height = 4.5, family = FONT_FAMILY
  )
  print(E2)
  dev.off()
}


######## LOUVAIN CLUSTERING ########
reduction = ifelse(
  analyses$denovo_run_harmony, 
  "harmonyRNA", 
  "pca"
  )

for (resolution in config$RESOLUTIONS) {
  SeuratObj <- GEX_louvain(
    SeuratObj, resolution = resolution, reduction = reduction
  )
}

if (USE_ADT & ObjName != "RNA") {
  # ADT
  for (resolution in config$RESOLUTIONS) {
    SeuratObj <- ADT_louvain(SeuratObj, resolution)
  }
  # WNN
  for (resolution in config$RESOLUTIONS) {
    SeuratObj <- WNN_louvain(SeuratObj, resolution)
  }
}
# no default resolution!
SeuratObj$seurat_clusters <- NULL


######## RUN UMAP ########
SeuratObj <- RunUMAP(
  SeuratObj, 
  reduction = reduction, 
  dims = c(1:analyses$denovo_pcaGEX_dims),
  reduction.name = "umapRNA",
  reduction.key = "umapRNA_"
)

if (USE_ADT & ObjName != "RNA") {
  # ADT
  SeuratObj <- RunUMAP(
    SeuratObj, 
    reduction = "pcaADT", 
    dims = c(1:analyses$pcaADT_dims),
    reduction.name = "umapADT",
    reduction.key = "umapADT_"
  )
  
  # WNN
  SeuratObj <- RunUMAP(
    SeuratObj, 
    nn.name = "weighted.nn",
    reduction.name = "umapWNN",
    reduction.key = "UMAPWNN_"
  )
}


SeuratObj@meta.data[, ] <- lapply(
  SeuratObj@meta.data, 
  function(x) type.convert(as.character(x), as.is = TRUE)
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
  paste0(ConfigDirectory, ObjName, "_", Subset, "_sessionInfo.txt")
)
file.copy(
  from = here("config.json"), 
  to = paste0(ConfigDirectory, ObjName, "_", Subset, "_config_params.json")
)
file.copy(
  from = here("analysis.json"), 
  to = paste0(ConfigDirectory, ObjName, "_", Subset, "_analysis_params.json")
)


