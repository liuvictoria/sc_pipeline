################# LOAD UTILS  ##############
source("~/Documents/victoria_liu/matching_patients/R_Code/utils.R")

# capture session info, versions, etc.
write_experimental_configs()

############# SUBSET CELLS FROM PARENT ############
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
  if (ObjName == "GEX" & ! USE_ADT) {
    SeuratObj <- subset(
      SeuratObj, subset = ClusterRNA != cluster_remove
    )
  }
  if (USE_ADT & ObjName == "ADT") {
    SeuratObj <- subset(
      SeuratObj, subset = ClusterADT != cluster_remove
    )
  }
}

######## NORMALIZATION ########
SeuratObj <- GEX_normalization(SeuratObj)

# includes (r)PCA for ADT
if (USE_ADT & ObjName == "ADT") {
  SeuratObj <- ADT_integrate(SeuratObj)
}

#################  CELL POPULATION: SINGLER ###################
if (! exists("refs_singler")) {
  define_refs_singler()
}
for (reference_name in names(refs_singler)) {
  SeuratObj <- SingleR_wrapper(
    SeuratObj, refs_singler[[reference_name]],
    reference_name, celltype = Subset
  )
}



################# (T CELLS ONLY) PROJECTILS ASSIGNMENT ####################
if (analyses$denovo_lineage == "T" & grepl("T", Subset, fixed = T)) {
  ################# CELL POPULATION: PROJECTILS WITH GATING ###################
  # assuming ref projectils has already been preprocessed; 
  # otherwise need to pass in a SeuratObj
  ref_projectils <- load_ref_projectils()
  SeuratObj <- projecTILs_wrapper(SeuratObj, ref_projectils, gating = TRUE)
  
  ##################### REMOVE NK CELLS ####################
  # remove NK cells based on projecTILs and SingleR
  SeuratObj <- subset_Tcells(SeuratObj)
  
  
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
############## SAVE SEURAT ################
saveRDS(
  SeuratObj, 
  file = paste0(
    RobjDirectory, ObjName, Subset, 
    "_resAll.rds"
    )
)

