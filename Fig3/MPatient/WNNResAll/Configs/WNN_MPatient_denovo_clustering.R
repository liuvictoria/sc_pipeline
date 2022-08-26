################# LOAD UTILS  ##############
source("~/Documents/victoria_liu/matching_patients/R_Code/utils.R")

# capture session info, versions, etc.
write_experimental_configs(code_file = "denovo_clustering")

############# SUBSET CELLS FROM PARENT ############
# this object is fully pre-processed for GEX
SeuratObj_filename <- paste0(
  RobjDir, analyses$denovo_superset, "/",
  ObjName, analyses$denovo_superset, 
  "_res", analyses$denovo_superset_resolution, ".rds"
)

print(paste0("loading SeuratObj from file: ", SeuratObj_filename))
# change the loading info in the config file
SeuratObj <- readRDS(SeuratObj_filename)

# subset cells / remove unwanted clusters
# unwanted clusters can be visually identified in parent UMAPs
stopifnot(
  length(analyses$denovo_subset_cols) == length(analyses$denovo_subset)
  )
for (idx in 1:length(analyses$denovo_subset_cols)) {
  denovo_subset_col <- analyses$denovo_subset_cols[idx]
  Idents(SeuratObj) <- denovo_subset_col
  SeuratObj <- subset(
    SeuratObj, idents = analyses$denovo_subset[[idx]]
  )
}

for (cluster_remove in analyses$denovo_clusters_remove) {
  if (ObjName == "GEX") {
    SeuratObj <- subset(
      SeuratObj, subset = ClusterRNA != cluster_remove
    )
  }
  if (USE_ADT & ObjName == "ADT") {
    SeuratObj <- subset(
      SeuratObj, 
      subset = ClusterADT != cluster_remove
    )
  }
  if (USE_ADT & ObjName == "WNN") {
    SeuratObj <- subset(
      SeuratObj, 
      subset = ClusterWNN != cluster_remove
    )
  }
}

######## NORMALIZATION ########
SeuratObj <- GEX_normalization(SeuratObj)

# includes (r)PCA for ADT
if (USE_ADT & ObjName != "GEX") {
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
E1 <- ElbowPlot(SeuratObj, ndims = 50)
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


######## (GEX + ADT) RUN UMAP ########
GEX_reduction = ifelse(
  analyses$denovo_run_harmony,
  "harmonyRNA",
  "pca"
)
SeuratObj <- RunUMAP(
  SeuratObj, 
  reduction = GEX_reduction, 
  dims = c(1:analyses$harmonyRNA_dims),
  reduction.name = "umapRNA",
  reduction.key = "umapRNA_"
)

if (USE_ADT & ObjName != "GEX") {
  # ADT
  SeuratObj <- RunUMAP(
    SeuratObj, 
    reduction = "pcaADT", 
    dims = c(1:analyses$pcaADT_dims),
    reduction.name = "umapADT",
    reduction.key = "umapADT_"
  )
}



######## LOUVAIN CLUSTERING ########
for (resolution in config$RESOLUTIONS) {
  SeuratObj <- GEX_louvain(
    SeuratObj, resolution = resolution
  )
}

if (USE_ADT & ObjName != "GEX") {
  # ADT
  for (resolution in config$RESOLUTIONS) {
    SeuratObj <- ADT_louvain(
      SeuratObj, resolution = resolution
      )
  }
  # WNN
  for (resolution in config$RESOLUTIONS) {
    SeuratObj <- WNN_louvain(SeuratObj, resolution = resolution)
  }
}
# no default resolution!
SeuratObj$seurat_clusters <- NULL



######## (WNN) RUN UMAP ######
# for WNN, we need to run louvain before running UMAP
if (USE_ADT & ObjName != "GEX") {
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


############### VALIDATION VISUALIZATION ##############
U0 <- plot_umap(
  seurat_object = SeuratObj, 
  group_by = "wsnn_res.0.3",
  reduction = paste0("umapWNN"),
  title = "Clusters", xlab = "UMAP1", ylab = "UMAP2",
  legend_position = "bottom",
  ncol_guide = 5,
  color_reverse = F, label_clusters = TRUE
)
U0


U1 <- plot_umap(
  seurat_object = SeuratObj, 
  group_by = "Assignment",
  reduction = paste0("umapWNN"),
  title = "Clusters", xlab = "UMAP1", ylab = "UMAP2",
  legend_position = "bottom",
  ncol_guide = 5,
  color_reverse = F, label_clusters = TRUE
)
U1

U2 <- plot_umap(
  seurat_object = SeuratObj, 
  group_by = "Sample",
  reduction = paste0("umapWNN"),
  title = "Clusters", xlab = "UMAP1", ylab = "UMAP2",
  legend_position = "bottom",
  ncol_guide = 5,
  color_reverse = F, label_clusters = TRUE,
  color_scheme = analyses$hue_palette,
  subtype_by = analyses$chromatose_subtype_by
)
U2




############## SAVE SEURAT ################
if (ObjName == "WNN") {
  ObjName <- "ADT"
}

saveRDS(
  SeuratObj, 
  file = paste0(
    RobjDirectory, ObjName, Subset, 
    "_resAll.rds"
    )
)


# need to save this so we can load smoothly during standard_viz
if (ObjName != "GEX") {
  saveRDS(
    SeuratObj, 
    file = paste0(
      RobjDirectory, "GEX", Subset, 
      "_resAll.rds"
    )
  )
}

SeuratObj <- subset(
  SeuratObj,
  subset = 
    wsnn_res.0.4 != 0 &
    wsnn_res.0.4 != 6 &
    wsnn_res.0.4 != 7 &
    wsnn_res.0.4 != 9 &
    wsnn_res.0.4 != 12
)
