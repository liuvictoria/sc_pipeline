################# LOAD UTILS  ##############
source("~/Box/Yun lab projects/victoria_liu/matching_patients/R_Code/utils.R")

# laod visualizing parameters
analyses <- fromJSON(file = here("analysis.json"))

#
############# SUBSET CELLS (GEX) ############
# this object is fully pre-processed for GEX
# change the loading info in the config file
SeuratObj <- readRDS(
  paste0(
    RobjDir, analyses$denovo_superset, "/",
    "GEX", analyses$denovo_superset, 
    "_res", analyses$denovo_superset_resolution, ".rds"
  )
)
Idents(SeuratObj) <- "Assignment"
SeuratObj <- subset(SeuratObj, idents = analyses$denovo_subset)


######## NORMALIZATION (GEX) ########
SeuratObj <- GEX_normalization(SeuratObj)


######## PCA (GEX) ########
SeuratObj <- GEX_pca(SeuratObj, paste0("Denovo ", config$Subset))
E1 <- ElbowPlot(SeuratObj)
pdf(paste0(
  ElbowDirectory, ObjName, Subset, 
  " pca.pdf"
), width = 12, height = 4.5, family = FONT_FAMILY
)
E1
dev.off()


if (analyses$denovo_run_harmony) {
  SeuratObj@reductions$harmonyRNA <- NULL
  SeuratObj <- SeuratObj %>% 
    RunHarmony(
      group.by.vars = config$batch_norm_by,
      reduction.save = "harmonyRNA"
    )
  E2 <- ElbowPlot(SeuratObj, ndims = 15, reduction = "harmonyRNA")

  pdf(paste0(
    ElbowDirectory, ObjName, Subset,
    " elbow plot after harmony.pdf"
  ), width = 12, height = 4.5, family = FONT_FAMILY
  )
  E2
  dev.off()
}


######## LOUVAIN CLUSTERING (GEX) ########
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
# no default resolution!
SeuratObj$seurat_clusters <- NULL



######## RUN UMAP (GEX) ########
SeuratObj <- RunUMAP(
  SeuratObj, 
  reduction = reduction, 
  dims = c(1:analyses$denovo_pcaGEX_dims),
  reduction.name = "umapRNA",
  reduction.key = "umapRNA_"
)

################# (CUSTOM) CELL POPULATION: PROJECTILS ###################
# load murine TIL reference
ref_projectils_filename <- paste0(
  dataDirectory, "projecTILs/processed_mouse2human_reference.rds"
)
# ideally, it has already been pre-processed from murine to human
if (file.exists(ref_projectils_filename)) {
  ref_projectils <- readRDS(ref_projectils_filename)
} else {
  # we will preprocess from murine to human
  # also includes recomp of umap / pca according to projecTILs spec
  # first read in the original reference
  ref_projectils <- readRDS(paste0(
    dataDirectory, "projecTILs/ref_TILAtlas_mouse_v1.rds"
  ))
  ref_projectils <- projectils_ref_full_translation (
    ref_projectils, SeuratObj
  )
}

# plot ref umap
U_TILs1 <-  plot_umap(
  ref_projectils, group_by = "functional.cluster",
  title = "projecTILs reference with filtered Seurat genes", 
  xlab = "UMAP_1", ylab = "UMAP2",
  legend_position = "right", reduction = "umap", 
  label_clusters = TRUE
)

pdf(paste0(
  dataDirectory, "projecTILs/", 
  " projecTILs recalculated with human genes and prcmp_umap UMAP.pdf"
), width = 12, height = 7, family = FONT_FAMILY
)
print(U_TILs1)
dev.off()


# TIME TO ACTUALLY ASTRAL PROJECT 
# project per sample. Don't project on integrated samples
DefaultAssay(SeuratObj) <- analyses$projecTILs_assay
manual_colors <- get_colors(ref_projectils, color_by = "functional.cluster")
projectils_metadata <- data.frame()
for (sample_name in unique(SeuratObj$Sample)) {
  
  # subset based on sample and assignment
  query.data <- subset(
    SeuratObj, Sample == sample_name
    ) 
  
  # proJECT
  query.projected <- make.projection(query.data, ref = ref_projectils)
  query.projected <- cellstate.predict(
    ref <- ref_projectils, query = query.projected
  )
  
  # save metadata to add back to SeuratObj later
  sample_metadata <- as.data.frame(
    query.projected[[c('functional.cluster', 'functional.cluster.conf')]]
    )
  projectils_metadata <- rbind(projectils_metadata, sample_metadata)

  ### PLOTS ###
  # draw projection
  U_TILs2 <- plot.projection(
    ref_projectils, 
    query.projected, 
    cols = manual_colors
  )
  pdf(paste0(
    projTILDirectory, 
    sample_name, " ", analyses$projecTILs_assay, "assay projecTILs",  
    " query projected onto ref UMAP.pdf"
  ), width = 12, height = 7, family = FONT_FAMILY
  )
  print(U_TILs2)
  dev.off()
  
  # cell type frequency table
  write.csv(
    table(query.projected$functional.cluster), 
    paste0(
      projTILDirectory, 
      sample_name, " ", analyses$projecTILs_assay, "assay projecTILs", 
      " cell type distribution.csv"
    ), 
    row.names = FALSE
  )
}

# add functional.cluster and functional.cluster.conf info
SeuratObj <- AddMetaData(
  SeuratObj, 
  metadata = projectils_metadata, 
  col.name = c("projecTILs", "projecTILs_conf")
  )

# add back NK cell info 
# (BTW, for TTumor, do clustifyr FIRST, to get NK designations)
SeuratObj$projecTILs <- ifelse(
  is.na(SeuratObj$projecTILs) & 
    grepl("NK", SeuratObj$Assignment, fixed = TRUE),
  "NK-like", SeuratObj$projecTILs
)


# plot total distribution
P1 <- plot_bargraph (
  seurat_object = SeuratObj, aesX = "Sample", fill = "projecTILs",
  y_label = "Composition (percentage of cells)", x_label = NULL,
  title = "projecTILs T cell assignment; clustifyr NK-like assignment",
  y_lower_limit = 0, y_break = 1000, position = "fill"
)

pdf(paste0(
  projTILDirectory, 
  "ALL_SAMPLES ", analyses$projecTILs_assay, "assay projecTILs", 
  " percentage of cells per sample and Assignment barplot.pdf"
), width = 6, height = 5.5, family = FONT_FAMILY
)
print(P1)
dev.off()


################# (CUSTOM) CELL POPULATION: SINGLER ###################
ref_singler <- celldex::BlueprintEncodeData()
singler_predictions <-SingleR(
  test = SeuratObj[[analyses$SingleR_assay]]@data, 
  ref = ref_singler, 
  labels = ref_singler$label.main
  )

# add metadata back to Seurat
SeuratObj$SingleR <- singler_predictions[["pruned.labels"]] 

# cell type frequency table
write.csv(
  table(SeuratObj$SingleR), 
  paste0(
    singleRDirectory, 
    analyses$SingleR_assay, "assay SingleR", 
    " cell type distribution.csv"
  ), 
  row.names = FALSE
)

# plot total distribution
P2 <- plot_bargraph (
  seurat_object = SeuratObj, aesX = "Sample", fill = "SingleR",
  y_label = "Composition (percentage of cells)", x_label = NULL,
  title = "SingleR T cell assignment, pruned labels",
  y_lower_limit = 0, y_break = 1000, position = "fill"
)

pdf(paste0(
  singleRDirectory, 
  "ALL_SAMPLES ", analyses$SingleR_assay, "assay SingleR", 
  " percentage of cells per sample and Assignment barplot.pdf"
), width = 6, height = 5.5, family = FONT_FAMILY
)
print(P2)
dev.off()

# QC pruning
P3 <- plotDeltaDistribution(singler_predictions)
pdf(paste0(
  singleRDirectory, 
  "ALL_SAMPLES ", analyses$SingleR_assay, "assay SingleR", 
  " pruning distribution.pdf"
), width = 6, height = 5.5, family = FONT_FAMILY
)
print(P3)
dev.off()

# QC heatmap
P4 <- plotScoreHeatmap(singler_predictions)
pdf(paste0(
  singleRDirectory, 
  "ALL_SAMPLES ", analyses$SingleR_assay, "assay SingleR", 
  " score heatmap.pdf"
), width = 6, height = 5.5, family = FONT_FAMILY
)
print(P4)
dev.off()

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
  paste0(ConfigDirectory, ObjName, "_", Subset, "_sessionInfo_RNA.txt")
)
file.copy(
  from = here("config.json"), 
  to = paste0(ConfigDirectory, ObjName, "_", Subset, "_config_params_RNA.json")
)
file.copy(
  from = here("analysis.json"), 
  to = paste0(ConfigDirectory, ObjName, "_", Subset, "_analysis_params_RNA.json")
)

