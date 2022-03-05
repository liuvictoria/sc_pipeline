################# LOAD UTILS  ##############
source("~/Box/Yun lab projects/victoria_liu/matching_patients/R_Code/utils.R")

# laod visualizing parameters
analyses <- fromJSON(file = here("analysis.json"))


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
SeuratObj <- GEX_louvain(
  SeuratObj, reduction = reduction, cluster_prefix = "TC"
  )

######## RUN UMAP (GEX) ########
SeuratObj <- RunUMAP(
  SeuratObj, 
  reduction = reduction, 
  dims = c(1:analyses$denovo_pcaGEX_dims),
  reduction.name = "umapRNA",
  reduction.key = "umapRNA_"
)

################ FIND CLUSTER BIOMARKERS (GEX) ################
markers <- GEX_cluster_markers(SeuratObj)

# save, because it takes a little time to calculate
write.csv(
  markers, 
  paste0(
    OutputDirectory, ObjName, Subset, 
    " RNA cluster markers (by RNA)", "res", RESOLUTION, ".csv"
  )
)

########## SAVE LOUPE PROJECTIONS (GEX) ##########
write.csv(
  Embeddings(SeuratObj, "umapRNA"), 
  paste0(LoupeDirectory, ObjName, Subset, "_umapRNA.csv")
)

# sample and Seurat cluster by cell barcode
write.csv(
  select(SeuratObj@meta.data, ClusterRNA), 
  paste0(LoupeDirectory, ObjName, Subset, "_clusters.csv")
)

write.csv(
  select(SeuratObj@meta.data, Sample), 
  paste0(LoupeDirectory, ObjName, Subset, "_samples.csv")
)


############## SAVE SEURAT AND SESSION INFO, LOOSE ENDS ################
saveRDS(
  SeuratObj, 
  file = paste0(
    RobjDirectory, ObjName, Subset, 
    "_res", config$RESOLUTION, ".rds"
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

