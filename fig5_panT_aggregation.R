# this module contains code that is more customized / less standard
# than the analysis from standard_viz
# none of this code is particularly reusable, except for its specific use case
################# LOAD UTILS ##############
source("~/hpc/Sumner/R_Code/utils_HPC.R")

# capture session info, versions, etc.
write_experimental_configs(code_file = "panT_aggregation")

oldwarning <- getOption("warn")
options(warn = -1)


#################### AGGREGATE PAN T CANCER OBJECTS #####################
cancer_objects <- list()
for (idx in 1:length(FILES)) {
  print(FILES[idx])
  cancer_objects[[idx]] <- readRDS(paste0(RobjDirectory, FILES[idx], ".rds"))
}

# merge into one Seurat object
external_pan <- merge(
  x = cancer_objects[[1]], y = cancer_objects[-1],
  add.cell.ids = config$FILES
  )

# beautify metadata
add_tumor_suffix <- function(word) {
  return (paste0(word, "_pan"))
}
external_pan$Sample <- unlist(lapply(external_pan$libraryID, add_tumor_suffix))
external_pan$Group <- "pan"
external_pan$Patient <- external_pan$patient
external_pan$Type <- external_pan$cancerType
external_pan$Assignment <- external_pan$meta.cluster

# external_pan$libraryID <- NULL
# external_pan$patient <- NULL
# external_pan$cancerType <- NULL
# external_pan$meta.cluster <- NULL


# normalize data
external_pan <- GEX_normalization(external_pan)

# # save before CC regression
saveRDS(
  external_pan,
  paste0(RobjDirectory, ObjName, Subset, "_resAll.rds")
)

# CC scaling regression
external_pan <- GEX_cc_regression(external_pan)


# PCA and harmony batch correciton
DefaultAssay(external_pan) <- "RNA"

external_pan <- external_pan %>% 
  GEX_pca(file = NULL) %>% 
  RunHarmony(
    group.by.vars = config$batch_norm_by,
    reduction.save = "harmonyRNA"
  )

E1 <- ElbowPlot(external_pan, ndims = 50, reduction = "harmonyRNA")
pdf(paste0(
  ElbowDirectory, ObjName, Subset, 
  " elbow plot after CC regression and harmony.pdf"
), width = 12, height = 4.5, family = FONT_FAMILY
)
print(E1)
dev.off()


# RUN UMAP (GEX)
external_pan <- RunUMAP(
  external_pan, 
  reduction = "harmonyRNA", 
  dims = c(1:analyses$harmonyRNA_dims),
  reduction.name = "umapRNA",
  reduction.key = "umapRNA_"
)

# LOUVAIN CLUSTERING (GEX) 
for (resolution in config$RESOLUTIONS) {
  external_pan <- GEX_louvain(
    external_pan, resolution = resolution
  )
}
# make sure there are no factors
external_pan@meta.data[, ] <- lapply(
  external_pan@meta.data, 
  function(x) type.convert(as.character(x), as.is = TRUE)
)

# no default resolution!
external_pan$seurat_clusters <- NULL


saveRDS(
  external_pan,
  paste0(RobjDirectory, ObjName, Subset, "_resAll.rds")
)

options(warn = oldwarning)
