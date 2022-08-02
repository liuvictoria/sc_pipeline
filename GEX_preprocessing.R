################# LOAD UTILS  ##############
source("~/Documents/victoria_liu/matching_patients/R_Code/utils.R")

# capture session info, versions, etc.
write_experimental_configs(
  suffix = "_RNA", code_file = "GEX_preprocessing"
  )



######### QC + Doublet removal ########
SeuratSamples <- list()
if (!config$aggr_cells & !config$preprocess_existing_atlas) {
  
  for(i in 1:length(FILES)){
    print(FILES[i])
    # read in GEX (and also ADT, if applicable)
    MYSC <- Read10X(
      data.dir = paste0(
        Directory,
        "/Data/",
        FILES[i], "/",
        config$count_matrices_dir, "/"
      )
    ) 
    
    # if there is also ADT data
    both_assays <- ADT_PRESENT[[FILES[i]]]
    
    
    # create Seurat obj
    # GEX
    SeuratObjMYSC <- CreateSeuratObject(
      counts = if (both_assays) MYSC[["Gene Expression"]] else MYSC, 
      min.features = 100, project = FILES[i]
    )
    
    if (USE_ADT) {  
      # ADT data, for removing bad cells QC
      SeuratObjMYSC[["ADT"]] <- CreateAssayObject(
        MYSC[["Antibody Capture"]][, colnames(x = SeuratObjMYSC)])
    }
    
    # subset ADT, if present
    if (USE_ADT) {  
      SeuratObjMYSC <- subset(
        SeuratObjMYSC, 
        subset = nFeature_ADT > config$nFeature_ADT
      )
    }
    
    SeuratObjMYSC <- GEX_QC(SeuratObjMYSC, FILES[i])
    
    # save object
    saveRDS(SeuratObjMYSC, file = paste0(RobjDirectory, FILES[i], ".rds"))
    SeuratSamples[[i]] <- SeuratObjMYSC
  }
  
} else if (config$preprocess_existing_atlas & !config$aggr_cells) {
SeuratSamples <- atlas_QC()
}


######## AGGREGATION ########## 
if (config$aggr_cells) {
  stopifnot(length(config$aggr_FILES) > 0)
  for (idx in 1:length(config$aggr_FILES)) {
    aggr_file <- config$aggr_FILES[idx]
    aggr_filename <- paste0(
      RobjDir, aggr_file, "/",
      "GEX", aggr_file, "_res", 
      config$aggr_resolutions[idx], ".rds"
    )
    print(paste0("reading RDS file ", aggr_filename))
    SeuratSample <- readRDS(aggr_filename)
    DefaultAssay(SeuratSample) <- "RNApreSCT"
    
    SeuratSample <- DietSeurat(
        SeuratSample, assays = c("RNApreSCT", "ADT")
      )
    
    SeuratSamples[[idx]] <- SeuratSample
  }
  
  SeuratObj <- merge(
    x = SeuratSamples[[1]], y = SeuratSamples[-1], 
    add.cell.ids = config$aggr_FILES
  )
  
  # subset immune cells as indicated
  # default identity is assignment
  Idents(SeuratObj) <- SeuratObj@meta.data[[config$aggr_category_col]]
  aggr_categories <- intersect(
    unique(SeuratObj@meta.data[[config$aggr_category_col]]), 
    config$aggr_categories
    )
  SeuratObj <- subset(
    SeuratObj, 
    cells = WhichCells(SeuratObj, idents = aggr_categories)
    )
  
} else {
  stopifnot(length(SeuratSamples) > 0)
  if (length(SeuratSamples) == 1) {
    SeuratObj <- SeuratSamples[[1]]
  } else {
    # merge samples into one object
    SeuratObj <- merge(
      x = SeuratSamples[[1]], y = SeuratSamples[-1], 
      add.cell.ids = config$FILES
    )
  }
}


########### QC PLOTS ########
# plot doublet status
PDS <- plot_bargraph (
  seurat_object = SeuratObj, aesX = "DoubletStatus", fill = "orig.ident",
  y_label = "Composition (Number of cells)", x_label = NULL,
  y_lower_limit = 0, y_break = 2000
)
# save plot
pdf(paste0(
  QCDirectory, ObjName, Subset,
  "Barplot number of Doublets per sample calc by doubletfinder",
  "at sevenPointFive percent freq.pdf"
), width = 8, height = 5.5, family = FONT_FAMILY)
PDS
dev.off()

if (! config$aggr_cells) {
  # default identity is double status (for now)
  Idents(SeuratObj) <- SeuratObj$DoubletStatus
  
  # remove doublets from Seurat data structure
  SeuratObj <- subset(
    SeuratObj, cells = WhichCells(SeuratObj, idents = "singlet")
  )
}

### ADD METADATA
if (! config$preprocess_existing_atlas) {
  SeuratObj <- add_all_master_metadata(SeuratObj)
}

# add categories for ribosomal expression (for dimplot)
percentage_to_category <- function(percentage) {
  if (percentage < 0.1) {
    return ("10 percent or less")
  } else if (0.1 <= percentage & percentage < 0.2) {
    return ("10 to 20 percent")
  } else if (0.2 <= percentage & percentage < 0.3) {
    return ("20 to 30 percent")
  } else if (0.3 <= percentage & percentage < 0.4) {
    return ("30 to 40 percent")
  } else {
    return ("40 percent or more")
  }
}

SeuratObj$ribo_category <- map_chr(
  SeuratObj$ribo_ratio, percentage_to_category
)

# plot
P1 <- plot_bargraph (
  seurat_object = SeuratObj, aesX = "Sample", fill = "Sample",
  y_label = "Composition (Number of cells)", x_label = NULL,
  y_lower_limit = 0, y_break = 1000
)

pdf(paste0(
  QCDirectory, ObjName, Subset,
  " Barplot number of cells per sample.pdf"),
  width = 8, height = 5.5, family = FONT_FAMILY
)
P1
dev.off()


# plot
x1 <- plot_densitygraph (
  seurat_object = SeuratObj, aesX = "nCount_RNA", fill = "Sample",
  color_by = "Sample", scale_x_log10 = TRUE
)

x2 <- plot_densitygraph (
  seurat_object = SeuratObj, aesX = "nFeature_RNA", fill = "Sample",
  color_by = "Sample", xintercept = config$nFeature_RNA, scale_x_log10 = TRUE
)

x3 <- plot_densitygraph (
  seurat_object = SeuratObj, aesX = "mito_ratio", fill = "Sample",
  color_by = "Sample", xintercept = config$mito_ratio
)

x4 <- plot_densitygraph (
  seurat_object = SeuratObj, aesX = "log10GenesPerUMI", fill = "Sample",
  color_by = "Sample", scale_x_log10 = TRUE
)

x5 <- plot_densitygraph (
  seurat_object = SeuratObj, aesX = "ribo_ratio", fill = "Sample",
  color_by = "Sample",  xintercept = config$ribo_ratio
)


XX <- grid_arrange_shared_legend(
  x1 + theme_min(),
  x2 + theme_min(),
  x3 + theme_min(),
  x4 + theme_min(),
  x5 + theme_min(),
  position = "right",
  ncol = 2, nrow = 3
)

#violin plots
X <- VlnPlot(
  SeuratObj,
  features = c("nFeature_RNA", "nCount_RNA", "mito_ratio", "ribo_ratio"),
  ncol = 4, pt.size = 0,
  cols = get_colors(seurat_object = SeuratObj, color_by = "Sample"),
  group.by = "Sample"
) +
  RotatedAxis()





x6 = ggarrange(
  X[[1]] + theme_min() + NoLegend() + RotatedAxis(),
  X[[2]] + theme_min() + NoLegend() + RotatedAxis(),
  X[[3]] + theme_min() + NoLegend() + RotatedAxis(),
  X[[4]] + theme_min() + NoLegend() + RotatedAxis(),
  nrow = 2
)

pdf(paste0(
  QCDirectory, ObjName, Subset,
  " GEX feature and count.pdf"
), width = 12, height = 16, family = FONT_FAMILY
)
grid.arrange(XX, x6, heights = c(0.8, 0.5))
dev.off()



######## NORMALIZATION (GEX) ########
# function placed in utils; to avoid code duplication
# when denovo clustering
SeuratObj <- GEX_normalization(SeuratObj)
########### CELL CYCLE REGRESSION (GEX) #############
SeuratObj <- GEX_cc_regression(SeuratObj)


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



########### PCA & HARMONY BATCH CORRECTION (GEX) ##########
SeuratObj <- SeuratObj %>% 
  GEX_pca("all samples", specific_PCA_features = T) %>% 
  RunHarmony(
    group.by.vars = config$batch_norm_by,
    reduction.save = "harmonyRNA"
  )

E1 <- ElbowPlot(SeuratObj, ndims = 50, reduction = "harmonyRNA")
pdf(paste0(
  ElbowDirectory, ObjName, Subset, 
  " elbow plot after CC regression and harmony.pdf"
), width = 12, height = 4.5, family = FONT_FAMILY
)
E1
dev.off()


########### RUN UMAP (GEX) ########
SeuratObj <- RunUMAP(
  SeuratObj,
  reduction = "harmonyRNA",
  dims = c(1:analyses$harmonyRNA_dims),
  reduction.name = "umapRNA",
  reduction.key = "umapRNA_"
)


########### LOUVAIN CLUSTERING (GEX) ########
for (resolution in config$RESOLUTIONS) {
  SeuratObj <- GEX_louvain(
    SeuratObj, resolution = resolution
  )
}
# make sure there are no factors
SeuratObj@meta.data[, ] <- lapply(
  SeuratObj@meta.data,
  function(x) type.convert(as.character(x), as.is = TRUE)
  )

# no default resolution!
SeuratObj$seurat_clusters <- NULL




############### VALIDATION VISUALIZATION ##############
U1 <- plot_umap(
  seurat_object = SeuratObj, 
  group_by = "RNA_snn_res.0.4",
  reduction = paste0("umapRNA"),
  title = "Clusters", xlab = "UMAP1", ylab = "UMAP2",
  legend_position = "bottom",
  ncol_guide = 5,
  color_reverse = F, label_clusters = TRUE
)
U1





############## SAVE SEURAT AND SESSION INFO, LOOSE ENDS ################
saveRDS(
  SeuratObj, 
  file = paste0(
    RobjDirectory, ObjName, Subset, 
    "_resAll.rds"
    )
)