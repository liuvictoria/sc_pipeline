################# LOAD UTILS  ##############
source("~/Box/Yun lab projects/victoria_liu/matching_patients/R_Code/utils.R")
SeuratSamples <- list()

######### QC + Doublet removal ########
if (!config$aggr_cells) {
  
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
  regexp <- "[[:digit:]]+"
  both_assays <- ADT_PRESENT[[str_extract(FILES[i], regexp)]]
  
  
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
  
  # determine doublets, need to change to sce
  SeuratObjMYSC.sce <- as.SingleCellExperiment(SeuratObjMYSC)
  sce <- scDblFinder(SeuratObjMYSC.sce)
  
  # Converting back to Seurat Obj
  SeuratObjMYSC <- as.Seurat(sce)
  CellInfo <- SeuratObjMYSC@meta.data
  Idents(SeuratObjMYSC) <- SeuratObjMYSC$orig.ident
  
  
  # SCT
  if (SCT) {
    SeuratObjMYSC <- RenameAssays(SeuratObjMYSC, RNA = "RNApreSCT")
    SeuratObjMYSC <- SeuratObjMYSC %>%
      SCTransform(
        assay = "RNApreSCT", method = "glmGamPoi", 
        new.assay.name = "RNA",
        vst.flavor = "v2", verbose = FALSE
      )
  } else {
    SeuratObjMYSC <- SeuratObjMYSC %>% 
      NormalizeData(assay = "RNA") %>%
      FindVariableFeatures(assay = "RNA") %>%
      ScaleData(assay = "RNA")
  }
  
  #Run PCA (code in utils)
  SeuratObjMYSC <- GEX_pca(SeuratObjMYSC, FILES[i])
  
  
  # QC for mitochondrial RNA
  SeuratObjMYSC$log10GenesPerUMI <- log10(SeuratObjMYSC$nFeature_RNA) / 
    log10(SeuratObjMYSC$nCount_RNA)
  
  SeuratObjMYSC$percent.mt <- PercentageFeatureSet(
    SeuratObjMYSC, pattern = "^MT-"
  )
  SeuratObjMYSC$mito_ratio <- SeuratObjMYSC$percent.mt / 100
  
  SeuratObjMYSC$ribo_ratio <- PercentageFeatureSet(
    SeuratObjMYSC, pattern = "^RP[LS]"
  )
  SeuratObjMYSC$ribo_ratio <- SeuratObjMYSC$ribo_ratio / 100
  
  # miQC
  SeuratObjMYSC <- RunMiQC(
    SeuratObjMYSC, 
    percent.mt = "percent.mt", 
    nFeature_RNA = "nFeature_RNA",
    backup.option = "percent",
    backup.percent = config$mito_ratio * 100
  )
  
  if (max(SeuratObjMYSC$miQC.probability) > 0) {
    # plots
    QC1 <- PlotMiQC(
      SeuratObjMYSC, color.by = "miQC.probability"
    ) + 
      ggplot2::scale_color_gradient(low = "grey", high = "purple")
    
    pdf(paste0(
      QCDirectory, FILES[i],
      " Scatterplot - Determining Probability of Low Quality Cells.pdf"
    ), width = 8, height = 5.5, family = "ArialMT")
    print(QC1)
    dev.off()
    
    QC2 <- PlotMiQC(SeuratObjMYSC, color.by = "miQC.keep")
    pdf(paste0(
      QCDirectory, FILES[i],
      " Scatterplot - Cells Passing QC.pdf"
    ), width = 8, height = 5.5, family = "ArialMT")
    print(QC2)
    dev.off()
  }
  
  # save plot
  D <- DimPlot(
    SeuratObjMYSC, split.by = "scDblFinder.class",
    order = TRUE, shuffle = TRUE
  )
  pdf(paste0(
    QCDirectory, FILES[i],
    "Doublet status dimplot.pdf"
  ), width = 8, height = 5.5, family = FONT_FAMILY
  )
  print(D)
  dev.off()
  
  # subset based on miQC
  SeuratObjMYSC <- subset(SeuratObjMYSC, miQC.keep == "keep")
  
  # rename columns
  SeuratObjMYSC$DoubletStatus <- SeuratObjMYSC$scDblFinder.class
  SeuratObjMYSC$scDblFinder.class <- NULL
  SeuratObjMYSC$scDblFinder.score <- NULL
  SeuratObjMYSC$scDblFinder.cxds_score <- NULL
  SeuratObjMYSC$scDblFinder.weighted <- NULL
  SeuratObjMYSC$miQC.probability <- NULL
  
  # save object
  saveRDS(SeuratObjMYSC, file = paste0(RobjDirectory, FILES[i], ".rds"))
  SeuratSamples[[i]] <- SeuratObjMYSC
}

}
######## AGGREGATION ########## 
if (config$aggr_cells) {
  stopifnot(length(config$aggr_FILES) > 0)
  for (idx in 1:length(config$aggr_FILES)) {
    aggr_file <- config$aggr_FILES[idx]
    aggr_filename <- paste0(
      RobjDir, aggr_file, "/",
      "GEX", aggr_file, "_res0.9.rds"
    )
    SeuratSamples[[idx]] <- readRDS(aggr_filename)
  }
  
  SeuratObj <- merge(
    x = SeuratSamples[[1]], y = SeuratSamples[-1], 
    add.cell.ids = config$aggr_FILES
  )
  
  # subset immune cells as indicated
  # default identity is double status (for now)
  Idents(SeuratObj) <- SeuratObj$Assignment
  
  # remove doublets from Seurat data structure 
  SeuratObj <- subset(
    SeuratObj, 
    cells = WhichCells(SeuratObj, idents = config$aggr_categories)
  )
  
} else {
  # merge samples into one object
  SeuratObj <- merge(
    x = SeuratSamples[[1]], y = SeuratSamples[-1], 
    add.cell.ids = config$FILES
  )
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



# extract current metadata
CellInfo <- SeuratObj@meta.data

# add sample information
sample_info <- config$sample_info

CellInfo$Sample <- NA
for (i in 1:length(sample_info)) {
  CellInfo$Sample[
    which(str_detect(CellInfo$orig.ident, names(sample_info[i])))
  ] <- sample_info[[i]]
}

Sample <- CellInfo$Sample
names(Sample) <- colnames(x = SeuratObj)
SeuratObj <- AddMetaData(
  object = SeuratObj, metadata = Sample, col.name = "Sample"
)

# add patient information (anonymized, deidentified)
patient_info <- config$patient_info

CellInfo$Patient <- NA
for (i in 1:length(patient_info)) {
  CellInfo$Patient[
    which(str_detect(CellInfo$orig.ident, names(patient_info[i])))
  ] <- patient_info[[i]]
}

Patient <- CellInfo$Patient
names(Patient) <- colnames(x = SeuratObj)
SeuratObj <- AddMetaData(
  object = SeuratObj, metadata = Patient, col.name = "Patient"
)

#add group information
group_info <- config$group_info

CellInfo$Group <- NA
for (i in 1:length(group_info)) {
  CellInfo$Group[
    which(str_detect(CellInfo$orig.ident, names(group_info[i])))
  ] <- group_info[[i]]
}

Group <- CellInfo$Group
names(Group) <- colnames(x = SeuratObj)
SeuratObj <- AddMetaData(
  object = SeuratObj, metadata = Group, col.name = "Group"
)


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
DefaultAssay(SeuratObj) <- "RNA"
# more general set of CC genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# score cells based on what part of cell cycle they are in
SeuratObj <- CellCycleScoring(
  SeuratObj, 
  s.features = s.genes, 
  g2m.features = g2m.genes, 
  set.ident = TRUE
)

# regress out cell cycle
# there is also an option to score out G2 vs S, while keeping
# the dichotomy of dividing vs not dividing
# Reference: https://satijalab.org/seurat/articles/cell_cycle_vignette.html

# GEX
if (SCT) {
  # sct / cellcyclescoring discussion:
  # https://github.com/satijalab/seurat/issues/1679
  # use pre-SCT RNA when normalizing; don't want to normalize twice
  # RNA assay already exists, but we are overwriting it w scaled CC
  SeuratObj <- SCTransform(
    SeuratObj, method = "glmGamPoi", 
    assay = "RNApreSCT", new.assay.name = "RNA",
    vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE
  )
} else {
  SeuratObj <- ScaleData(
    SeuratObj,
    vars.to.regress = c("S.Score", "G2M.Score"), 
    features = rownames(SeuratObj)
  )
}

# regressing is a really long process, so save for future
saveRDS(
  SeuratObj, 
  file = paste0(    
    RobjDirectory, ObjName, Subset, 
    "_res", config$RESOLUTION, ".rds"
    )
)
 
########### PCA & HARMONY BATCH CORRECTION (GEX) ##########
DefaultAssay(SeuratObj) <- "RNA"

SeuratObj <- SeuratObj %>% 
  RunPCA(
    features = VariableFeatures(object = SeuratObj)
  ) %>% 
  RunHarmony(
    group.by.vars = config$batch_norm_by,
    reduction.save = "harmonyRNA"
  )
E1 <- ElbowPlot(SeuratObj, ndims = 15, reduction = "harmonyRNA")

pdf(paste0(
  ElbowDirectory, ObjName, Subset, 
  " elbow plot after CC regression and harmony.pdf"
), width = 12, height = 4.5, family = FONT_FAMILY
)
E1
dev.off()


######## LOUVAIN CLUSTERING (GEX) ########
analyses <- fromJSON(file = here("analysis.json"))

for (resolution in config$RESOLUTIONS) {
  SeuratObj <- GEX_louvain(
    SeuratObj, resolution = resolution, reduction = "harmonyRNA"
  )
}
# no default resolution!
SeuratObj$seurat_clusters <- NULL


######## RUN UMAP (GEX) ########
SeuratObj <- RunUMAP(
  SeuratObj, 
  reduction = "harmonyRNA", 
  dims = c(1:analyses$harmonyRNA_dims),
  reduction.name = "umapRNA",
  reduction.key = "umapRNA_"
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
