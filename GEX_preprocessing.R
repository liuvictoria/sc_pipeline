################# LOAD UTILS  ##############
source("~/Box/Yun lab projects/victoria_liu/matching_patients/R_Code/utils.R")
######### QC + Doublet removal ########
SeuratSamples <- list()

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
  
  
  # QC for mitochondrial RNA
  SeuratObjMYSC$log10GenesPerUMI <- log10(SeuratObjMYSC$nFeature_RNA) / 
    log10(SeuratObjMYSC$nCount_RNA)
  
  SeuratObjMYSC$mito_ratio <- PercentageFeatureSet(
    SeuratObjMYSC, pattern = "^MT-"
  )
  SeuratObjMYSC$mito_ratio <- SeuratObjMYSC@meta.data$mito_ratio / 100
  
  #filter by removing cells that didn't pass qc
  SeuratObjMYSC <- subset(
    SeuratObjMYSC, 
    subset = nFeature_RNA > config$nFeature_RNA & 
      mito_ratio < config$mito_ratio
  )
  # subset ADT, if present
  if (USE_ADT) {  
    SeuratObjMYSC <- subset(
      SeuratObjMYSC, 
      subset = nFeature_ADT > config$nFeature_ADT
    )
  }
  
  
  #prep work for DoubletFinder
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
  
  #Run PCA
  SeuratObjMYSC <- RunPCA(
    SeuratObjMYSC,
    features = VariableFeatures(object = SeuratObjMYSC)
  )
  # most representative (by absolute value) genes for PCs
  VS = VizDimLoadings(SeuratObjMYSC, dims = 1:2, reduction = "pca")
  VSD = DimPlot(SeuratObjMYSC, reduction = "pca") + theme_min()

  # save plot
  pdf(paste0(
    QCDirectory, FILES[i], " GEX PCA features and plot.pdf"
  ),
  width = 12, height = 4.5, family = FONT_FAMILY
  )
  grid.arrange(
    VS[[1]] + theme_min(),
    VS[[2]] + theme_min(),
    VSD,
    nrow = 1, widths = c(0.6, 0.6, 1)
  )
  dev.off()

  # remove doublets
  # documentation: https://github.com/chris-mcginnis-ucsf/DoubletFinder
  # pK Identification (no ground-truth)
  # if using sctransform, mark sct as TRUE
  sweep.res.list <- paramSweep_v3(SeuratObjMYSC, PCs = 1:10, sct = SCT)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)

  # to plot max bcmvn for identifying pK (refer to graph on wiki page)
  df = as.data.frame(bcmvn)
  pK_graph = df %>%
    ggplot(aes( x = ParamID, y = BCmetric)) +
    geom_point(color = "royalblue") +
    geom_line(color = "royalblue") +
    theme_classic() +
    # plot pK on x-axis
    geom_vline(xintercept = df$ParamID[df$BCmetric == max(df$BCmetric)])

  pdf(paste0(
    QCDirectory, FILES[i],
    " doubletfinder BCmvn distribution.pdf"),
    width = 8, height = 5.5, family = FONT_FAMILY
  )
  print(pK_graph)
  dev.off()

  # Homotypic Doublet Proportion Estimate
  # I think $Clusters may be empty / just a placeholder?
  homotypic.prop <- modelHomotypic(SeuratObjMYSC@meta.data$Clusters)
  # poisson distribution
  nExp_poi <- round(DOUBLET_FORMATION_RATE * nrow(SeuratObjMYSC@meta.data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

  ## Run DoubletFinder with varying classification stringencies
  # use 10 principle components, doublet count as 0.25 (not too important)
  # pK and pANN (neighborhood size) calculated using parameter sweep
  SeuratObjMYSC <- doubletFinder_v3(
    SeuratObjMYSC, PCs = 1:10,
    pN = 0.25, pK = as.numeric(df$pK[df$BCmetric == max(df$BCmetric)]),
    nExp = nExp_poi, reuse.pANN = FALSE, sct = SCT
  )

  SeuratObjMYSC <- doubletFinder_v3(
    SeuratObjMYSC, PCs = 1:10,
    pN = 0.25, pK = as.numeric(df$pK[df$BCmetric == max(df$BCmetric)]),
    nExp = nExp_poi.adj, reuse.pANN = colnames(SeuratObjMYSC@meta.data)[5],
    sct = SCT
  )

  # rename columns
  colnames(SeuratObjMYSC@meta.data)[
    length(SeuratObjMYSC@meta.data) - 1
  ] <- "pANN"
  colnames(SeuratObjMYSC@meta.data)[
    length(SeuratObjMYSC@meta.data)
  ] <- "DoubletStatus"

  # save plot
  pdf(paste0(
    QCDirectory, FILES[i],
    "Doublet status dimplot.pdf"
  ), width = 8, height = 5.5, family = FONT_FAMILY
  )
  DimPlot(
    SeuratObjMYSC, split.by = "DoubletStatus",
    order = TRUE, shuffle = TRUE
  )
  dev.off()

  # save object
  saveRDS(SeuratObjMYSC, file = paste0(RobjDirectory, FILES[i], ".rds"))
  SeuratSamples[[i]] <- SeuratObjMYSC
}


######## AGGREGATION + QC PLOTS ########

# merge samples into one object
SeuratObj <- merge(
  x = SeuratSamples[[1]], y = SeuratSamples[-1], 
  add.cell.ids = config$FILES
)


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

# default identity is double status (for now)
Idents(SeuratObj) <- SeuratObj$DoubletStatus

# remove doublets from Seurat data structure 
SeuratObj <- subset(
  SeuratObj, cells = WhichCells(SeuratObj, idents = "Singlet")
)


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
  color_by = "Sample", xintercept = 1000, scale_x_log10 = TRUE
)

x2 <- plot_densitygraph (
  seurat_object = SeuratObj, aesX = "nFeature_RNA", fill = "Sample",
  color_by = "Sample", xintercept = 500, scale_x_log10 = TRUE
)

x3 <- plot_densitygraph (
  seurat_object = SeuratObj, aesX = "mito_ratio", fill = "Sample",
  color_by = "Sample", xintercept = 0.015
)

x4 <- plot_densitygraph (
  seurat_object = SeuratObj, aesX = "log10GenesPerUMI", fill = "Sample",
  color_by = "Sample", xintercept = 0.8, scale_x_log10 = TRUE
)

XX <- grid_arrange_shared_legend(
  x1 + theme_min(), 
  x2 + theme_min(), 
  x3 + theme_min(), 
  x4 + theme_min(), 
  position = "right", 
  ncol = 2, nrow = 2
)

#violin plots
X <- VlnPlot(
  SeuratObj, 
  features = c("nFeature_RNA", "nCount_RNA", "mito_ratio"), 
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
  " GEX feature and count.pdf"
), width = 12, height = 10, family = FONT_FAMILY
)
grid.arrange(XX, x5, heights = c(0.8, 0.5))
dev.off()



######## NORMALIZATION (GEX) ########
DefaultAssay(SeuratObj) <- "RNA"

# SCT for GEX
if (SCT) {
  # rename for ease of use; no matter if SCT or NormalizeData,
  # we would like to use "RNA" assay from now on
  SeuratObj <- SCTransform(
    SeuratObj,
    assay = "RNApreSCT",
    new.assay.name = "RNA",
    method = "glmGamPoi", verbose = FALSE
  )

} else {
  SeuratObj <- SeuratObj %>% 
    NormalizeData(assay = "RNA") %>%
    FindVariableFeatures(assay = "RNA") %>%
    ScaleData(assay = "RNA")
}


# plot variable features
top10 <- head(VariableFeatures(SeuratObj), 10, assay = "RNA")
plot1 <- VariableFeaturePlot(SeuratObj)
plot2 <- LabelPoints(plot = plot1, points = top10, size = 3)
pdf(paste0(
  QCDirectory, ObjName, Subset, 
  " GEX variable freatures.pdf"
), width = 8, height = 5.5, family = FONT_FAMILY
)
plot2 + theme_min()
dev.off()



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
  file = paste0(RobjDirectory, ObjName, Subset, "_GEX", ".rds")
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
  " elbow plot after CC regression and harmony (GEX).pdf"
), width = 12, height = 4.5, family = FONT_FAMILY
)
E1
dev.off()


######## LOUVAIN CLUSTERING (GEX) ########
DefaultAssay(SeuratObj) <- "RNA"
analyses <- fromJSON(file = here("analysis.json"))
  
SeuratObj <- SeuratObj %>%
  FindNeighbors(
    reduction = "harmonyRNA", 
    dims = c(1:analyses$harmonyRNA_dims)
  ) %>%
  FindClusters(resolution = RESOLUTION)


CellInfo <- SeuratObj@meta.data
# Rename Clusters
cluster_count <- length(
  levels(as.factor(SeuratObj@meta.data$seurat_clusters))
)
for(j in 1 : cluster_count){
  if (j < 10){
    CellInfo$ClusterRNA[CellInfo$seurat_clusters == j - 1] <- paste0("C0", j)
  }
  else {
    CellInfo$ClusterRNA[CellInfo$seurat_clusters == j - 1] <- paste0("C", j)
  }
}
SeuratObj@meta.data <- CellInfo
Idents(SeuratObj) <- CellInfo$ClusterRNA

# Get number of cells per cluster and per Sample
write.csv(
  as.matrix(table(SeuratObj@meta.data$ClusterRNA, SeuratObj@meta.data$Sample)),
  file = paste0(
    OutputDirectory, ObjName, Subset, 
    " number of cells per cluster and sample.csv"
  )
)

######## RUN UMAP (GEX) ########
SeuratObj <- RunUMAP(
  SeuratObj, 
  reduction = "harmonyRNA", 
  dims = c(1:analyses$harmonyRNA_dims),
  reduction.name = "umapRNA",
  reduction.key = "umapRNA_"
)


################ FIND CLUSTER BIOMARKERS (GEX) ################
# set default assay and identity
DefaultAssay(SeuratObj) <- "RNA"
Idents(SeuratObj) <- CellInfo$ClusterRNA
if (SCT) {
  SeuratObj <- PrepSCTFindMarkers(SeuratObj, assay = "RNA")
}
markers <- FindAllMarkers(
  SeuratObj,
  assay = "RNA",
  logfc.threshold = 0.5, 
  test.use = "wilcox", 
  only.pos = TRUE
)

# save, because it takes a little time to calculate
write.csv(
  markers, 
  paste0(
    OutputDirectory, ObjName, Subset, 
    " RNA cluster markers (by GEX)", "res", RESOLUTION, ".csv"
  )
)



############## SAVE SEURAT AND SESSION INFO, LOOSE ENDS ################
saveRDS(
  SeuratObj, 
  file = paste0(RobjDirectory, ObjName, Subset, "_GEX", ".rds")
)


# capture session info, versions, etc.
writeLines(
  capture.output(sessionInfo()), 
  paste0(ConfigDirectory, ObjName, Subset, "sessionInfo_GEX.txt")
)
file.copy(
  from = here("config.json"), 
  to = paste0(ConfigDirectory, ObjName, Subset, "config_params_GEX.json")
)
file.copy(
  from = here("analysis.json"), 
  to = paste0(ConfigDirectory, ObjName, Subset, "analysis_params_GEX.json")
)
