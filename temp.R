################# LOAD UTILS  ##############
source("~/Box/Yun lab projects/victoria_liu/matching_patients/R_Code/utils.R")

# laod visualizing parameters
analyses <- fromJSON(file = here("analysis.json"))

# this object is fully pre-processed for GEX
RDS_filename <- paste0(
  RobjDirectory, ObjName, Subset, 
  "_res", RESOLUTION, ".rds"
)
if (! file.exists(RDS_filename)) {
  RDS_filename <- paste0(
    RobjDirectory, ObjName, Subset, 
    "_resAll.rds"
  )
}
SeuratObj <- readRDS(RDS_filename)

# if this is a denovo object, it may contain type info
# inherited from its superset parent
if (analyses$denovo) {
  SeuratObj$type <- NULL
  SeuratObj$r <- NULL
}


############### ADD CLUSTER METADATA  #################
CellInfo <- SeuratObj@meta.data
# Rename Clusters
cluster_count <- length(levels(as.factor(
  SeuratObj[[paste0("RNA_snn_res.", RESOLUTION)]][, 1]
)))

for(j in 1 : cluster_count){
  if (j < 10){
    CellInfo$ClusterRNA[
      CellInfo[[paste0("RNA_snn_res.", RESOLUTION)]] == j - 1
    ] <- paste0(analyses$cluster_prefix, "0", j)
  }
  else {
    CellInfo$ClusterRNA[
      CellInfo[[paste0("RNA_snn_res.", RESOLUTION)]] == j - 1
    ] <- paste0(analyses$cluster_prefix, j)
  }
}
SeuratObj@meta.data <- CellInfo
Idents(SeuratObj) <- CellInfo$ClusterRNA


# remove unnecessary metadata; focusing on one resolution at a time
for (metadata_col in colnames(SeuratObj@meta.data)) {
  if (
    grepl("RNA_snn_res.", metadata_col, fixed = TRUE) 
    & !grepl(RESOLUTION, metadata_col, fixed = TRUE)
  ) {
    print (paste0("removing ", metadata_col, " from metadata"))
    SeuratObj[[metadata_col]] <- NULL
  }
}


# Get number of cells per cluster and per Sample
write.csv(
  as.matrix(table(SeuratObj@meta.data$ClusterRNA, SeuratObj@meta.data$Sample)),
  file = paste0(
    OutputDirectory, ObjName, Subset, 
    " number of cells per cluster and sample.csv"
  )
)

########## SAVE LOUPE PROJECTIONS ##########
write.csv(
  Embeddings(SeuratObj, paste0("umap", analyses$viz_clustering)), 
  paste0(
    LoupeDirectory, ObjName, Subset, "_umap", analyses$viz_clustering, ".csv"
  )
)


write.csv(
  select(SeuratObj@meta.data, Sample), 
  paste0(LoupeDirectory, ObjName, Subset, "_samples.csv")
)

write.csv(
  select(SeuratObj@meta.data, ClusterRNA), 
  paste0(LoupeDirectory, ObjName, Subset, "_clusters.csv")
)

############## FIND CLUSTER BIOMARKERS (GEX) ################
# set default assay and identity
markers <- GEX_cluster_markers(SeuratObj)
# save, because it takes a little time to calculate
write.csv(
  markers, 
  paste0(
    OutputDirectory, ObjName, Subset, 
    " RNA cluster markers (by RNA)", "res", RESOLUTION, ".csv"
  )
)

############### CLUSTER / SAMPLE HEATMAP ##################

HM_object <- plot_heatmap (
  seurat_object = SeuratObj, downsample_n = 5000,
  markers = markers, top_n = 20, label_n = 2, 
  cluster = paste0("Cluster", analyses$viz_clustering),
  data_type = "logcounts",
  use_raster = TRUE
)

pdf(paste0(
  heatDirectory, "heatmap", ObjName, Subset, 
  "_res", RESOLUTION, "_top20 genes per ", 
  analyses$viz_clustering, "Cluster.pdf"
), width = 7, height = 6
)
draw(
  HM_object[[1]], 
  heatmap_legend_list = list(HM_object[[2]], HM_object[[3]]),
  heatmap_legend_side = "right"
)
dev.off()



################# CLUSTER / SAMPLE BARPLOTS ###################
P2 <- plot_bargraph (
  seurat_object = SeuratObj, aesX = "Sample", 
  fill = paste0("Cluster", analyses$viz_clustering), 
  y_label = "Composition (Number of cells)", x_label = NULL, 
  y_lower_limit = 0, y_break = 1000,
  color_reverse = TRUE
)

pdf(paste0(
  densityplotDirectory, ObjName, Subset, "res", RESOLUTION, 
  "_number of cells per sample and ", analyses$viz_clustering,
  " Cluster barplot.pdf"
),
width = 6, height = 5.5, family = FONT_FAMILY
)
P2
dev.off()



P3 <- plot_bargraph (
  seurat_object = SeuratObj, 
  aesX = paste0("Cluster", analyses$viz_clustering), fill = "Sample",
  y_label = "Composition (percentage of cells)", x_label = NULL, 
  y_lower_limit = 0, y_break = 0.2,
  position = "fill",
  plot_margin = unit(c(0.2, 0.5, 0.2, 0.5), "cm")
)

pdf(paste0(
  densityplotDirectory, ObjName, Subset, 
  "res", RESOLUTION, 
  "_percent of cells per sample and ", analyses$viz_clustering,
  " Cluster barplot.pdf"
), width = 6, height = 5.5, family = FONT_FAMILY
)
P3
dev.off()



P4 <- plot_bargraph (
  seurat_object = SeuratObj, 
  aesX = paste0("Cluster", analyses$viz_clustering), fill = "Sample", 
  y_label = "Composition (Number of cells)", x_label = NULL, 
  y_lower_limit = 0, y_break = 1000
)

pdf(paste0(
  densityplotDirectory, ObjName, Subset, 
  "res", RESOLUTION, "_number of cells per ", analyses$viz_clustering,
  " Cluster and sample barplot.pdf"
), width = 6, height = 5.5, family = FONT_FAMILY
)
P4
dev.off()


################# CLUSTER / SAMPLE DIMPLOTS ###################
U1 <- plot_umap(
  seurat_object = SeuratObj, 
  group_by = paste0("Cluster", analyses$viz_clustering),
  reduction = paste0("umap", analyses$viz_clustering),
  title = "Clusters", xlab = "UMAP1", ylab = "UMAP2",
  legend_position = "bottom",
  ncol_guide = 5,
  color_reverse = TRUE, label_clusters = TRUE
)

pdf(paste0(
  UMAPDirectory, ObjName, Subset,
  "_res", RESOLUTION, 
  "_", analyses$viz_clustering, "Clusters_UMAP.pdf"
), width = 5.5, height = 6, family = FONT_FAMILY
)
U1
dev.off()



U2 <- plot_umap(
  seurat_object = SeuratObj, 
  group_by = "Sample",
  reduction = paste0("umap", analyses$viz_clustering),
  title = "Samples", xlab = "UMAP1", ylab = "UMAP2",
  legend_position = "bottom",
  ncol_guide = 3
)


pdf(paste0(
  UMAPDirectory, ObjName, Subset,
  "_res", RESOLUTION, 
  "_Samples_UMAP.pdf"
), width = 5.5, height = 6, family = FONT_FAMILY
)
U2
dev.off()



U3 <- plot_umap (
  seurat_object = SeuratObj, 
  group_by = paste0("Cluster", analyses$viz_clustering),
  reduction = paste0("umap", analyses$viz_clustering),
  title = "Clusters", xlab = "UMAP1", ylab = "UMAP2",
  legend_position = "right",
  title_font_size = 16, x_font_size = 16, y_font_size = 16, 
  pt_size = 0.2, split_by = "Sample", ncol_dimplot = 3,
  color_reverse = TRUE
)

pdf(paste0(
  UMAPDirectory, ObjName, Subset,
  "_res", RESOLUTION, "_",
  analyses$viz_clustering, "Clusters_UMAP_Iteration_by_sample.pdf"
), width = 12, height = 8, family = FONT_FAMILY
)
U3
dev.off()


################# (CUSTOM) CELL POPULATION: CLUSTIFYR ###################
REF_MATRIX = cbmc_ref

# clustify using 
correlation_matrix <- clustify(
  input = SeuratObj[[analyses$clustifyr_assay]]@data, 
  metadata = SeuratObj@meta.data,
  cluster_col = paste0("Cluster", analyses$viz_clustering), 
  ref_mat = REF_MATRIX,
  # query_genes = FindVariableFeatures(
  #     SeuratObj, assay = "RNApreSCT"
  #   )[["RNApreSCT"]]@var.features[1:500]
)

# predicted type with correlation coefficients
correlation_coefficients <- cor_to_call(
  cor_mat = correlation_matrix,
  cluster_col = paste0("Cluster", analyses$viz_clustering)
)

# plot heatmap
pdf(paste0(
  heatDirectory, "heatmap", ObjName, Subset, 
  "_res", RESOLUTION, "_cellIdentities_", 
  analyses$viz_clustering, "Cluster.pdf"
), width = 7, height = 6
)
heatmap.2(
  correlation_matrix, 
  col=viridis, trace = "none", 
  dendrogram = "none", 
  offsetRow= -29, margins = c(8, 5)
)
dev.off()

# add metadata to SeuratObj
SeuratObj@meta.data <- call_to_metadata(
  res = correlation_coefficients, 
  metadata = SeuratObj@meta.data,
  cluster_col = "ClusterRNA"
)

SeuratObj$Assignment <- SeuratObj$type
SeuratObj$type <- NULL


################# (CUSTOM) CELL POPULATIONS: ASSIGNMENT ###############


## MANUAL (i.e. when clustifyr is not sufficient)
## NEEDS A COMPLETE REVISION


# add loupe projection
write.csv(
  select(SeuratObj@meta.data, Assignment), 
  paste0(LoupeDirectory, ObjName, Subset, "_assignments.csv")
)

################# ASSIGNMENT BARGRAPHS ###################
P5 <- plot_bargraph (
  seurat_object = SeuratObj, aesX = "Sample", fill = "Assignment",
  y_label = "Composition (Number of cells)", x_label = NULL,
  y_lower_limit = 0, y_break = 1000,
)

pdf(paste0(
  densityplotDirectory, ObjName, Subset, 
  "res", RESOLUTION,
  "_number of cells per sample and Assignment barplot.pdf"
), width = 6, height = 5.5, family = FONT_FAMILY
)
P5
dev.off()



P6 <- plot_bargraph (
  seurat_object = SeuratObj, aesX = "Assignment", fill = "Sample",
  y_label = "Composition (percentage of cells)", x_label = NULL, 
  y_lower_limit = 0, y_break = 0.2,
  position = "fill",
  plot_margin = unit(c(0.2, 0.5, 0.2, 0.5), "cm")
)
pdf(paste0(
  densityplotDirectory, ObjName, Subset, 
  "res", RESOLUTION, 
  "percent of cells per sample and Assignment barplot.pdf"
), width = 6, height = 5.5, family = FONT_FAMILY
)
P6
dev.off()



P7 <- plot_bargraph (
  seurat_object = SeuratObj, aesX = "Assignment", fill = "Sample",
  y_label = "Composition (Number of cells)", x_label = NULL, 
  y_lower_limit = 0, y_break = 1000
)

pdf(paste0(
  densityplotDirectory, ObjName, Subset,
  "res", RESOLUTION, 
  "_percentage of cells per major population and Sample barplot.pdf"
), width = 7, height = 5.5, family = FONT_FAMILY
)
P7
dev.off()

P8 <- plot_bargraph (
  seurat_object = SeuratObj, aesX = "Sample", fill = "Assignment",
  y_label = "Composition (percentage of cells)", x_label = NULL,
  y_lower_limit = 0, y_break = 1000, position = "fill"
)

pdf(paste0(
  densityplotDirectory, ObjName, Subset, 
  "res", RESOLUTION,
  "_percentage of cells per sample and Assignment barplot.pdf"
), width = 6, height = 5.5, family = FONT_FAMILY
)
P8
dev.off()

plots = ggarrange(P2, P5, P8, P3, P4, P6, P7, ncol = 2)
pdf(paste0(
  densityplotDirectory, ObjName, Subset, 
  "res", RESOLUTION, "_all barplots.pdf"
), width = 18, height = 22, family = FONT_FAMILY
)
print(plots)
dev.off()


################# ASSIGNMENT DIMPLOTS ###################
U4 <- plot_umap(
  seurat_object = SeuratObj, group_by = "Assignment",
  reduction = paste0("umap", analyses$viz_clustering),
  title = "Assignment", xlab = "UMAP1", ylab = "UMAP2",
  legend_position = "bottom",
  ncol_guide = 4,
  label_clusters = TRUE
)

pdf(paste0(
  UMAPDirectory, ObjName, Subset, 
  "_res", RESOLUTION, 
  "_", analyses$viz_clustering,
  "Clusters_Assignment UMAP.pdf"
), width = 5.5, height = 6, family = FONT_FAMILY
)
U4
dev.off()

plots2 <- ggarrange(U1, U2, U4, ncol = 3)

pdf(paste0(
  UMAPDirectory, ObjName, Subset, 
  "_res", RESOLUTION, "_all UMAPs.pdf"
), width = 15, height = 6, family = FONT_FAMILY
)
print(plots2)
dev.off()



U5 <- plot_umap(
  seurat_object = SeuratObj, group_by = "Assignment",
  reduction = paste0("umap", analyses$viz_clustering),
  title = "Assignment", xlab = "UMAP1", ylab = "UMAP2",
  legend_position = "bottom",
  title_font_size = 16, x_font_size = 16, y_font_size = 16, 
  pt_size = 0.2, split_by = "Sample", ncol_dimplot = 2
)

pdf(paste0(
  UMAPDirectory, ObjName, Subset, 
  "_res", RESOLUTION, 
  "_", analyses$viz_clustering,
  "Clusters_Assignment UMAP Iteration by sample.pdf"
), width = 12, height = 12, family = FONT_FAMILY
)
U5
dev.off()


######### RIBO RATIO QC ##########
# umap
U0 <- plot_umap (
  seurat_object = SeuratObj, 
  group_by = paste0("Cluster", analyses$viz_clustering),
  reduction = paste0("umap", analyses$viz_clustering),
  title = "Ribosomal Ratios", xlab = "UMAP1", ylab = "UMAP2",
  legend_position = "right",
  title_font_size = 16, x_font_size = 16, y_font_size = 16, 
  pt_size = 0.2, split_by = "ribo_category", ncol_dimplot = 3,
  color_reverse = TRUE, label_clusters = TRUE
)

pdf(paste0(
  RiboQCDirectory, ObjName, Subset,
  "_res", RESOLUTION, 
  "_", analyses$viz_clustering, "clusters_ribo_ratio_UMAP.pdf"
), width = 12, height = 6, family = FONT_FAMILY
)
U0
dev.off()

# violin
#violin plots
V1 <- VlnPlot(
  SeuratObj, 
  features = "ribo_ratio", 
  pt.size = 0.01, 
  cols = get_colors(
    seurat_object = SeuratObj, 
    color_by = paste0("Cluster", analyses$viz_clustering)
  ), 
  group.by = paste0("Cluster", analyses$viz_clustering)
) +
  labs(title = paste0("Cluster", analyses$viz_clustering))

V2 <- VlnPlot(
  SeuratObj, 
  features = "ribo_ratio", 
  pt.size = 0.01, 
  cols = get_colors(
    seurat_object = SeuratObj, 
    color_by = "Assignment"
  ), 
  group.by = "Assignment"
) +
  labs(title = "Assignment")

V3 <- VlnPlot(
  SeuratObj, 
  features = "ribo_ratio", 
  pt.size = 0.01, 
  cols = get_colors(
    seurat_object = SeuratObj, 
    color_by = "Sample"
  ), 
  group.by = "Sample"
) +
  labs(title = "Sample")

V_combined <- ggarrange(
  V1 + theme_min() + NoLegend() + RotatedAxis(), 
  V2 + theme_min() + NoLegend() + RotatedAxis(), 
  V3 + theme_min() + NoLegend() + RotatedAxis(), 
  nrow = 1
)

pdf(paste0(
  RiboQCDirectory, ObjName, Subset,
  "_res", RESOLUTION, 
  "_", analyses$viz_clustering, "clusters_ribo_ratio_violin_all.pdf"
), width = 14, height = 6, family = FONT_FAMILY
)
V_combined
dev.off()



x1 <- plot_densitygraph (
  seurat_object = SeuratObj, aesX = "ribo_ratio", 
  fill = paste0("Cluster", analyses$viz_clustering),
  color_by = paste0("Cluster", analyses$viz_clustering), 
  facet_category = paste0("Cluster", analyses$viz_clustering),
  xintercept = config$ribo_ratio,
  title = paste0("Ribosome ratio by ", "Cluster", analyses$viz_clustering)
)

x2 <- plot_densitygraph (
  seurat_object = SeuratObj, aesX = "ribo_ratio", 
  fill = "Assignment",
  color_by = "Assignment", 
  facet_category = "Assignment",
  xintercept = config$ribo_ratio,
  title = "Ribosome ratio by Assignment"
)

x3 <- plot_densitygraph (
  seurat_object = SeuratObj, aesX = "ribo_ratio", 
  fill = "Sample",
  color_by = "Sample", 
  facet_category = "Sample",
  xintercept = config$ribo_ratio,
  title = "Ribosome ratio by Sample"
)

x_combined <- ggarrange(
  x1 + theme_min(), 
  x2 + theme_min(), 
  x3 + theme_min(), 
  nrow = 3
)

pdf(paste0(
  RiboQCDirectory, ObjName, Subset,
  "_res", RESOLUTION, 
  "_", analyses$viz_clustering, 
  "clusters ribosome ratios split by cluster assignment sample.pdf"
), width = 20, height = 20, family = FONT_FAMILY
)
x_combined
dev.off()



################# CELL CYCLE SORTING ###################

U6 <- plot_umap (
  seurat_object = SeuratObj, group_by = "Phase",
  reduction = paste0("umap", analyses$viz_clustering),
  title = "Seurat Cell Cycle scoring", xlab = "UMAP1", ylab = "UMAP2",
  legend_position = "right",
  title_font_size = 16, x_font_size = 16, y_font_size = 16, 
  pt_size = 0.2, split_by = "Sample", ncol_dimplot = 2,
  label_clusters = TRUE, repel_labels = TRUE, label_size = 3
)

pdf(paste0(
  UMAPDirectory, ObjName, Subset, 
  "CellCycle UMAP by Iteration Sample.pdf"
), width = 12, height = 12, family = FONT_FAMILY
)
U6
dev.off()


############## SAVE SEURAT AND SESSION INFO, LOOSE ENDS ################
saveRDS(
  SeuratObj, 
  file = paste0(
    RobjDirectory, ObjName, Subset, 
    "_res", RESOLUTION, ".rds"
  )
)
