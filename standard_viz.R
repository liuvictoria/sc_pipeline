################# LOAD UTILS  ##############
source("~/Box/Yun lab projects/victoria_liu/matching_patients/R_Code/utils.R")

# laod visualizing parameters
analyses <- fromJSON(file = here("analysis.json"))

# this object is fully pre-processed for GEX
SeuratObj <- readRDS(
  paste0(
    RobjDirectory, ObjName, Subset, 
    "_res", config$RESOLUTION, ".rds"
    )
)

# contains results of FindAllMarkers
markers <- read.csv(
  paste0(
    OutputDirectory, ObjName, Subset, " ",
    analyses$viz_assay," cluster markers (by ", analyses$viz_clustering, ")", 
    "res", RESOLUTION, ".csv"
  )
)
markers <- markers[order(markers$cluster), , drop = F]


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
REF_MARKERS = cbmc_m
# determine which reference matrix / gene list to use

# clustify using predefined gene lists
# only supports using jaccard metric
# r values are weird; try to avoid using lists if possible
list_res <- clustify_lists(
  input = SeuratObj[[analyses$clustifyr_assay]]@data,
  metadata = SeuratObj@meta.data,
  cluster_col = paste0("Cluster", analyses$viz_clustering),
  marker = REF_MARKERS,
  metric = "jaccard"
)

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
SeuratObj$type <- NA

################# (CUSTOM) CELL POPULATIONS: DOTPLOTS ###################
# loads Clusterspecificgenes
load(paste0(here(), "/", "Clusterspecificgenes"))

dotgraphs <- list()
for (cluster_name in names(Clusterspecificgenes)) {
  cluster_genes <- Clusterspecificgenes[[cluster_name]]
  D1 <- plot_dotgraph(
    seurat_object = SeuratObj,
    group_by = paste0("Cluster", analyses$viz_clustering),
    features = cluster_genes, title = cluster_name
  )
  dotgraphs[[cluster_name]] <- D1
}

plots <- ggarrange(plots = dotgraphs, nrow = length(Clusterspecificgenes))


pdf(paste0(
  OutputDirectory, "dotplot ", ObjName, Subset,
  "by predefined Cluster.pdf"
), width = 18, height = 18, family = FONT_FAMILY
)
print(plots)
dev.off()



D2 <- plot_dotgraph(
  seurat_object = SeuratObj, 
  paste0(group_by = "Cluster", analyses$viz_clustering),
  features = unique(get_top_cluster_markers(markers, 2)$gene),
  title = "Top 2 genes by cluster"
)


pdf(paste0(
  dotDirectory, "dotplot ", ObjName, Subset, "top 2 genes by Cluster.pdf"
), width = 10, height = 5
)
D2
dev.off()

################# (CUSTOM) CELL POPULATIONS: ASSIGNMENT ###############


## MANUAL (i.e. when clustifyr is not sufficient)
## NEEDS A COMPLETE REVISION


# add loupe projection
write.csv(
  select(SeuratObj@meta.data, Assignment), 
  paste0(LoupeDirectory, Subset, "_assignments.csv")
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

plots = ggarrange(P2, P3, P4, P5, P6, P7, ncol = 2)
pdf(paste0(
  densityplotDirectory, ObjName, Subset, 
  "res", RESOLUTION, "_all barplots.pdf"
), width = 18, height = 18, family = FONT_FAMILY
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


########### PREDEFINED CLUSTER FEATURE PLOTS #############
for (cluster_name in names(Clusterspecificgenes)) {
  predefined_cluster_plots <- list()
  cluster_genes <- Clusterspecificgenes[[cluster_name]]
  
  for (gene in cluster_genes){
    if (gene %in% rownames(SeuratObj)) {
      F1 <- plot_featureplot (
        seurat_object = SeuratObj, feature_gene = gene, split_by = "Group"
      )
      predefined_cluster_plots[[gene]] <- F1
    }
  }
  
  plots <- ggarrange(
    plots = predefined_cluster_plots, 
    ncol = 2,
    nrow = length(cluster_genes) %/% 2 + length(cluster_genes) %% 2
  )
  
  pdf(paste0(
    featuremapDirectory, cluster_name, 
    " featuremap ", ObjName, " ", Subset, ".pdf"
  ), width = 26, height = 46
  )
  print(plots)
  dev.off()
  
}


################# CELL CYCLE SORTING ###################
RidgePlot(SeuratObj, features = c("Pcna", "Top2a", "Mcm6", "Mki67"), ncol = 2)
DimPlot(SeuratObj, reduction = "umap", group.by = "Phase")


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
    "_res", config$RESOLUTION, ".rds"
    )
)
