################# LOAD UTILS ##############
source("~/Box/Yun lab projects/victoria_liu/matching_patients/R_Code/utils.R")

# load visualizing parameters
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
  if (! file.exists(RDS_filename)) {
    print ("no Seurat object to analyze")
  }
}
print(paste0("reading SeuratObj from file ", RDS_filename))
SeuratObj <- readRDS(RDS_filename)

# if this is a denovo object, it may contain type info
# inherited from its superset parent
# remove, so that clustifyr has space to add new info, when called
if (analyses$denovo) {
  SeuratObj$r <- NULL
  SeuratObj$type <- NULL
}

# loads Clusterspecificgenes
Clusterspecificgenes <- analyses[[
  paste0(analyses[["denovo_lineage"]], "_lineage_markers")
]]

############### ADD CLUSTER METADATA  #################
# i.e. cluster_res_column = "RNA_snn_res.0.5"
cluster_res_column <- ifelse(
  analyses$viz_clustering == "RNA",
  paste0("RNA_snn_res.", RESOLUTION),
  paste0("wsnn_res.", RESOLUTION)
)
# i.e. cluster_name = "ClusterRNA"
cluster_name <- ifelse(
  analyses$viz_clustering == "RNA",
  "ClusterRNA",
  "ClusterWNN"
)
CellInfo <- SeuratObj@meta.data

# Rename Clusters
cluster_count <- length(levels(as.factor(
  SeuratObj[[cluster_res_column]][, 1]
  )))

for(j in 1 : cluster_count){
  if (j < 10){
    CellInfo[[cluster_name]][
      CellInfo[[cluster_res_column]] == j - 1
      ] <- paste0(analyses$cluster_prefix, "0", j)
  }
  else {
    CellInfo[[cluster_name]][
      CellInfo[[cluster_res_column]] == j - 1
      ] <- paste0(analyses$cluster_prefix, j)
  }
}
SeuratObj@meta.data <- CellInfo
Idents(SeuratObj) <- CellInfo[[cluster_name]]


# remove unnecessary metadata; focusing on one resolution at a time
for (metadata_col in colnames(SeuratObj@meta.data)) {
  if (
    grepl("snn_res.", metadata_col, fixed = TRUE) & 
    ! grepl(RESOLUTION, metadata_col, fixed = TRUE)
    ) {
    print (paste0("removing ", metadata_col, " from metadata"))
    SeuratObj[[metadata_col]] <- NULL
  }
}


# Get number of cells per cluster and per Sample
write.csv(
  as.matrix(table(
    SeuratObj@meta.data[[cluster_name]], SeuratObj@meta.data$Sample
  )),
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
  SeuratObj$Sample, 
  paste0(LoupeDirectory, ObjName, Subset, "_samples.csv")
)

write.csv(
  SeuratObj[[cluster_name]], 
  paste0(LoupeDirectory, ObjName, Subset, "_clusters.csv")
)

############## FIND CLUSTER BIOMARKERS ################
if (analyses$viz_clustering == "RNA") {
  markers_filename <- paste0(
    OutputDirectory, ObjName, Subset, 
    " RNA cluster markers (by RNA)", "res", RESOLUTION, ".csv"
    )
  markers <- cluster_markers(
    SeuratObj, assay = "RNA", default_ident = "ClusterRNA"
    )
  # save, because it takes a little time to calculate
  write.csv(markers, markers_filename)

  
  
} else if (analyses$viz_clustering == "WNN") {
  markersRNA_filename <- paste0(
    OutputDirectory, ObjName, Subset, 
    " RNA cluster markers (by WNN)", "res", RESOLUTION, ".csv"
  )
  markersRNA <- cluster_markers(
    SeuratObj, assay = "RNA", default_ident = "ClusterWNN"
  )
  # save, because it takes a little time to calculate
  write.csv(markersRNA, markersRNA_filename)
  
  
  markersADT_filename <- paste0(
    OutputDirectory, ObjName, Subset, 
    " ADT cluster markers (by WNN)", "res", RESOLUTION, ".csv"
  )
  markersADT <- cluster_markers(
    SeuratObj, assay = "ADTpreInt", default_ident = "ClusterWNN"
  )
  # save, because it takes a little time to calculate
  write.csv(markersADT, markersADT_filename)
  markers <- markersRNA
}

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
  "_res", RESOLUTION, "_top20 ", analyses$markers_assay, " features per ", 
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
), width = 8, height = 5.5, family = FONT_FAMILY
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
), width = 8, height = 5.5, family = FONT_FAMILY
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
), width = 8, height = 7, family = FONT_FAMILY
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
  pt_size = 0.2, split_by = "Sample", ncol_dimplot = 2,
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
# determine which reference matrix / gene list to use
REF_MATRIX = cbmc_ref

# refSeuratObj <- readRDS(
#   paste0(
#     RobjDir,
#     "GBMAtlas/",
#     "Allhuman-11-3-21.rds"
#   )
# )
# 
# # if column with assignment has space in it
# REF_MATRIX <- seurat_ref(
#   seurat_object = refSeuratObj,
#   cluster_col = "Assignment"
# )


# clustify using
correlation_matrix <- clustify(
  input = SeuratObj[[analyses$clustifyr_assay]]@data,
  metadata = SeuratObj@meta.data,
  cluster_col = paste0("Cluster", analyses$viz_clustering),
  ref_mat = REF_MATRIX,
  query_genes = FindVariableFeatures(
      SeuratObj, assay = "RNApreSCT"
    )[["RNApreSCT"]]@var.features
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
  cluster_col = paste0("Cluster", analyses$viz_clustering)
)

SeuratObj$clustifyr <- SeuratObj$type
SeuratObj$clustifyr_r <- SeuratObj$r
SeuratObj$type <- NULL
SeuratObj$r <- NULL



################# CELL POPULATIONS: ASSIGNMENT ###############
if (
  ! is.na(analyses$which_assignment) &
  analyses$which_assignment %in% colnames(SeuratObj@meta.data)
) {
  SeuratObj$Assignment <- SeuratObj[[analyses$which_assignment]]
}

# for manual hand correction
for (new_assignment in names(analyses$cluster_to_assignment)) {
  cluster = paste0(analyses$cluster_prefix, new_assignment)
  SeuratObj$Assignment[
    which(str_detect(
      SeuratObj[[paste0("Cluster", analyses$viz_clustering)]][, 1], cluster
      ))
    ] <- analyses$cluster_to_assignment[[new_assignment]]
}


# add loupe projection
write.csv(
  dplyr::select(SeuratObj@meta.data, Assignment),
  paste0(
    LoupeDirectory, ObjName, Subset,
    "_assignments (", analyses$which_assignment, ") .csv")
)

################# ASSIGNMENT BARPLOTS ###################
P5 <- plot_bargraph (
  seurat_object = SeuratObj, aesX = "Sample", fill = "Assignment",
  y_label = "Composition (Number of cells)", x_label = NULL,
  y_lower_limit = 0, y_break = 1000, 
  title = paste0(analyses$which_assignment, " Assignment")
)

pdf(paste0(
  densityplotDirectory, ObjName, Subset, 
  "res", RESOLUTION,
  "_number of cells per sample and Assignment (", 
  analyses$which_assignment,
  ") barplot.pdf"
), width = 6, height = 5.5, family = FONT_FAMILY
)
P5
dev.off()



P6 <- plot_bargraph (
  seurat_object = SeuratObj, aesX = "Assignment", fill = "Sample",
  y_label = "Composition (percentage of cells)", x_label = NULL, 
  y_lower_limit = 0, y_break = 0.2,
  position = "fill",
  plot_margin = unit(c(0.2, 0.5, 0.2, 0.5), "cm"),
  title = paste0(analyses$which_assignment, " Assignment")
)
pdf(paste0(
  densityplotDirectory, ObjName, Subset, 
  "res", RESOLUTION, 
  "percent of cells per sample and Assignment (", 
  analyses$which_assignment,
  ") barplot.pdf"
), width = 6, height = 5.5, family = FONT_FAMILY
)
P6
dev.off()


P7 <- plot_bargraph (
  seurat_object = SeuratObj, aesX = "Assignment", fill = "Sample",
  y_label = "Composition (Number of cells)", x_label = NULL, 
  y_lower_limit = 0, y_break = 1000,
  title = paste0(analyses$which_assignment, " Assignment")
)

pdf(paste0(
  densityplotDirectory, ObjName, Subset,
  "res", RESOLUTION, 
  "_percentage of cells per major population and Sample (", 
  analyses$which_assignment,
  ") barplot.pdf"
), width = 7, height = 5.5, family = FONT_FAMILY
)
P7
dev.off()

P8 <- plot_bargraph (
  seurat_object = SeuratObj, aesX = "Sample", fill = "Assignment",
  y_label = "Composition (percentage of cells)", x_label = NULL,
  y_lower_limit = 0, y_break = 1000, position = "fill",
  title = paste0(analyses$which_assignment, " Assignment")
)

pdf(paste0(
  densityplotDirectory, ObjName, Subset, 
  "res", RESOLUTION,
  "_percentage of cells per sample and Assignment (", 
  analyses$which_assignment,
  ") barplot.pdf"
), width = 6, height = 5.5, family = FONT_FAMILY
)
P8
dev.off()




# multiple
bar_plots <- ggarrange(P2, P5, P8, P3, P4, P6, P7, ncol = 2)
pdf(paste0(
  densityplotDirectory, ObjName, Subset, 
  "res", RESOLUTION, "_all barplots (Assignment ", 
  analyses$which_assignment,
  ").pdf"
), width = 18, height = 22, family = FONT_FAMILY
)
print(bar_plots)
dev.off()


################# ASSIGNMENT DIMPLOTS ###################
U4 <- plot_umap(
  seurat_object = SeuratObj, group_by = "Assignment",
  reduction = paste0("umap", analyses$viz_clustering),
  title = paste0(analyses$which_assignment, " Assignment"), 
  xlab = "UMAP1", ylab = "UMAP2",
  legend_position = "bottom",
  ncol_guide = 4,
  label_clusters = TRUE,
  label_size = 3,
  color_reverse = FALSE
)

pdf(paste0(
  UMAPDirectory, ObjName, Subset, 
  "_res", RESOLUTION, 
  "_", analyses$viz_clustering,
  "Clusters_Assignment (", 
  analyses$which_assignment,
  ") UMAP.pdf"
), width = 8, height = 6, family = FONT_FAMILY
)
U4
dev.off()

U5 <- plot_umap(
  seurat_object = SeuratObj, group_by = "Assignment",
  reduction = paste0("umap", analyses$viz_clustering),
  title = paste0(analyses$which_assignment, " Assignment"),
  xlab = "UMAP1", ylab = "UMAP2",
  legend_position = "bottom",
  title_font_size = 16, x_font_size = 16, y_font_size = 16, 
  pt_size = 0.2, split_by = "Sample", ncol_dimplot = 2
)

pdf(paste0(
  UMAPDirectory, ObjName, Subset, 
  "_res", RESOLUTION, 
  "_", analyses$viz_clustering,
  "Clusters_Assignment (", 
  analyses$which_assignment,
  ") UMAP Iteration by sample.pdf"
), width = 12, height = 12, family = FONT_FAMILY
)
U5
dev.off()


U6 <- plot_umap(
  seurat_object = SeuratObj, group_by = "Assignment",
  reduction = paste0("umap", analyses$viz_clustering),
  title = paste0(analyses$which_assignment, " Assignment"),
  xlab = "UMAP1", ylab = "UMAP2",
  legend_position = "bottom",
  title_font_size = 16, x_font_size = 16, y_font_size = 16, 
  pt_size = 0.5, split_by = "Assignment", ncol_dimplot = 2
)

pdf(paste0(
  UMAPDirectory, ObjName, Subset, 
  "_res", RESOLUTION, 
  "_", analyses$viz_clustering,
  "Clusters_Assignment (", 
  analyses$which_assignment,
  ") UMAP Iteration by assignment.pdf"
), width = 12, height = 30, family = FONT_FAMILY
)
U6
dev.off()




# multiple
plots2 <- ggarrange(U1, U4, U2, ncol = 3)

pdf(paste0(
  UMAPDirectory, ObjName, Subset, 
  "_res", RESOLUTION, "(Assignment ", 
  analyses$which_assignment,
  ") all UMAPs.pdf"
), width = 20, height = 7, family = FONT_FAMILY
)
print(plots2)
dev.off()



################# CELL POPULATIONS: DOTPLOTS ###################
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

dot_plots <- ggarrange(plots = dotgraphs, nrow = length(Clusterspecificgenes))


pdf(paste0(
  dotDirectory, "dotplot ", ObjName, Subset,
  "by predefined Cluster.pdf"
), width = 10, height = 56, family = FONT_FAMILY
)
dot_plots
dev.off()



D2 <- plot_dotgraph(
  seurat_object = SeuratObj,
  group_by = paste0("Cluster", analyses$viz_clustering),
  features = unique(get_top_cluster_markers(markers, 5)$gene),
  title = "Top 5 genes by cluster"
)


pdf(paste0(
  dotDirectory, "dotplot ", ObjName, Subset, "top 5 genes by Cluster.pdf"
), width = 22, height = 8
)
D2
dev.off()

########### (GEX) PREDEFINED CLUSTER FEATURE PLOTS #############
# loop each assignment
for (cluster_name in names(Clusterspecificgenes)) {
  predefined_cluster_plots_GEX <- list()
  cluster_genes <- Clusterspecificgenes[[cluster_name]]
  # loop each marker gene in assignment
  for (gene in cluster_genes){
    print(gene)
    F1 <- plot_featureplot (
      seurat_object = SeuratObj,
      feature_gene = gene,
      split_by = NULL,
      reduction = paste0("umap", analyses$viz_clustering),
      pt_size = 0.6
    )
    predefined_cluster_plots_GEX[[gene]] <- F1[[2]]
  }


  feature_plots <- ggpubr::ggarrange(
    plotlist = predefined_cluster_plots_GEX,
    ncol = 2,
    nrow = length(cluster_genes) %/% 2 + length(cluster_genes) %% 2
  )

  pdf(paste0(
    featuremapDirectory, cluster_name,
    " featuremap umap_", analyses$viz_clustering, " ", 
    ObjName, " ", Subset, ".pdf"
  ), width = 26, height = ceiling(length(cluster_genes) / 2) * 10
  )
  print(feature_plots)
  dev.off()

}

#################### (ADT) FEATURE PLOTS #######################
# loop each assignment
if (USE_ADT & analyses$viz_clustering == "WNN") {
# loop each assignment
for (cluster_name in names(analyses$ADT_markers)) {
  predefined_cluster_plots_ADT_GEX <- list()
  cluster_features <- analyses$ADT_markers[[cluster_name]]

  # loop each marker gene in assignment
  for (ADT_feature in cluster_features){
    print(ADT_feature)
    F1 <- plot_featureplot (
      seurat_object = SeuratObj,
      feature_gene = ADT_feature,
      assay = "ADT",
      split_by = NULL,
      reduction = paste0("umap", analyses$viz_clustering),
      pt_size = 0.6
    )

    RNA_feature <- analyses$ADT_to_RNA[[ADT_feature]]
    print(RNA_feature)
    F2 <- plot_featureplot (
      seurat_object = SeuratObj,
      feature_gene = RNA_feature,
      assay = "RNA",
      split_by = NULL,
      reduction = paste0("umap", analyses$viz_clustering),
      pt_size = 0.6
    )
    predefined_cluster_plots_ADT_GEX[[ADT_feature]] <- F1[[2]]
    predefined_cluster_plots_ADT_GEX[[RNA_feature]] <- F2[[2]]
  }


  feature_plots <- ggpubr::ggarrange(
    plotlist = predefined_cluster_plots_ADT_GEX,
    ncol = 2,
    nrow = length(cluster_genes) %/% 2 + length(cluster_genes) %% 2
  )

  pdf(paste0(
    ADTDirectory, cluster_name,
    " featuremap umap_", analyses$viz_clustering, " ",
    ObjName, " ", Subset, ".pdf"
  ), width = 26, height = ceiling(length(cluster_genes) / 2) * 10
  )
  print(feature_plots)
  dev.off()

}

}
#################### (ADT) CORRELATION PLOTS ###################
if (USE_ADT & analyses$viz_clustering == "WNN") {
# loop each assignment
for (cluster_name in names(analyses$ADT_markers)) {
  predefined_cluster_plots_correlation <- list()
  cluster_features <- analyses$ADT_markers[[cluster_name]]

  # loop each marker gene in assignment
  for (ADT_feature in cluster_features){
    RNA_feature <- analyses$ADT_to_RNA[[ADT_feature]]
    print(paste0(RNA_feature, ", ", ADT_feature))

    # make sure x is RNA_feature, y is ADT_feature
    sp1 <- plot_correlation (
      SeuratObj,
      x = RNA_feature, y = ADT_feature, split_by = "Assignment",
      lab_x = RNA_feature, lab_y = ADT_feature,
      method = "pearson"
    )

    predefined_cluster_plots_correlation[[ADT_feature]] <- sp1
  }


  feature_plots <- ggpubr::ggarrange(
    plotlist = predefined_cluster_plots_correlation,
    ncol = 2,
    nrow = length(cluster_genes) %/% 2 + length(cluster_genes) %% 2
  )

  pdf(paste0(
    ADTDirectory, cluster_name,
    " ADT RNA correlation ",
    ObjName, " ", Subset, ".pdf"
  ), width = 26, height = ceiling(length(cluster_genes) / 2) * 10
  )
  print(feature_plots)
  dev.off()

}

}

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
  labs(title = paste0("Assignment", analyses$which_assignment))

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
  "_", "Assignment (", analyses$which_assignment, ") ",
   analyses$viz_clustering, "clusters_ribo_ratio_violin_all.pdf"
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
  title = paste0(
    "Ribosome ratio by Assignment (", analyses$which_assignment,
  ")")
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
  "clusters ribosome ratios split by cluster assignment (", 
  analyses$which_assignment,
  ") sample.pdf"
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

# capture session info, versions, etc.
writeLines(
  capture.output(sessionInfo()), 
  paste0(ConfigDirectory, ObjName, "_", Subset, "_sessionInfo.txt")
)
file.copy(
  from = here("config.json"), 
  to = paste0(ConfigDirectory, ObjName, "_", Subset, "_config_params.json")
)
file.copy(
  from = here("analysis.json"), 
  to = paste0(ConfigDirectory, ObjName, "_", Subset, "_analysis_params.json")
)