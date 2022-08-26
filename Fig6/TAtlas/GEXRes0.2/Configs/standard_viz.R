################# LOAD UTILS + SEURATOBJ  ##############
source("~/Documents/victoria_liu/matching_patients/R_Code/utils.R")

# capture session info, versions, etc.
write_experimental_configs(code_file = "standard_viz")

# this object is fully pre-processed for GEX and/or ADT
if (ObjName != "WNN") {
  tmpObjName <- ObjName
} else {
  # ADT and WNN are preprocessed in the same object, called ADT
  tmpObjName <- "ADT"
}

RDS_filename <- paste0(
  RobjDirectory, tmpObjName, Subset, 
  "_resAll.rds"
)
if (! file.exists(RDS_filename)) {
  print ("no Seurat object to analyze")
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
Clusterspecificgenes <- master[[
  paste0(analyses[["denovo_lineage"]], "_lineage_markers")
]]

############### ADD CLUSTER METADATA  #################
SeuratObj <- add_cluster_metadata(SeuratObj)

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


########## SAVE LOUPE PROJECTIONS ##########
cluster_name <- get_cluster_name()

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
  markersRNA <- cluster_markers(
    SeuratObj, assay = "RNA", default_ident = "ClusterRNA"
    )
  # save, because it takes a little time to calculate
  write.csv(markersRNA, markers_filename)

  
} else if (analyses$viz_clustering == "ADT") {
  markers_filename <- paste0(
    OutputDirectory, ObjName, Subset, 
    " ADT cluster markers (by ADT)", "res", RESOLUTION, ".csv"
  )
  markersADT <- cluster_markers(
    SeuratObj, assay = "ADT", default_ident = "ClusterADT"
  )
  # save, because it takes a little time to calculate
  write.csv(markersADT, markers_filename)
  

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
}

############### CLUSTER / SAMPLE HEATMAP ##################
# RNA only
if (analyses$viz_clustering != "ADT") {
  heatmap_wrapper(SeuratObj, markersRNA, "RNA")
} 



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
width = 6, height = 6.5, family = FONT_FAMILY
)
print(P2)
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
print(P3)
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
print(P4)
dev.off()


################# CLUSTER / SAMPLE DIMPLOTS ##################
U1 <- plot_umap(
  seurat_object = SeuratObj,
  group_by = paste0("Cluster", analyses$viz_clustering),
  reduction = paste0("umap", analyses$viz_clustering),
  title = "Clusters", xlab = "UMAP1", ylab = "UMAP2",
  legend_position = "bottom",
  label_size = 5,
  ncol_guide = 5,
  color_reverse = TRUE, label_clusters = TRUE,
  shuffle = TRUE
)

pdf(paste0(
  UMAPDirectory, ObjName, Subset,
  "_res", RESOLUTION,
  "_", analyses$viz_clustering, "Clusters_UMAP.pdf"
), width = 5.5, height = 6, family = FONT_FAMILY
)
print(U1)
dev.off()

U2 <- plot_umap(
  seurat_object = SeuratObj,
  group_by = "Sample",
  color_scheme = analyses$hue_palette,
  subtype_by = analyses$chromatose_subtype_by,
  reduction = paste0("umap", analyses$viz_clustering),
  title = "Samples", xlab = "UMAP1", ylab = "UMAP2",
  legend_position = "bottom",
  ncol_guide = 3, shuffle = TRUE
)

pdf(paste0(
  UMAPDirectory, ObjName, Subset,
  "_res", RESOLUTION,
  "_Samples_UMAP.pdf"
), width = 8, height = 7, family = FONT_FAMILY
)
print(U2)
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
), width = 12, height = 2.5 * length(unique(SeuratObj$Sample)),
  family = FONT_FAMILY
)
print(U3)
dev.off()


################# (CUSTOM) CELL POPULATION: CLUSTIFYR ###################
# refSeuratObj <- readRDS(
#   paste0(
#     RobjDir,
#     "GBMAtlas/Allhuman-11-3-21.rds"
#   )
# )

# gam_git_gib <- list(
#   Myeloid = "GAM",
#   TCells = "GIT",
#   BCells = "GIB"
# )
#
# CellInfo <- refSeuratObj@meta.data
# for (i in 1:length(gam_git_gib)) {
#   CellInfo$Assignment[
#     which(str_detect(refSeuratObj$Assignment, names(gam_git_gib[i])))
#   ] <- gam_git_gib[[i]]
# }
#
# refSeuratObj$Assignment <- CellInfo$Assignment
#
# saveRDS(
#   refSeuratObj,
#   file = paste0(
#     RobjDir,
#     "GBMAtlas/",
#     "Allhuman-11-3-21.rds"
#   )
# )


# seurat_ref_matrix <- seurat_ref(
#   seurat_object = refSeuratObj,
#   cluster_col = "Assignment"
# )

REF_MATRICES <- list(
  cbmc_ref
  # seurat_ref_matrix
)
REF_MATRICES_NAMES <- list(
  "cbmc"
  # "nour_all"
)
stopifnot(length(REF_MATRICES) == length(REF_MATRICES_NAMES))

for (idx in 1:length(REF_MATRICES)) {
  REF_MATRIX <- REF_MATRICES[[idx]]
  REF_MATRIX_NAME <- REF_MATRICES_NAMES[idx]
  clustifyr_colname <- paste0("clustifyr_", REF_MATRIX_NAME)

  SeuratObj <- clustifyr_wrapper(SeuratObj, REF_MATRIX, clustifyr_colname)

}


################# (CUSTOM) MANUAL CELL ASSIGNMENT ###############
# add desired automated assignment to Assignment column
if (
  ! is.na(analyses$which_assignment) &
  analyses$which_assignment %in% colnames(SeuratObj@meta.data)
) {
  SeuratObj$Assignment <- SeuratObj[[analyses$which_assignment]]
} else (
  print ("did not update Assignment; could not find column")
)
if (analyses$which_assignment == "Assignment") {
  analyses$which_assignment <- "Manual"
}

# for manual hand correction
if (length(analyses$cluster_to_assignment) > 0) {
  # if Assignment isn't a column, make sure to add it as a column
  if (! "Assignment" %in% SeuratObj@meta.data) {
    SeuratObj$Assignment <- NA
  }
  
  # manual correction
  for (new_assignment in names(analyses$cluster_to_assignment)) {
    cluster = paste0(analyses$cluster_prefix, new_assignment)
    SeuratObj$Assignment[
      which(str_detect(
        SeuratObj[[paste0("Cluster", analyses$viz_clustering)]][, 1], cluster
      ))
    ] <- analyses$cluster_to_assignment[[new_assignment]]
  }
  # label plots as manual assignment
  analyses$which_assignment <- "Manual"
} 



write.csv(
  SeuratObj[["Assignment"]],
  paste0(
    LoupeDirectory, ObjName, Subset,
    "_assignments (", analyses$which_assignment, ") .csv"
  )
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
print(P5)
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
print(P6)
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
), width = 8, height = 6, family = FONT_FAMILY
)
print(P7)
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
print(P8)
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


################# (CUSTOM) SINGLER ASSIGNMENT DIMPLOTS ###################
if (! exists("refs_singler")) {
  define_refs_singler()
}
# will have three singler plots
singler_assignment_plots <- list()
for (ref_name in names(refs_singler)) {
  singler_colname <- paste0("SingleR_", ref_name)
  stopifnot(singler_colname %in% colnames(SeuratObj@meta.data))

  # plot list
  singler_assignment_plots[[ref_name]] <-
    plot_umap(
      seurat_object = SeuratObj, group_by = singler_colname,
      reduction = paste0("umap", analyses$viz_clustering),
      title = paste0(ref_name, " SingleR Assignment"),
      xlab = "UMAP1", ylab = "UMAP2",
      legend_position = "bottom",
      ncol_guide = 4,
      label_clusters = TRUE,
      label_size = 5,
      color_reverse = FALSE,
      repel_labels = TRUE
    )
}

# egg package
singler_plots <- ggarrange(
  singler_assignment_plots[[1]],
  singler_assignment_plots[[2]],
  singler_assignment_plots[[3]],
  ncol = 3
  )

pdf(paste0(
  singleRDirectory, ObjName, Subset,
  "_resAll SingleR comparisons.pdf"
), width = 20, height = 7, family = FONT_FAMILY
)
print(singler_plots)
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
  label_size = 5,
  color_reverse = FALSE,
  repel_labels = TRUE
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
print(U4)
dev.off()

U5 <- plot_umap(
  seurat_object = SeuratObj, group_by = "Assignment",
  reduction = paste0("umap", analyses$viz_clustering),
  title = paste0(analyses$which_assignment, " Assignment"),
  xlab = "UMAP1", ylab = "UMAP2",
  legend_position = "bottom",
  title_font_size = 16, x_font_size = 16, y_font_size = 16, 
  pt_size = 0.2, split_by = "Sample", ncol_dimplot = 2,
  repel_labels = TRUE
)

pdf(paste0(
  UMAPDirectory, ObjName, Subset, 
  "_res", RESOLUTION, 
  "_", analyses$viz_clustering,
  "Clusters_Assignment (", 
  analyses$which_assignment,
  ") UMAP Iteration by sample.pdf"
), width = 12, height = 2.5 * length(unique(SeuratObj$Sample)), 
  family = FONT_FAMILY
)
print(U5)
dev.off()


U6 <- plot_umap(
  seurat_object = SeuratObj, 
  group_by = "Assignment",
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
), width = 12, height = 2.5 * length(unique(SeuratObj$Sample)), family = FONT_FAMILY
)
print(U6)
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



################# (GEX) DOTPLOTS ###################
# we want GEX dotplots for GEX or WNN umaps
if (analyses$viz_clustering != "ADT") {
dotgraphs_GEX <- list()
for (cluster_name in names(Clusterspecificgenes)) {
  cluster_genes <- Clusterspecificgenes[[cluster_name]]
  D1 <- plot_dotgraph(
    seurat_object = SeuratObj,
    group_by = paste0("Cluster", analyses$viz_clustering),
    features = cluster_genes, title = cluster_name
  )
  
  dotgraphs_GEX[[cluster_name]] <- D1

}

dot_plots_GEX <- ggpubr::ggarrange(
  plotlist = dotgraphs_GEX,
  ncol = 1,
  nrow = 1
)

pdf(paste0(
  dotDirectory, "GEX dotplot ", ObjName, Subset,
  " by predefined Cluster ",
  analyses[["denovo_lineage"]],
  ".pdf"
),
width = 12,
height = 8,
family = FONT_FAMILY
)
print(dot_plots_GEX)
dev.off()


D2 <- plot_dotgraph(
  seurat_object = SeuratObj,
  group_by = paste0("Cluster", analyses$viz_clustering),
  features = unique(get_top_cluster_markers(markersRNA, 5)$marker),
  title = "Top 5 genes by cluster",
  features_sorted = TRUE
)
pdf(paste0(
  dotDirectory, "dotplot ", ObjName, Subset, "top 5 genes by Cluster.pdf"
), width = 22, height = 8
)
print(D2)
dev.off()
}

################# (ADT) DOTPLOTS ###############
if (USE_ADT & analyses$viz_clustering != "RNA") {
dotgraphs_ADT <- list()
# loop each assignment
for (cluster_name in names(Clusterspecificgenes)) {
  print(cluster_name)
  # gene names, case sensitive to match RNA_to_ADT dict
  cluster_genes <- Clusterspecificgenes[[cluster_name]]
  cluster_genes <- case_sensitive_features(
    features = cluster_genes, reference = names(master$RNA_to_ADT)
  )
  # translate to features; some CSP might not exist
  cluster_features <- list()
  for (RNA_feature in cluster_genes) {
    for (ADT_feature in master$RNA_to_ADT[[RNA_feature]]) {
      cluster_features[[ADT_feature]] <- RNA_feature
    }
  }
  if (length(cluster_features) == 0) {
    print(paste0("skipping cluster ", cluster_name))
    next
  }
  
  D1 <- plot_dotgraph(
    seurat_object = SeuratObj,
    group_by = paste0("Cluster", analyses$viz_clustering),
    features = names(cluster_features), title = cluster_name,
    assay = "ADT"
  )
  dotgraphs_ADT[[cluster_name]] <- D1
  
  
}


dot_plots_ADT <- ggpubr::ggarrange(
  plotlist = dotgraphs_ADT,
  ncol = 1,
  nrow = 1
)

pdf(paste0(
  dotDirectory, "CSP dotplot ", ObjName, Subset,
  "by predefined Cluster ", analyses[["denovo_lineage"]],
  ".pdf"
), 
width = 10, 
height = 8, 
family = FONT_FAMILY
)
print(dot_plots_ADT)
dev.off()



D3 <- plot_dotgraph(
  seurat_object = SeuratObj,
  group_by = paste0("Cluster", analyses$viz_clustering),
  features = unique(get_top_cluster_markers(markersADT, 5)$marker),
  title = "Top 5 CSP by cluster",
  features_sorted = TRUE,
  assay = "ADT"
)
pdf(paste0(
  dotDirectory, "dotplot ", ObjName, Subset, 
  "top 5 CSP by ", analyses$viz_clustering, " Cluster.pdf"
), width = 22, height = 8
)
print(D3)
dev.off()
  
}
########### (GEX) PREDEFINED CLUSTER FEATURE PLOTS #############
if (analyses$viz_clustering == "RNA") {
# loop each assignment
for (cluster_name in names(Clusterspecificgenes)) {
  predefined_cluster_plots_GEX <- list()
  cluster_genes <- Clusterspecificgenes[[cluster_name]]
  # loop each marker gene in assignment
  for (gene in cluster_genes){
    gene <- case_sensitive_features(
      SeuratObj,
      c(gene),
      assay = "RNA"
    )
    
    if (length(gene) == 1) {
      print(gene)
      gene <- gene[[1]]
      F1 <- plot_featureplot (
        seurat_object = SeuratObj,
        feature_gene = gene,
        split_by = NULL,
        reduction = paste0("umap", analyses$viz_clustering),
        pt_size = 0.6
      )
      predefined_cluster_plots_GEX[[gene]] <- F1[[2]]
    }
  }
  
  plot_count <- length(predefined_cluster_plots_GEX)

  if (plot_count > 0) {
    feature_plots <- ggpubr::ggarrange(
      plotlist = predefined_cluster_plots_GEX,
      ncol = 2,
      nrow = plot_count %/% 2 + plot_count %% 2
    )
  
    pdf(paste0(
      featuremapDirectory, cluster_name,
      " featuremap umap_", analyses$viz_clustering, " ", 
      ObjName, " ", Subset, ".pdf"
    ), width = 26, height = ceiling(plot_count / 2) * 10
    )
    print(feature_plots)
    dev.off()
  }

}
}

#################### (ADT) FEATURE PLOTS #######################
# only plot ADT if we're using WNN UMAP
if (USE_ADT & analyses$viz_clustering != "RNA") {
# loop each assignment
for (cluster_name in names(Clusterspecificgenes)) {
  print(cluster_name)
  # gene names, case sensitive to match RNA_to_ADT dict
  cluster_genes <- Clusterspecificgenes[[cluster_name]]
  cluster_genes <- case_sensitive_features(
    features = cluster_genes, reference = names(master$RNA_to_ADT)
  )
  # translate to features; some CSP might not exist
  cluster_features <- list()
  for (RNA_feature in cluster_genes) {
    for (ADT_feature in master$RNA_to_ADT[[RNA_feature]]) {
        cluster_features[[ADT_feature]] <- RNA_feature
        }
    }
  if (length(cluster_features) == 0) {
    print(paste0("skipping cluster ", cluster_name))
    next
  }
  
  
  predefined_cluster_plots_ADT_GEX <- list()  
  # loop each marker gene in assignment
  for (idx in 1:length(names(cluster_features))){
    ADT_feature <- names(cluster_features)[idx]
    F1 <- plot_featureplot (
      seurat_object = SeuratObj,
      feature_gene = ADT_feature,
      assay = "ADT",
      split_by = NULL,
      reduction = paste0("umap", analyses$viz_clustering),
      pt_size = 0.6
    )

    RNA_feature <- cluster_features[[ADT_feature]]
    F2 <- plot_featureplot (
      seurat_object = SeuratObj,
      feature_gene = RNA_feature,
      assay = "RNA",
      split_by = NULL,
      reduction = paste0("umap", analyses$viz_clustering),
      pt_size = 0.6
    )
    predefined_cluster_plots_ADT_GEX[[2 * idx - 1]] <- F1[[2]]
    predefined_cluster_plots_ADT_GEX[[2 * idx]] <- F2[[2]]
  }


  feature_plots <- ggpubr::ggarrange(
    plotlist = predefined_cluster_plots_ADT_GEX,
    ncol = 2,
    nrow = 1
  )

  pdf(paste0(
    featuremapDirectory, cluster_name,
    " featuremap umap_", analyses$viz_clustering, " ",
    ObjName, " ", Subset, ".pdf"
  ), width = 26, height = 10
  )
  print(feature_plots)
  dev.off()

}

}
#################### (ADT) CORRELATION PLOTS ###################
if (USE_ADT & analyses$viz_clustering != "RNA") {
# loop each assignment
for (cluster_name in names(Clusterspecificgenes)) {
  print(cluster_name)
  # gene names, case sensitive to match RNA_to_ADT dict
  cluster_genes <- Clusterspecificgenes[[cluster_name]]
  cluster_genes <- case_sensitive_features(
    features = cluster_genes, reference = names(master$RNA_to_ADT)
  )
  # translate to features; some CSP might not exist
  cluster_features <- list()
  for (RNA_feature in cluster_genes) {
    for (ADT_feature in master$RNA_to_ADT[[RNA_feature]]) {
      cluster_features[[ADT_feature]] <- RNA_feature
    }
  }
  if (length(cluster_features) == 0) {
    print(paste0("skipping cluster ", cluster_name))
    next
  }
  
  predefined_cluster_plots_correlation <- list()

  # loop each marker gene in assignment
  for (ADT_feature in names(cluster_features)){
    RNA_feature <- cluster_features[[ADT_feature]]
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
    ncol = 1,
    nrow = 1
  )

  pdf(paste0(
    ADTDirectory, cluster_name,
    " ADT RNA correlation ",
    ObjName, " ", Subset, ".pdf"
  ), width = 26, height = 15
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
print(U0)
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
print(V_combined)
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
print(x_combined)
dev.off()



################# CELL CYCLE SORTING ###################

U7 <- plot_umap (
  seurat_object = SeuratObj, group_by = "Phase",
  reduction = paste0("umap", analyses$viz_clustering),
  title = "Seurat Cell Cycle scoring", xlab = "UMAP1", ylab = "UMAP2",
  legend_position = "bottom",
  ncol_guide = 3,
  title_font_size = 16, x_font_size = 16, y_font_size = 16, 
  pt_size = 0.2,
  label_clusters = TRUE, repel_labels = TRUE, label_size = 3,
  shuffle = TRUE
)


pdf(paste0(
  UMAPDirectory, ObjName, Subset, 
  "CellCycle UMAP.pdf"
), width = 8, height = 7, family = FONT_FAMILY
)
print(U7)
dev.off()

U8 <- plot_umap (
  seurat_object = SeuratObj, group_by = "Phase",
  reduction = paste0("umap", analyses$viz_clustering),
  title = "Seurat Cell Cycle scoring", xlab = "UMAP1", ylab = "UMAP2",
  legend_position = "right", split_by = "Sample",
  title_font_size = 16, x_font_size = 16, y_font_size = 16, 
  pt_size = 0.2, ncol_dimplot = 2,
  label_clusters = TRUE, repel_labels = TRUE, label_size = 3
)


pdf(paste0(
  UMAPDirectory, ObjName, Subset, 
  "CellCycle UMAP by Iteration Sample.pdf"
), width = 12, height = 2.5 * length(unique(SeuratObj$Sample)), family = FONT_FAMILY
)
print(U8)
dev.off()

# ################### (T CELLS ONLY) PSEUDOTIME SLINGSHOT ##################
# # figure out directory
# SlingDir <- paste0(PseudoDirectory, "/slingshot/")
# if (! dir.exists(SlingDir)) {
#   dir.create(SlingDir)
# }
# 
# # figure out column for colors (numbered by int, rather than char)
# SeuratObj <- meta_category_to_int(SeuratObj, "Assignment", "Assignment_ColNo")
# SeuratObj <- meta_category_to_int(SeuratObj, "ClusterRNA", "Cluster_ColNo")
# 
# 
# # color palettes
# ClusterColors <- Nour_pal("all", reverse = T)(
#   length(unique(SeuratObj$ClusterRNA))
# )
# AssignmentColors <- Nour_pal("all")(
#   length(unique(SeuratObj$Assignment))
# )
# 
# DefaultAssay(SeuratObj) <- "RNA"
# SeuratObj.sce = as.SingleCellExperiment(SeuratObj)
# 
# start_clusters <- c("CD4_NaiveLike", "CD8_NaiveLike")
# for (start_cluster in start_clusters) {
# 
#   # remove slingshot on second round
# 
#   if ("slingPseudotime_1" %in% colnames(colData(SeuratObj.sce))) {
#     SeuratObj.sce$slingPseudotime_1 <- NULL
#   }
# 
#   SeuratObj.sce <- slingshot(
#     SeuratObj.sce,
#     clusterLabels = "Assignment", reducedDim = "UMAPRNA",
#     reassign = F, maxit = 10, allow.breaks = T,
#     start.clus = start_cluster
#   )
# 
#   # add metadata back to SeuratObj
#   SeuratObj <- AddMetaData(
#     object = SeuratObj,
#     metadata = SeuratObj.sce$slingPseudotime_1,
#     col.name = 'slingPseudotime_1'
#   )
# 
#   SlingDirectory <- paste0(SlingDir, "/", start_cluster, "/")
#   if (! dir.exists(SlingDirectory)) {
#     dir.create(SlingDirectory)
#   }
# 
# 
#   pdf(paste0(
#     SlingDirectory, ObjName, Subset,
#     "res", RESOLUTION, "slingshot umap.pdf"
#   ), width = 5, height = 5, family = "ArialMT", onefile = T
#   )
# 
#   # colored by cluster
#   plot_slingshot_umap(
#     SeuratObj.sce,
#     ClusterColors, "Cluster_ColNo",
#     "Slingshot pseudotime, colored by cluster"
#   )
# 
#   # colored by assignment
#   plot_slingshot_umap(
#     SeuratObj.sce,
#     AssignmentColors, "Assignment_ColNo",
#     "Slingshot pseudotime, colored by assignment"
#   )
#   dev.off()
# 
# 
# 
# 
#   # by assignment 1
#   G1 <- plot_slingshot_pseudotime(
#     SeuratObj, "slingPseudotime_1", "Assignment", AssignmentColors
#   )
# 
#   # by cluster 1
#   G2 <- plot_slingshot_pseudotime(
#     SeuratObj, "slingPseudotime_1", "ClusterRNA", ClusterColors
#   )
# 
# 
#   pdf(paste0(
#     SlingDirectory, ObjName, Subset,
#     "res", RESOLUTION, "slingshot time.pdf"
#   ), width = 5, height = 5, family = "ArialMT", onefile = T
#   )
#   print(G1)
#   print(G2)
#   dev.off()
# 
# }
# ################### (T CELLS ONLY) PSEUDOTIME MONOCLE3 ##################
# # figure out directories
# MonocleDir <- paste0(PseudoDirectory, "/monocle3/")
# if (! dir.exists(MonocleDir)) {
#   dir.create(MonocleDir)
# }
# 
# 
# # rename umap so that monocle3 recognizes it
# SeuratObj@reductions$umap <- SeuratObj@reductions$umapRNA
# # run monocle
# SeuratObj.cds <- as.cell_data_set(SeuratObj)
# SeuratObj.cds <- cluster_cells(
#   cds = SeuratObj.cds, reduction_method = "UMAP"
# )
# # learn mst
# SeuratObj.cds <- learn_graph(
#   SeuratObj.cds, verbose = T
# )
# # find order of graph
# get_earliest_principal_node <- function(cds, early_bin = start_cluster){
#   cell_ids <- which(colData(cds)[, "Assignment"] == early_bin)
#   closest_vertex <-
#     cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
#   closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
#   root_pr_nodes <-
#     igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(
#       which.max(table(closest_vertex[cell_ids,]))
#     ))]
# 
#   root_pr_nodes
# }
# 
# for (start_cluster in start_clusters) {
#   MonocleDirectory <- paste0(MonocleDir, "/", start_cluster, "/")
#   if (! dir.exists(MonocleDirectory)) {
#     dir.create(MonocleDirectory)
#   }
# 
#   MonocleFeatureDirectory <- paste0(MonocleDirectory, "/FeatureMaps/")
#   if (! dir.exists(MonocleFeatureDirectory)) {
#     dir.create(MonocleFeatureDirectory)
#   }
# 
# 
#   SeuratObj.cds <- order_cells(
#     SeuratObj.cds,
#     root_pr_nodes = get_earliest_principal_node(
#       SeuratObj.cds, early_bin = start_cluster
#       )
#   )
# 
# 
#   xx = fData(SeuratObj.cds)
#   xx$gene_short_name = row.names(xx)
#   fData(SeuratObj.cds) = xx
# 
#   # feature maps
#   for (cluster_name in names(Clusterspecificgenes)) {
#     cluster_genes <- Clusterspecificgenes[[cluster_name]]
#     cluster_genes <- case_sensitive_features(
#       SeuratObj,
#       c(cluster_genes),
#       assay = "RNA"
#     )
#     PP1 <- plot_cells(
#       cds = SeuratObj.cds, scale_to_range = T, cell_size = 0.7,
#       label_cell_groups = FALSE,
#       show_trajectory_graph = TRUE, trajectory_graph_color = "black",
#       trajectory_graph_segment_size = 0.7,
#       label_leaves = TRUE, min_expr = 1,
#       label_branch_points = TRUE, genes = cluster_genes
#     ) + theme_min()
# 
# 
#     # plot cells
#     pdf(paste0(
#       MonocleFeatureDirectory, ObjName, Subset,
#       " monocle 3 ", cluster_name, " markers.pdf"
#     ),
#     width = length(cluster_genes) * 1.5,
#     height = length(cluster_genes),
#     family = "ArialMT"
#     )
#     print(PP1)
#     dev.off()
#   }
# 
# 
#   # by cluster
#   PP2 <- plot_cells(
#     cds = SeuratObj.cds,
#     label_cell_groups = FALSE,
#     color_cells_by = "cluster",
#     show_trajectory_graph = TRUE,
#     label_leaves = TRUE,
#     cell_size = 1,
#     label_branch_points = TRUE,
#     trajectory_graph_color = "white",
#     trajectory_graph_segment_size = 1
#   ) + scale_color_Nour(palette="all",reverse = F,discrete = T) +
#     theme_gray() +
#     scale_y_continuous(breaks=NULL) +
#     scale_x_continuous(breaks=NULL) +
#     xlab("UMAP1") +
#     ylab("UMAP2") +
#     labs(title = "Monocle3 Pseudotime") +
#     FontSize(x.title = 16, y.title = 16, main = 16) +
#     theme(legend.position = "bottom",text =element_text(size=15)) +
#     guides(
#       colour = guide_legend(
#         title = "Cluster",
#         override.aes = list(size = 5),
#         title.theme = element_text(size = 15, face = "bold"),
#         title.position = "top",
#         label.theme = element_text(size=15)
#       )
#     )
# 
#   pdf(paste0(
#     MonocleDirectory, ObjName, Subset,
#     " monocle 3 colored by cluster.pdf"
#   ),width = 5, height = 6, family = "ArialMT"
#   )
#   print(PP2)
#   dev.off()
# 
#   # by time
#   PP3 <- plot_cells(
#     SeuratObj.cds,
#     color_cells_by = "pseudotime",
#     label_cell_groups = FALSE,
#     label_leaves = TRUE,
#     label_branch_points = TRUE
#   )
#   pdf(paste0(
#     MonocleDirectory, ObjName, Subset,
#     " monocle 3 colored by pseudotime.pdf"
#   ),width = 7, height = 6, family = "ArialMT"
#   )
#   print(PP3)
#   dev.off()
# }
# 

############## SAVE SEURAT ################
saveRDS(
  SeuratObj,
  file = paste0(
    RobjDirectory, ObjName, Subset,
    "_res", RESOLUTION, ".rds"
    )
)