################# LOAD UTILS  ##############
source("~/Box/Yun lab projects/victoria_liu/matching_patients/R_Code/utils.R")

# this object is fully pre-processed for GEX
SeuratObj <- readRDS(
  paste0(RobjDirectory, ObjName, Subset, "_GEX", ".rds")
)

############### CLUSTER / SAMPLE HEATMAP ##################
HM_object <- plot_heatmap (
  seurat_object = SeuratObj, downsample_n = 5000,
  markers = markers, top_n = 20, label_n = 2, 
  cluster = "ClusterRNA",
  data_type = "logcounts",
  use_raster = TRUE
)

pdf(paste0(
  heatDirectory, "heatmap", ObjName, Subset, 
  "res", RESOLUTION, "top20 genes per RNAcluster.pdf"
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
  seurat_object = SeuratObj, aesX = "Sample", fill = "ClusterRNA", 
  y_label = "Composition (Number of cells)", x_label = NULL, 
  y_lower_limit = 0, y_break = 1000,
  color_reverse = TRUE
)

pdf(paste0(
  densityplotDirectory, ObjName, Subset, "res", RESOLUTION, 
  "number of cells per sample and RNAclusters barplot.pdf"
),
width = 6, height = 5.5, family = FONT_FAMILY
)
P2
dev.off()



P3 <- plot_bargraph (
  seurat_object = SeuratObj, aesX = "ClusterRNA", fill = "Sample",
  y_label = "Composition (percentage of cells)", x_label = NULL, 
  y_lower_limit = 0, y_break = 0.2,
  position = "fill",
  plot_margin = unit(c(0.2, 0.5, 0.2, 0.5), "cm")
)

pdf(paste0(
  densityplotDirectory, ObjName, Subset, 
  "res", RESOLUTION, 
  "percent of cells per sample and RNAclusters barplot.pdf"
), width = 6, height = 5.5, family = FONT_FAMILY
)
P3
dev.off()



P4 <- plot_bargraph (
  seurat_object = SeuratObj, aesX = "ClusterRNA", fill = "Sample", 
  y_label = "Composition (Number of cells)", x_label = NULL, 
  y_lower_limit = 0, y_break = 1000
)

pdf(paste0(
  densityplotDirectory, ObjName, Subset, 
  "res", RESOLUTION, "number of cells per RNAcluster and sample barplot.pdf"
), width = 6, height = 5.5, family = FONT_FAMILY
)
P4
dev.off()


################# CLUSTER / SAMPLE DIMPLOTS ###################
U1 <- plot_umap(
  seurat_object = SeuratObj, group_by = "ClusterRNA",
  title = "Clusters", xlab = "UMAP1", ylab = "UMAP2",
  legend_position = "bottom",
  ncol_guide = 5, # this might be unnecessary?
  color_reverse = TRUE, label_clusters = TRUE
)

pdf(paste0(
  UMAPDirectory, ObjName, Subset, 
  "res", RESOLUTION, 
  " RNAClusters UMAP.pdf"
), width = 5.5, height = 6, family = FONT_FAMILY
)
U1
dev.off()



U2 <- plot_umap(
  seurat_object = SeuratObj, group_by = "Sample",
  title = "Samples", xlab = "UMAP1", ylab = "UMAP2",
  legend_position = "bottom",
  ncol_guide = 3
)


pdf(paste0(
  UMAPDirectory, ObjName, Subset, 
  "Samples UMAP.pdf"
), width = 5.5, height = 6, family = FONT_FAMILY
)
U2
dev.off()



U3 <- plot_umap (
  seurat_object = SeuratObj, group_by = "ClusterRNA",
  title = "Clusters", xlab = "UMAP1", ylab = "UMAP2",
  legend_position = "right",
  title_font_size = 16, x_font_size = 16, y_font_size = 16, 
  pt_size = 0.2, split_by = "Sample", ncol_dimplot = 3,
  color_reverse = TRUE
)

pdf(paste0(
  UMAPDirectory, ObjName, Subset, 
  "RNACluster UMAP Iteration by sample.pdf"
), width = 12, height = 8, family = FONT_FAMILY
)
U3
dev.off()


################# PREDEFINED POPULATION GENES ###################
# loads Glusterspecificgenes
load(paste0(
  "~/Box/Yun\ lab\ projects/scRNAseq\ Project/2808combined\ old\ and\ new/",
  "Victori2808aAnalysis/R_Objects/ClusterSpecificGenes.rds"
))

dotgraphs <- list()
for (cluster_name in names(Glusterspecificgenes)) {
  cluster_genes <- Glusterspecificgenes[[cluster_name]]
  
  D1 <- plot_dotgraph(
    seurat_object = SeuratObj, group_by = "Cluster",
    features = cluster_genes, title = cluster_name
  )
  dotgraphs[[cluster_name]] <- D1
}

plots <- ggarrange(plots = dotgraphs, nrow = length(Glusterspecificgenes))


pdf(paste0(
  OutputDirectory, "dotplot ", ObjName, Subset,
  "by predefined Cluster.pdf"
), width = 18, height = 18, family = FONT_FAMILY
)
print(plots)
dev.off()



D2 <- plot_dotgraph(
  seurat_object = SeuratObj, group_by = "Cluster",
  features = unique(get_top_cluster_markers(markers, 2)$gene),
  title = "Top 2 genes by cluster"
)


pdf(paste0(
  OutputDirectory, "dotplot ", ObjName, Subset, "top 2 genes by Cluster.pdf"
), width = 10, height = 3.5
)
D2
dev.off()

################# CLUSTER IDENTITY ASSIGNMENT ###############
## assign clusters (work with Dr. Yun)
Assign = list()
Assign[["Tcells"]] = c("C04", "C10")
Assign[["Microglia"]] = c("C01", "C08")
Assign[["Glioma"]] = c("C03", "C05", "C07", "C11", "C12")
Assign[["BCells"]] = c("C09")
Assign[["Myeloid"]] = c("C02", "C06")
Assign[["Pericytes and Vasc"]] = c("C12")

Assignment <- NA
for(i in 1 : length(Assign)){
  Assignment[SeuratObj@meta.data$Cluster %in% Assign[[i]]] <- names(Assign[i])
}
SeuratObj <- AddMetaData(
  object = SeuratObj, metadata = Assignment, col.name = "Assignment"
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
  "number of cells per sample and Assignment barplot.pdf"
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
  "percentage of cells per major population and Sample barplot.pdf"
), width = 7, height = 5.5, family = FONT_FAMILY
)
P7
dev.off()

plots = ggarrange(P1, P2, P3, P4, P5, P6, P7, ncol = 2)
pdf(paste0(
  densityplotDirectory, ObjName, Subset, 
  "res", RESOLUTION, "all barplots.pdf"
), width = 18, height = 18, family = FONT_FAMILY
)
print(plots)
dev.off()


################# ASSIGNMENT DIMPLOTS ###################
U4 <- plot_umap(
  seurat_object = SeuratObj, group_by = "Assignment",
  title = "Assignment", xlab = "UMAP1", ylab = "UMAP2",
  legend_position = "bottom",
  ncol_guide = 4,
  label_clusters = TRUE
)

pdf(paste0(
  UMAPDirectory, ObjName, Subset, 
  "res", RESOLUTION, "Assignment UMAP.pdf"
), width = 5.5, height = 6, family = FONT_FAMILY
)
U4
dev.off()

plots2 <- ggarrange(U1, U2, U4, ncol = 3)

pdf(paste0(
  UMAPDirectory, ObjName, Subset, 
  "res", RESOLUTION, "all UMAPs.pdf"
), width = 15, height = 6, family = FONT_FAMILY
)
print(plots2)
dev.off()



U5 <- plot_umap(
  seurat_object = SeuratObj, group_by = "Assignment",
  title = "Assignment", xlab = "UMAP1", ylab = "UMAP2",
  legend_position = "top",
  title_font_size = 16, x_font_size = 16, y_font_size = 16, 
  pt_size = 0.2, split_by = "Sample", ncol_dimplot = 2
)

pdf(paste0(
  UMAPDirectory, ObjName, Subset, 
  "Assignment UMAP Iteration by sample.pdf"
), width = 12, height = 12, family = FONT_FAMILY
)
U5
dev.off()


########### PREDEFINED CLUSTER FEATURE PLOTS #############
for (cluster_name in names(Glusterspecificgenes)) {
  predefined_cluster_plots <- list()
  cluster_genes <- Glusterspecificgenes[[cluster_name]]
  
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
saveRDS(SeuratObj, file = paste0(RobjDirectory, "AllClusters.rds"))

# example for how to read it into current R session
SeuratObj = readRDS(
  paste0(RobjDirectory, "AllClusters.rds")
)

# capture session info, versions, etc.
writeLines(
  capture.output(sessionInfo()), paste0(ConfigDirectory, "sessionInfo.txt")
)
file.copy(
  from = here("config.json"), 
  to = paste0(ConfigDirectory, "config_params.json")
)
file.copy(
  from = here("analysis.json"), 
  to = paste0(ConfigDirectory, "analysis_params.json")
)
sink(type = "output")

# send first email to authorize
send_bot_email(
  traceback_file, output_file,
  subject = "pipeline finished",
  message_body = "congratulations."
)
