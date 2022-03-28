split_by = "Assignment"
method = "pearson"

x = "rna_PTPRC"
y = "adt_CD45RO"

lab_x = "PTPRC (RNA)"
lab_y = "CD45RO (ADT)"


sp1 <- plot_correlation (
  SeuratObj, 
  x = "rna_PTPRC", y = "adt_CD45RO", split_by = "Assignment",
  lab_x = "PTPRC (RNA)", lab_y = "CD45RO (ADT)",
  method = "pearson"
)

# loop each assignment
for (cluster_name in names(analyses$ADT_markers)) {
  predefined_cluster_plots <- list()
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
    predefined_cluster_plots[[ADT_feature]] <- F1[[2]]
    predefined_cluster_plots[[RNA_feature]] <- F2[[2]]
  }
  
  
  feature_plots <- ggpubr::ggarrange(
    plotlist = predefined_cluster_plots,
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
