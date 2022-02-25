get_top_cluster_markers <- function(
  markers, n
) {
  top_n_markers <- markers %>% 
    group_by(cluster) %>% 
    top_n(n = n, wt = avg_log2FC)
  
  top_n_markers$cluster <- as.character(top_n_markers$cluster)
  top_n_markers <- top_n_markers[order(top_n_markers$cluster), ]
  return (top_n_markers)
}




seurat_object = SeuratObj
downsample_n = 1000
top_n = 20
label_n = 2
cluster = "ClusterRNA"
data_type = "logcounts"
color_map = c("#007dd1", "white", "#ab3000")
use_raster = TRUE


top_n_markers <- get_top_cluster_markers(markers, top_n)
label_n_markers <- get_top_cluster_markers(markers, label_n)

subset <- subset(seurat_object, downsample = downsample_n)
sorted_barcodes <- names(sort(subset$ClusterRNA))

# plot data rows are genes, cols are cells
plot_data <- as.data.frame(subset@assays$RNA@scale.data)
plot_data <- plot_data[top_n_markers$gene, ]
plot_data <- na.omit(plot_data)
plot_data <- plot_data - rowMeans(plot_data)



# col_anno_df rows are cells, 
# cols are "ClusterRNA" / "ClusterWNN", "Sample"
col_anno_df <- subset@meta.data[, c(cluster, "Sample"), drop = F] 
# heatmap is ordered in cluster (primary) and then sample (secondary)
col_anno_df <- col_anno_df[
  order(col_anno_df[[cluster]], col_anno_df$Sample), 
  , 
  drop = F
]

# order cells by cluster (primary) and sample (secondary)
plot_data <- plot_data[rownames(col_anno_df)]

# sample and cluster colors are reverse of each other
sample_colors <- get_colors(
  seurat_object = SeuratObj,
  color_by = "Sample"
)
cluster_colors <- get_colors(
  seurat_object = SeuratObj,
  color_by = cluster,
  color_reverse = TRUE    
)
column_colors = list()
column_colors[["Sample"]] <- sample_colors
column_colors[[cluster]] <- cluster_colors



col_anno <- columnAnnotation(
  df = col_anno_df,
  show_annotation_name = TRUE,
  show_legend = TRUE,
  col = column_colors
)

row_anno <- rowAnnotation(
  # only label select genes (defined by input args)
  sel = anno_mark(
    at = match(label_n_markers$gene, row.names(plot_data)),
    labels = label_n_markers$gene,
    labels_gp = gpar(col = "black", fontsize = 9)
  )
)

data_colors <- circlize::colorRamp2(
  breaks = c(-2, 0, 2),
  colors = color_map
)

HM <- Heatmap(
  name = "logcounts",
  as.matrix(plot_data),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  top_annotation = col_anno,
  right_annotation = row_anno,
  row_names_gp = gpar(fontsize = 5),
  col = data_colors,
  show_column_names = FALSE,
  show_row_names = FALSE,
  border = FALSE,
  show_heatmap_legend = TRUE,
  use_raster = use_raster
)
HM

