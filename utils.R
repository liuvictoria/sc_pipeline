################# INSTALL / LOAD PACKAGES  ##############
# # install them first if necessary. Only need to install once
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# BiocManager::install(version = "3.14")
# 
# install.packages("devtools")

library(devtools)
library(easypackages)

# only need to install once
# devtools::install_github("satijalab/sctransform", ref = "develop")
# devtools::install_github("jokergoo/ComplexHeatmap")
# BiocManager::install("ComplexHeatmap")
# BiocManager::install("SingleCellExperiment")
# BiocManager::install('limma')
# BiocManager::install("clustifyr")
# # gamma poisson generalized linear model
# BiocManager::install("glmGamPoi")

# remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
# remotes::install_github("mojaveazure/seurat-disk")
# devtools::install_github("sjessa/ggmin")
# install.packages("magick")
# install.packages("gmailr", repos="http://cran.r-project.org")
# install.packages("rjson")

MyPackages <- c(
  "dplyr", "ggplot2", "ggpubr", "gridExtra", "viridis", "egg",
  "grid", "lattice", "gtools", "Biobase", "RColorBrewer",
  "Seurat", "cowplot", "patchwork", "stringr", "ComplexHeatmap", 
  "SingleCellExperiment", "ggmin", "Nourpal", "Cairo",
  "DoubletFinder", "harmony", "magick", "viridis", "limma", 
  "glmGamPoi", "gmailr", "rjson", "here", "SeuratDisk",
  "gplots", "clustifyr", "fgsea"
)

# similar to libaries, but will install package as well
packages(MyPackages)

# get experiment parameters
config <- fromJSON(file = here("config.json"))

# Nour's palette. Install and load
NourpalDirectory <- config$NourpalDirectory
setwd(NourpalDirectory)
devtools::load_all("Nourpal.Rproj")
library(Nourpal)



################ SET PARAMETERS ###############
DOUBLET_FORMATION_RATE <- config$DOUBLET_FORMATION_RATE
ADT_PRESENT <- config$ADT_PRESENT
USE_ADT <- config$USE_ADT
if (USE_ADT) {
  for (presence in ADT_PRESENT) {
    if (! presence) {
      stop("Trying to use ADT, but ADT is not present for all samples")
    }
  }
}
# resolution for cluster-finding after PCA
RESOLUTION <- config$RESOLUTION
# if we are doing scTransform
SCT <- config$SCT
# font family for output plots
FONT_FAMILY <- config$FONT_FAMILY
# files to analyze
FILES = config$FILES
# saving files
ObjName <- config$ObjName
Subset <- config$Subset
################# CONFIG DIRECTORIES #################
# Set working directory and create input / output folders
Directory <- config$Directory
if (! dir.exists(Directory)) dir.create(Directory)
setwd(Directory)

RobjDir <- paste0(Directory, "/R_Objects/")
if (! dir.exists(RobjDir)) dir.create(RobjDir)
RobjDirectory <- paste0(RobjDir, Subset, "/")
if (! dir.exists(RobjDirectory)) dir.create(RobjDirectory)
RcodeDirectory <- paste0(Directory, "/R_Code/")
if (! dir.exists(RcodeDirectory)) dir.create(RcodeDirectory)
dataDirectory <- paste0(Directory, "/Data/")
if (! dir.exists(dataDirectory)) dir.create(dataDirectory)
outdir <- paste0(Directory, "/Output/")
if (! dir.exists(outdir)) dir.create(outdir)


OutDirectory <- paste0(Directory, "/Output/", Subset, "/")
if (! dir.exists(OutDirectory)) dir.create(OutDirectory)
OutputDirectory <- paste0(
  OutDirectory, "/", ObjName, "Res", config$RESOLUTION, "/"
  )
if (! dir.exists(OutputDirectory)) dir.create(OutputDirectory)
ConfigDirectory <- paste0(OutputDirectory, "/Configs/")
if (! dir.exists(ConfigDirectory)) dir.create(ConfigDirectory)


QCDirectory <- paste0(OutputDirectory, "/QualityControl/")
if (! dir.exists(QCDirectory)) dir.create(QCDirectory)
ElbowDirectory <- paste0(OutputDirectory, "/ElbowPlots/")
if (! dir.exists(ElbowDirectory)) dir.create(ElbowDirectory)
densityplotDirectory <- paste0(OutputDirectory, "/DensityPlots/")
if (! dir.exists(densityplotDirectory)) dir.create(densityplotDirectory)

heatDirectory <- paste0(OutputDirectory, "/HeatMaps/")
if (! dir.exists(heatDirectory)) dir.create(heatDirectory)
featuremapDirectory <- paste0(OutputDirectory, "/FeatureMaps/")
if (! dir.exists(featuremapDirectory)) dir.create(featuremapDirectory)
UMAPDirectory <- paste0(OutputDirectory, "/UMAPs/")
if (! dir.exists(UMAPDirectory)) dir.create(UMAPDirectory)

dotDirectory <- paste0(OutputDirectory, "/DotPlots/")
if (! dir.exists(dotDirectory)) dir.create(dotDirectory)
LoupeDirectory <- paste0(OutputDirectory, "/LoupeProjections/")
if (! dir.exists(LoupeDirectory)) dir.create(LoupeDirectory)



############################# HELPER FUNCTIONS ###############################

{
###### PLOT HEIGHT ######
roundUpNice <- function(x, nice = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) {
  if (length(x) != 1) stop("'x' must be of length 1")
  return(
    10 ^ floor(log10(x)) * 
      nice[[which(x <= 10 ^ floor(log10(x)) * nice)[[1]]]]
  )
}

get_height <- function(category) {
  max_value = max(as.matrix(table(SeuratObj[[category]])))
  return(roundUpNice(max_value))
}

###### GET NOURPAL COLORS #######
get_colors <- function(
  seurat_object,
  color_by,
  color_reverse = FALSE
) {
  colors = Nour_pal("all", reverse = color_reverse)(
    length(unique(as.vector(as.matrix(SeuratObj[[color_by]]))))
  )
  names(colors) = sort(unique(as.vector(as.matrix(SeuratObj[[color_by]]))))
  return (colors)
}

############### GET CASE CORRECT FEATURES ##############
case_sensitive_features <- function(
  seurat_object, features
) {
  final_features <- c()
  for (feature in features) {
    idx_match <- match(
      tolower(feature), 
      tolower(rownames(SeuratObj)),
      nomatch = -1
    )
    if (idx_match != -1) {
      final_features <- c(
        final_features, rownames(SeuratObj)[idx_match])
    }
  }
  return (unlist(final_features))
}

###### HEATMAP ANNOTATIONS ########
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


################ TEMPLATE GGPLOT BARPLOT ###################

plot_bargraph <- function (
  seurat_object, aesX, fill,
  y_label, x_label, y_lower_limit, y_break,
  color_reverse = FALSE,
  position = "stack", geom_bar_color = "black", geom_bar_width = 0.7,
  x_text_angle = 45, plot_margin = unit(c(0.5, 0.5, 0.2, 0.5), "cm")
) {
  # get colors
  color_by <- fill
  manual_colors <- get_colors(seurat_object, color_by, color_reverse)
  # # determine upper limit
  if (position == "fill") {
    y_upper_limit = 1
  }
  else {
    y_upper_limit = get_height(aesX)
  }
  
  # get cell info metadata
  data_frame <- seurat_object@meta.data
  
  # plot object
  plot_object <- data_frame %>%   
    ggplot(aes_string(x = aesX, fill = fill)) +
    scale_fill_manual(values = manual_colors) +
    geom_bar(
      position = position, color = geom_bar_color, width = geom_bar_width
    ) +
    theme(
      # CONSTANT EVERYWHERE START
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      legend.text = element_text(color = "black", size = 13, face = "bold"),
      legend.title = element_text(color = "black", size = 13, face = "bold"),
      axis.line = element_line(colour = "black"),
      axis.text.y = element_text(color = "black", size = 12),
      # CONSTANT EVERYWHERE END
      axis.text.x = element_text(
        color = "black", angle = x_text_angle, hjust = 1, size = 15
      ),
      axis.text = element_text(size = 15, face = "bold"),
      axis.title = element_text(size = 15, face = "bold"),
      plot.margin = plot_margin
    ) +
    labs(y = y_label, x = x_label) +
    scale_y_continuous(
      expand = c(0, 0), limits = c(y_lower_limit, y_upper_limit),
      breaks = seq(y_lower_limit, y_upper_limit, by = y_break)
    )
  return (plot_object)
}

################ TEMPLATE DENSITY PLOT #################
plot_densitygraph <- function (
  seurat_object, aesX, fill, color_by, 
  xintercept = Inf,
  scale_x_log10 = FALSE,
  alpha = 0.2, color_reverse = FALSE
) {
  # get colors
  manual_colors <- get_colors(seurat_object, color_by, color_reverse)
  # cell info
  data_frame <- seurat_object@meta.data
  
  
  plot_object <- data_frame %>%
    ggplot(aes_string(x = aesX, fill = fill, color = color_by)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    scale_fill_manual(values = manual_colors) +
    scale_color_manual(values = manual_colors)
  
  if (xintercept != Inf) {
    plot_object <- plot_object +
      geom_vline(xintercept = xintercept)
  }
  
  if (scale_x_log10) {
    plot_object <- plot_object +
      scale_x_log10()
  }
  
  return (plot_object)
}


################ UMAP PLOT TEMPLATE ##################
plot_umap <- function(
  seurat_object, group_by,
  title, xlab, ylab,
  legend_position, reduction,
  color_reverse = FALSE,
  title_font_size = 20, x_font_size = 20, y_font_size = 20, 
  pt_size = NULL, split_by = NULL, ncol_dimplot = 1, ncol_guide = 1,
  label_clusters = FALSE, repel_labels = FALSE, label_size = 4
) {
  color_by <- group_by
  # get colors
  manual_colors <- get_colors(seurat_object, color_by, color_reverse)
  
  plot_object <- DimPlot(
    seurat_object, reduction = reduction, group.by = group_by,
    cols = manual_colors,
    pt.size = pt_size, split.by = split_by, ncol = ncol_dimplot,
    label = label_clusters, repel = repel_labels, label.size = label_size
  ) +
    ggmin::theme_min() +
    scale_y_continuous(breaks = NULL) +
    scale_x_continuous(breaks = NULL) +
    xlab(xlab) +
    ylab(ylab) +
    FontSize(
      main = title_font_size,
      x.title = x_font_size,
      y.title = y_font_size
    ) +
    theme(legend.position = legend_position) +
    labs(title = title) +
    guides(colour = guide_legend(
      override.aes = list(size = 5),
      title.theme = element_text(size = 15, face = "bold"),
      title.position = "top",
      label.theme = element_text(size = 15),
      ncol = ncol_guide
    ))
  return (plot_object)
}


#################### TEMPLATE DOTPLOT #####################
plot_dotgraph <- function (
  seurat_object, group_by,
  features, title,
  x_text_angle = 45,
  dot_scale = 5, scale = TRUE,
  color_palette_option = "plasma"
) {
  features <- case_sensitive_features(seurat_object, features)
  plot_object <- DotPlot(
    seurat_object, group.by = group_by, dot.scale = dot_scale,
    features = features, scale = scale
  ) + 
    scale_colour_viridis_c(option = color_palette_option) + 
    labs(title = title) +
    ggmin::theme_min() +
    RotatedAxis() +
    theme(
      # CONSTANT EVERYWHERE START
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      title = element_text(size = 15, face = "bold"),
      legend.text = element_text(color = "black", size = 13, face = "bold"),
      legend.title = element_text(color = "black", size = 13, face = "bold"),
      axis.line = element_line(colour = "black"),
      axis.text.y = element_text(color = "black", size = 12),
      # CONSTANT EVERYWHERE END
      axis.text.x = element_text(
        color = "black", angle = x_text_angle, hjust = 1, size = 15
      ),
      axis.text = element_text(size = 15, face = "bold"),
      axis.title = element_text(size = 15, face = "bold"),
      plot.margin = unit(c(0.5, 0.5, 0.2, 0.5), "cm")
    )
  
  return (plot_object)
}


#################### TEMPLATE FEATURE PLOT ###################
plot_featureplot <- function(
  seurat_object, feature_gene, split_by, reduction,
  order = TRUE, label = TRUE, label_size = 2,
  cols = Nour_cols(c("darkpurple", "lightorange")), pt_size = 0.1,
  legend_title_position = "top", legend_title_size = 10,
  legend_location = "right",
  ncol = 3, nrow = 1, widths = c(0.03, 1, 1)
) {
  feature_gene <- case_sensitive_features(
    seurat_object,
    c(feature_gene)
  )
  if (length(feature_gene) == 0) {
    stop(paste0("could not find gene feature: ", feature_gene))
  } else {
    feature_gene <- feature_gene[[1]]
  }
  
  F1 <- FeaturePlot(
    seurat_object, features = feature_gene, split.by = split_by,
    order = order, label = label, label.size = label_size,
    cols = cols, pt.size = pt_size, reduction = reduction
  )
  legend <- get_legend(
    F1[[2]] +
      theme_min() +
      NoAxes() +
      guides(colour = guide_colourbar(
        title = feature_gene, title.position = legend_title_position,
        title.theme = element_text(size = legend_title_size)))
  )
  
  plot_object <- ggpubr::ggarrange(
    ggparagraph(text = " ",  size = 0),
    F1[[1]] + theme_min2() + NoLegend() + NoAxes(),
    F1[[2]] + theme_min2() + NoLegend() + NoAxes(),
    ncol = ncol, nrow = nrow, widths = widths,
    legend.grob = legend, legend = legend_location
  )
  
  return (plot_object)
}


################### TEMPLATE HEATMAP #################
plot_heatmap <- function(
  seurat_object, downsample_n,
  markers, top_n, label_n, 
  cluster,
  data_type = "logcounts",
  color_map = c("#007dd1", "white", "#ab3000"),
  use_raster = TRUE
) {
  top_n_markers <- get_top_cluster_markers(markers, top_n)
  label_n_markers <- get_top_cluster_markers(markers, label_n)
  
  subset <- subset(seurat_object, downsample = downsample_n)
  sorted_barcodes <- names(sort(subset$ClusterRNA))
  
  # plot data rows are genes, cols are cells
  plot_data <- as.data.frame(subset@assays$RNA@scale.data)
  plot_data <- plot_data[top_n_markers$gene, ]
  plot_data <- na.omit(plot_data)
  plot_data <- plot_data - rowMeans(plot_data)
  
  # col_anno_df rows are cells
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
    show_legend = FALSE,
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
  
  lgd1 <- Legend(
    labels = unique(as.vector(as.matrix(col_anno_df[[cluster]]))),
    title = "Cluster",
    legend_gp = gpar(fill = cluster_colors, fontsize = 25)
  )
  lgd2 <- Legend(
    labels = levels(as.factor(col_anno_df$Sample)),
    title = "Sample", 
    legend_gp = gpar(fill = sample_colors, fontsize = 25)
  )
  
  return (list(HM, lgd1, lgd2))
}


######### TEMPLATE ARRANGE PLOTS W COMMON LEGEND ####
grid_arrange_shared_legend <- function(
  ...,
  ncol = length(list(...)),
  nrow = 1,
  position = c("bottom", "right")
) {
  plots <- list(...)
  position <- match.arg(position)
  g <- 
    ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x)
    x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x)
    x + theme(legend.position = "none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(
    position,
    "bottom" = arrangeGrob(
      do.call(arrangeGrob, gl),
      legend,
      ncol = 1,
      heights = unit.c(unit(1, "npc") - lheight, lheight)
    ),
    "right" = arrangeGrob(
      do.call(arrangeGrob, gl),
      legend,
      ncol = 2,
      widths = unit.c(unit(1, "npc") - lwidth, lwidth)
    )
  )
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
}



######### TEMPLATE THEME2  #########
theme_min2 <- function(base_size = 11, base_family = "") {
  theme_light(base_size = 11, base_family = "") +
    theme(     
      plot.title = element_text(size = rel(0.9), hjust = 0.5, vjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, colour = "grey90", size = 1),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(colour = "black", size = rel(1.2)),
      strip.text.y = element_text(colour = "black", size = rel(1.2)),
      title = element_text(size = rel(0.9), hjust = 0.5, vjust = 0.5),
      axis.text = element_text(colour = "black", size = rel(0.8)),
      axis.title = element_blank(),
      legend.title = element_text(colour = "black", size = rel(0.9), hjust = 0.5),
      legend.key.size = unit(0.9, "lines"),
      legend.text = element_text(size = rel(0.7), colour = "black"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA),
      plot.margin = unit(c(0.1, 0, 0, -0.2), "lines")
    )
}
}

