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
# BiocManager::install("celldex")
# BiocManager::install("SingleR")
# BiocManager::install("biomaRt")
# BiocManager::install("scDblFinder")
# BiocManager::install("miQC")

# install.packages("tidyverse")
# install.packages('clustree') 
# install.packages("pheatmap")
# install.packages("corrplot")


# remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
# # remotes::install_github("mojaveazure/seurat-disk")
# remotes::install_github("carmonalab/UCell", ref="v1.3")
# remotes::install_github("carmonalab/scGate")
# remotes::install_github("carmonalab/ProjecTILs")
# remotes::install_github('satijalab/seurat-wrappers')
# devtools::install_github("sjessa/ggmin")
# install.packages("magick")
# install.packages("gmailr", repos="http://cran.r-project.org")
# install.packages("rjson")
# install.packages("PerformanceAnalytics")

MyPackages <- c(
  "dplyr", "ggplot2", "ggpubr", "gridExtra", "viridis", "egg",
  "grid", "lattice", "gtools", "Biobase", "RColorBrewer",
  "Seurat", "cowplot", "patchwork", "stringr", "ComplexHeatmap", 
  "SingleCellExperiment", "ggmin", "Nourpal", "Cairo",
  "DoubletFinder", "harmony", "magick", "viridis", "limma", 
  "glmGamPoi", "gmailr", "rjson", "here", "SeuratDisk",
  "gplots", "clustifyr", "fgsea", "purrr", "clustree",
  "SingleR", "celldex", "ProjecTILs", "biomaRt", "data.table",
  "umap", "pheatmap", "scDblFinder", "miQC", "SeuratWrappers",
  "PerformanceAnalytics", "corrplot"
)

# similar to libaries, but will install package as well
packages(MyPackages)

# get experiment parameters
config <- fromJSON(file = here("config.json"))
master <- fromJSON(file = here("master.json"))

# Nour's palette. Install and load
NourpalDirectory <- config$NourpalDirectory
setwd(NourpalDirectory)
devtools::load_all("Nourpal.Rproj")
library(Nourpal)



################ SET PARAMETERS ###############

# resolution for cluster-finding after PCA
RESOLUTIONS <- config$RESOLUTIONS
RESOLUTION <- config$RESOLUTION

# if we are doing scTransform
SCT <- config$SCT
# font family for output plots
FONT_FAMILY <- config$FONT_FAMILY
# files to analyze
FILES = config$FILES

ADT_PRESENT <- master$ADT_PRESENT
USE_ADT <- config$USE_ADT
if (USE_ADT) {
  regexp <- "[[:digit:]]+"
  for (file in FILES) {
    if (! ADT_PRESENT[[str_extract(file, regexp)]]) {
      stop("Trying to use ADT, but ADT is not present for all samples")
    }
  }
}

# saving files
ObjName <- config$ObjName
Subset <- config$Subset
################# CONFIG DIRECTORIES #################
# Set working directory and create input / output folders
Directory <- config$Directory
if (! dir.exists(Directory)) dir.create(Directory)
setwd(Directory)

# DO NOT CHANGE RobjDir (needed for reading in .rds)
# kind of like extra level of specification
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
utils::browseURL(OutputDirectory)
if (! dir.exists(OutputDirectory)) dir.create(OutputDirectory)
ConfigDirectory <- paste0(OutputDirectory, "/Configs/")
if (! dir.exists(ConfigDirectory)) dir.create(ConfigDirectory)


QCDirectory <- paste0(OutputDirectory, "/QualityControl/")
if (! dir.exists(QCDirectory)) dir.create(QCDirectory)
RiboQCDirectory <- paste0(QCDirectory, "RiboQC/")
if (! dir.exists(RiboQCDirectory)) dir.create(RiboQCDirectory)
ElbowDirectory <- paste0(OutputDirectory, "/ElbowPlots/")
if (! dir.exists(ElbowDirectory)) dir.create(ElbowDirectory)
densityplotDirectory <- paste0(OutputDirectory, "/BarPlots/")
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

projTILDirectory <- paste0(OutputDirectory, "/projecTILs/")
if (! dir.exists(projTILDirectory) & RESOLUTION == "All") {
  dir.create(projTILDirectory)
} 
singleRDirectory <- paste0(OutputDirectory, "/singleR/")
if (! dir.exists(singleRDirectory) & RESOLUTION == "All") {
  dir.create(singleRDirectory)
} 
clustifyrDirectory <- paste0(OutputDirectory, "/clustifyr/")
if (! dir.exists(clustifyrDirectory) & RESOLUTION != "All") {
  dir.create(clustifyrDirectory)
} 

ADTDirectory <- paste0(OutputDirectory, "/ADT_plots/")
if (! dir.exists(ADTDirectory) & config$USE_ADT == TRUE) {
  dir.create(ADTDirectory)
}

SubtypeCorrDirectory <- paste0(OutputDirectory, "/SubtypeCorrelations/")
if (config$aggr_cells & config$preprocess_existing_RDS) {
  dir.create(SubtypeCorrDirectory)
}

##################### GENERIC HELPER FUNCTIONS ######################

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
  subtype_by,
  color_scheme = "nourpal",
  color_reverse = FALSE
) {
  if (color_scheme == "nourpal") {
    colors = Nour_pal("all", reverse = color_reverse)(
      length(unique(as.vector(as.matrix(seurat_object[[color_by]]))))
    )
    names(colors) = sort(unique(
      as.vector(as.matrix(seurat_object[[color_by]]))
    ))
  }
  else if (color_scheme == "chromatose") {
    colors <- list()
    samples <- unique(as.vector(as.matrix(seurat_object[[color_by]])))
    for (idx_subtype in 1:length(subtype_by)) {
      subtype <- subtype_by[[idx_subtype]]
      color_hex <- master$chromatose[[
        master$chromatose_keys[[idx_subtype]]
      ]]
      if (color_reverse) {
        color_hex <- rev(color_hex)
      }
      color_sample_names <- c()
      for (sample in samples){
        if (grepl(subtype, sample, fixed = TRUE)) {
          color_sample_names <- c(color_sample_names, sample)
        }
      }
      for (idx_sample in 1:length(color_sample_names)) {
        sample <- color_sample_names[idx_sample]
        idx_color <- length(color_hex) %/% (length(color_sample_names) + 1) * 
          (idx_sample) + 1
        
        colors[[sample]] <- color_hex[idx_color]
      }
    }
    colors <- colors[order(match(colors, samples))]
  }
  
  return (colors)
}

###### HEATMAP ANNOTATIONS ########
get_top_cluster_markers <- function(
  markers, n
) {
  top_n_markers <- markers %>% 
    group_by(cluster) %>% 
    top_n(n = n, wt = avg_log2FC)
  
  top_n_markers <- top_n_markers[order(top_n_markers$cluster),]
    
  
  top_n_markers$cluster <- as.character(top_n_markers$cluster)
  top_n_markers <- top_n_markers[order(top_n_markers$cluster), ]
  return (top_n_markers)
}



############### GET CASE CORRECT FEATURES ##############
case_sensitive_features <- function(
  seurat_object = NULL, features, reference = NULL, assay = "RNA"
) {
  final_features <- list()
  
  if (! is.null(seurat_object)) {
    DefaultAssay(seurat_object) <- assay
    reference <- rownames(seurat_object)
  }
  

  for (seurat_gene in reference) {
    idx_match <- match(
      tolower(seurat_gene),
      tolower(features),
      nomatch = -1
    )
    if (idx_match != -1) {
      final_features[[features[[idx_match]]]] <- seurat_gene
    }
  }
  
  return (final_features)
}
  
################ RENAME SEURAT FEATURE ROWS #############
# renaming function
# https://github.com/satijalab/seurat/issues/1049
# Replace gene names in different slots of a Seurat object. 
# It only changes obj@assays$integrated@counts, @data and @scale.data.
# also changes var.features of integrated assay
RenameGenesSeurat <- function(
  obj, newnames, var_features, assay_name = "integrated"
) { 
  my_assay <- obj[[assay_name]]
  
  if (nrow(my_assay) == length(newnames)) {
    
    if (length(my_assay@counts)) {
      rownames(my_assay@counts) <- newnames
    }
    if (length(my_assay@data)) {
      rownames(my_assay@data) <- newnames
    }
    if (length(my_assay@scale.data)) {
      rownames(my_assay@scale.data) <- newnames
    }
    
    # add new var features
    my_assay@var.features <- unlist(var_features, use.names=FALSE)
  } else {
    "Unequal gene sets: nrow(my_assay) != nrow(newnames)"
  }
  obj[[assay_name]] <- my_assay
  return(obj)
}


}


##################### PROJECTILS PREPROCESSING ######################
{
################ MOUSE TO HUMAN GENE BIOMART CONVERSION ###########
convertMouseGeneList <- function(x){
  if (!exists("MOUSE_BIOMART") | !exists("HUMAN_BIOMART")) {
    HUMAN_BIOMART <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    MOUSE_BIOMART <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  }
  genesV2 = getLDS(
    attributes = c("mgi_symbol"), 
    filters = "mgi_symbol", 
    values = x, 
    mart = MOUSE_BIOMART, 
    attributesL = c("hgnc_symbol"), 
    martL = HUMAN_BIOMART, 
    uniqueRows = T
  )
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}


########## PROJECTILS CUSTOM REF MAP RECOMP ###########
# https://carmonalab.github.io/ProjecTILs.demo/build_ref_atlas.html
projectils_ref_spec_dimreducs <- function(ref_projectils) {
  # recalc pca / umap using projecTILs prcomp and umap
  ### PCA ###
  varfeat <- ref_projectils@assays[["integrated"]]@var.features
  
  refdata <- data.frame(t(
    ref_projectils@assays[["integrated"]]@data[varfeat,]
  ))
  refdata <- refdata[, sort(colnames(refdata))]
  
  ref.pca <- prcomp(
    refdata, rank. = 50, scale. = TRUE, center = TRUE, retx = TRUE)
  
  ref.pca$rotation[1:5,1:5]
  
  
  ### UMAP ###
  seed <- 1234
  n.neighbors <- 30
  min.dist <- 0.3
  metric <- "cosine"
  ndim <- 16
  
  umap.config <- umap.defaults
  umap.config$n_neighbors <- n.neighbors
  umap.config$min_dist <- min.dist
  umap.config$metric <- metric
  umap.config$n_components <- 2
  umap.config$random_state <- seed
  umap.config$transform_state <- seed
  
  ref.umap <- umap(ref.pca$x[,1:ndim], config = umap.config)
  colnames(ref.umap$layout) <- c("UMAP_1","UMAP_2")
  
  ### add pca and umap info back into ref ### 
  ref_projectils@reductions$umap@cell.embeddings <- ref.umap$layout
  ref_projectils@reductions$pca@cell.embeddings <- ref.pca$x
  ref_projectils@reductions$pca@feature.loadings <- ref.pca$rotation
  colnames(ref_projectils@reductions$pca@cell.embeddings) <- 
    gsub("PC(\\d+)", "PC_\\1", colnames(ref.pca$x), perl = TRUE)
  colnames(ref_projectils@reductions$pca@feature.loadings) <- 
    gsub("PC(\\d+)", "PC_\\1", colnames(ref.pca$rotation), perl = TRUE)
  
  #Store the complete PCA and UMAP object in @misc
  ref_projectils@misc$pca_object <- ref.pca
  ref_projectils@misc$umap_object <- ref.umap
  ref_projectils@misc$projecTILs <- "custom_atlas"
  
  return (ref_projectils)
}

############### PROJECTILS MURINE TO HUMAN #############
projectils_ref_full_translation <- function (
  ref_projectils, 
  SeuratObj,
  ref_genes_human_filename = paste0(
    dataDirectory, "projecTILs/", analyses$projecTILs_ref,
    "_mouse_to_human_genes.csv"
  )
) {
  # we will be recalculating almost everything except integrated assay
  DefaultAssay(ref_projectils) <- "integrated"
  ref_projectils <- DietSeurat(  
    ref_projectils,
    assays = "integrated"
  )
  ref_projectils@misc <- list()
  
  
  
  # start translation
  ref_genes_mouse <- rownames(ref_projectils)
  
  # translate all mouse genes to human genes
  if (file.exists(ref_genes_human_filename)) {
    # read the precalculated ref_genes_human via csv
    ref_genes_human <- as.list(read.csv(
      ref_genes_human_filename
    ))
  } else {
    # translate reference genes to human
    # (may already be human, but might still be mouse...)
    ref_genes_human <- list()
    for (mouse_gene in ref_genes_mouse) {
      human_gene <- convertMouseGeneList(mouse_gene)
      if (length(human_gene) != 0) {
        # biomart human >> ref mouse
        ref_genes_human[[human_gene[1]]] <- mouse_gene
      }
    }
    # this takes a long time, so write to csv:
    write.table(
      ref_genes_human,
      file = ref_genes_human_filename,
      quote = FALSE, sep = ",", row.names = TRUE
    )
  }
  
  # ref_genes_human: biomart human >> mouse gene
  
  # translate biomart human to Seurat human identifiers
  # genes_untranslated: biomart human >> Seurat human
  genes_untranslated <- case_sensitive_features(
    SeuratObj, names(ref_genes_human)
  )
  # we ultimately want 
  # genes_translated: mouse gene >> Seurat human
  genes_translated <- list()
  for (human_gene_ref in names(genes_untranslated)) {
    genes_translated[[ref_genes_human[[human_gene_ref]]]] <-
      genes_untranslated[[human_gene_ref]]
  }

  
  # translate reference gene names, in order
  ref_to_Seurat <- list()
  for (ref_gene in ref_genes_mouse) {
    if (ref_gene %in% names(genes_translated)) {
      ref_to_Seurat <- append(ref_to_Seurat, genes_translated[[ref_gene]])
    } else {
      # no translation, so use original name
      # we will exclud these genes later via var.features
      ref_to_Seurat <- append(ref_to_Seurat, ref_gene)
    }
  }
  
  # rename feature rows across entire rej_projectils object
  ref_projectils <- RenameGenesSeurat(
    obj = ref_projectils, 
    newnames = ref_to_Seurat, 
    var_features = genes_translated
  )
  
  # scale data
  ref_projectils <- ScaleData(ref_projectils, assay = "integrated")
  
  # compute pca and umaps, using Seurat methods
  ref_projectils <- GEX_pca(ref_projectils, " projecTILs pca rerun ")
  E1 <- ElbowPlot(ref_projectils)
  pdf(paste0(
    dataDirectory, "projecTILs/ref", analyses$projecTILs_ref,
    " projecTILs rerun pca.pdf"
  ), width = 12, height = 4.5, family = FONT_FAMILY
  )
  E1
  dev.off()
  
  # umap
  ref_projectils <- RunUMAP(
    ref_projectils, 
    dims = c(1:20)
  )
  
  # plot
  UprojTILs <-  plot_umap(
    ref_projectils, group_by = "functional.cluster",
    title = "projecTILs reference with filtered Seurat genes", 
    xlab = "UMAP_1", ylab = "UMAP2",
    legend_position = "right", reduction = "umap", 
    label_clusters = TRUE
  )
  
  pdf(paste0(
    dataDirectory, "projecTILs/ref", analyses$projecTILs_ref,
    " projecTILs recalculated with human genes UMAP.pdf"
  ), width = 12, height = 7, family = FONT_FAMILY
  )
  UprojTILs
  dev.off()
  
  # custom dimreducs, following projecTILs specs
  ref_projectils <- projectils_ref_spec_dimreducs(ref_projectils)
  
  
  saveRDS(
    ref_projectils, 
    file = paste0(
      dataDirectory, "projecTILs/", analyses$projecTILs_ref,
      "_processed_mouse2human_reference.rds"
    )
  )
  return (ref_projectils)
}


}

############################# PLOTTING FUNCTIONS ###############################
{
################ TEMPLATE GGPLOT BARPLOT ###################

plot_bargraph <-  function (
  seurat_object, aesX, fill,
  y_label, x_label, y_lower_limit, y_break,
  title = NULL,
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
    labs(y = y_label, x = x_label, title = title) +
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
  alpha = 0.2, color_reverse = FALSE,
  facet_category = NULL,
  title = NULL
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
  
  if (!is.null(facet_category)) {
    plot_object <- plot_object + 
      facet_wrap(as.formula(paste("~", facet_category)))
  }
  
  if (!is.null(title)) {
    plot_object <- plot_object + 
      labs(title = title)
  }
  
  return (plot_object)
}



################ UMAP PLOT TEMPLATE ##################
plot_umap <- function(
  seurat_object, group_by,
  title, xlab, ylab,
  legend_position, reduction,
  subtype_by,
  color_scheme = "nourpal",
  color_reverse = FALSE,
  title_font_size = 20, x_font_size = 20, y_font_size = 20, 
  pt_size = NULL, split_by = NULL, ncol_dimplot = 1, ncol_guide = 1,
  label_clusters = FALSE, repel_labels = FALSE, label_size = 4
) {
  color_by <- group_by
  # get colors
  manual_colors <- get_colors(
    seurat_object = seurat_object, 
    color_by = color_by, 
    color_scheme = color_scheme,
    subtype_by = subtype_by,
    color_reverse = color_reverse
    )

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
  color_palette_option = "plasma",
  assay = "RNA", features_sorted = FALSE
) {
  DefaultAssay(seurat_object) <- assay
  if (features_sorted == FALSE) {
    features <- case_sensitive_features(
      seurat_object, features, assay = assay
      )
    features <- unlist(unname(features))
  }
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
  order = TRUE, label = FALSE, label_size = 7, repel = FALSE,
  cols = Nour_cols(c("darkpurple", "lightorange")), pt_size = 0.3,
  legend_title_position = "top", legend_title_size = 10,
  legend_location = "right",
  ncol = 1, nrow = 1, widths = c(0.03, 1, 1),
  assay = "RNA"
) {
  DefaultAssay(seurat_object) <- assay
  
  feature_gene <- case_sensitive_features(
    seurat_object,
    c(feature_gene),
    assay = assay
  )

  feature_gene <- unlist(unname(feature_gene))
  if (length(feature_gene) == 0) {
    stop(paste0("could not find gene feature: ", feature_gene))
  } else {
    feature_gene <- feature_gene[[1]]
  }
  
  F1 <- FeaturePlot(
    seurat_object, features = feature_gene, split.by = split_by,
    order = order, cols = cols, pt.size = pt_size, reduction = reduction,
    label = label, repel = repel, label.size = label_size
  )
  
  legend <- get_legend(
    F1 + NoAxes() + theme_min() +
      guides(colour = guide_colourbar(
        title = feature_gene, title.position = legend_title_position,
        title.theme = element_text(size = legend_title_size)))
  )
  
  plot_object <- ggpubr::ggarrange(
    ggparagraph(text = " ",  size = 0),
    F1 + NoLegend() + NoAxes()
    # + theme_min2()
    ,
    ncol = ncol, nrow = nrow, widths = widths,
    legend.grob = legend, legend = legend_location
  )
  
  return (plot_object)
}


##################### TEMPLATE ADT CORRELATION ########################
plot_correlation <- function(
  SeuratObj, 
  x, y, split_by, lab_x, lab_y,
  rna_threshold = 0,
  adt_threshold = 0,
  slot = "scale.data",
  method = "pearson"
) {
  # translate to features
  x <- paste0(
    "rna_",
    case_sensitive_features(SeuratObj, features = x, assay = "RNA")
  )
  y <- paste0(
    "adt_",
    case_sensitive_features(SeuratObj, features = y, assay = "ADT")
  )
  
  df <- FetchData(
    SeuratObj, vars = c(x, y, split_by), slot = slot
  )
  df <- df[df[x] > rna_threshold & df[y] > adt_threshold,]
  
  plot <- ggplot(
    df, 
    aes_string(
      x = paste0("`", x, "`"),
      y = paste0("`", y, "`")
    )
  ) + 
    geom_point()+
    geom_smooth(method = lm) +
    facet_wrap(
      as.formula(paste("~", split_by))
    ) + stat_cor(method = method) + 
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.text.y = element_text(color = "black", size = 12),
      axis.text.x = element_text(
        color = "black", hjust = 1, size = 15
      ),
      title = element_text(size = 15, face = "bold"),
      axis.text = element_text(size = 15),
      axis.title = element_text(size = 15, face = "bold"),
      plot.margin = unit(c(0.5, 0.5, 0.2, 0.5), "cm"),
      strip.text.x = element_text(size = 15),
      strip.background = element_rect(fill = "white", colour = "black", size = 1)
    ) +
    labs(x = lab_x, y = lab_y, title = "ADT vs RNA")
  
  return(plot)
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
  if (grepl("RNA", cluster, fixed = TRUE)) {
    sorted_barcodes <- names(sort(subset$ClusterRNA))
  } else if (grepl("ADT", cluster, fixed = TRUE)) {
    sorted_barcodes <- names(sort(subset$ClusterADT))
  } else if (grepl("WNN", cluster, fixed = TRUE)) {
    sorted_barcodes <- names(sort(subset$ClusterWNN))
  }
  
  # plot data rows are genes, cols are cells
  plot_data <- as.data.frame(subset@assays$RNA@scale.data)
  plot_data <- plot_data[top_n_markers$marker, ]
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
      at = match(label_n_markers$marker, row.names(plot_data)),
      labels = label_n_markers$marker,
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
###################### GEX PREPROCESSING FUNCTIONS ###################
# created here, so as to avoid code duplication
{
############ SCDBLFINDER & QC ##############
GEX_QC <- function(SeuratObjMYSC, file) {
  # determine doublets, need to change to sce
  SeuratObjMYSC.sce <- as.SingleCellExperiment(SeuratObjMYSC)
  sce <- scDblFinder(SeuratObjMYSC.sce)
  
  # Converting back to Seurat Obj
  SeuratObjMYSC <- as.Seurat(sce)
  Idents(SeuratObjMYSC) <- SeuratObjMYSC$orig.ident
  DefaultAssay(SeuratObjMYSC) <- "RNA"
  
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
  SeuratObjMYSC <- GEX_pca(SeuratObjMYSC, file)
  
  
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
      QCDirectory, file,
      " Scatterplot - Determining Probability of Low Quality Cells.pdf"
    ), width = 8, height = 5.5, family = "ArialMT")
    print(QC1)
    dev.off()
    
    QC2 <- PlotMiQC(SeuratObjMYSC, color.by = "miQC.keep")
    pdf(paste0(
      QCDirectory, file,
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
    QCDirectory, file,
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
  
  return (SeuratObjMYSC)
}
  
############ NORMALIZATION ##############
GEX_normalization <- function(SeuratObj){
  
  # SCT for GEX
  if (SCT) {
    # rename for ease of use; no matter if SCT or NormalizeData,
    # we would like to use "RNA" assay from now on
    SeuratObj <- SCTransform(
      SeuratObj,
      assay = "RNApreSCT",
      new.assay.name = "RNA",
      method = "glmGamPoi", verbose = FALSE,
      return.only.var.genes = FALSE
    )
    
  } else {
    SeuratObj <- SeuratObj %>% 
      NormalizeData(assay = "RNA") %>%
      FindVariableFeatures(assay = "RNA") %>%
      ScaleData(assay = "RNA")
  }
  
  DefaultAssay(SeuratObj) <- "RNA"
  
  # plot variable features
  top10 <- head(VariableFeatures(SeuratObj), 10, assay = "RNA")
  plot1 <- VariableFeaturePlot(SeuratObj)
  plot2 <- LabelPoints(plot = plot1, points = top10, size = 3)
  pdf(paste0(
    QCDirectory, ObjName, Subset, 
    " GEX variable freatures.pdf"
  ), width = 8, height = 5.5, family = FONT_FAMILY
  )
  print(plot2 + theme_min())
  dev.off()
  return (SeuratObj)
}

############ PCA ##############
GEX_pca <- function(SeuratObjMYSC, file) {
  DefaultAssay(SeuratObjMYSC) <- "RNA"
  SeuratObjMYSC <- RunPCA(
    SeuratObjMYSC,
    features = VariableFeatures(object = SeuratObjMYSC)
  )
  # most representative (by absolute value) genes for PCs
  VS = VizDimLoadings(SeuratObjMYSC, dims = 1:2, reduction = "pca")
  VSD = DimPlot(SeuratObjMYSC, reduction = "pca") + theme_min()
  
  # save plot
  pdf(paste0(
    QCDirectory, file, " GEX PCA features and plot.pdf"
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
  return (SeuratObjMYSC)
}

######## LOUVAIN CLUSTERING (GEX) ########
GEX_louvain <- function(
  SeuratObj, resolution,
  reduction = "harmonyRNA"
) {
  DefaultAssay(SeuratObj) <- "RNA"
  
  SeuratObj <- SeuratObj %>%
    FindNeighbors(
      reduction = reduction, 
      dims = c(1:analyses$harmonyRNA_dims)
    ) %>%
    FindClusters(resolution = resolution)
  
  return (SeuratObj)
}


}


###################### ADT + WNN PREPROCESSING FUNCTIONS ###################
{
#################### ADT BATCH INTEGRATION #################
  
ADT_integrate <- function(SeuratObj) {  
  DefaultAssay(SeuratObj) <- "ADTpreInt"
  # normalize per sample
  # using CLR normalization based on: 
  # https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html
  samples_list <- SplitObject(SeuratObj, split.by = config$batch_norm_by) %>%
    lapply(
      FUN = function(x) {
        x <- NormalizeData(
          x, normalization.method = "CLR", margin = 2, assay = "ADTpreInt"
        )
      }
    )
  # select features that are repeatedly variable across data sets
  features <- SelectIntegrationFeatures(
    object.list = samples_list
  )
  
  samples_list <- lapply(
    X = samples_list, 
    FUN = function(x) {
      x <- ScaleData(x, features = features, verbose = FALSE)
      x <- RunPCA(x, features = features, verbose = FALSE)
    })
  
  # find anchor cells (this step stakes a while)
  # use rpca; reference:
  # https://satijalab.org/seurat/archive/v3.2/integration.html
  cell_anchors <- FindIntegrationAnchors(
    object.list = samples_list, anchor.features = features, reduction = "rpca"
  )
  # combine into integrated Seurat object
  SeuratObj_ADT <- IntegrateData(
    anchorset = cell_anchors, 
    normalization.method = "LogNormalize",
    new.assay.name = "ADT"
  )
  
  # scale data & dimreduc PCA, in prep for WNN
  DefaultAssay(SeuratObj_ADT) <- "ADT"
  SeuratObj_ADT <- SeuratObj_ADT %>%
    ScaleData() %>%
    RunPCA(reduction.name = "pcaADT")
  
  # add ADT data back into original SeuratObj
  # necessary bc IntegrateData destroys scale.data for all assays
  # ignore any warnings about offending keys
  # ref: https://github.com/satijalab/seurat/issues/3843
  SeuratObj[["ADT"]] <- SeuratObj_ADT[["ADT"]]
  SeuratObj[["pcaADT"]] <- SeuratObj_ADT[["pcaADT"]]
  
  E2 <- ElbowPlot(SeuratObj, ndims = 15, reduction = "pcaADT")
  pdf(paste0(
    ElbowDirectory, ObjName, Subset, 
    " elbow plot after CC scaling integration and pca (ADT).pdf"
  ), width = 12, height = 4.5, family = FONT_FAMILY
  )
  print(E2)
  dev.off()
  
  return (SeuratObj)
}
  
  
####################### ADT LOUVAIN #####################
ADT_louvain <- function(
    SeuratObj, resolution,
    reduction = "pcaADT"
) {
  DefaultAssay(SeuratObj) <- "ADT"
  
  SeuratObj <- SeuratObj %>%
    FindNeighbors(
      reduction = reduction, 
      dims = c(1:analyses$pcaADT_dims)
    ) %>%
    FindClusters(resolution = resolution)
  
  return (SeuratObj)
}


####################### WNN LOUVAIN #####################
WNN_louvain <- function(
  SeuratObj, resolution, reduction_list = list("harmonyRNA", "pcaADT"),
  algorithm = 3, verbose = FALSE
) {
  print(paste0("Louvain for res", resolution))
  SeuratObj <- SeuratObj %>%
    FindMultiModalNeighbors(
      reduction.list = reduction_list, 
      dims.list = list(
        1 : analyses$harmonyRNA_dims, 1 : analyses$pcaADT_dims
      )
    ) %>%
    FindClusters(
      graph.name = "wsnn",
      algorithm = algorithm,
      resolution = resolution,
      verbose = verbose
    )
  return (SeuratObj)
}

}
###################### VIZ FUNCTIONS ###################
{
################ FIND CLUSTER BIOMARKERS (GEX + ADT) ################
cluster_markers <- function(
  SeuratObj, assay = "RNA", default_ident = "ClusterRNA"
  ) {
  DefaultAssay(SeuratObj) <- assay
  Idents(SeuratObj) <- default_ident
  if (SCT & assay == "RNA") {
    SeuratObj <- PrepSCTFindMarkers(SeuratObj, assay = assay)
  }
  markers <- FindAllMarkers(
    SeuratObj,
    assay = assay,
    logfc.threshold = 0.5, 
    test.use = "wilcox", 
    only.pos = TRUE
  )
  markers$marker <- markers$gene
  markers$gene <- NULL
  return (markers)
  }
}

####################### ATLAS FUNCTIONS ######################
add_atlas_TMB <- function(SeuratObj) {
  CellInfo <- SeuratObj@meta.data
  for (sample in unique(SeuratObj$Sample)) {
    tmb_file_contents <- readLines(paste0(
      config$TMBDirectory, sample, "-C_1.fq.gz_snp_indel.TMB_results.txt"
      ))
                                   
    tmb <- strsplit(
     tmb_file_contents[grepl("TMB= ", tmb_file_contents)], "TMB= "
    )[[1]][2]
    tmb <- as.double(tmb)
    
    CellInfo[["TMB_SampleWide"]][
      CellInfo$Sample == sample
    ] <- tmb
  }

}

atlas_QC <- function() {
  # make sure we are indeed doing QC on a predefined atlas
  stopifnot(config$preprocess_existing_RDS)
  
  # to contain preprocessed samples
  SeuratSamples <- list()
  # load atlas
  seurat_object <- readRDS(paste0(
    RobjDir, config$existing_RDS_name
  ))
  # split by sample
  seurat_list <- SplitObject(seurat_object, split.by = "Sample")
  
  for (idx in 1:length(seurat_list)) {
    # get sample
    seurat_split <- seurat_list[idx][[1]]
    sample_name <- unique(seurat_split$Sample)[1]
    print(paste0("Sample", sample_name))
    # make sure there are enough cells and it's not a sample to be removed
    if (
      ncol(seurat_split) > 55 & 
      ! sample_name %in% config$existing_RDS_remove_cluster
    ) {
      # QC
      seurat_split <- GEX_QC(seurat_split, sample_name)

      # save, add to list
      saveRDS(
        seurat_split,
        file = paste0(RobjDirectory, sample_name, ".rds")
      )
      SeuratSamples[[length(SeuratSamples) + 1]] <- seurat_split
      config$FILES[[length(SeuratSamples)]] <<- sample_name
    } else {
      print ("not enough cells, moving onto next sample")
    }
  }
  
  return (SeuratSamples)
}
