################# INSTALL / LOAD PACKAGES  ##############
if (F) {
  # R -f test.R
  # install them first if necessary. Only need to install once
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install(version = "3.14")
  
  install.packages("devtools")
  devtools::install_github("jakesherman/easypackages")
  
  
  library(devtools)
  library(easypackages)

  # only need to install once
  remotes::install_github("satijalab/sctransform", ref = "develop")
  BiocManager::install(c(
    'BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
    'limma', 'S4Vectors', 'SingleCellExperiment',
    'SummarizedExperiment', 'batchelor', 'Matrix.utils'
  ))

  
  devtools::install_github("jokergoo/ComplexHeatmap")
  devtools::install_github('cole-trapnell-lab/leidenbase')
  devtools::install_github('cole-trapnell-lab/monocle3')
  BiocManager::install("ComplexHeatmap")
  BiocManager::install("clustifyr")
  BiocManager::install("slingshot")
  # gamma poisson generalized linear model
  BiocManager::install("glmGamPoi")
  BiocManager::install("celldex")
  BiocManager::install("SingleR")
  BiocManager::install("biomaRt")
  BiocManager::install("scDblFinder")
  BiocManager::install("miQC")
  BiocManager::install("GEOquery")
  
  BiocManager::install(c(
    "clustifyr", "slingshot", "glmGamPoi", "celldex", "SingleR",
    "biomaRt", "scDblFinder", "miQC", "GEOquery"
    ))
  
  BiocManager::install(c("celldex", "biomaRt"))
  
  install.packages("ggpubr")
  install.packages('Seurat')
  install.packages("tidyverse")
  install.packages('clustree')
  install.packages("pheatmap")
  install.packages("corrplot")
  install.packages("extrafont")
  
  
  remotes::install_github("mojaveazure/seurat-disk")
  remotes::install_github("carmonalab/UCell", ref = "v1.3")
  remotes::install_github("carmonalab/scGate")
  remotes::install_github("carmonalab/ProjecTILs")
  devtools::install_github('satijalab/seurat-wrappers')
  devtools::install_github("sjessa/ggmin")
  install.packages("magick")
  install.packages("rjson")
  install.packages("PerformanceAnalytics")
  MyPackages <- c(
    "dplyr", "ggplot2", 
    # "ggpubr", 
    "gridExtra", "viridis", 
    "egg",
    "grid", "lattice", "gtools", "Biobase", "RColorBrewer", "tibble",
    "Seurat", "cowplot", "patchwork", "stringr", "ComplexHeatmap", 
    "SingleCellExperiment", 
    # "ggmin", 
    # "Nourpal", 
    # "Cairo",
    "harmony", 
    # "magick", 
    "viridis", "limma", 
    "glmGamPoi", "rjson", "here", 
    # "SeuratDisk",
    "gplots", "clustifyr", "purrr", "clustree",
    "SingleR", 
    # "celldex", 
    # "ProjecTILs", 
    # "biomaRt", 
    "data.table",
    "umap", 
    "pheatmap", "scDblFinder", "miQC", 
    # "SeuratWrappers",
    "PerformanceAnalytics", "corrplot", "GEOquery",
    "slingshot", "ggbeeswarm", 
    # "monocle3", 
    "extrafont"
    # "jpeg"
  )
  
}

MyPackages <- c(
  "dplyr", "ggplot2", "ggpubr", "gridExtra", "viridis", "egg",
  "grid", "lattice", "gtools", "Biobase", "RColorBrewer", "tibble",
  "Seurat", "cowplot", "patchwork", "stringr", "ComplexHeatmap", 
  "SingleCellExperiment", "ggmin", "Nourpal", "Cairo",
  "harmony", "magick", "viridis", "limma", 
  "glmGamPoi", "rjson", "here", "SeuratDisk",
  "gplots", "clustifyr", "fgsea", "purrr", "clustree",
  "SingleR", "celldex", "ProjecTILs", "biomaRt", "data.table",
  "umap", "pheatmap", "scDblFinder", "miQC", "SeuratWrappers",
  "PerformanceAnalytics", "corrplot", "GEOquery",
  "slingshot", "ggbeeswarm", "monocle3", "extrafont",
  "jpeg"
)

library(devtools)
library(easypackages)

# similar to libaries, but will install package as well
packages(MyPackages)
loadfonts()

# get experiment parameters
config <- fromJSON(file = here("config.json"))
master <- fromJSON(file = here("master.json"))
analyses <- fromJSON(file = here("analysis.json"))

# Nour's palette. Install and load
NourpalDirectory <- config$NourpalDirectory
setwd(NourpalDirectory)
devtools::load_all("Nourpal.Rproj")
library(Nourpal)



################ CONFIRM + SET + SAVE PARAMETERS ###############

# make sure viz_clustering is good
stopifnot (analyses$viz_clustering %in% c("RNA", "ADT", "WNN"))

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
  for (file in FILES) {
    if (! ADT_PRESENT[[file]]) {
      stop("Trying to use ADT, but ADT is not present for all samples")
    }
  }
}

# saving files
ObjName <- config$ObjName
Subset <- config$Subset

write_experimental_configs <- function(suffix = "") {
  writeLines(
    capture.output(sessionInfo()), 
    paste0(ConfigDirectory, ObjName, "_", Subset, "_sessionInfo.txt")
  )
  
  
  config_params_filename <- paste0(
    ConfigDirectory, ObjName, "_", Subset, "_config_params", suffix, ".json"
  )
  if (file.exists(config_params_filename)) {
    unlink(config_params_filename)
  }
  file.copy(
    from = here("config.json"), to = config_params_filename
  )
  
  
  analysis_params_filename <- paste0(
    ConfigDirectory, ObjName, "_", Subset, "_analysis_params", suffix, ".json"
  )
  if (file.exists(analysis_params_filename)) {
    unlink(analysis_params_filename)
  }
  file.copy(
    from = here("analysis.json"), to = analysis_params_filename
  )
}

################# CONFIG DIRECTORIES #################
# Set working directory and create input / output folders
Directory <- config$Directory
if (! dir.exists(Directory)) dir.create(Directory)
setwd(Directory)

GEOdir <- paste0(Directory, "/", analyses$GEO, "/")
if (! dir.exists(GEOdir)) dir.create(GEOdir)

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
# open folder to ensure it's correct
utils::browseURL(
  paste0("Output/", Subset, "/", ObjName, "Res", config$RESOLUTION, "/")
  )
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
if (RESOLUTION != "All") {
  projTILDirectory <- str_replace(
    projTILDirectory, paste0("Res", RESOLUTION), paste0("Res", "All")
  )
}
if (! dir.exists(projTILDirectory) & RESOLUTION == "All") {
  dir.create(projTILDirectory)
} 
singleRDirectory <- paste0(OutputDirectory, "/singleR/")
if (! dir.exists(singleRDirectory)) {
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
if (! dir.exists(SubtypeCorrDirectory)) {
  dir.create(SubtypeCorrDirectory)
}

PseudoDirectory <- paste0(OutputDirectory, "/Pseudotime/")
if (! dir.exists(PseudoDirectory)) {
  dir.create(PseudoDirectory)
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
  subtype_by = c("blood", "tumor"),
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
  } else if (color_scheme == "nourpal_hot_cool") {
    stopifnot(length(subtype_by) == 2)
    hues <- c("cool", "hot")
    samples <- unique(as.vector(as.matrix(seurat_object[[color_by]])))
    colors <- list()
    for (idx_subtype in 1:length(subtype_by)) {
      subtype <- subtype_by[[idx_subtype]]
      subtype_samples <- c()
      for (sample in samples) {
        if (grepl(subtype, sample, fixed = TRUE)) {
          subtype_samples <- c(subtype_samples, sample)
        }
        subtype_colors <- Nour_pal(hues[idx_subtype])(length(subtype_samples))
        names(subtype_colors) <- subtype_samples
      }
      colors <- c(colors, subtype_colors)
      }
    } else if (color_scheme == "chromatose") {
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


################ RENAME CLUSTERRNA OR CLUSTERADT OR CLUSTERWNN #######
add_cluster_metadata <- function (SeuratObj) {
  # i.e. cluster_res_column = "RNA_snn_res.0.5"
  if (analyses$viz_clustering == "RNA") {
    cluster_res_column = paste0("RNA_snn_res.", RESOLUTION)
  } else if (analyses$viz_clustering == "ADT") {
    cluster_res_column = paste0("ADT_snn_res.", RESOLUTION)
  } else if (analyses$viz_clustering == "WNN") {
    cluster_res_column = paste0("wsnn_res.", RESOLUTION)
  }
  
  cluster_name <- get_cluster_name()
  
  CellInfo <- SeuratObj@meta.data
  
  # Rename Clusters; i.e. WNNC01, WNNC02
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
  
  return (SeuratObj)
}

############ ADD MASTER METADATA ##########
add_master_metadata <- function(SeuratObj, info, colname) {
  CellInfo <- SeuratObj@meta.data
  CellInfo[[colname]] <- NA
  for (i in 1:length(info)) {
    CellInfo[[colname]][
      which(str_detect(CellInfo$orig.ident, names(info[i])))
    ] <- info[[i]]
  }
  
  NewData <- CellInfo[[colname]]
  names(NewData) <- colnames(x = SeuratObj)
  SeuratObj <- AddMetaData(
    object = SeuratObj, metadata = NewData, col.name = colname
  )
}
}

add_all_master_metadata <- function(SeuratObj) {
  # this func adds metadata as defined by master.json file
  
  SeuratObj <- add_master_metadata(SeuratObj, master$sample_info, "Sample")
  SeuratObj <- add_master_metadata(SeuratObj, master$patient_info, "Patient")
  SeuratObj <- add_master_metadata(SeuratObj, master$group_info, "Group")
  SeuratObj <- add_master_metadata(SeuratObj, master$grade_info, "Grade")
  SeuratObj <- add_master_metadata(SeuratObj, master$type_info, "Type")
  SeuratObj <- add_master_metadata(SeuratObj, master$sex_info, "sex")
  
  # Fragment / IDFrag
  add_single_fragment_postfix <- function(orig_ident) {
    return (paste0(orig_ident, "_1"))
  }
  Fragment <- unlist(lapply(SeuratObj$Patient, add_single_fragment_postfix))
  SeuratObj <- AddMetaData(
    object = SeuratObj, metadata = Fragment, col.name = "Fragment"
  )
  IDFrag <- unlist(lapply(SeuratObj$Sample, add_single_fragment_postfix))
  SeuratObj <- AddMetaData(
    object = SeuratObj, metadata = IDFrag, col.name = "IDFrag"
  )
}

add_all_master_metadata_multiple_objs <- function(resolutions, subsets) {
  # resolutions = c("All", "0.3")
  # subsets = 
  for (subset in subsets) {
    for (res in resolutions) {
      filename <- paste0(
        RobjDir, subset, "/",
        "GEX", subset, "_res", res, ".rds"
      )
      print(filename)
      SeuratObj <- readRDS(filename)
      add_all_master_metadata(SeuratObj)
      saveRDS(SeuratObj, filename)
    }
  }
}



################## PROJECTILS AND SINGLER AND CLUSTIFYR ###################
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
  UprojTILs <- plot_umap(
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


projecTILs_wrapper <- function(SeuratObj, ref_projectils, gating = TRUE) {
  
  # TIME TO ACTUALLY ASTRAL PROJECT
  # project per sample. Don't project on integrated samples
  DefaultAssay(SeuratObj) <- analyses$projecTILs_assay
  manual_colors <- get_colors(ref_projectils, color_by = "functional.cluster")
  projectils_metadata <- data.frame()
  for (sample_name in unique(SeuratObj$Sample)) {
    print(sample_name)
    # subset based on sample and assignment
    query.data <- subset(
      SeuratObj, Sample == sample_name
    )
    
    # proJECT
    query.projected <- make.projection(
      query.data, 
      ref = ref_projectils,
      filter.cells = gating
    )
    query.projected <- cellstate.predict(
      ref <- ref_projectils, query = query.projected
    )
    
    # save metadata to add back to SeuratObj later
    sample_metadata <- as.data.frame(
      query.projected[[c('functional.cluster', 'functional.cluster.conf')]]
    )
    projectils_metadata <- rbind(projectils_metadata, sample_metadata)
    
  }
  
  # add functional.cluster and functional.cluster.conf info
  SeuratObj <- AddMetaData(
    SeuratObj,
    metadata = projectils_metadata,
    col.name = c("projecTILs", "projecTILs_conf")
  )
  
  if (
    ! is.na(analyses$which_assignment) &
    analyses$which_assignment %in% colnames(SeuratObj@meta.data)
  ) {
    SeuratObj$Assignment <- SeuratObj[[analyses$which_assignment]]
  }
  return (SeuratObj)
}

################# READ PROJECTILS REFERENCE ##################
load_ref_projectils <- function(SeuratObj = NULL) {
  # load murine TIL reference
  ref_projectils_filename <- paste0(
    dataDirectory, "projecTILs/", analyses$projecTILs_ref,
    "_processed_mouse2human_reference.rds"
  )
  print(paste0(
    "trying to read projecTILs reference: ",
    ref_projectils_filename
  ))
  # ideally, it has already been pre-processed from murine to human
  if (file.exists(ref_projectils_filename)) {
    ref_projectils <- readRDS(ref_projectils_filename)
  } else {
    print("projecTILs reference has not been processed, processing now")
    # we will preprocess from murine to human
    # also includes recomp of umap / pca according to projecTILs spec
    # first read in the original reference
    ref_projectils <- readRDS(paste0(
      dataDirectory, "projecTILs/", analyses$projecTILs_ref,
      "_mouse_atlas.rds"
    ))
    ref_projectils <- projectils_ref_full_translation (
      ref_projectils, SeuratObj
    )
  }
  return (ref_projectils)
}


############### SINGLER WRAPPER #############
SingleR_wrapper <- function(
    SeuratObj, ref_singler, 
    ref_name = "BlueprintEncode", celltype = "T cell"
    ) {
  singler_colname <- paste0("SingleR_", ref_name)
  
  singler_predictions <- SingleR(
    test = SeuratObj[[analyses$SingleR_assay]]@data,
    ref = ref_singler,
    labels = ref_singler$label.main
  )
  
  # add metadata back to Seurat
  SeuratObj@meta.data[[singler_colname]] <- singler_predictions[["pruned.labels"]]
  
  # cell type frequency table
  write.csv(
    table(SeuratObj@meta.data[[singler_colname]]),
    paste0(
      singleRDirectory,
      analyses$SingleR_assay, "assay SingleR",
      " cell type distribution.csv"
    ),
    row.names = FALSE
  )
  
  # plot total distribution
  P2 <- plot_bargraph (
    seurat_object = SeuratObj, aesX = "Sample", fill = singler_colname,
    y_label = "Composition (percentage of cells)", x_label = NULL,
    title = paste0(
      "SingleR ", celltype, " assignment using ", ref_name, " reference"
      ),
    y_lower_limit = 0, y_break = 1000, position = "fill"
  )
  
  pdf(paste0(
    singleRDirectory,
    celltype, " assignment using ", ref_name, " reference_", 
    analyses$SingleR_assay, " assay_",
    "percentage of cells per sample and Assignment barplot.pdf"
  ), width = 6, height = 5.5, family = FONT_FAMILY
  )
  print(P2)
  dev.off()
  
  # QC pruning
  P3 <- plotDeltaDistribution(singler_predictions)
  pdf(paste0(
    singleRDirectory,
    celltype, " assignment using ", ref_name, " reference_",
    analyses$SingleR_assay, "assay_",
    " pruning distribution.pdf"
  ), width = 6, height = 5.5, family = FONT_FAMILY
  )
  print(P3)
  dev.off()
  
  # QC heatmap
  P4 <- plotScoreHeatmap(singler_predictions)
  pdf(paste0(
    singleRDirectory,
    celltype, " assignment using ", ref_name, " reference_",
    analyses$SingleR_assay, "assay_",
    " score heatmap.pdf"
  ), width = 6, height = 5.5, family = FONT_FAMILY
  )
  print(P4)
  dev.off()
  
  
  if (
    ! is.na(analyses$which_assignment) &
    analyses$which_assignment %in% colnames(SeuratObj@meta.data)
  ) {
    SeuratObj$Assignment <- SeuratObj[[analyses$which_assignment]]
    if (! grepl("SingleR", analyses$which_assignment, fix = F)) {
      print(paste0(" not using SingleR for cell 'Assignment' ", ref_name))
    }
  }
  

  return (SeuratObj)
}

################ DEFINE SINGLER REFERENCES ##################
define_refs_singler <- function() {
  refs_singler <<- list(
    "BlueprintEncode" = celldex::BlueprintEncodeData(),
    "HPCA" = celldex::HumanPrimaryCellAtlasData(),
    "Novershtern" = celldex::NovershternHematopoieticData()
  )
}

################ CLUSTIFYR WRAPPER ################
clustifyr_wrapper <- function(SeuratObj, REF_MATRIX, clustifyr_colname) {
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
  
  SeuratObj@meta.data[clustifyr_colname] <- SeuratObj$type
  SeuratObj$type <- NULL
  SeuratObj$r <- NULL
  
  
  clustifyr_plot <- plot_umap(
    seurat_object = SeuratObj, group_by = clustifyr_colname,
    reduction = paste0("umap", analyses$viz_clustering),
    title = paste0(clustifyr_colname, " Assignment"),
    xlab = "UMAP1", ylab = "UMAP2",
    legend_position = "bottom",
    ncol_guide = 4,
    label_clusters = TRUE,
    label_size = 5,
    color_reverse = FALSE,
    repel_labels = TRUE
  )
  
  pdf(paste0(
    clustifyrDirectory, ObjName, Subset,
    "_res", RESOLUTION, "_", clustifyr_colname, ".pdf"
  ), width = 10, height = 7, family = FONT_FAMILY
  )
  print(clustifyr_plot)
  dev.off()
  
  
  return (SeuratObj)
  
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
  pt_size = 0.6, split_by = NULL, ncol_dimplot = 1, ncol_guide = 1,
  label_clusters = FALSE, repel_labels = FALSE, label_size = 4,
  shuffle = F, order = NULL
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
    label = label_clusters, repel = repel_labels, label.size = label_size,
    shuffle = shuffle, order = order
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
  Idents(seurat_object) <- group_by
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
  assay = "RNA", min_cutoff = NA
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
    label = label, repel = repel, label.size = label_size,
    min.cutoff = min_cutoff
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

heatmap_wrapper <- function (
    SeuratObj, markers, markers_assay
) {
  HM_object <- plot_heatmap (
    seurat_object = SeuratObj, downsample_n = 5000,
    markers = markers, top_n = 20, label_n = 2, 
    cluster = paste0("Cluster", analyses$viz_clustering),
    data_type = "logcounts",
    use_raster = TRUE
  )
  
  pdf(paste0(
    heatDirectory, "heatmap", ObjName, Subset, 
    "_res", RESOLUTION, "_top20 ", markers_assay, " features per ", 
    analyses$viz_clustering, "Cluster.pdf"
  ), width = 7, height = 6
  )
  draw(
    HM_object[[1]], 
    heatmap_legend_list = list(HM_object[[2]], HM_object[[3]]),
    heatmap_legend_side = "right"
  )
  dev.off()
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

############ GEX FIRST SCT ############
GEX_first_SCT <- function(SeuratObjMYSC) {
  DefaultAssay(SeuratObjMYSC) <- "RNA"
  
  # SCT
  if (SCT) {
    SeuratObjMYSC <- RenameAssays(SeuratObjMYSC, RNA = "RNApreSCT")
    SeuratObjMYSC <- SeuratObjMYSC %>%
      SCTransform(
        assay = "RNApreSCT", method = "glmGamPoi", 
        new.assay.name = "RNA",
        vst.flavor = "v2", verbose = TRUE
      )
  } else {
    SeuratObjMYSC <- SeuratObjMYSC %>% 
      NormalizeData(assay = "RNA") %>%
      FindVariableFeatures(assay = "RNA") %>%
      ScaleData(assay = "RNA")
  }
  return (SeuratObjMYSC)
}  
############ SCDBLFINDER & QC ##############
GEX_QC <- function(SeuratObjMYSC, file) {
  # determine doublets, need to change to sce
  SeuratObjMYSC.sce <- as.SingleCellExperiment(SeuratObjMYSC)
  sce <- scDblFinder(SeuratObjMYSC.sce)
  
  # Converting back to Seurat Obj
  SeuratObjMYSC <- as.Seurat(sce)
  Idents(SeuratObjMYSC) <- SeuratObjMYSC$orig.ident
  
  SeuratObjMYSC <- GEX_first_SCT(SeuratObjMYSC)

  
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
  


############ AGGREGATION TIDYING ##############
tidy_metadata <- function(SeuratObj, remove_cols = NA, replacements = NA) {
  for (col_name in remove_cols) {
    if (col_name %in% colnames(SeuratObj@meta.data)) {
      SeuratObj@meta.data[col_name] <- NULL
    }
  }
  
  for (original in names(replacements)) {
    if (original %in% colnames(SeuratObj@meta.data)) {
      replacement <- replacements[[original]]
      SeuratObj@meta.data[replacement] <- 
        SeuratObj@meta.data[original]
      SeuratObj@meta.data[original] <- NULL
    }
  }
  
  return (SeuratObj)
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
      method = "glmGamPoi", verbose = TRUE,
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
GEX_pca <- function(SeuratObjMYSC, file, specific_PCA_features = F) {
  DefaultAssay(SeuratObjMYSC) <- "RNA"
  
  PCA_features <- VariableFeatures(object = SeuratObjMYSC)[
    1:config$pca_variable_features
    ]
  if (specific_PCA_features & config$pca_features != "variable") {
    stopifnot (config$pca_features == "M" | config$pca_features == "T")
    
    # features include VariableFeatures of :
    # 1. SeuratObj
    # 2. Yun lab curated T (or M) cell markers
    # 3. projecTILs reference

    PCA_features <- union(
      VariableFeatures(object = SeuratObjMYSC)[1:config$pca_variable_features],
      unlist(case_sensitive_features(
        SeuratObjMYSC,
        unlist(flatten(
          master[[paste0(config$pca_features, "_lineage_markers")]]
        ))
        ))
    )
    
    # add projectils features
    if (config$pca_features == "T") {
      if (! exists("ref_projectils")) {
        ref_projectils <- load_ref_projectils()
      }
      
      PCA_features <- union(
        PCA_features,
        VariableFeatures(ref_projectils)
      )
    }
    
  }
  
  SeuratObjMYSC <- RunPCA(
    SeuratObjMYSC,
    features = PCA_features
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

############ CELL CYCLE REGRESSION #########
GEX_cc_regression <- function(SeuratObj) {
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
      vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE,
      return.only.var.genes = FALSE
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
    file = paste0(    
      RobjDirectory, ObjName, Subset, 
      "_res", config$RESOLUTION, ".rds"
    )
  )
  return (SeuratObj)
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
  
  E2 <- ElbowPlot(SeuratObj, ndims = 50, reduction = "pcaADT")
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
################ GET CLUSTER NAME ###############
  get_cluster_name <- function() {
    # i.e. cluster_name = "ClusterRNA"
    if (analyses$viz_clustering == "RNA") {
      cluster_name = "ClusterRNA"
    } else if (analyses$viz_clustering == "ADT") {
      cluster_name = "ClusterADT"
    } else if (analyses$viz_clustering == "WNN") {
      cluster_name = "ClusterWNN"
    } else {
      stop("invalid cluster_name")
    }
    return (cluster_name)
  }
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
{
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
  stopifnot(config$preprocess_existing_atlas)
  
  # to contain preprocessed samples
  SeuratSamples <- list()
  # load atlas
  seurat_object <- readRDS(paste0(
    RobjDir, config$existing_atlas_name
  ))
  # want to plot atlas Patient, not Sample
  # so switch them out. Also, orig.ident can be changed to Patient
  seurat_object$orig.ident <- seurat_object$Patient
  seurat_object$Patient <- seurat_object$Sample
  seurat_object$Sample <- seurat_object$orig.ident
  
  # split by patient
  seurat_list <- SplitObject(
    seurat_object, split.by = config$existing_atlas_split_by
    )
  
  for (idx in 1:length(seurat_list)) {
    # get sample
    seurat_split <- seurat_list[idx][[1]]
    sample_name <- unique(
      seurat_split@meta.data[[config$existing_atlas_split_by]]
      )[1]
    print(paste0(config$existing_atlas_split_by, " ", sample_name))
    # make sure there are enough cells and it's not a sample to be removed
    if (
      ncol(seurat_split) > 55 & 
      sample_name %in% config$existing_atlas_use_cluster
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
      print ("not enough cells or bad sample, moving onto next sample")
    }
  }
  
  return (SeuratSamples)
}

get_GEO_unzip <- function() {
  # get GEO files
  getGEOSuppFiles(analyses$GEO)
  files_to_untar <- list.files(GEOdir, pattern = ".gz")
  for (file in files_to_untar) {
    gunzip(paste0(GEOdir, file))
  }
}

################ CUSTOM TIDY ATLAS FOR MERGE ###################
{
  determine_orig_ident <- function(idx, cellnames, libraries, IDfrags) {
    library <- libraries[[idx]]
    cell_name <- cellnames[[idx]]
    IDfrag <- IDfrags[[idx]]
    return_value <- str_split(cell_name, "_")[[1]][1]
    if (! grepl("MYSC", return_value, fixed = F)) {
      if (! is.na(library)) {
        return_value <- library
      } else {
        return_value <- IDfrag
      }
      
    }
    return (return_value)
  }
  
  
  set_orig_ident <- function(refSeuratObj) {
    orig.idents <- lapply(
      seq_along(orig_library), 
      determine_orig_ident, 
      rownames(refSeuratObj@meta.data),
      refSeuratObj$Library,
      refSeuratObj$IDfrag
    )
    refSeuratObj$orig.ident <- unlist(orig.idents)
    refSeuratObj$Library <- refSeuratObj$orig.ident
    return (refSeuratObj)
  }
  
  remove_unnecessary_metadata <- function(
    refSeuratObj, keep_cols = NA, replacements = NA
  ) {
    for (col_name in colnames(refSeuratObj@meta.data)) {
      if (! col_name %in% keep_cols) {
        refSeuratObj@meta.data[col_name] <- NULL
      }
    }
    
    for (original in names(replacements)) {
      if (original %in% colnames(refSeuratObj@meta.data)) {
        replacement <- replacements[[original]]
        refSeuratObj@meta.data[replacement] <- 
          refSeuratObj@meta.data[original]
        refSeuratObj@meta.data[original] <- NULL
      }
    }
    
    # call refSeuratObj <- tidy_metadata(refSeuratObj, config$)
    return (refSeuratObj)
  }
  
  master_tidy_atlas_metadata <- function(refSeuratObj) {
    # remove annoying cols
    refSeuratObj <- remove_unnecessary_metadata(
      refSeuratObj, keep_cols = master$metadata_merger_keep
    )
    # set orig.ident and Library cols
    refSeuratObj <- set_orig_ident(refSeuratObj)
    
    # # remove sample ndGBM-05 (repeat of MYSC70)
    # refSeuratObj <- subset(
    #   refSeuratObj, Patient != "ndGBM-05"
    #   )
    return (refSeuratObj)
  }
  
  
}

############### CUSTOM TIDY MATCHING PATIENTS FOR MERGE #################
{
  
}

}


#################### DOWNSTREAM  #################
{
#################### SUBSET T CELLS ONLY #################
subset_Tcells <- function(seurat_object) {
  seurat_object <- subset(
    x = seurat_object,
    subset = 
      # (
      #   SingleR_BlueprintEncode == "CD8+ T-cells" |
      #   SingleR_BlueprintEncode == "CD4+ T-cells"
      # ) &
      # (
      #   SingleR_HPCA == "T_cells"
      # ) &
      (
        SingleR_Novershtern == "NK T cells" |
        SingleR_Novershtern == "CD4+ T cells" |
        SingleR_Novershtern == "CD8+ T cells"
      ) |
      (
        projecTILs == "CD8_EffectorMemory" |
        projecTILs == "Th1" |
        projecTILs == "CD8_Tex" |
        projecTILs == "Treg" |
        projecTILs == "CD8_Tpex" |
        projecTILs == "Tfh" |
        projecTILs == "CD8_EarlyActiv" |
        projecTILs == "CD8_NaiveLike" |
        projecTILs == "CD4_NaiveLike"
      )
  )
  
  return (seurat_object)
}
#################### WRITE ATLAS NORMALIZATION COUNTS ##################
preprocess_atlas_objects_corr <- function(
    RDS_filename, population, flip_patient_sample = T
    ) {
  seurat_all <- readRDS(file = RDS_filename)
  seurat_all <- subset(
    seurat_all, 
    subset = Patient != "LGG-04" &
      Patient != "LGG-03"
  )
  
  if (flip_patient_sample) {
    seurat_all$orig.ident <- seurat_all$Patient
    seurat_all$Patient <- seurat_all$Sample
    seurat_all$Sample <- seurat_all$orig.ident
  }


  seurat_all$TumorType <- ifelse(
    grepl("ndGBM", seurat_all$Sample, fixed = TRUE),
    "ndGBM",
    "rGBM"
  )
  
  write.table(
    table(seurat_all$Fragment),
    paste0(
      dataDirectory, "GBMAtlas/", 
      population, "_normalization_counts.csv"
    ),
    col.names = c("Fragment", "Freq"),
    sep = ",",
    row.names = F
  )
  return (seurat_all)
}
################## READ ATLAS NORMALIZATION COUNTS #####################
read_normalization_counts <- function(csv_name) {
  counts <- read.csv(
    paste0(dataDirectory, "GBMAtlas/", csv_name)
  )
  normalization_counts <- list()
  for (idx in 1:length(counts[[1]])) {
    normalization_counts[[counts$Fragment[idx]]] <- counts$Freq[idx]
  }
  return (normalization_counts)
}

#################### SUBTYPE COUNTS PER FRAGMENT ######################
write_subtype_counts_per_fragment <- function (
    seurat_object, fragment_col, subtype_col, 
    patient_population, cell_population
) {
  write.table(
    table(
      seurat_object@meta.data[[fragment_col]], 
      seurat_object@meta.data[[subtype_col]]
    ),
    paste0(
      dataDirectory, "GBMAtlas/", 
      patient_population, "_", cell_population, "_", 
      fragment_col, "_", subtype_col, "_populations.csv"
    ),
    sep = ",",
    row.names = T,
    col.names = NA
  )
}

#################### SLINGSHOT COLOR COLUMN #######################
# modify SeuratObj for color scheme
private_category_to_int <- function(category, translation) {
  return (translation[[category]])
}
meta_category_to_int <- function(SeuratObj, column, translation_column) {
  translation <- list()
  for (idx_category in 1:length(unique(SeuratObj@meta.data[[column]]))) {
    category <- sort(unique(SeuratObj@meta.data[[column]]))[idx_category]
    translation[[category]] <- idx_category
  }
  SeuratObj <- AddMetaData(
    object = SeuratObj, 
    metadata = map_int(
      SeuratObj@meta.data[[column]], private_category_to_int, translation
      ),
    col.name = translation_column
  )
  
  return (SeuratObj)
}

#################### SLINGSHOT PLOT UMAP #######################
plot_slingshot_umap <- function (SeuratObj.sce, Colors, ColNo, title) {
  plot (
    reducedDims(SeuratObj.sce)$UMAPRNA, 
    col = Colors[SeuratObj@meta.data[[ColNo]]], 
    pch = 16, asp = 1, cex = 0.5, main = title
  )
  lines(SlingshotDataSet(SeuratObj.sce), lwd = 2, col = "black")
}

#################### SLINGSHOT PLOT PSEUDOTIME #######################
plot_slingshot_pseudotime <- function (
    SeuratObj,
    pseudotime_version,
    cluster_by,
    color_scheme,
    x_lab = "Slingshot pseudotime",
    y_lab = "Subtype"
) {
  G1 <- ggplot(
    as.data.frame(SeuratObj@meta.data), 
    aes_string(
      x = SeuratObj@meta.data[[pseudotime_version]], 
      y = cluster_by, 
      colour = cluster_by
    )
  ) + 
    ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
    theme_classic() +  
    xlab(x_lab) + 
    ylab(y_lab) + 
    ggmin::theme_min() + 
    scale_color_manual(values = color_scheme)
  return (G1)
}

}

