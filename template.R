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
# # gamma poisson generalized linear model
# BiocManager::install("glmGamPoi")

# remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
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
  "glmGamPoi", "gmailr", "rjson", "here"
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



################# CONFIG DIRECTORIES #################
# Set working directory and create input / output folders
Directory <- config$Directory
#dir.create(Directory)
setwd(Directory)
RobjDirectory <- paste0(Directory, "/R_Objects/")
# dir.create(RobjDirectory)
RcodeDirectory <- paste0(Directory, "/R_Code/")
# dir.create(RcodeDirectory)
dataDirectory <- paste0(Directory, "/Data/")
# dir.create(dataDirectory)

outdir <- paste0(Directory, "/Output/")
# dir.create(outdir)
OutputDirectory <- paste0(Directory, "/Output/AllClusters/")
# dir.create(OutputDirectory)

QCDirectory <- paste0(OutputDirectory, "/QualityControl/")
# dir.create(QCDirectory)
ElbowDirectory <- paste0(OutputDirectory, "/ElbowPlots/")
# dir.create(ElbowDirectory)
featuremapDirectory <- paste0(OutputDirectory, "/FeatureMaps/")
# dir.create(featuremapDirectory)
densityplotDirectory <- paste0(OutputDirectory, "/DensityPlots/")
# dir.create(densityplotDirectory)
UMAPDirectory <- paste0(OutputDirectory, "/UMAPs/")
# dir.create(UMAPDirectory)
ConfigDirectory <- paste0(OutputDirectory, "/Configs/")
# dir.create(ConfigDirectory)
QCDirectory <- paste0(OutputDirectory, "/QualityControl/")
# dir.create(QCDirectory)

# NOTE: tempdir() directory is destroyed upon closing R session
errorDir <- paste0(Directory, "/temp/")
# dir.create(errorDir)

################# CONFIG ERROR LOG & EMAIL NOTIFS #################
PERSON_TO = config$PERSON_TO
EMAIL_TO = config$EMAIL_TO
CATCH_ERRORS = config$CATCH_ERRORS
SEND_WARNINGS = config$SEND_WARNINGS

if (SEND_WARNINGS) {
  warn_number <- 2
} else {warn_number <- 0}

# create output, input, and traceback files
# initialize and populate files upon startup
traceback_file <- paste0(errorDir, "traceback.txt")
output_file <- paste0(errorDir, "output.txt")
sink(file = output_file, append = FALSE, type = "output", split = TRUE)
print("STARTING PIPELINE")
write(
  capture.output(traceback(2)), traceback_file,
  append = FALSE
)

# setup emailing service
gm_auth_configure(
  path = paste0(
    "~/Box/Yun lab projects/victoria_liu/Notifications/Rbotcredentials.json"
  )
)
# username: friendlyrbot@gmail.com
# password: 1lovemedicine@yay

send_bot_email <- function(
  traceback_file, output_file,
  person_to = PERSON_TO,
  email_to = EMAIL_TO,
  signature = paste0(
    "<br><br>Cheers,",
    "<br><strong>R notification bot</strong>",
    "<br>maintained by Victoria Liu (victoria.liu@jax.org)"
  ),
  subject = "R pipeline suspended",
  message_body = paste0(
    "Dear ", person_to, ",", "<br>",
    "Your pipeline has been aborted due to a warning or failure. ",
    "Please refer to the attached trace and output history.<br>"
  )
) {
  email <- gm_mime() %>%
    gm_to(paste0(person_to, " <", email_to, ">")) %>%
    gm_from("Friendly R Bot <friendlyrbot@gmail.com>") %>%
    gm_subject(subject) %>%
    gm_html_body(
      paste0(
        "<p><font face = 'courier' size = '+1'>",
        message_body,
        signature,
        "</font></p>"
      )
    ) %>%
    gm_attach_file(traceback_file) %>%
    gm_attach_file(output_file)

  gm_send_message(email)
}

options(
  warn = warn_number,
  # error = NULL,
  error = function () {
    write(
      capture.output(traceback(2)), traceback_file,
      append = FALSE
    )

    if (CATCH_ERRORS) {
      send_bot_email(traceback_file, output_file)
      system(
        "say -v Samantha Pipeline aborted, sent a warning email. Going to sleep"
      )
      Sys.sleep(Inf)
    } else {
      system("say -v Samantha Error detected, continuing anyway.")
    }
  }
)
options(error = NULL)

# send first email to authorize
send_bot_email(
  traceback_file, output_file,
  subject = "pipeline sarting",
  message_body = "EOM"
)

################ SET PARAMETERS ###############
DOUBLET_FORMATION_RATE <- config$DOUBLET_FORMATION_RATE
# resolution for cluster-finding after PCA
RESOLUTION <- config$RESOLUTION
# if we are doing scTransform
SCT <- config$SCT
# font family for output plots
FONT_FAMILY <- config$FONT_FAMILY
# files to analyze
FILES = config$FILES


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
    xintercept,
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
      geom_vline(xintercept = xintercept) + 
      scale_fill_manual(values = manual_colors) +
      scale_color_manual(values = manual_colors)
    
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
    legend_position,
    color_reverse = FALSE, reduction = "umap",
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
    seurat_object, feature_gene, split_by,
    order = TRUE, label = TRUE, label_size = 2,
    cols = Nour_cols(c("darkpurple", "lightorange")), pt_size = 0.1,
    legend_title_position = "top", legend_title_size = 10,
    legend_location = "right",
    ncol = 3, nrow = 1, widths = c(0.03, 1, 1)
  ) {
    F1 <- FeaturePlot(
      seurat_object, features = feature_gene, split.by = split_by,
      order = order, label = label, label.size = label_size,
      cols = cols, pt.size = pt_size
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
    data_type = "logcounts",
    color_map = c("#007dd1", "white", "#ab3000"),
    use_raster = TRUE
  ) {
    top_n_markers <- get_top_cluster_markers(markers, top_n)
    label_n_markers <- get_top_cluster_markers(markers, label_n)
    
    subset <- subset(seurat_object, downsample = downsample_n)
    # seurat_sce <- as.SingleCellExperiment(subset)
    
    # plot data rows are genes, cols are cells
    plot_data <- as.data.frame(SeuratObj@assays$RNA@scale.data)
    # plot_data <- as.data.frame(assay(seurat_sce, data_type))
    plot_data <- plot_data[top_n_markers$gene, ]
    plot_data <- na.omit(plot_data)
    plot_data <- plot_data - rowMeans(plot_data)
    
    # column_annot rows are cells, cols are "Cluster", "Sample"
    col_anno_df <- subset@meta.data[, c("Cluster", "Sample"), drop = F] 
    col_anno_df$Sample = as.factor(col_anno_df$Sample)
    # heatmap is ordered in cluster (primary) and then sample (secondary)
    col_anno_df <- with(col_anno_df, col_anno_df[order(Sample), , drop = F])
    col_anno_df <- with(col_anno_df, col_anno_df[order(Cluster), , drop = F])
    # order cells by cluster (primary) and sample (secondary)
    plot_data <- plot_data[, row.names(col_anno_df)]
    
    # sample and cluster colors are reverse of each other
    sample_colors <- get_colors(
      seurat_object = SeuratObj,
      color_by = "Sample"
    )
    cluster_colors <- get_colors(
      seurat_object = SeuratObj,
      color_by = "Cluster",
      color_reverse = TRUE    
    )
    column_colors = list()
    column_colors[["Sample"]] <- sample_colors
    column_colors[["Cluster"]] <- cluster_colors
    
    Sample = as.matrix(col_anno_df[, c("Sample"), drop = F])
    Cluster = as.matrix(col_anno_df[, c("Cluster"), drop = F])
    
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
      labels = levels(as.factor(col_anno_df$Cluster)),
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

######### QC + Doublet removal ########
SeuratSamples <- list()



########### remove start #########
for(i in 1:length(FILES)){
  print(FILES[i])
  # read in GEX and ADT
  MYSC <- Read10X(
    data.dir = paste0(
      Directory,
      "/Data/",
      FILES[i], 
      "/filtered_feature_bc_matrix/"
    )
  )
  
  # create Seurat obj
  # GEX
  SeuratObjMYSC <- CreateSeuratObject(
    counts = MYSC[["Gene Expression"]], 
    min.features = 100, project = FILES[i]
  )
  # ADT
  SeuratObjMYSC[["ADT"]] <- CreateAssayObject(
    MYSC[["Antibody Capture"]][, colnames(x = SeuratObjMYSC)])
  
  SeuratSamples[[i]] <- SeuratObjMYSC

}
############ remove end ###########


for(i in 1:length(FILES)){
  print(FILES[i])
  # read in GEX and ADT
  MYSC <- Read10X(
    data.dir = paste0(
      Directory,
      "/Data/",
      FILES[i], 
      "/filtered_feature_bc_matrix/"
    )
  )
  
  # create Seurat obj
  # GEX
  SeuratObjMYSC <- CreateSeuratObject(
    counts = MYSC[["Gene Expression"]], 
    min.features = 100, project = FILES[i]
  )
  # ADT
  SeuratObjMYSC[["ADT"]] <- CreateAssayObject(
    MYSC[["Antibody Capture"]][, colnames(x = SeuratObjMYSC)])
  
  # QC for mitochondrial RNA
  SeuratObjMYSC$log10GenesPerUMI <- log10(SeuratObjMYSC$nFeature_RNA) / 
    log10(SeuratObjMYSC$nCount_RNA)
  
  SeuratObjMYSC$mito_ratio <- PercentageFeatureSet(
    SeuratObjMYSC, pattern = "^mt-"
  )
  SeuratObjMYSC$mito_ratio <- SeuratObjMYSC@meta.data$mito_ratio / 100
  
  #filter by removing cells that didn't pass qc
  SeuratObjMYSC <- subset(
    SeuratObjMYSC, 
    subset = nFeature_RNA > config$nFeature_RNA & 
      mito_ratio < config$mito_ratio & 
      nFeature_ADT > config$nFeature_ADT
  )
  
  #prep work for DoubletFinder
  if (SCT) {
    SeuratObjMYSC <- SeuratObjMYSC %>%
      SCTransform(
        assay = "RNA", method = "glmGamPoi", 
        vst.flavor = "v2", verbose = FALSE
      )
  } else {
    SeuratObjMYSC <- SeuratObjMYSC %>% 
      NormalizeData(assay = "RNA") %>%
      FindVariableFeatures(assay = "RNA") %>%
      ScaleData(assay = "RNA")
  }
  
  
  #Run PCA
  SeuratObjMYSC <- RunPCA(
    SeuratObjMYSC, features = VariableFeatures(object = SeuratObjMYSC)
  )
  print(SeuratObjMYSC[["pca"]], dims = 1:5, nfeatures = 5)
  # most representative (by absolute value) genes for PCs
  VS = VizDimLoadings(SeuratObjMYSC, dims = 1:2, reduction = "pca")
  VSD = DimPlot(SeuratObjMYSC, reduction = "pca") + theme_min()
  
  # save plot
  pdf(paste0(
    QCDirectory, FILES[i], " GEX PCA features and plot.pdf"
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
  
  # remove doublets
  # documentation: https://github.com/chris-mcginnis-ucsf/DoubletFinder
  # pK Identification (no ground-truth)
  # if using sctransform, mark sct as TRUE
  sweep.res.list <- paramSweep_v3(SeuratObjMYSC, PCs = 1:10, sct = SCT)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  # to plot max bcmvn for identifying pK (refer to graph on wiki page)
  df = as.data.frame(bcmvn)
  pK_graph = df %>%
    ggplot(aes( x = ParamID, y = BCmetric)) +
    geom_point(color = "royalblue") + 
    geom_line(color = "royalblue") +
    theme_classic() + 
    # plot pK on x-axis
    geom_vline(xintercept = df$ParamID[df$BCmetric == max(df$BCmetric)])
  
  pdf(paste0(
    QCDirectory, FILES[i],
    " doubletfinder BCmvn distribution.pdf"),
    width = 8, height = 5.5, family = FONT_FAMILY
  )
  print(pK_graph)
  dev.off()
  
  # Homotypic Doublet Proportion Estimate
  # I think $Clusters may be empty / just a placeholder? 
  homotypic.prop <- modelHomotypic(SeuratObjMYSC@meta.data$Clusters)    
  # possion distribution
  nExp_poi <- round(DOUBLET_FORMATION_RATE * nrow(SeuratObjMYSC@meta.data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies
  # use 10 principle components, doublet count as 0.25 (not too important)
  # pK and pANN (neighborhood size) calculated using parameter sweep
  SeuratObjMYSC <- doubletFinder_v3(
    SeuratObjMYSC, PCs = 1:10, 
    pN = 0.25, pK = as.numeric(df$pK[df$BCmetric == max(df$BCmetric)]),
    nExp = nExp_poi, reuse.pANN = FALSE, sct = SCT
  )
  
  SeuratObjMYSC <- doubletFinder_v3(
    SeuratObjMYSC, PCs = 1:10, 
    pN = 0.25, pK = as.numeric(df$pK[df$BCmetric == max(df$BCmetric)]),
    nExp = nExp_poi.adj, reuse.pANN = colnames(SeuratObjMYSC@meta.data)[5], 
    sct = SCT
  )
  
  # rename columns
  colnames(SeuratObjMYSC@meta.data)[
    length(SeuratObj@meta.data) - 1
  ] <- "pANN"
  colnames(SeuratObjMYSC@meta.data)[
    length(SeuratObj@meta.data)
  ] <- "DoubletStatus"
  
  # save plot
  pdf(paste0(
    QCDirectory, FILES[i],
    "Doublet status dimplot.pdf"
  ), width = 8, height = 5.5, family = FONT_FAMILY
  )
  DimPlot(
    SeuratObjMYSC, split.by = "DoubletStatus", 
    order = TRUE, shuffle = TRUE)
  dev.off()
  
  # save object
  saveRDS(SeuratObjMYSC, file = paste0(RobjDirectory, FILES[i], ".rds"))
  
  SeuratSamples[[i]] <- SeuratObjMYSC
  # save space
  rm(SeuratObjMYSC)
}

######## AGGREGATION ########

# names for saving plots
ObjName <- config$ObjName
Subset <- config$Subset

# merge samples into one object
SeuratObj <- merge(
  x = SeuratSamples[[1]], y = SeuratSamples[[-1]], 
  add.cell.ids = config$FILES
)


# plot doublet status
PDS <- plot_bargraph (
  seurat_object = SeuratObj, aesX = "DoubletStatus", fill = "orig.ident",
  y_label = "Composition (Number of cells)", x_label = NULL,
  y_lower_limit = 0, y_break = 2000
)
# save plot
pdf(paste0(
  QCDirectory, ObjName, Subset,
  "Barplot number of Doublets per sample calc by doubletfinder",
  "at sevenPointFive percent freq.pdf"
), width = 8, height = 5.5, family = FONT_FAMILY)
PDS
dev.off()

# default identity is double status (for now)
Idents(SeuratObj) <- SeuratObj$DoubletStatus

# remove doublets from Seurat data structure (including for ADT)
SeuratObj <- subset(
  SeuratObj, cells = WhichCells(SeuratObj, idents = "Singlet")
)


# extract current metadata
CellInfo <- SeuratObj@meta.data

# add sample information
sample_info <- config$sample_info

CellInfo$Sample <- NA
for (i in 1:length(sample_info)) {
  CellInfo$Sample[
    which(str_detect(CellInfo$orig.ident, names(sample_info[i])))
  ] <- sample_info[[i]]
}

Sample <- CellInfo$Sample
names(Sample) <- colnames(x = SeuratObj)
SeuratObj <- AddMetaData(
  object = SeuratObj, metadata = Sample, col.name = "Sample"
)

# add patient information (anonymized, deidentified)
patient_info <- config$patient_info

CellInfo$Patient <- NA
for (i in 1:length(patient_info)) {
  CellInfo$Patient[
    which(str_detect(CellInfo$orig.ident, names(patient_info[i])))
  ] <- patient_info[[i]]
}

Patient <- CellInfo$Patient
names(Patient) <- colnames(x = SeuratObj)
SeuratObj <- AddMetaData(
  object = SeuratObj, metadata = Patient, col.name = "Patient"
)

#add group information
group_info <- config$group_info

CellInfo$Group <- NA
for (i in 1:length(group_info)) {
  CellInfo$Group[
    which(str_detect(CellInfo$orig.ident, names(group_info[i])))
  ] <- group_info[[i]]
}

Group <- CellInfo$Group
names(Group) <- colnames(x = SeuratObj)
SeuratObj <- AddMetaData(
  object = SeuratObj, metadata = Group, col.name = "Group"
)


# plot
P1 <- plot_bargraph (
  seurat_object = SeuratObj, aesX = "Sample", fill = "Sample",
  y_label = "Composition (Number of cells)", x_label = NULL, 
  y_lower_limit = 0, y_break = 1000
)

pdf(paste0(
  densityplotDirectory, ObjName, Subset, 
  " Barplot number of cells per sample.pdf"), 
  width = 8, height = 5.5, family = FONT_FAMILY
)
P1
dev.off()


# plot
x1 <- plot_densitygraph (
  seurat_object = SeuratObj, aesX = "nCount_RNA", fill = "Sample",
  color_by = "Sample", xintercept = 1000, scale_x_log10 = TRUE
)

x2 <- plot_densitygraph (
  seurat_object = SeuratObj, aesX = "nFeature_RNA", fill = "Sample",
  color_by = "Sample", xintercept = 500, scale_x_log10 = TRUE
)

x3 <- plot_densitygraph (
  seurat_object = SeuratObj, aesX = "mito_ratio", fill = "Sample",
  color_by = "Sample", xintercept = 0.015
)

x4 <- plot_densitygraph (
  seurat_object = SeuratObj, aesX = "log10GenesPerUMI", fill = "Sample",
  color_by = "Sample", xintercept = 0.8, scale_x_log10 = TRUE
)

x5 <- plot_densitygraph (
  seurat_object = SeuratObj, aesX = "nCount_ADT", fill = "Sample",
  color_by = "Sample", xintercept = 0.8, scale_x_log10 = TRUE
)

x6 <- plot_densitygraph (
  seurat_object = SeuratObj, aesX = "nFeature_ADT", fill = "Sample",
  color_by = "Sample", xintercept = 0.8, scale_x_log10 = TRUE
)

XX <- grid_arrange_shared_legend(
  x1 + theme_min(), 
  x2 + theme_min(), 
  x3 + theme_min(), 
  x4 + theme_min(), 
  position = "right", 
  ncol = 2, nrow = 3
)

#violin plots
X <- VlnPlot(
  SeuratObj, 
  features = c("nFeature_RNA", "nCount_RNA", "mito_ratio"), 
  ncol = 3, pt.size = 0, 
  cols = get_colors(seurat_object = SeuratObj, color_by = "Sample"), 
  group.by = "Sample"
) +
  RotatedAxis()
x5 = ggarrange(
  X[[1]] + theme_min() + NoLegend() + RotatedAxis(), 
  X[[2]] + theme_min() + NoLegend() + RotatedAxis(), 
  X[[3]] + theme_min() + NoLegend() + RotatedAxis(), 
  nrow = 1
)

pdf(paste0(
  QCDirectory, ObjName, Subset, 
  " qualitycontrol.pdf"
), width = 12, height = 10, family = FONT_FAMILY
)
grid.arrange(XX, x5, heights = c(0.8, 0.5))
dev.off()



######## NORMALIZATION ((GEX + ADT) ########
DefaultAssay(SeuratObj) <- "RNA"

# SCT for GEX
if (SCT) {
  # rename for ease of use; no matter if SCT or NormalizeData,
  # we would like to use "RNA" assay from now on
  SeuratObj <- RenameAssays(SeuratObj, RNA = "RNApreSCT")
  SeuratObj <- SCTransform(
    SeuratObj,
    assay = "RNApreSCT",
    new.assay.name = "RNA",
    method = "glmGamPoi", verbose = FALSE
  )

} else {
  SeuratObj <- SeuratObj %>% 
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData()
}

# Normalize for ADT 
# not a lot of features, so use all antibody rows
# scaling happens after CC regression and integration
DefaultAssay(SeuratObj) <- "ADT"
SeuratObj <- NormalizeData(
  SeuratObj, normalization.method = "CLR", margin = 2
)
VariableFeatures(SeuratObj) <- rownames(SeuratObj[["ADT"]])
# return to GEX assay for plotting
DefaultAssay(SeuratObj) <- "RNA"

# plot variable features
top10 <- head(VariableFeatures(SeuratObj), 10)
plot1 <- VariableFeaturePlot(SeuratObj)
plot2 <- LabelPoints(plot = plot1, points = top10, size = 3)
pdf(paste0(
  OutputDirectory, ObjName, Subset, 
  " GEX variable freatures.pdf"
), width = 8, height = 5.5, family = FONT_FAMILY
)
plot2 + theme_min()
dev.off()

if (config$run_Harmony_CC_controls & FALSE) {
  ############ FOR DEMO PURPOSES START #############
  #Run PCA
  
  SeuratObj <- RunPCA(
    SeuratObj, features = VariableFeatures(object = SeuratObj)
  )
  print(SeuratObj[["pca"]], dims = 1:5, nfeatures = 5)
  
  # plot important genes of first two PCs
  # plot dim map of first two PCs, grouped by sample
  # ideally, shouldn't see much difference in clustering amongst samples
  # batch effects can be removed with Harmony
  VS1 <- VizDimLoadings(SeuratObj, dims = 1:2, reduction = "pca")
  VSD1 <- DimPlot(
    SeuratObj, reduction = "pca", group.by = "Sample", 
    cols = get_colors(seurat_object = SeuratObj, color_by = "Sample")
  ) + theme_min()
  
  pdf(paste0(
    OutputDirectory, ObjName, Subset, 
    " demo PCA features and plot (no harmony or CC regression).pdf"
  ), width = 12, height = 4.5, family = FONT_FAMILY
  )
  grid.arrange(
    VS1[[1]] + theme_min(), 
    VS1[[2]] + theme_min(), 
    VSD1, 
    nrow = 1, 
    widths = c(0.6, 0.6, 1)
  )
  dev.off()
  
  # this takes a long time
  #SeuratObj <- JackStraw(SeuratObj, num.replicate = 100)
  #SeuratObj <- ScoreJackStraw(SeuratObj, dims = 1:20)
  #JackStrawPlot(SeuratObj, dims = 1:15)
  
  # elbow plot is quick and dirty, to get number of dims
  elbowplot <- ElbowPlot(SeuratObj, ndims = 50)
  
  # UMAP: based on GEX PCA only
  dims <- 1:20
  SeuratObj <- SeuratObj %>%
    # get "euclidian" distance based on PC feature vectors
    FindNeighbors(dims = dims) %>%
    # find clusters based on resolution and neighbor graph
    FindClusters(resolution = RESOLUTION) %>%
    #initial clustering
    RunUMAP(dims = dims)
  # save reduction
  SeuratObj@reductions$GEX_UMAP_noharmony_noCC <- SeuratObj@reductions$umap
  
  # harmony batch correction
  SeuratObj <- SeuratObj %>% 
    RunHarmony(
      group.by.vars = config$batch_norm_by,
      reduction.save = "harmonyRNA"
    ) %>%
    FindNeighbors(reduction = "harmonyRNA", dims = dims) %>%
    FindClusters(resolution = RESOLUTION) %>%
    RunUMAP(reduction = "harmonyRNA", dims = dims)
  # save info
  SeuratObj@reductions$GET_UMAP_harmony_noCC <- SeuratObj@reductions$umap
  ############ FOR DEMO PURPOSES END #############
}

########### CELL CYCLE REGRESSION (GEX + ADT) #############
DefaultAssay(SeruatObj) <- "RNA"
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
    vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE
  )
} else {
  SeuratObj <- ScaleData(
    SeuratObj,
    vars.to.regress = c("S.Score", "G2M.Score"), 
    features = rownames(SeuratObj)
  )
}

# ADT
DefaultAssay(SeuratObj) <- "ADT"
SeuratObj <- ScaleData(
  SeuratObj,
  vars.to.regress = c("S.Score", "G2M.Score"),
  features = rownames(SeuratObj)
)
DefaultAssay(SeuratObj) <- "RNA"

# regressing is a really long process, so save for future
saveRDS(SeuratObj, file = paste0(RobjDirectory, "AllClusters.rds"))

########### PCA & HARMONY BATCH CORRECTION (GEX) ##########
DefaultAssay(SeruatObj) <- "RNA"

SeuratObj <- SeuratObj %>% 
  RunPCA(
    features = VariableFeatures(object = SeuratObj)
  ) %>% 
  RunHarmony(
    group.by.vars = config$batch_norm_by,
    reduction.save = "harmonyRNA"
  )
E1 <- ElbowPlot(SeuratObj, ndims = 15, reduction = "harmonyRNA")

pdf(paste0(
  ElbowDirectory, ObjName, Subset, 
  " elbow plot after CC regression and harmony (GEX).pdf"
), width = 12, height = 4.5, family = FONT_FAMILY
)
E1
dev.off()


########### BATCH INTEGRATION (ADT) ###########
DefaultAssay(SeuratObj) <- "ADT"
# normalize per sample
# using CLR normalization based on: 
# https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html
samples_list <- SplitObject(SeuratObj, split.by = config$batch_norm_by) %>%
  lapply(
    FUN = function(x) {
      x <- NormalizeData(
        x, normalization.method = "CLR", margin = 2, assay = "ADT"
      )
    }
  )
# select features that are repeatedly variable across data sets
features <- SelectIntegrationFeatures(
  object.list = samples_list
)

# find anchor cells (this step stakes a while)
cell_anchors <- FindIntegrationAnchors(
  object.list = samples_list, anchor.features = features
)
# combine into integrated Seurat object
SeuratObj_ADT <- IntegrateData(
  anchorset = cell_anchors, 
  normalization.method = "LogNormalize",
  new.assay.name = "integratedADT"
)

# scale data & dimreduc PCA, in prep for WNN
DefaultAssay(SeuratObj_ADT) <- "integratedADT"
SeuratObj_ADT <- SeuratObj_ADT %>%
  ScaleData() %>%
  RunPCA(reduction.name = "pcaADT")

# add ADT data back into original SeuratObj
# necessary bc IntegrateData destroys scale.data for all assays
# ignore any warnings about offending keys
# ref: https://github.com/satijalab/seurat/issues/3843
SeuratObj[["integratedADT"]] <- SeuratObj_ADT[["integratedADT"]]
SeuratObj[["pcaADT"]] <- SeuratObj_ADT[["pcaADT"]]

E2 <- ElbowPlot(SeuratObj, ndims = 15, reduction = "pcaADT")
pdf(paste0(
  ElbowDirectory, ObjName, Subset, 
  " elbow plot after CC scaling integration and pca (ADT).pdf"
), width = 12, height = 4.5, family = FONT_FAMILY
)
E2
dev.off()

DefaultAssay(SeuratObj) <- "RNA"

######## LOUVAIN CLUSTERING ########
# elbow plot
# get downstream experiment parameters
analyses <- fromJSON(file = here("analysis.json"))
if (config$clustering == "WNN") {
  SeuratObj <- SeuratObj %>%
    FindMultiModalNeighbors(
      reduction.list = list("harmonyRNA", "pcaADT"), 
      dims.list = list(
        1 : analyses$harmonyRNA_dims, 1 : analyses$pcaADT_dims
      ), 
      modality.weight.name = "RNA.weight"
  )
} else {
  DefaultAssay(SeuratObj) <- "RNA"
  SeuratObj <- SeuratObj %>%
    FindNeighbors(reduction = "harmonyRNA", dims = c(1:8)) %>%
    FindClusters(resolution = RESOLUTION)
}

######## RUN UMAP ########
SeuratObj <- RunUMAP(SeuratObj, reduction = "harmonyRNA", dims = c(1:8))

# batch effects should be minimal at this point
DimPlot(SeuratObj, split.by = "Sample")
FeaturePlot(SeuratObj, features = "Ptprc", order = T)


CellInfo <- SeuratObj@meta.data
# Rename Clusters
clusters <- levels(as.factor(SeuratObj@meta.data$seurat_clusters))
for(j in 1 : length(clusters)){
  if (j < 10){
    CellInfo$Cluster[CellInfo$seurat_clusters == j - 1] <- paste0("C0", j)
  }
  else {
    CellInfo$Cluster[CellInfo$seurat_clusters == j - 1] <- paste0("C", j)
  }
}
SeuratObj@meta.data <- CellInfo
Idents(SeuratObj) <- CellInfo$Cluster


DimPlot(SeuratObj, group.by = "Cluster")

# Get number of cells per cluster and per Sample
write.csv(
  as.matrix(table(SeuratObj@meta.data$Cluster, SeuratObj@meta.data$Sample)),
  file = paste0(
    OutputDirectory, ObjName, Subset, 
    " number of cells per cluster and sample.csv"
  )
)

################ FIND CLUSTER BIOMARKERS ################

markers <- FindAllMarkers(
  SeuratObj,
  logfc.threshold = 0.5, test.use = "wilcox", only.pos = TRUE
)

# save, because it takes a little time to calculate
write.csv(
  markers, 
  paste0(
    OutputDirectory, ObjName, Subset, 
    " cluster markers ", "res", RESOLUTION, ".csv"
  )
)
# read in, if it's been saved previously
# markers = read.csv(paste0(
#   OutputDirectory, ObjName, Subset, 
#   " cluster markers ", "res", RESOLUTION, ".csv"
#     )
#   )


############### CLUSTER / SAMPLE HEATMAP ##################
HM_object <- plot_heatmap (
  seurat_object = SeuratObj, downsample_n = 5000,
  markers = markers, top_n = 20, label_n = 2, 
  data_type = "logcounts",
  use_raster = TRUE
)

pdf(paste0(
  OutputDirectory, "heatmap", ObjName, Subset, 
  "res", RESOLUTION, "top20 genes per cluster.pdf"
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
  seurat_object = SeuratObj, aesX = "Sample", fill = "Cluster", 
  y_label = "Composition (Number of cells)", x_label = NULL, 
  y_lower_limit = 0, y_break = 1000,
  color_reverse = TRUE
)

pdf(paste0(
  densityplotDirectory, ObjName, Subset, 
  "res", RESOLUTION, "number of cells per sample and clusters barplot.pdf"
),
width = 6, height = 5.5, family = FONT_FAMILY
)
P2
dev.off()



P3 <- plot_bargraph (
  seurat_object = SeuratObj, aesX = "Cluster", fill = "Sample",
  y_label = "Composition (percentage of cells)", x_label = NULL, 
  y_lower_limit = 0, y_break = 0.2,
  position = "fill",
  plot_margin = unit(c(0.2, 0.5, 0.2, 0.5), "cm")
)

pdf(paste0(
  densityplotDirectory, ObjName, Subset, 
  "res", RESOLUTION, 
  "percent of cells per sample and clusters barplot.pdf"
), width = 6, height = 5.5, family = FONT_FAMILY
)
P3
dev.off()



P4 <- plot_bargraph (
  seurat_object = SeuratObj, aesX = "Cluster", fill = "Sample", 
  y_label = "Composition (Number of cells)", x_label = NULL, 
  y_lower_limit = 0, y_break = 1000
)

pdf(paste0(
  densityplotDirectory, ObjName, Subset, 
  "res", RESOLUTION, "number of cells per Cluster and Sample barplot.pdf"
), width = 6, height = 5.5, family = FONT_FAMILY
)
P4
dev.off()


################# CLUSTER / SAMPLE DIMPLOTS ###################
U1 <- plot_umap(
  seurat_object = SeuratObj, group_by = "Cluster",
  title = "Clusters", xlab = "UMAP1", ylab = "UMAP2",
  legend_position = "bottom",
  ncol_guide = 5, # this might be unnecessary?
  color_reverse = TRUE, label_clusters = TRUE
)

pdf(paste0(
  UMAPDirectory, ObjName, Subset, 
  "res", RESOLUTION, 
  "Clusters UMAP.pdf"
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
  seurat_object = SeuratObj, group_by = "Cluster",
  title = "Clusters", xlab = "UMAP1", ylab = "UMAP2",
  legend_position = "right",
  title_font_size = 16, x_font_size = 16, y_font_size = 16, 
  pt_size = 0.2, split_by = "Sample", ncol_dimplot = 3,
  color_reverse = TRUE
)

pdf(paste0(
  UMAPDirectory, ObjName, Subset, 
  "Cluster UMAP Iteration by sample.pdf"
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
