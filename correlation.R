################# LOAD UTILS ##############
source("~/Documents/victoria_liu/matching_patients/R_Code/utils.R")

# capture session info, versions, etc.
write_experimental_configs(
  code_file = "correlation"
)

################# RUN ONCE ##################
# seurat_atlas_all <- preprocess_atlas_objects_corr(
#   paste0(RobjDir, "GBMAtlas/", "Allhuman-11-3-21.rds"),
#   "all_fragment"
# )
normalization_counts_fragment <- read_normalization_counts(
  "all_fragment_normalization_counts.csv"
)

seurat_glioma_all <- preprocess_atlas_objects_corr(
  paste0(RobjDir, "GBMAtlas/", "GliomaClusters-11-3-21 patients renamed.rds"),
  "all_fragment_glioma"
)

seurat_glioma_all$Assignment <- seurat_glioma_all$Cluster

normalization_counts_glioma <- read_normalization_counts(
  "all_fragment_glioma_normalization_counts.csv"
)

seurat_myeloid_all <- preprocess_atlas_objects_corr(
  paste0(RobjDir, "GBMAtlas/", "MyeloidClusters-11-4-21-patients renamed.rds"),
  "all_fragment_myeloid"
)
normalization_counts_myeloid <- read_normalization_counts(
  "all_fragment_myeloid_normalization_counts.csv"
)

# group = "ndGBM" OR "rGBM" or "Atlas"
# group <- str_sub(strsplit(Subset, "-")[[1]][1], 2, -1)
group <- "Atlas"
stopifnot (group == "ndGBM" | group == "rGBM" | group == "Atlas")
seurat_T_all <- preprocess_atlas_objects_corr(
  paste0(
    RobjDir, "T", group, "/", 
    "GEXT", group, "_res", analyses$correlation_T_resolution, ".rds"
    ),
  paste0(group, "_fragment_T"),
  flip_patient_sample = F
)
normalization_counts_T <- read_normalization_counts(
  paste0(group, "_fragment_T_normalization_counts.csv")
)


############## CORRELATION BETWEEN CELL TYPES BULK NORMALIZED ###############
# some variables
population <- str_sub(Subset, 2, -1)
assignment_types <- c("Assignment", "Cluster")
T_correlation_with <- c("myeloid", "glioma")
normalizations <- c("fragment_all", "fragment_subtype", "none")

# T cell subset
seurat_T <- subset(
  seurat_T_all,
  subset = TumorType == population
)
seurat_T$Cluster <- seurat_T$Assignment


# do myeloid or glioma, one at a time
for (correlation_with in T_correlation_with) {
  
  # create directories for outputs
  SubDirectory <- paste0(SubtypeCorrDirectory, correlation_with, "/")
  if (! dir.exists(SubDirectory)) {dir.create(SubDirectory)}
  
  # hard code which correlation type it is (I mean, there's only two, sooo)
  if (correlation_with == "myeloid") {
    seurat_correlation_with <- seurat_myeloid_all
  } else if (correlation_with == "glioma") {
    seurat_correlation_with <- seurat_glioma_all
  } else {
    stop("invalid correlation Seurat object")
  }
  
  # figure out patient subset
  seurat_correlation_with_subset <- subset(
    seurat_correlation_with,
    subset = TumorType == population
  )
  
  # merge objects
  seurat_sample <- merge(
    x = seurat_correlation_with_subset, y = seurat_T, 
    add.cell.ids = c(correlation_with, "T")
  )
  
  for (assignment_type in assignment_types) {
    if (correlation_with == "glioma" & assignment_type == "Assignment") {
      next
      }
    # create directories for outputs
    SubtypeDirectory <- paste0(SubDirectory, assignment_type, "/")
    if (! dir.exists(SubtypeDirectory)) {dir.create(SubtypeDirectory)}
    
    my_data <- table(
      seurat_sample$Fragment,
      seurat_sample@meta.data[[assignment_type]]
    )
    
    # take away fragments that don't contain both T and myeloid cells
    rows <- intersect(
      unique(seurat_T$Fragment), unique(seurat_correlation_with_subset$Fragment)
    )
    my_data <- my_data[rownames(my_data) %in% rows, ]
    
    # create custom column order
    col_order <- c(
      sort(unique(seurat_correlation_with_subset@meta.data[[assignment_type]])), 
      sort(unique(seurat_T@meta.data[[assignment_type]]))
    )
    # # CD4 naive cells have too few counts
    # col_order <- col_order[!is.na(col_order) & col_order != "CD4_NaiveLike"]
    my_data <- my_data[, col_order]
    
    for (normalization in normalizations) {
      # normalize
      my_data_normalized <- my_data
      for (subtype in colnames(my_data_normalized)) {
        for (fragment in rownames(my_data_normalized)) {
          if (normalization == "none") {
            normalize_denominator <- 1
          } else if (normalization == "fragment_all") {
            normalize_denominator <- normalization_counts_fragment[[
              fragment
            ]]
          } else if (normalization == "fragment_subtype") {
            if (
              subtype %in% unique(seurat_myeloid_all@meta.data[[assignment_type]])
              ) {
              normalize_denominator <- normalization_counts_myeloid[[fragment]]
            } else if (
              subtype %in% unique(seurat_glioma_all@meta.data[[assignment_type]])
              ) {
              normalize_denominator <- normalization_counts_glioma[[fragment]]
            } else if (
              subtype %in% unique(seurat_T@meta.data[[assignment_type]])
              ) {
              normalize_denominator <- normalization_counts_T[[fragment]]
            } else {
              print(subtype)
              print(fragment)
              stop("normalization denominator could not be found")
            }
          }
          stopifnot(! is.null(normalize_denominator))
          
          my_data_normalized[fragment, subtype] <- 
            my_data_normalized[fragment, subtype] / normalize_denominator
        }
      }
      
      res <- cor(my_data_normalized)
      res_with_conf <- cor.mtest(my_data_normalized, conf.level = 0.95)
      
      pdf(paste0(
        SubtypeDirectory, ObjName, Subset, 
        "__Sample ", population,
        "__Normalization ", normalization,
        "__Subtype Correlation heatmap.pdf"
      ), width = 15, height = 15, 
      family = FONT_FAMILY
      )

      corrplot(
        res, 
        p.mat = res_with_conf$p,
        addCoef.col = "black",
        number.cex = 0.8,
        number.font = 1,
        type = "upper", order = "original",
        insig = "pch",
        tl.col = "black", tl.srt = 45, 
        tl.cex = 1.5,
        cl.cex = 1.5,
        mar = c(0, 0, 10, 0),
        col = rev(COL2("RdBu")),
        title = paste0(
          "Sample: ", population,
          " __ Normalization: ", normalization,
          " __ Subtype Correlation heatmap"
        )
      )
      dev.off()
    }
  }
}
