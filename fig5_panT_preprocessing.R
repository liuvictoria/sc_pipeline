# this module contains code that is more customized / less standard
# than the analysis from standard_viz
# none of this code is particularly reusable, except for its specific use case
################# LOAD UTILS ##############
# this is run locally and requires SeuratWrapper package
source("~/hpc/Sumner/R_Code/utils_HPC.R")

# capture session info, versions, etc.
write_experimental_configs()
#
#################### PAN CANCER T CELLS #####################
# # get GEO files and unzip; only need to do once
# get_GEO_unzip()

# get external counts and metadata
counts_files <- list.files(GEOdir, pattern = "counts")
meta_files <- list.files(GEOdir, pattern = "metadata")
meta_file = meta_files[1]

# create seurat objects of external samples
requested_count_files <- list()
for (i in 1: length(counts_files)){
  count_file = counts_files[i]

  cancer_subtype <- paste0(
    strsplit(count_file, "_")[[1]][2], "_",
    strsplit(strsplit(count_file, "CD")[[1]][2], ".counts")[[1]][1]
  )
  
  if (cancer_subtype %in% config$FILES) {
    requested_count_files[[cancer_subtype]] <- count_file
  }
}
stopifnot(length(requested_count_files) == length(config$FILES))

# doublet removal
for (cancer_subtype in names(requested_count_files)) {
  print(cancer_subtype)
  count_file <- requested_count_files[[cancer_subtype]]
  
  seurat_data <- read.csv(
    paste0(GEOdir, count_file), sep = "\t", row.names = 1
  )
  meta_data <- read.csv(
    paste0(GEOdir, meta_file), sep = "\t", row.names = 1
  )
  
  seurat_obj <- CreateSeuratObject(
    counts = seurat_data,
    min.features = 100,
    project = cancer_subtype,
    meta.data = meta_data
  )
  # only want tumor cells, not normal cells
  seurat_obj <- subset(
    seurat_obj,
    subset = loc == "T"
  )
  seurat_object <- GEX_QC(seurat_obj, cancer_subtype)
  
  # only want tumor cells, not normal cells
  seurat_object <- subset(
    seurat_object,
    subset = DoubletStatus == "singlet"
  )
  
  # save object
  saveRDS(seurat_object, file = paste0(RobjDirectory, cancer_subtype, ".rds"))
}
