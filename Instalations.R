library(pacman)
library(devtools)
library(BiocManager,quietly = TRUE)


p_load(BiocGenerics, DelayedArray, DelayedMatrixStats, limma, S4Vectors, SingleCellExperiment,
       SummarizedExperiment, batchelor, Matrix.utils,clustifyr,slingshot,glmGamPoi,
       SingleR,scDblFinder,miQC,GEOquery,Seurat,clustree,pheatmap,corrplot,extrafont,rjson,PerformanceAnalytics)
install.packages("cli")
BiocManager::install(c("celldex", "biomaRt"),force = TRUE)
install.packages("tidyverse")
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}


library(celldex)
library(biomaRt)


p_load_gh(c("jokergoo/ComplexHeatmap","cole-trapnell-lab/leidenbase","cole-trapnell-lab/monocle3","sjessa/ggmin"))
      
# I had an error installing the line below which was solved by following this:
#(https://programs.wiki/wiki/error-failed-to-install-unknown-package-from-github-http-error-403.html)
remotes::install_github("satijalab/sctransform", ref="develop")
remotes::install_github("mojaveazure/seurat-disk")
remotes::install_github("carmonalab/UCell", ref="v1.3")
remotes::install_github("carmonalab/scGate")# loading this will raise error so do not use p_load_gh for loading
remotes::install_github('satijalab/seurat-wrappers')


install.packages("magick")
install.packages("ggpubr")

remotes::install_github("carmonalab/ProjecTILs")# Entered (None)

source("/projects/compsci/USERS/alizae/GBM/matching_patients/R_Code/sc_pipeline/utils.R")#ggmin and umap were installed from here

