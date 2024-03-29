# this module contains code that is more customized / less standard
# than the analysis from standard_viz
# none of this code is particularly reusable, except for its specific use case
################# LOAD UTILS ##############
source("~/Documents/victoria_liu/matching_patients/R_Code/utils.R")

# capture session info, versions, etc.
# capture session info, versions, etc.
write_experimental_configs(code_file = "downstream")
#
#################### PAN CANCER T CELLS #####################
# # get GEO files
# get_GEO_unzip()

# get external counts and metadata
counts_files <- list.files(GEOdir, pattern = "counts")
meta_files <- list.files(GEOdir, pattern = "metadata")
meta_file = meta_files[1]

# create seurat objects of external samples
cancer_names <- c()
cancer_objects <- list()
for (i in 1: length(counts_files)){
  count_file = counts_files[i]

  cancer_subtype <- paste0(
    strsplit(count_file, "_")[[1]][2], "_",
    strsplit(strsplit(count_file, "CD")[[1]][2], ".counts")[[1]][1]
  )
  
  if (cancer_subtype %in% config$FILES) {
    print(cancer_subtype)

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
    cancer_names[length(cancer_names) + 1] <- cancer_subtype
    cancer_objects[[length(cancer_objects) + 1]] <- seurat_obj
  }
}
stopifnot(length(cancer_objects) == length(config$FILES))

# doublet removal
for (idx in 1:length(cancer_names)) {
  print(cancer_names[idx])
  cancer_object <- GEX_QC(cancer_objects[[idx]], FILES[idx])
  
  # save object
  saveRDS(cancer_object, file = paste0(RobjDirectory, FILES[idx], ".rds"))
  cancer_objects[[idx]] <- cancer_object
}


# merge into one Seurat object
external_pan <- merge(
  x = cancer_objects[[1]], y = cancer_objects[-1], 
  add.cell.ids = config$FILES
  )
# beautify metadata
add_tumor_suffix <- function(word) {
  return (paste0(word, "_tumor"))
}
external_pan$Sample <- unlist(lapply(external_pan$libraryID, add_tumor_suffix))
external_pan$Group <- "tumor"
external_pan$Patient <- external_pan$patient
external_pan$Type <- external_pan$cancerType
external_pan$Assignment <- external_pan$meta.cluster

# external_pan$libraryID <- NULL
# external_pan$patient <- NULL
# external_pan$cancerType <- NULL
# external_pan$meta.cluster <- NULL



# saving takes a long time
saveRDS(
  external_pan,
  paste0(RobjDirectory, ObjName, Subset, "_external.rds")
)
external_pan <- readRDS(paste0(RobjDirectory, ObjName, Subset, "_external.rds"))


# normalize data
external_pan <- GEX_normalization(external_pan)
# CC scaling regression
external_pan <- GEX_cc_regression(external_pan)

# PCA and harmony batch correciton
DefaultAssay(external_pan) <- "RNA"

external_pan <- external_pan %>% 
  RunPCA(
    features = VariableFeatures(object = external_pan)
  ) %>% 
  RunHarmony(
    group.by.vars = config$batch_norm_by,
    reduction.save = "harmonyRNA"
  )

refObj <- readRDS(paste0(
  RobjDir, "TAtlas/", 
  "GEXTAtlas_res0.3.rds"
  ))

pan_T <- merge(x = external_pan, y = refObj)
genes_atlas <- row.names(Loadings(refObj, reduction = "pca"))

genes_use <- intersect(genes_atlas, rownames(external_pan))
refObj_loadings <- Loadings(refObj, reduction = "pca")[genes_use, ]
refObj_embeddings <- Embeddings(refObj, reduction = "pca")
refObj_umap <- Embeddings(refObj, reduction = "umapRNA")

merged_seurat2 <- NormalizeData(merged_seurat2, normalization.method = "LogNormalize", scale.factor = 10000)
merged_seurat2 <- FindVariableFeatures(merged_seurat2, selection.method = "vst", nfeatures = 2000)

#scale data
merged_seurat2 <- ScaleData(merged_seurat2,features = genes_use)

#Run PCA
merged_seurat2 <- RunPCA(merged_seurat2, features = VariableFeatures(object = merged_seurat2))

merged_seurat2@meta.data$ID <- merged_seurat2@meta.data$orig.ident
data <- as.matrix(merged_seurat2@assays$RNA@scale.data[genes_use,])
SeuratObj.embeddings <- apply(refSeuratObj_loadings, 2 , FUN = function(i) colSums(i*data))
library(caret)
train.umap1 <- knnreg(refSeuratObj_embeddings, refSeuratObj_umap[, 1], k = 100)
train.umap2 <- knnreg(refSeuratObj_embeddings, refSeuratObj_umap[, 2], k = 100)
test.umap1 <- predict(train.umap1, SeuratObj.embeddings)
test.umap2 <- predict(train.umap2, SeuratObj.embeddings)
test.umap <- cbind(test.umap1, test.umap2)
rownames(test.umap) <- rownames(SeuratObj.embeddings)
colnames(test.umap) <- paste0('umap_', 1:2)
merged_seurat2@reductions$refSeuratObjumap <- CreateDimReducObject(key='umap_', embeddings = test.umap)
DimPlot(merged_seurat2,reduction = "refSeuratObjumap")
DimPlot(merged_seurat2,reduction = "refSeuratObjumap",split.by = "library_id")
#merged_seurat2$
metaall=merged_seurat2@meta.data
metaall$TissueType=metaall$orig.ident
metaall$TissueType[str_detect(metaall$TissueType,pattern = "CNSTM")]="Glioma"
metaall$TissueType[str_detect(metaall$TissueType,pattern = "MDAG")]="Glioma"

#metaall2=merged_seurat2@meta.data
merged_seurat2@meta.data=metaall
DimPlot(merged_seurat2,reduction = "refSeuratObjumap",group.by = "TissueType")
DimPlot(merged_seurat2,reduction = "refSeuratObjumap",split.by = "TissueType",ncol = 4,pt.size = 2,group.by = "Assignment")
DimPlot(merged_seurat2,reduction = "refSeuratObjumap",split.by = "tissue",ncol = 4,pt.size = 2,group.by = "Assignment")

saveRDS(merged_seurat2,"~/Library/CloudStorage/Box-Box/Yun lab projects/scRNAseq Project/Human Glioma/pan myeloid/panmyeloidobjmergedwithourmyeloid2.rds")
merged_seurat2=readRDS("~/Library/CloudStorage/Box-Box/Yun lab projects/scRNAseq Project/Human Glioma/pan myeloid/panmyeloidobjmergedwithourmyeloid2.rds")

library(Nourpal)
TissueColors= Nour_pal("all")(length(levels(as.factor(merged_seurat2$TissueType))))
names(TissueColors)=levels(as.factor(merged_seurat2$TissueType))

U1=DimPlot(merged_seurat2, reduction = "refSeuratObjumap",shuffle = T, raster = T,pt.size = 0.5,
           group.by = "TissueType" ,cols =  TissueColors)+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+
  FontSize(x.title = 20, y.title = 20,main = 20)+ theme(legend.position = "bottom")+labs(title="TissueType")+
  guides(colour = guide_legend(override.aes = list(size=5),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=15),ncol = 5))
AssignmentColors=c(`a-microglia` = "#003F5C", `s-mac 1` = "#2D5187", `AP-microglia` = "#74599E", 
                   DCs = "#B863A5", `h-microglia` = "#E46388", `i-microglia` = "#E25159", 
                   `MDSC` = "#F77B38", Proliferating = "#F6CC4B", `s-mac 2` = "#8DB032"
)

Info=merged_seurat2@meta.data
Assign= list()
Assign[["NonCancer"]]=c("N")
Assign[["PeripheralBlood"]]=c("P")
Assign[["OtherTumors"]]=c("T")


Info$Origin=NA
for(i in 1:length(Assign)){
  Info$Origin[Info$tissue %in% Assign[[i]]]=names(Assign[i])
}
Info$Origin[is.na(Info$Origin)]="GBM"
merged_seurat2@meta.data=Info

U2=DimPlot(merged_seurat2, reduction = "refSeuratObjumap",raster = T,pt.size = 0.5,
           group.by = "Assignment" ,cols =  AssignmentColors,split.by = "Origin")+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+
  FontSize(x.title = 20, y.title = 20,main = 20)+ theme(legend.position = "bottom")+labs(title="TissueOrigin")+
  guides(colour = guide_legend(override.aes = list(size=5),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=15),ncol = 5))

U3=DimPlot(merged_seurat2, reduction = "refSeuratObjumap",raster = T,pt.size = 1,
           group.by = "Assignment" ,cols =  AssignmentColors,split.by = "TissueType")+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+
  FontSize(x.title = 20, y.title = 20,main = 20)+ theme(legend.position = "bottom")+labs(title="TissueOrigin")+
  guides(colour = guide_legend(override.aes = list(size=5),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=15),ncol = 11))

mSeuratObj.sce=as.SingleCellExperiment(merged_seurat2)

reducedDimNames(mSeuratObj.sce)

ElbowPlot(merged_seurat2,ndims = 25)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

merged_seurat2 <- CellCycleScoring(merged_seurat2, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
DimPlot(merged_seurat2,group.by = "Phase")
merged_seurat2 <- ScaleData(merged_seurat2, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(merged_seurat2))
merged_seurat2@meta.data$Patient[is.na(merged_seurat2@meta.data$Patient)]=merged_seurat2@meta.data$patient
PatientColors= Nour_pal("all")(length(levels(as.factor(merged_seurat2$Patient))))
names(PatientColors)=levels(as.factor(merged_seurat2$Patient))

#initial clustering
library(harmony)
merged_seurat2 <- RunHarmony(merged_seurat2, c("TissueType"))
merged_seurat2 <- FindNeighbors(merged_seurat2,reduction = "harmony", dims = 1:20)

merged_seurat2 <- RunUMAP(merged_seurat2, reduction = "harmony",dims = 1:20 )

#merged_seurat2 <- FindNeighbors(merged_seurat2, dims = dimz)
resolution=0.25
merged_seurat2 <- FindClusters(merged_seurat2, resolution = 0.5)
DimPlot(merged_seurat2,reduction = "umap")
DimPlot(merged_seurat2,split.by = "TissueType",label = T,reduction = "umap")
DimPlot(merged_seurat2,split.by = "TissueType",label = T)

#head(Idents(merged_seurat2), 5)
#merged_seurat2 <- RunUMAP(merged_seurat2, dims = dimz)

DimPlot(merged_seurat2)

#supplementary fig 5a
Info=merged_seurat2@meta.data
#Info$oldclusters=Info$Cluster
#Rename Clusters
clusters=levels(as.factor(merged_seurat2@meta.data$seurat_clusters))
abb="PMC"
for(j in 1:length(clusters)){
  if(j<10){
    Info$Cluster[Info$seurat_clusters==j-1]=paste0(abb,"0",j)
  }else{
    Info$Cluster[Info$seurat_clusters==j-1]=paste0(abb,j)
  }
}
merged_seurat2@meta.data=Info
# Cluster color 
set.seed(001) # just to make it reproducible

ClusterColors= sample(Nour_pal("all",reverse = T)(length(levels(as.factor(merged_seurat2@meta.data$Cluster)))))

names(ClusterColors)=levels(as.factor(merged_seurat2@meta.data$Cluster))

DimPlot(merged_seurat2, reduction = "refSeuratObjumap",raster = T,pt.size = 1,
        group.by = "Cluster" ,cols =  ClusterColors,split.by = "TissueType")+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+
  FontSize(x.title = 20, y.title = 20,main = 20)+ theme(legend.position = "bottom")+labs(title="Cluster")+
  guides(colour = guide_legend(override.aes = list(size=5),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=15),ncol = 11))
DimPlot(merged_seurat2, reduction = "umap",raster = T,pt.size = 1,
        group.by = "Cluster" ,cols =  ClusterColors,split.by = "TissueType")+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+
  FontSize(x.title = 20, y.title = 20,main = 20)+ theme(legend.position = "bottom")+labs(title="Cluster")+
  guides(colour = guide_legend(override.aes = list(size=5),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=15),ncol = 11))

DimPlot(merged_seurat2, reduction = "refSeuratObjumap",raster = T,pt.size = 1,
        group.by = "Cluster" ,cols =  ClusterColors,split.by = "Cluster",ncol = 11)+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+
  FontSize(x.title = 20, y.title = 20,main = 20)+ theme(legend.position = "bottom")+labs(title="Cluster")+
  guides(colour = guide_legend(override.aes = list(size=5),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=15),ncol = 11))


DimPlot(merged_seurat2, reduction = "umap",group.by = "Cluster",cols =ClusterColors, pt.size = 0.2,label = T)+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+
  FontSize(x.title = 16, y.title = 16,main = 16)+ theme(legend.position = "right",text =element_text(size=15) )+
  guides(colour = guide_legend(title="Clusters", override.aes = list(size=5),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=15)))

#Info= Info %>% filter(Cluster!="MC07")
#merged_seurat2= subset(merged_seurat2, cells=row.names(Info))

# Get number of cells per cluster and per Sample
write.csv(as.matrix(table(merged_seurat2@meta.data$Cluster,merged_seurat2@meta.data$Patient)),
          file = paste0(OutputDirectory,ObjName,Subset," number of cells per cluster and sample.csv"))
write.csv(as.matrix(table(merged_seurat2@meta.data$Cluster,merged_seurat2@meta.data$Type)),
          file = paste0(OutputDirectory,ObjName,Subset," number of cells per cluster and Type.csv"))
write.csv(as.matrix(table(merged_seurat2@meta.data$Cluster,merged_seurat2@meta.data$sex)),
          file = paste0(OutputDirectory,ObjName,Subset," number of cells per cluster and Sex.csv"))
write.csv(as.matrix(table(merged_seurat2@meta.data$Cluster,merged_seurat2@meta.data$Fragment)),
          file = paste0(OutputDirectory,ObjName,Subset," number of cells per cluster and Fragment.csv"))
write.csv(as.matrix(table(merged_seurat2@meta.data$Cluster,merged_seurat2@meta.data$Grade)),
          file = paste0(OutputDirectory,ObjName,Subset," number of cells per cluster and Grade.csv"))
#make heatmap
resolution=0.2
OutputDirectory="/Users/nourhan/Library/CloudStorage/Box-Box/Yun lab projects/scRNAseq Project/Human Glioma/pan myeloid/Output/"
ObjName="Pan"
Subset=" Myeloid_harmony with tissue type"


Idents(merged_seurat2)=Info$Cluster
markers = FindAllMarkers(merged_seurat2,logfc.threshold = 0.25,test.use = "wilcox",only.pos = T )
write.csv(markers,paste0(OutputDirectory,ObjName,Subset," Cluster markers ","res",resolution,".csv"))
#markers = read.table(paste0(OutputDirectory,ObjName,Subset," Cluster markers ","res",resolution,".csv"))
markers=markers[order(markers$cluster),]
top20.markers <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
top20.markers$cluster=as.character(top20.markers$cluster)
top20.markers=top20.markers[order(top20.markers$cluster),]
set.seed(12)
Idents(merged_seurat2)=Info$Patient
subset=subset(merged_seurat2, downsample=500)
#DimPlot(subset,split.by = "Sample",group.by = "Cluster")
merged_seurat2.sce=as.SingleCellExperiment(subset)
plot.data<-as.data.frame(assay(merged_seurat2.sce, "logcounts"))
plot.data<-plot.data[top20.markers$gene,]
plot.data <- plot.data - rowMeans(plot.data)
plot.data=na.omit(plot.data)
InfoS=subset@meta.data
column_annot <-InfoS[,c("Cluster","Patient","TissueType"),drop=F]
column_annot$Patient = as.factor(as.character(column_annot$Patient))
column_annot=with(column_annot, column_annot[order(Patient), , drop=F])
column_annot=with(column_annot, column_annot[order(TissueType), , drop=F])
column_annot=with(column_annot, column_annot[order(Cluster), , drop=F])

plot.data<-plot.data[,row.names(column_annot)]
column.col= PatientColors


column.colors=list()
column.colors[["Patient"]]<-column.col
column.colors[["Cluster"]]<-ClusterColors
column.colors[["TissueType"]]<-TissueColors

Patient=as.matrix(column_annot[,c("Patient"),drop=F])
Cluster=as.matrix(column_annot[,c("Cluster"),drop=F])
TissueType=as.matrix(column_annot[,c("TissueType"),drop=F])

colanno <- columnAnnotation(df=column_annot,
                            show_annotation_name =T,show_legend = F,col=column.colors)
top2 <- markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
genes= top2$gene
genes= c("CCL4","CCL3","IL1B","CD83","IL1A","FN1","MIF","BNIP3","S100A8","S100A9","S100A4","MALAT1","SPRY1","MKI67", "P2RY12","TMEM119","HLA-DRA","CD1C","BATF3",
         "HLA-DRB5",  "CCL3L1", "CCL4L2" , "ISG15", 
         "CXCL10",  "STMN1", "RNASE1", 
         "SELENOP",  "SPRY1", "AREG", "HLA-DQA1")
rows=rowAnnotation(sel = anno_mark(at = match(genes,row.names(plot.data)), labels = genes,labels_gp =gpar(col = "black", fontsize = 9)))

col<- circlize::colorRamp2(breaks = c(-2, 0, 2), colors = c("#007dd1", "white", "#ab3000"))
HM=Heatmap(name="logcounts",as.matrix(plot.data),cluster_rows =F,cluster_columns = F,top_annotation = colanno, right_annotation = rows,row_names_gp = gpar(fontsize=1),
           col = col,show_column_names= F,show_row_names=F,border = F,show_heatmap_legend = F,use_raster = T)        
lgd=Legend(title = "logcounts", at=  c(-2,0, 2),col_fun = col)
lgd1=Legend(labels = levels(as.factor(column_annot$Cluster)),title="Cluster",legend_gp = gpar(fill=ClusterColors,fontsize=5))
lgd2=Legend(labels = levels(as.factor(column_annot$Patient)),title="Patient",legend_gp = gpar(fill=column.col,fontsize=5))
lgd4=Legend(labels = levels(as.factor(column_annot$TissueType)),title="TissueType",legend_gp = gpar(fill=TissueColors,fontsize=5))

pdf(paste0(OutputDirectory,"heatmap",ObjName,Subset,"res",resolution,"top20 genes per cluster .pdf"),width=20,height=15)
draw(HM,heatmap_legend_list = list( lgd,lgd1,lgd4, lgd2), heatmap_legend_side = "right")
dev.off()
Genes=c('IL1B','IL6','IL12A','IL23A','CXCL9'
        ,'NOS2','CD80','CD86','TNF','MRC1','IL1R1',
        'CCL17','TGFB1','IGF1','FN1','CCL1','IL10','TNFSF14',
        'MERTK','CXCL13','CD163','VEGFA',"FN1", "S100A4", "IBSP",
        "S100A8","CCL3", "CCL4", "CD83", "IL1B", "IL1A","HLA-DRA",
        "CD74", "CD14", "APOE","B2M", "SPRY1", "BHLHE41", "SORL1",
        "IFNGR1", "MALAT1", "CX3CR1", "BNIP3", "MIF","MKI67", "IER3",
        "NFKBIZ","TMEM119","P2RY12","CX3CR1","SLC2A5","CST3","P2RY13","CCL2","EGR2","CCL4","CD83","EGR3")

pdf(paste0(OutputDirectory,"dotplot ",ObjName,Subset, "by Cluster.pdf"),width=14,height=5)
DotPlot(merged_seurat2,group.by = "Cluster",dot.scale = 5 ,features = unique(Genes) ,scale = T)+scale_colour_viridis_c(option = "plasma")+
  ggmin::theme_min()+ RotatedAxis()+theme(axis.title = element_blank())
DotPlot(merged_seurat2,group.by = "Cluster" ,features = c('IL1B','IL6','IL12A','IL23A','CXCL9',"NFKBIZ"
                                                          ,'NOS2','CD80','CD86','TNF','MRC1','IL1R1',
                                                          'CCL17','TGFB1','IGF1','FN1','CCL1','IL10','TNFSF14',
                                                          'MERTK','CXCL13','CD163','VEGFA'),
        scale = T)+RotatedAxis()+RotatedAxis()+scale_colour_viridis_c(option = "plasma")
DotPlot(merged_seurat2,group.by = "Cluster" ,features = c('ITGAM','ITGAX','CD14','CD33','CD68',"MSR1"
                                                          ,'HLA-DRB1','IFNGR1','CD1C','S100A4','MIF','LYZ',
                                                          'CX3CR1','TMEM119','P2RY12','BHLHE41','SPRY1','NFKBIZ','CCL4',
                                                          'CCL3','ITGA4'),
        scale = T)+RotatedAxis()+RotatedAxis()+scale_colour_viridis_c(option = "plasma")
DotPlot(merged_seurat2,group.by = "Cluster" ,features = c('ITGAM','ITGAX','CD14','CD33','CD68',"MSR1"
                                                          ,'HLA-DRB1',"CD74","B2M",'IFNGR1','CD1C','S100A4','MIF','LYZ',
                                                          'CX3CR1','TMEM119','P2RY12','BHLHE41','SPRY1','CCL4',
                                                          'CCL3','ITGA4','IL1B','IL6','IL12A','IL23A','CXCL9',"NFKBIZ"
                                                          ,'NOS2','CD80','CD86','TNF','MRC1','IL1R1',
                                                          'CCL17','TGFB1','IGF1','FN1','CCL1','IL10','TNFSF14',
                                                          'MERTK','CXCL13','CD163','VEGFA',"MKI67"),
        scale = T)+RotatedAxis()+RotatedAxis()+scale_colour_viridis_c(option = "plasma")

DotPlot(merged_seurat2,group.by = "Cluster" ,features = unique(genes),
        scale = T)+RotatedAxis()+RotatedAxis()+scale_colour_viridis_c(option = "plasma")

dev.off()

DotPlot(merged_seurat2,group.by = "Cluster" ,features = c('KIT','TPSAB1','CPA3','LILRA4','GZMB',"IL3RA"
                                                          ,'CLEC9A',"FLT3","IDO1",'CD1C','FCER1A','HLA-DQA1','LAMP3','CCR7',
                                                          'FSCN1','FCN1','S100A9','S100A8','FCGR3A','LST1',
                                                          'LILRB2','INHBA','IL1RN','CCL4','NLRP3','EREG','IL1B',"LYVE1"
                                                          ,'PLTP','SEPP1','C1QC','C1QA','APOE',"CD14","CD16","TREM2"),
        scale = T)+RotatedAxis()+RotatedAxis()+scale_colour_viridis_c(option = "plasma")
FeaturePlot(merged_seurat2,features = "CLEC9A",reduction = "umap")








#################### AGGR CELLS VARIOUS CLUSTIFYR REFS #################
{
  ################# ASSIGNMENT BARPLOTS ###################
P5 <- plot_bargraph (
  seurat_object = SeuratObj, aesX = "Sample", fill = "Assignment",
  y_label = "Composition (Number of cells)", x_label = NULL,
  y_lower_limit = 0, y_break = 1000, 
  title = paste0(analyses$which_assignment, " Assignment")
)

pdf(paste0(
  densityplotDirectory, ObjName, Subset, 
  "res", RESOLUTION,
  "_number of cells per sample and Assignment (", 
  analyses$which_assignment,
  ") barplot cbmc reference.pdf"
), width = 6, height = 5.5, family = FONT_FAMILY
)
print(P5)
dev.off()



P6 <- plot_bargraph (
  seurat_object = SeuratObj, aesX = "Assignment", fill = "Sample",
  y_label = "Composition (percentage of cells)", x_label = NULL, 
  y_lower_limit = 0, y_break = 0.2,
  position = "fill",
  plot_margin = unit(c(0.2, 0.5, 0.2, 0.5), "cm"),
  title = paste0(analyses$which_assignment, " Assignment")
)
pdf(paste0(
  densityplotDirectory, ObjName, Subset, 
  "res", RESOLUTION, 
  "percent of cells per sample and Assignment (", 
  analyses$which_assignment,
  ") barplot cbmc reference.pdf"
), width = 6, height = 5.5, family = FONT_FAMILY
)
print(P6)
dev.off()


P7 <- plot_bargraph (
  seurat_object = SeuratObj, aesX = "Assignment", fill = "Sample",
  y_label = "Composition (Number of cells)", x_label = NULL, 
  y_lower_limit = 0, y_break = 1000,
  title = paste0(analyses$which_assignment, " Assignment")
)

pdf(paste0(
  densityplotDirectory, ObjName, Subset,
  "res", RESOLUTION, 
  "_percentage of cells per major population and Sample (", 
  analyses$which_assignment,
  ") barplot cbmc reference.pdf"
), width = 8, height = 6, family = FONT_FAMILY
)
print(P7)
dev.off()

P8 <- plot_bargraph (
  seurat_object = SeuratObj, aesX = "Sample", fill = "Assignment",
  y_label = "Composition (percentage of cells)", x_label = NULL,
  y_lower_limit = 0, y_break = 1000, position = "fill",
  title = paste0(analyses$which_assignment, " Assignment")
)

pdf(paste0(
  densityplotDirectory, ObjName, Subset, 
  "res", RESOLUTION,
  "_percentage of cells per sample and Assignment (", 
  analyses$which_assignment,
  ") barplot cbmc reference.pdf"
), width = 6, height = 5.5, family = FONT_FAMILY
)
print(P8)
dev.off()




# multiple
bar_plots <- ggarrange(P2, P5, P8, P3, P4, P6, P7, ncol = 2)
pdf(paste0(
  densityplotDirectory, ObjName, Subset, 
  "res", RESOLUTION, "_all barplots (Assignment ", 
  analyses$which_assignment,
  ") cbmc reference.pdf"
), width = 18, height = 22, family = FONT_FAMILY
)
print(bar_plots)
dev.off()


################# ASSIGNMENT DIMPLOTS ###################
U4 <- plot_umap(
  seurat_object = SeuratObj, group_by = "Assignment",
  reduction = paste0("umap", analyses$viz_clustering),
  title = paste0(analyses$which_assignment, " Assignment"), 
  xlab = "UMAP1", ylab = "UMAP2",
  legend_position = "bottom",
  ncol_guide = 4,
  label_clusters = TRUE, 
  label_size = 5,
  color_reverse = FALSE
)

pdf(paste0(
  UMAPDirectory, ObjName, Subset, 
  "_res", RESOLUTION, 
  "_", analyses$viz_clustering,
  "Clusters_Assignment (", 
  analyses$which_assignment,
  ") UMAP cbmc reference.pdf"
), width = 8, height = 6, family = FONT_FAMILY
)
print(U4)
dev.off()

U5 <- plot_umap(
  seurat_object = SeuratObj, group_by = "Assignment",
  reduction = paste0("umap", analyses$viz_clustering),
  title = paste0(analyses$which_assignment, " Assignment"),
  xlab = "UMAP1", ylab = "UMAP2",
  legend_position = "bottom",
  title_font_size = 16, x_font_size = 16, y_font_size = 16, 
  pt_size = 0.2, split_by = "Sample", ncol_dimplot = 2
)

pdf(paste0(
  UMAPDirectory, ObjName, Subset, 
  "_res", RESOLUTION, 
  "_", analyses$viz_clustering,
  "Clusters_Assignment (", 
  analyses$which_assignment,
  ") UMAP Iteration by sample cbmc reference.pdf"
), width = 12, height = 2.5 * length(unique(SeuratObj$Sample)), 
family = FONT_FAMILY
)
print(U5)
dev.off()


U6 <- plot_umap(
  seurat_object = SeuratObj, 
  group_by = "Assignment",
  reduction = paste0("umap", analyses$viz_clustering),
  title = paste0(analyses$which_assignment, " Assignment"),
  xlab = "UMAP1", ylab = "UMAP2",
  legend_position = "bottom",
  title_font_size = 16, x_font_size = 16, y_font_size = 16, 
  pt_size = 0.5, split_by = "Assignment", ncol_dimplot = 2
)

pdf(paste0(
  UMAPDirectory, ObjName, Subset, 
  "_res", RESOLUTION, 
  "_", analyses$viz_clustering,
  "Clusters_Assignment (", 
  analyses$which_assignment,
  ") UMAP Iteration by assignment cbmc reference.pdf"
), width = 12, height = 2.5 * length(unique(SeuratObj$Sample)), family = FONT_FAMILY
)
print(U6)
dev.off()




# multiple
plots2 <- ggarrange(U1, U4, U2, ncol = 3)

pdf(paste0(
  UMAPDirectory, ObjName, Subset, 
  "_res", RESOLUTION, "(Assignment ", 
  analyses$which_assignment,
  ") all UMAPs cbmc reference.pdf"
), width = 20, height = 7, family = FONT_FAMILY
)
print(plots2)
dev.off()



}