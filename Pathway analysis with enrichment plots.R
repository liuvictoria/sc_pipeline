library(presto)
setwd('/Users/nourhan/Box/Yun lab projects/scRNAseq Project/2808combined old and new/Output/AllClusters/Pathway per assignment per sample/')
setwd("~/Box/Yun lab projects/scRNAseq Project/2808combined old and new/Output/AllClusters/pathway2")
Idents(SeuratObj)=SeuratObj@meta.data$Assignment
CellInfo=SeuratObj@meta.data
objname="x2808"
Cluster= levels(as.factor(CellInfo$Group))
Iterate= levels(as.factor(CellInfo$Assignment))
#Iterate=c("Myeloid","Glioma","Bcells",  "Microglia",  "NKCells", "Tcells")
load("~/Box/Yun lab projects/ESC Analysis October 2020/r objects/C8mousegenesconverted.RData")
library(presto)
library(msigdbr)
library(dplyr)
library(fgsea)
library(tibble)
library(pheatmap)
Z=""
pathways = "CP:KEGG"
names= c("KEGG")
category=c("C2")
for (j in 1:length(Iterate)){
  Z=Iterate[j]
  print(Z)
  Annot.subset <-subset(SeuratObj, cells=WhichCells(SeuratObj,idents =Z ))
  Annot=Annot.subset@meta.data
  #2 Run Test
  Clusters.genes <- wilcoxauc(Annot.subset , 'Group')
  write.csv(Clusters.genes,paste(Z,"DE genes between Samples in", Z,"-presto output.csv"))
  #check if test was successful
  dplyr::count(Clusters.genes, group)
  #pathways = c("H","C7","C8","CP:KEGG","CP:BIOCARTA","C6")
  #names= c("Hallmark","C7","C8","KEGG","Biocarta","C6")
  #category=c("H","C7","C8","C2","C2","C6")
  pathways = c("CP:BIOCARTA")
  names= c("BIOCARTA")
  category=c("C2")
  for (i in 1: length(pathways)) {
    cat=pathways[i]
    print(cat)
    name=names[i]
    
    if (category[i]=="C8") {
      fgsea_sets<- mousemygo
      Annot.pathway2=as.data.frame(names(mousemygo))
      names(Annot.pathway2)="pathway"
    } else {
      if (category[i]=="C2") {
        m_df_H<- msigdbr(species = "Mus musculus", category = "C2",subcategory = cat)
      } else {
        m_df_H<- msigdbr(species = "Mus musculus", category = cat)
      }
      Annot.pathway2=as.data.frame(levels(as.factor(m_df_H$gs_name)))
      names(Annot.pathway2)="pathway"
      fgsea_sets<- m_df_H %>% split(x = .$gene_symbol, f = .$gs_name)
    }
    Clusters=levels(as.factor(Annot$Group))
    
    for (i in 1: length(Clusters)){
      X=  Clusters[i]
      print(X)
      Genes<-Clusters.genes %>%
        dplyr::filter(group == Clusters[[i]]) %>%
        arrange(desc(auc)) %>% 
        dplyr::select(feature, auc)
      ranks<- deframe(Genes)
      head(ranks)
      
      fgseaRes<- fgsea(fgsea_sets, stats = ranks,minSize=10)
      fgseaResTidy <- fgseaRes %>%
        as_tibble() %>%
        arrange(desc(NES))
      fgseaResTidy %>% 
        dplyr::select(-leadingEdge, -ES, -log2err) %>% 
        arrange(padj) %>% 
        head()
      fgseaResTidy$Enrichment = ifelse(fgseaResTidy$NES > 0, "Up-regulated", "Down-regulated")
      Filtidy<-fgseaResTidy %>% filter(padj < 0.05) 
      filtRes = rbind(Filtidy %>% filter(Enrichment == "Up-regulated") %>% head(n= 30),
                      Filtidy %>% filter(Enrichment == "Down-regulated") %>% head(n= 15))
      Y<-fgseaResTidy
      hmc=as.data.frame(Y)
      hmc=apply(Y,2,as.character) 
      write.csv(hmc,paste(name, X,Z,".csv"))
      names(Y)[names(Y)=="NES"]=paste(X)
      Annot.pathway<-Y[,c("pathway",paste(X))]
      Annot.pathway2<-merge(Annot.pathway, Annot.pathway2, by.x="pathway", by.y="pathway")
      
      plot=ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
        geom_segment( aes(reorder(pathway, NES), xend=pathway, y=0, yend=NES)) +
        geom_point( size=3, aes( fill = Enrichment), shape=21, stroke=1) +
        scale_fill_manual(values = c("Down-regulated" = "dodgerblue", "Up-regulated" = "firebrick") ) +
        scale_y_continuous(limits = NULL,expand = expansion(mult = 0.2, add = 0))+
        theme(axis.text= element_text(size=8, face="bold"))+theme_min()+
        coord_flip() +labs(x="Pathway", y="Normalized Enrichment Score",title=paste(name,"-",X,Z))
      ggsave(plot, file=paste(name,"-",objname, X,Z,".pdf"), scale=1,width=10)
      
    }
    ##make a heatmap
    rownames(Annot.pathway2)=Annot.pathway2$pathway
    Annot.pathway2=Annot.pathway2[,-1]
    positions=c("B6","S100a4KI")

    Annot.pathway2=Annot.pathway2[,positions]
    Annot.pathway2[is.na(Annot.pathway2)]=0
    Annot.pathway3=rbind(Annot.pathway2%>%filter_all(any_vars(.<= 2)),Annot.pathway2%>%filter_all(any_vars(.<= -2)))
    redblackblu=CustomPalette(low = "blue", high = "red", mid = "white", k = 1000)
    pheatmap(Annot.pathway3,scale = "none",cluster_rows = T,cluster_cols = F,border_color = "black",treeheight_row = 0,legend_breaks = c(-5,-2,0,2,5),
             fontsize =8, cellwidth = 15,cellheight =10,color= plasma(10),angle_col = 90,labels_col = c("B6",expression(paste("S100a4"^"-/-"))),
             main = paste("                                          ",name, Z,objname," (NES)"),filename = paste("GSEA",name,Z," (NES)-Heatmap.pdf"))
    
    
    pathnames=c(row.names(Annot.pathway2 %>% top_n(n = 2, wt = eval(parse(text=names(Annot.pathway2)[1])))),row.names(Annot.pathway2 %>% top_n(n = -2, wt = eval(parse(text=names(Annot.pathway2)[1])))))
    for (i in 2: length(names(Annot.pathway2))){
      pathnames1=c(row.names(Annot.pathway2 %>% top_n(n = 2, wt = eval(parse(text=names(Annot.pathway2)[i])))),row.names(Annot.pathway2 %>% top_n(n = -2, wt = eval(parse(text=names(Annot.pathway2)[i])))))
      pathnames=c(pathnames,pathnames1)
    }
    pathnames =pathnames[!duplicated(pathnames)]
    Annot.pathway4=Annot.pathway2[pathnames,]
    Annot.pathway4[is.na(Annot.pathway4)]=0
    redblackblu=CustomPalette(low = "blue", high = "red", mid = "white", k = 1000)
    pheatmap(Annot.pathway4,scale = "none",cluster_rows = T,cluster_cols = F,border_color = "black",treeheight_row = 0,legend_breaks = c(-5,-2,0,2,5),
             fontsize =8, cellwidth = 15,cellheight =10,color= plasma(10),angle_col = 90,labels_col = c("B6",expression(paste("S100a4"^"-/-"))),
             main = paste("                                          ",name,Z,objname," (NES)"),filename = paste("GSEA",name,Z," (NES) top and bottom 2each-Heatmap.pdf"))
    pathnames=c(row.names(Annot.pathway2 %>% top_n(n = 5, wt = eval(parse(text=names(Annot.pathway2)[1])))),row.names(Annot.pathway2 %>% top_n(n = -5, wt = eval(parse(text=names(Annot.pathway2)[1])))))
    for (i in 2: length(names(Annot.pathway2))){
      pathnames1=c(row.names(Annot.pathway2 %>% top_n(n = 5, wt = eval(parse(text=names(Annot.pathway2)[i])))),row.names(Annot.pathway2 %>% top_n(n = -5, wt = eval(parse(text=names(Annot.pathway2)[i])))))
      pathnames=c(pathnames,pathnames1)
    }
    pathnames =pathnames[!duplicated(pathnames)]
    Annot.pathway4=Annot.pathway2[pathnames,]
    Annot.pathway4[is.na(Annot.pathway4)]=0
    redblackblu=CustomPalette(low = "blue", high = "red", mid = "white", k = 1000)
    pheatmap(Annot.pathway4,scale = "none",cluster_rows = T,cluster_cols = F,border_color = "black",treeheight_row = 0,legend_breaks = c(-5,-2,0,2,5),
             fontsize =8, cellwidth = 15,cellheight =10,color= plasma(10),angle_col = 90,labels_col = c("B6",expression(paste("S100a4"^"-/-"))),
             main = paste("                                          ",name,Z,objname," (NES)"),filename = paste("GSEA",name,Z," (NES) top 5 and bottom 5each-Heatmap.pdf"))
    pathnames=c(row.names(Annot.pathway2 %>% top_n(n = 7, wt = eval(parse(text=names(Annot.pathway2)[1])))),row.names(Annot.pathway2 %>% top_n(n = -7, wt = eval(parse(text=names(Annot.pathway2)[1])))))
    for (i in 2: length(names(Annot.pathway2))){
      pathnames1=c(row.names(Annot.pathway2 %>% top_n(n = 7, wt = eval(parse(text=names(Annot.pathway2)[i])))),row.names(Annot.pathway2 %>% top_n(n = -7, wt = eval(parse(text=names(Annot.pathway2)[i])))))
      pathnames=c(pathnames,pathnames1)
    }
    pathnames =pathnames[!duplicated(pathnames)]
    Annot.pathway4=Annot.pathway2[pathnames,]
    Annot.pathway4[is.na(Annot.pathway4)]=0
    redblackblu=CustomPalette(low = "blue", high = "red", mid = "white", k = 1000)
    pheatmap(Annot.pathway4,scale = "none",cluster_rows = T,cluster_cols = F,border_color = "black",treeheight_row = 0,legend_breaks = c(-5,-2,0,2,5),
             fontsize =8, cellwidth = 15,cellheight =10,color= plasma(10),angle_col = 90,labels_col = c("B6",expression(paste("S100a4"^"-/-"))),
             main = paste("                                          ",name,Z,objname," (NES)"),filename = paste("GSEA",name,Z," (NES) top 7 and bottom 5each-Heatmap.pdf"))
    pathnames=c(row.names(Annot.pathway2 %>% top_n(n = 10, wt = eval(parse(text=names(Annot.pathway2)[1])))),row.names(Annot.pathway2 %>% top_n(n = -10, wt = eval(parse(text=names(Annot.pathway2)[1])))))
    for (i in 2: length(names(Annot.pathway2))){
      pathnames1=c(row.names(Annot.pathway2 %>% top_n(n = 10, wt = eval(parse(text=names(Annot.pathway2)[i])))),row.names(Annot.pathway2 %>% top_n(n = -10, wt = eval(parse(text=names(Annot.pathway2)[i])))))
      pathnames=c(pathnames,pathnames1)
    }
    pathnames =pathnames[!duplicated(pathnames)]
    Annot.pathway4=Annot.pathway2[pathnames,]
    Annot.pathway4[is.na(Annot.pathway4)]=0
    redblackblu=CustomPalette(low = "blue", high = "red", mid = "white", k = 1000)
    pheatmap(Annot.pathway4,scale = "none",cluster_rows = T,cluster_cols = F,border_color = "black",treeheight_row = 0,legend_breaks = c(-5,-2,0,2,5),
             fontsize =8, cellwidth = 15,cellheight =10,color= plasma(10),angle_col = 90,labels_col = c("B6",expression(paste("S100a4"^"-/-"))),
             main = paste("                                          ",name,Z,objname," (NES)"),filename = paste("GSEA",name,Z," (NES) top 10 and bottom 5each-Heatmap.pdf"))
    
  } 
  
  
}
plotEnrichment(fgsea_sets[["KEGG_FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS"]],ranks)
name=c(expression(paste("S100a4"^"-/-")))
bquote(ring("S100a4")^"-/-")
fea <- plotEnrichment(fgsea_sets[["BIOCARTA_TCYTOTOXIC_PATHWAY"]],ranks)+ 
  labs(subtitle="BIOCARTA_TCYTOTOXIC_PATHWAY (NES=3.4)",title = expression(paste("S100a4"^"-/-"," vs B6")))
XX2 <- plotEnrichment(fgsea_sets[["GSE14415_INDUCED_TREG_VS_FAILED_INDUCED_TREG_UP"]],ranks)+ 
  labs(subtitle="INDUCED_TREG_VS_FAILED_TREG_UP (NES=9.4)",title = expression(paste("S100a4"^"-/-"," vs B6")))

XX=XX+theme(text = element_text(size=7))
pathwaynames=read.csv("../pathways.csv",header = F,as.is = T)
Annot.pathway2=as.data.frame(read.csv("../pathways.csv",header = F))
names(Annot.pathway2)="pathway"
Clusters=levels(as.factor(MB@meta.data$Group))
library("tidyverse") 
objname="X2808"
for (i in 1:length(Clusters)){
  ZZ=Clusters[i]
  filess <- list.files(pattern=paste(ZZ,"AllClusters .csv"))
  filessL <- as.list(filess)
  names(filessL) <- str_replace(filess,paste(ZZ,"AllClusters .csv"),"")
  filessL <- lapply(filessL, read.csv, row.names=1)
  filessDF=filessL %>% reduce(full_join) 
  filessDFS=filessDF %>% filter(filessDF$pathway %in%Annot.pathway2$pathway, )
  assign(paste0(ZZ,"pathways"), filessDFS)
  names(filessDFS)[names(filessDFS)=="NES"]=paste0(ZZ)
  Annot.pathway<-filessDFS[,c("pathway",paste0(ZZ))]
  Annot.pathway2<-merge(Annot.pathway, Annot.pathway2, by.x="pathway", by.y="pathway")  
}
Annot.pathway2=Annot.pathway2[!duplicated(Annot.pathway2$pathway),]
rownames(Annot.pathway2)=Annot.pathway2$pathway
Annot.pathway2=Annot.pathway2[,-1]
positions=Clusters
Annot.pathway2=Annot.pathway2[,positions]
Annot.pathway2[is.na(Annot.pathway2)]=0
write.csv(Annot.pathway2,"selectpathways.csv")

#Annot.pathway2=read.csv("selectpathways.csv",row.names = 1)
#Annot.pathway3=rbind(Annot.pathway2%>%filter_all(any_vars(.<= 2)),Annot.pathway2%>%filter_all(any_vars(.<= -2)))
redblackblu=CustomPalette(low = "blue", high = "red", mid = "white", k = 1000)
pheatmap(Annot.pathway2,scale = "none",cluster_rows = T,cluster_cols = F,border_color = "black",
         fontsize =8, cellwidth = 15,cellheight =10,color= plasma(10),angle_col = 90,
         main = paste("select pathways All Clusters (NES)"),filename = paste(objname, "GSEA select pathways allclusters (NES)-Heatmap.pdf"))
n=4
Annot.pathway3=filter(Annot.pathway2, across(, ~ .x >= n|.x <= -n))

pheatmap(Annot.pathway3,scale = "none",cluster_rows = T,cluster_cols = F,border_color = "black",
         fontsize =8, cellwidth = 15,cellheight =10,color= plasma(10),angle_col = 90,
         main = paste("select pathways All Clusters (NES)"),filename = paste0(objname, " GSEA select pathways allclusters (NES>",n,")-Heatmap.pdf"))

n=2
pathnames=c(row.names(Annot.pathway2 %>% top_n(n = n, wt = eval(parse(text=names(Annot.pathway2)[1])))),row.names(Annot.pathway2 %>% top_n(n = -n, wt = eval(parse(text=names(Annot.pathway2)[1])))))
for (i in 2: length(names(Annot.pathway2))){
  pathnames1=c(row.names(Annot.pathway2 %>% top_n(n = n, wt = eval(parse(text=names(Annot.pathway2)[i])))),row.names(Annot.pathway2 %>% top_n(n = -n, wt = eval(parse(text=names(Annot.pathway2)[i])))))
  pathnames=c(pathnames,pathnames1)
}
pathnames =pathnames[!duplicated(pathnames)]
Annot.pathway4=Annot.pathway2[pathnames,]
Annot.pathway4[is.na(Annot.pathway4)]=0
redblackblu=CustomPalette(low = "blue", high = "red", mid = "white", k = 1000)
pheatmap(Annot.pathway4,scale = "none",cluster_rows = T,cluster_cols = F,border_color = "black",
         fontsize =8, cellwidth = 15,cellheight =10,color= plasma(10),angle_col = 90,
         main = paste("select pathways All Clusters (NES)"),filename = paste(objname, "GSEA select top and bottom", n , "pathways all clusters (NES)-Heatmap.pdf"))


