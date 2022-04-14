library(easypackages)
MyPackages<-c("dplyr","ggplot2","ggpubr","gridExtra","viridis","egg",
              "grid","lattice","gtools","Biobase","RColorBrewer",
              "Seurat","cowplot","patchwork","stringr","ComplexHeatmap",
              "SingleCellExperiment", "ggmin", "clustifyr", 'harmony', 'devEMF',
              'presto', 'msigdbr', 'tibble', 'fgsea')
libraries(MyPackages)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("fgsea", dependencies = TRUE, INSTALL_opts = '--no-lock')
library(fgsea)

NourpalDirectory = "C:/Users/josem/Box/Jose Maldonado (Nourhan Abdelfattah)/ESC4-7 Single Cell Analysis/Nourpal"
NourpalDirectory = "C:/Users/TMHJAM38/Box/Jose Maldonado (Nourhan Abdelfattah)/ESC4-7 Single Cell Analysis/Nourpal"

setwd(NourpalDirectory)

devtools::load_all("Nourpal.Rproj")
library(Nourpal)



Directory="C:/Users/TMHJAM38/Box/Jose Maldonado (Nourhan Abdelfattah)/ESC4-7 Single Cell Analysis/"
Directory="C:/Users/josem/Box/Jose Maldonado (Nourhan Abdelfattah)/ESC4-7 Single Cell Analysis/"
RobjDirectory=paste0(Directory,"R_Objects/")
OutputDirectory=paste0(Directory,"Output/CancerCells/")
Pathwayoutput= paste0(OutputDirectory,"Pathway analysis by Sample in Cancer Clusters")
dir.create(Pathwayoutput)
setwd(Pathwayoutput)

SeuratObj <-(readRDS(paste0(RobjDirectory,"CancerCellsOnly.rds")))
SeuratObj <-(readRDS(paste0(RobjDirectory,"CancerCellsOnlyPostUpdates.rds")))


CellInfo=SeuratObj@meta.data
objname="Cancer Clusters"
Sample= levels(as.factor(CellInfo$Sample))


library (presto)
#2 Run Test
AnnotSamples.genes <- wilcoxauc(SeuratObj , 'Sample')
write.csv(AnnotSamples.genes,paste0("DE genes between ",objname, " Samples by presto.csv"))
#check if test was successful
dplyr::count(AnnotSamples.genes, group)
pathways = c("H","C7","C8","CP:KEGG","CP:BIOCARTA")
names= c("Hallmark","C7","C8","KEGG","Biocarta")
category=c("H","C7","C8","C2","C2")
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
  Samples=levels(as.factor(CellInfo$Sample))
  
  for (i in 1: length(Samples)){
    X=  Samples[i]
    print(X)
    Genes<-AnnotSamples.genes %>%
      dplyr::filter(group == Samples[[i]]) %>%
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
    filtRes = rbind(Filtidy %>% filter(Enrichment == "Up-regulated") %>% head(n= 10),
                    Filtidy %>% filter(Enrichment == "Down-regulated") %>% head(n= 10))
    Y<-fgseaResTidy
    hmc=as.data.frame(Y)
    hmc=apply(Y,2,as.character) 
    write.csv(hmc,paste(name, X,".csv"))
    names(Y)[names(Y)=="NES"]=paste(X)
    Annot.pathway<-Y[,c("pathway",paste(X))]
    Annot.pathway2<-merge(Annot.pathway, Annot.pathway2, by.x="pathway", by.y="pathway")
    
    plot=ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
      geom_segment( aes(reorder(pathway, NES), xend=pathway, y=0, yend=NES)) +
      geom_point( size=3, aes( fill = Enrichment), shape=21, stroke=1) +
      scale_fill_manual(values = c("Down-regulated" = "dodgerblue", "Up-regulated" = "firebrick") ) +
      scale_y_continuous(limits = NULL,expand = expansion(mult = 0.2, add = 0))+
      theme(axis.text= element_text(size=8, face="bold"))+
      coord_flip() +labs(x="Pathway", y="Normalized Enrichment Score",title=paste(name,"-", X))
    ggsave(plot, file=paste(name,"-",objname, X,".pdf"), scale=1,width=10)
    
  }
  ##make a heatmap
  rownames(Annot.pathway2)=Annot.pathway2$pathway
  Annot.pathway2=Annot.pathway2[,-1]
  positions=Sample
  Annot.pathway2=Annot.pathway2[,positions]
  Annot.pathway2[is.na(Annot.pathway2)]=0
  Annot.pathway3=rbind(Annot.pathway2%>%filter_all(any_vars(.<= 2)),Annot.pathway2%>%filter_all(any_vars(.<= -2)))
  redblackblu=CustomPalette(low = "blue", high = "red", mid = "white", k = 1000)
  pheatmap(Annot.pathway3,scale = "none",Sample_rows = T,Sample_cols = F,border_color = "black",
           fontsize =8, cellwidth = 15,cellheight =10,color= redblackblu,angle_col = 90,
           main = paste(name,objname,"Samples (NES)"),filename = paste("GSEA",name,"Samples (NES)-Heatmap.pdf"))
  colnames(Annot.pathway2)=str_replace_all(colnames(Annot.pathway2)," ",".")
  samplename=colnames(Annot.pathway2)
  pathnames=c(row.names(Annot.pathway2 %>% top_n(n = 2, wt = eval(parse(text=samplename[1])))),row.names(Annot.pathway2 %>% top_n(n = -2, wt = eval(parse(text=samplename[1])))))
  for (i in 2: length(samplename)){
    pathnames1=c(row.names(Annot.pathway2 %>% top_n(n = 2, wt = eval(parse(text=samplename[i])))),row.names(Annot.pathway2 %>% top_n(n = -2, wt = eval(parse(text=samplename[i])))))
    pathnames=c(pathnames,pathnames1)
  }
  pathnames =pathnames[!duplicated(pathnames)]
  Annot.pathway4=Annot.pathway2[pathnames,]
  Annot.pathway4[is.na(Annot.pathway4)]=0
  redblackblu=CustomPalette(low = "blue", high = "red", mid = "white", k = 1000)
  pheatmap(Annot.pathway4,scale = "none",Sample_rows = T,Sample_cols = F,border_color = "black",
           fontsize =8, cellwidth = 15,cellheight =10,color= redblackblu,angle_col = 90,
           main = paste(name,objname,"Samples (NES)"),filename = paste("GSEA",name,"Samples (NES) top and bottom 2each-Heatmap.pdf"))
} 


