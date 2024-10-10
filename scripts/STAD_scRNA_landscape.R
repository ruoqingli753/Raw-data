library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(dplyr)
library(SeuratData)
library(tidydr)
library(dplyr)
library(stringr)
#BiocManager::install("hdf5r")    #
library(hdf5r)
library(Seurat) 
library(ggplot2)
library(SingleCellExperiment)
library(cols4all)
#install.packages('harmony')
library(harmony)
library(Matrix)

#01#####
# 
base_path <- '00_origin_datas/GEO/scRNA/GSE150290_RAW' # 
patient_folders <- list.dirs(path = base_path, full.names = TRUE, recursive = FALSE)
samples_name = list.files(path =base_path )
samples_name=gsub('.raw_gene_bc_matrices',"",samples_name)
files <- list.files(base_path, pattern = "matrix.mtx", full.names =TRUE, recursive = TRUE)
files <- gsub('/matrix.mtx',"",files)
scRNAlist <- list()

# 
for(i in 1:length(files ) ) {
  #i=2
  counts <- Read10X(data.dir=files[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, project=samples_name[i], min.cells=3, min.features=200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id=samples_name[i])  
  
  
  #
  if(T){ 
    scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-")
  }
  #
  if(T){
    scRNAlist[[i]][["percent.rb"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^RP[SL]")
  }
  #
  if(T){
    HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
    HB.genes <- CaseMatch(HB.genes, rownames(scRNAlist[[i]]))
    scRNAlist[[i]][["percent.HB"]]=PercentageFeatureSet(scRNAlist[[i]],features=HB.genes)
    
  }
}

### 
names(scRNAlist) <- samples_name
# 
scRNA <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])
Samples=stringr::str_split_fixed(rownames(scRNA@meta.data),'_',3)[,1]
Patient=stringr::str_split_fixed(rownames(scRNA@meta.data),'_',3)[,2]
Type=stringr::str_split_fixed(Patient,'-',2)[,2]
# 
condition <- grepl('Pat15', Patient) | grepl('Pat21', Patient) | grepl('Pat03', Patient)
# 
HP <- ifelse(condition, 'negative', 'positive')
scRNA <- AddMetaData(scRNA , Samples,col.name = "Samples")
scRNA <- AddMetaData(scRNA , Patient,col.name = "Patient")
scRNA <- AddMetaData(scRNA , Type,col.name = "Type")
scRNA <- AddMetaData(scRNA , HP,col.name = "HP")
scRNA@meta.data[1:5,7:9]
colnames(scRNA@meta.data)
table(scRNA@meta.data[,c("Patient","Type","HP")])

#c4a_gui()
#colors <- sample(c4a('tableau.20',12),size = 12,replace = F)
VlnPlot(scRNA, features=c("nFeature_RNA","nCount_RNA",'percent.mt'), pt.size=0.0)
scRNA = subset(scRNA, subset=nFeature_RNA>200&nFeature_RNA<5000&percent.mt<10)#&  #
ggsave("singlecell/01_landscape/after_nFeature_count_mt.pdf",height = 6,width = 15)

#02################################################################
scRNA <- NormalizeData(scRNA)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
scRNA <- ScaleData(scRNA, features = rownames(scRNA))

# s.genes=Seurat::cc.genes.updated.2019$s.genes
# g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes
# scRNA <- CellCycleScoring(scRNA, s.features=s.genes, g2m.features=g2m.genes, set.ident=TRUE)
# scRNA = SCTransform(scRNA, vars.to.regress=c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"), verbose=FALSE)
scRNA = RunPCA(scRNA, verbose=FALSE)

#03跑harmony####
scRNA = RunHarmony(scRNA, group.by.vars="orig.ident", max.iter.harmony=50, lambda=0.5)
# my_harmony_embeddings <- HarmonyMatrix(
#   data_mat  = as.matrix(scRNA@reductions$pca@cell.embeddings),
#   meta_data = scRNA@meta.data,
#   vars_use  = 'orig.ident',
#   do_pca = FALSE)
# 
# rownames(my_harmony_embeddings) <- rownames(scRNA@reductions$pca@cell.embeddings)
# scRNA[["harmony"]] <- CreateDimReducObject(embeddings = my_harmony_embeddings, key = "harmony_", assay = DefaultAssay(scRNA))
ElbowPlot(scRNA,ndims = 50)
#04####
scRNA <- FindNeighbors(scRNA, dims=1:30, reduction="harmony")
scRNA <- RunUMAP(scRNA, dims=1:30, reduction="harmony")
DimPlot(scRNA, reduction="umap", group.by="orig.ident", pt.size=0.05)+
  theme(legend.position="right", plot.title=element_blank())+scale_color_manual(values=colors)
ggsave
##4.1################################ 
c4a_gui()
colors <- c4a('brewer.paired',	12)
mydata <- FindClusters(scRNA, resolution=0.1)
UMAPPlot(mydata, pt.size=1, label=T, label.size=5)+NoLegend()
ggsave("singlecell/01_landscape/cluster_umap.pdf",height = 10,width = 10)
markers <- FindAllMarkers(mydata,group_by='seurat_clusters',only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25,slot = 'scale.data')
write.table(markers, "singlecell/01_landscape/All_celltype_DEG_cell_type.txt", col.names=T, row.names=T, quote=F, sep="\t")
table(mydata$seurat_clusters)

##4.2 ##################

VlnPlot(mydata, features=c("CD34"), pt.size=0, cols=colors,group.by ='seurat_clusters')+
  NoLegend()+theme(axis.title.x=element_blank())
#0:B cells 1 :BANK1,MS4A1,IRF8,CD79B, CD79A
#1 :Plasma cells:SLAMF7,TNFRSF17,GNG7,FCRL5,CD79A
#2:NK/T cells:CD3G,NKG7,CD3E,CCL5
#3:Epithelial cells 1 :KRT7,EPCAM,KRT18,KRT19
#4 Macrophage:MS4A6A,CD68,
#5:Fibroblast cells :ACTA2，COL1A1, COL1A2, DCN, THY1
#6:Endothelial cells:CD34,CDH5,VWF
#7:Mast cell:HPGDS,TPSAB1,CPA3
#8:B cells 2 :CD79A,MS4A1,CD79B,CD19
#9:Epithelial cells 2:EPCAM,KRT18,KRT19
#10:Red blood cell (erythrocyte)	: MS4A1,CD79A
#11:Parietal cell:'GIF','ATP4A','ATP4B'



#mydata = subset(mydata, seurat_clusters %in% c(0,1,2,3,4,5,6,7,8,9,10))
#mydata@meta.data$seurat_clusters = droplevels(mydata@meta.data$seurat_clusters, exclude=c(11))
##4.3 ################
cell_label = c("B cells 1", "Plasma cells", 'NK/T cells',"Epithelial cells 1", "Macrophage", "Fibroblast cells",
               "Endothelial cells", "Mast cell", "B cells 2","Epithelial cells 2",'Red blood cell (erythrocyte)	','Parietal cells')
names(cell_label) <- levels(mydata)
mydata <- RenameIdents(mydata, cell_label)
mydata[["cell_type"]] = Idents(mydata)

pdf("singlecell/01_landscape/UMAP_subcluster.pdf",width=6, height=6)
UMAPPlot(mydata, pt.size=1, label=T, label.size=5)+NoLegend()+scale_color_manual(values=colors)
dev.off()
#ggsave("singlecell/01_landscape/UMAP_subcluster.pdf",  width=6, height=6)
genes = c('BANK1','MS4A1','IRF8','CD79B', 'CD79A',
          'SLAMF7','TNFRSF17','GNG7','FCRL5',
          'CD3G','NKG7','CD3E','CCL5',
          'KRT7','EPCAM','KRT18','KRT19',
          'MS4A6A','CD68',
          'ACTA2','COL1A1', 'COL1A2', 'DCN', 'THY1',
          'CD34','CDH5','VWF',
          'HPGDS','TPSAB1','CPA3',
         'HBA1','HBA2','HBB',#'CD19',
         'GIF','ATP4A','ATP4B'
         #'CLIC6'
)
DotPlot(mydata, features=genes)+coord_flip()+
  scale_color_gradientn(colors=c("#CCCCCC",  "firebrick1"))+
  theme_minimal()+
  theme(text=element_text(family="Times"), axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"), axis.text.y=element_text(face="bold", size=12), axis.title.x=element_blank(), axis.title.y=element_blank())
ggsave("singlecell/01_landscape/dotplot_gene_marker.pdf",height = 8,width = 12)

######
table(mydata $HP)
mydata.positive <- subset(mydata,subset=HP=="positive")
table(mydata.positive$HP)
#
# mydata@meta.data<- mydata@meta.data %>%
#   mutate(type = str_extract(orig.ident, "(?<=_).+"))
table(mydata.positive$Type)
mydata.positive$Type <-ifelse(mydata.positive$Type=="A",'Adjacent normal',"tumor" )
mydata.positive$Type <- factor(mydata.positive$Type ,levels = c('Adjacent normal','tumor'))

# ##############################################################################
# top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# pdf("Z:/users/liy/20240218_lymphoma_scRNA/01_landscape/Heatmap.pdf",height = 6,width = 8)
# DoHeatmap(mydata, features = genes,slot ='counts',group.colors=colors,label = F) + 
# scale_fill_gradientn(colors=c("white", "lightgray", "red"))
# dev.off()
# ggsave("Z:/users/liy/20240218_lymphoma_scRNA/01_landscape/Heatmap.pdf",height = 6,width = 8)
##############################################################
#' set =  c('CD3E','NKG7','CD3G','CD247',
#'           'COL3A1','COL6A1','ACTA2',
#'           'LYZ' ,'CD86','MS4A4A','CYBB',
#'           'PECAM1','VWF','CDH5','FLT1',
#'           'TP63','KRT5','KRT16',
#'          'FLT4','LYVE1',
#'           #'GNLY',
#'           'TAGLN','MYH11','CNN1',
#'           'CD79A','MS4A1','FCER2',
#'           'PMEL','MLANA','DCT'
#'           )
#' FeaturePlot(mydata, features=set, cols=c("snow", "purple"), ncol=4)
#' ggsave('Z:/users/liy/20240218_lymphoma_scRNA/01_landscape/FeaturePlot.pdf',height = 4,width = 11)
##4.4#####################################################
gene = c('MS4A1', 'CD79A',
          'TNFRSF17',
          'NKG7','CD3E',
          'EPCAM',
          'MS4A6A','CD68',
          'ACTA2','COL1A1', 
          'VWF',
          'CPA3',
          'HBA1',#'CD19',
          'GIF','ATP4A'
          #'CLIC6'
)
VlnPlot(mydata, features=gene, pt.size=0, cols=colors, ncol=5)+NoLegend()+theme(axis.title.x=element_blank())
ggsave('singlecell/01_landscape/violin_Plot.pdf',height = 10,width =15)

##4.5 ######################################################## 
bar <-  as.data.frame(with(mydata.positive@meta.data, table(Type, cell_type)))
ggplot(data=bar, aes(x=Type, y=Freq, fill=cell_type))+ 
  geom_bar(stat="identity", position=position_fill())+
  scale_fill_manual(values=colors)+theme_classic()+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(face="bold", size=12, angle=0, hjust=0.5), axis.text.y=element_text(face="bold", size=12), legend.text=element_text(face="bold", size=12), legend.title=element_blank(), legend.position="right")
ggsave("singlecell/01_landscape/cell_type_number.pdf",  width=6, height=6)
##4.6##########################################################  
# Type_label = c("PBMMC", "HHD")
# bar$type = factor(bar$type, levels=Type_label)
bar = bar %>% group_by(Type) %>% mutate(percent=100*Freq/sum(Freq))

ggplot(data=bar, aes(x=cell_type, y=percent, fill=Type,label = sprintf("%.2f", percent)))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_manual(values=c("#028391","#FEAE6F"))+theme_classic()+
  ggtitle("Percent(%)")+
  geom_text(position = position_dodge(width = 0.9), vjust = -0.5, size = 4)+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"), axis.text.y=element_text(face="bold", size=12), legend.text=element_text(face="bold", size=12), legend.title=element_blank())
ggsave("singlecell/01_landscape/barplot_pair_number.pdf",  width=10, height=6)


#5#########
#("B cells 1", "Plasma cells", 'NK/T cells',"Epithelial cells 1", "Macrophage", "Fibroblast cells",
# "Endothelial cells", "Mast cell", "B cells 2","Epithelial cells 2",'Red blood cell (erythrocyte)	','Parietal cells')
##5.1####
#2297
#Epithelial cells 1
mydata.positive.epith1 <- subset(mydata.positive,cell_type%in% c("Epithelial cells 1"))
"'EMB','CPVL','CTLA4','FAM241A','CXCR4'"
VlnPlot(mydata.positive.epith1,features = c('EMB','CPVL','CTLA4','FAM241A','CXCR4'),group.by ='Type' ,pt.size=0.1,)+ggtitle("Epithelial cells 1")
ggsave("singlecell/01_landscape/prognostic_genes_violin_epith1.pdf",  width=8, height=6)
DotPlot(mydata.positive.epith1,features = c('EMB','CPVL','CTLA4','FAM241A','CXCR4'),group.by ='Type' ,)+coord_flip()+ggtitle("Epithelial cells 1")
ggsave("singlecell/01_landscape/prognostic_genes_DotPlot_epith1.pdf",  width=6, height=6)

#####Epithelial cells 2
mydata.positive.epith2 <- subset(mydata.positive,cell_type%in% c("Epithelial cells 2"))
VlnPlot(mydata.positive.epith2,features = c('EMB','CPVL','CTLA4','FAM241A','CXCR4'),group.by ='Type' ,pt.size=0.1,)+labs(title = "Epithelial cells 2")
ggsave("singlecell/01_landscape/prognostic_genes_violin_epith2.pdf",  width=8, height=6)
DotPlot(mydata.positive.epith2,features = c('EMB','CPVL','CTLA4','FAM241A','CXCR4'),group.by ='Type' ,)+coord_flip()+ggtitle("Epithelial cells 2")
ggsave("singlecell/01_landscape/prognostic_genes_DotPlot_epith2.pdf",  width=6, height=6)

#####'Parietal cells'
mydata.positive.Parietal <- subset(mydata.positive,cell_type%in% c('Parietal cells'))
VlnPlot(mydata.positive.Parietal,features = c('EMB','CPVL','CTLA4','FAM241A','CXCR4'),group.by ='Type' ,pt.size=0.1,)+labs(title = 'Parietal cells')
ggsave("singlecell/01_landscape/prognostic_genes_violin_Parietal_cells.pdf",  width=8, height=6)
DotPlot(mydata.positive.Parietal,features = c('EMB','CPVL','CTLA4','FAM241A','CXCR4'),group.by ='Type' ,)+coord_flip()+ggtitle("Parietal cells")
ggsave("singlecell/01_landscape/prognostic_genes_DotPlot_Parietal.pdf",  width=6, height=6)

##5.2 ####
HP_genes <- readRDS("singlecell/HP_genes.RDS")
HP_genes <- unique(HP_genes)
DotPlot(mydata.positive.epith1,features =HP_genes,group.by ='Type' ,)+ggtitle("Epithelial cells 1")+theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("singlecell/01_landscape/HP_genes_DotPlot_epith1.pdf",  width=12, height=4)
##
DotPlot(mydata.positive.epith2,features =HP_genes,group.by ='Type' ,)+ggtitle("Epithelial cells 2")+theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("singlecell/01_landscape/HP_genes_DotPlot_epith2.pdf",  width=12, height=4)
##
DotPlot(mydata.positive.Parietal,features =HP_genes,group.by ='Type' ,)+ggtitle("Parietal cells")+theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("singlecell/01_landscape/HP_genes_DotPlot_Parietal.pdf",  width=12, height=4)
######
# ###########
# #2417############################################
# "0.263*RGS1+-0.535*CTLA4+0.492*FAM241A+0.259*BASP1"
# DotPlot(mydata.positive.epith, features =c('RGS1','CTLA4','FAM241A','BASP1'),group.by ='Type')+coord_flip()
# #2358 （epi1 ） ######################
# VlnPlot(mydata.positive.epith,features = c('EMB','CPVL','CTLA4','FAM241A','TNFRSF9','BASP1'),group.by ='Type' ,pt.size=0.1,)
# ####
# DoHeatmap(mydata.positive.epith, features =c('EMB','CPVL','CTLA4','FAM241A','TNFRSF9','BASP1'),slot ='counts',group.by ='Type',group.colors=colors,label = F) +
# scale_fill_gradientn(colors=c( "lightgray", "red"))
# ###
# DotPlot(mydata.positive.epith, features =c('EMB','CPVL','CTLA4','FAM241A','TNFRSF9','BASP1'),group.by ='Type')+coord_flip()
# ggsave("singlecell/01_landscape/DoHeatmap.pdf",height = 6,width = 5)
# #201################
# VlnPlot(mydata.positive.epith,features = c('CTLA4'   ,   'CIITA',    'FAM241A',      'BASP1',     'SLC2A3'),group.by ='Type' ,pt.size=0.1,)
# DotPlot(mydata.positive.epith,features = c('CTLA4'   ,   'CIITA',    'FAM241A',      'BASP1',     'SLC2A3'),group.by ='Type' ,)+coord_flip()

# #254#################
# "'CTLA4','FAM241A','BASP1','PKD2L1','SLC2A3'"
# VlnPlot(mydata.positive.epith,features = c("CTLA4","FAM241A","BASP1","PKD2L1","SLC2A3"),group.by ='Type' ,pt.size=0.1,)
# DotPlot(mydata.positive.epith,features = c("CTLA4","FAM241A","BASP1","PKD2L1","SLC2A3"),group.by ='Type' ,)+coord_flip()
# #380    (B CELL1: 'CXCR4'结果比较好)###
# "'DENND1C','STAT5A','FAM241A','CXCR4','BASP1','ZC3H12D'"
# VlnPlot(mydata.positive.epith,features = c('DENND1C','STAT5A','FAM241A','CXCR4','BASP1','ZC3H12D'),group.by ='Type' ,pt.size=0.1,)
# DotPlot(mydata.positive.epith,features = c('DENND1C','STAT5A','FAM241A','CXCR4','BASP1','ZC3H12D'),group.by ='Type' ,)+coord_flip()
# #425####
# "'CPVL','FAM241A','PTGDR','BASP1','PDCD1'"
# VlnPlot(mydata.positive.epith,features = c('CPVL','FAM241A','PTGDR','BASP1','PDCD1'),group.by ='Type' ,pt.size=0.1,)
# DotPlot(mydata.positive.epith,features = c('CPVL','FAM241A','PTGDR','BASP1','PDCD1'),group.by ='Type' ,)+coord_flip()
# #476  ()####
# "'TNFRSF1B','CARD11','FAM241A','CXCR4','PDCD1'"
# VlnPlot(mydata.positive.epith,features = c('TNFRSF1B','CARD11','FAM241A','CXCR4','PDCD1'),group.by ='Type' ,pt.size=0.5,)
# DotPlot(mydata.positive.epith,features = c('TNFRSF1B','CARD11','FAM241A','CXCR4','PDCD1'),group.by ='Type' ,)+coord_flip()
# #625 (B1：RGS1)####
# "'RGS1','CARD11','CTLA4','FAM241A','BASP1','PDCD1'"
# VlnPlot(mydata.positive.epith,features = c('CARD11','CTLA4','FAM241A','BASP1','PDCD1'),group.by ='Type' ,pt.size=0.1,)
# DotPlot(mydata.positive.epith,features = c('RGS1','CARD11','CTLA4','FAM241A','BASP1','PDCD1'),group.by ='Type' ,)+coord_flip()
# #1021#####
# "'CIITA','FAM241A','CXCR4','BASP1','TNFRSF18'"
# VlnPlot(mydata.positive.epith,features = c('CIITA','FAM241A','CXCR4','BASP1','TNFRSF18'),group.by ='Type' ,pt.size=0.1,)
# DotPlot(mydata.positive.epith,features = c('CIITA','FAM241A','CXCR4','BASP1','TNFRSF18'),group.by ='Type' ,)+coord_flip()
# #1026####
# "'FAM241A','CXCR4','BASP1','PDCD1','SLC2A3'"
# VlnPlot(mydata.positive.epith,features = c('FAM241A','CXCR4','BASP1','PDCD1','SLC2A3'),group.by ='Type' ,pt.size=0.1,)
# DotPlot(mydata.positive.epith,features = c('FAM241A','CXCR4','BASP1','PDCD1','SLC2A3'),group.by ='Type' ,)+coord_flip()
# #1594####
# "'FAM241A','CXCR4','BASP1','PDCD1'"
# VlnPlot(mydata.positive.epith,features = c('FAM241A','BASP1','PDCD1','CXCR4'),group.by ='Type' ,pt.size=0.1,)
# DotPlot(mydata.positive.epith,features = c('FAM241A','CXCR4','BASP1','PDCD1'),group.by ='Type' ,)+coord_flip()
# #1880####
# "'CTLA4','FAM241A','CXCR4','CHI3L2'"
# VlnPlot(mydata.positive.epith,features = c('CTLA4','FAM241A','CXCR4','CHI3L2'),group.by ='Type' ,pt.size=0.1,)
# DotPlot(mydata.positive.epith,features = c('CTLA4','FAM241A','CXCR4','CHI3L2'),group.by ='Type' ,)+coord_flip()
# #2074####
# "'EMB','CTLA4','FAM241A','CXCR4','LY96'"
# VlnPlot(mydata.positive.epith,features = c('EMB','CTLA4','FAM241A','LY96'),group.by ='Type' ,pt.size=0.1,)
# DotPlot(mydata.positive.epith,features = c('EMB','CTLA4','FAM241A','CXCR4','LY96'),group.by ='Type' ,)+coord_flip()
# #2099 (B1 CXCR4 可，其他变化不大)####
# "'CTLA4','CIITA','FAM241A','CXCR4'"
# VlnPlot(mydata.positive.epith,features = c('CTLA4','CIITA','FAM241A','CXCR4'),group.by ='Type' ,pt.size=0.1,)
# DotPlot(mydata.positive.epith,features = c('CTLA4','CIITA','FAM241A','CXCR4'),group.by ='Type' ,)+coord_flip()
# 
# #2322 （）####
# "'CPVL','CTLA4','FAM241A','SLC2A3'"
# VlnPlot(mydata.positive.epith,features = c('CPVL','CTLA4','FAM241A','SLC2A3'),group.by ='Type' ,pt.size=0.1,)
# DotPlot(mydata.positive.epith,features = c('CPVL','CTLA4','FAM241A','SLC2A3'),group.by ='Type' ,)+coord_flip()
# #2358####
# "'EMB','CPVL','CTLA4','FAM241A','TNFRSF9','BASP1'"
# VlnPlot(mydata.positive.epith,features = c('EMB','CPVL','CTLA4','FAM241A','TNFRSF9','BASP1'),group.by ='Type' ,pt.size=0.1,)
# DotPlot(mydata.positive.epith,features = c('EMB','CPVL','CTLA4','FAM241A','TNFRSF9','BASP1'),group.by ='Type' ,)+coord_flip()
# #2575####
# "'HCAR2','CTLA4','FAM241A','CXCR4'"
# VlnPlot(mydata.positive.epith,features = c('HCAR2','CTLA4','FAM241A','CXCR4'),group.by ='Type' ,pt.size=0.1,)
# DotPlot(mydata.positive.epith,features = c('HCAR2','CTLA4','FAM241A','CXCR4'),group.by ='Type' ,)+coord_flip()
# #2582 不好####
# "'IKBIP','CTLA4','FAM241A','CXCR4'"
# VlnPlot(mydata.positive.epith,features = c('IKBIP','CTLA4','FAM241A','CXCR4'),group.by ='Type' ,pt.size=0.1,)
# DotPlot(mydata.positive.epith,features = c('IKBIP','CTLA4','FAM241A','CXCR4'),group.by ='Type' ,)+coord_flip()
# #2774 （B1 CXCR4)####
# "'CPVL','CARD11','CXCR4','BASP1','PDCD1'"
# VlnPlot(mydata.positive.epith,features = c('CPVL','CARD11','BASP1','PDCD1'),group.by ='Type' ,pt.size=0.1,)
# DotPlot(mydata.positive.epith,features = c('CPVL','CARD11','CXCR4','BASP1','PDCD1'),group.by ='Type' ,)+coord_flip()
# #2877 （B1 CXCR4）####
# "'CTLA4','G0S2','FAM241A','CXCR4','BASP1'"
# VlnPlot(mydata.positive.epith,features = c('CTLA4','G0S2','FAM241A','BASP1'),group.by ='Type' ,pt.size=0.1,)
# DotPlot(mydata.positive.epith,features = c('CTLA4','G0S2','FAM241A','CXCR4','BASP1'),group.by ='Type' ,)+coord_flip()
