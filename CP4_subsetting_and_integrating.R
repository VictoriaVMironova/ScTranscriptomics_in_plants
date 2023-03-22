#The script for computer practice 3:
#Seurat: subsetting and integrating the data

#Step 0. Preparatory. Loading libraries.
#Use _packages.R to install all the necessary packages prior loading the libraries

library(Seurat)
library(tidyverse)
library(patchwork)

set.seed(42)

#Loading the leaf data set
leaf.dataset<-readRDS('Data/leaf.dataset.rds')
DimPlot(leaf.dataset, label = TRUE, pt.size = 1.5, label.size = 10) + NoLegend()

#Task 1: Exploring Seurat Object
str(leaf.dataset)
show(leaf.dataset)
show(leaf.dataset@assays)
show(leaf.dataset@assays$RNA)
show(leaf.dataset@assays$RNA@counts)
show(leaf.dataset@assays$RNA@counts@Dimnames[1])

head(leaf.dataset@meta.data)
head(leaf.dataset@meta.data$percent.mt)

#You can manipulate the data taken from Seurat Object externally
hist(leaf.dataset@meta.data$percent.mt)

MetaData <- as.data.frame(leaf.dataset@meta.data)

Status <- MetaData %>% 
  mutate(Status = case_when(
    nFeature_RNA <= 1500 ~ "old",
    (nFeature_RNA > 1500 & nFeature_RNA < 3000)  ~ "adult",
    nFeature_RNA >= 3000 ~ "young")) %>%
  select(Status)

#You can also add meta data into Seurat Objects
leaf.dataset <- AddMetaData(
  object = leaf.dataset,
  metadata = data.frame(Status),
  col.name = "Status"
)
show(leaf.dataset@meta.data)

#Lets explore what we can do with Seurat Object
utils::methods(class = 'Seurat')
?SeuratObject::Idents
Idents(leaf.dataset)
summary(Idents(leaf.dataset))

#Task 2: Adjusting the number of clusters
# What if we would like to see less clusters?
# here we decrease the resolution
leaf.dataset <- FindClusters(leaf.dataset, resolution = 0.2, verbose = FALSE)
DimPlot(leaf.dataset)
# here we decrease the resolution
leaf.dataset <- FindClusters(leaf.dataset, resolution = 1.5, verbose = FALSE)
DimPlot(leaf.dataset)
# The clusters taken with different resolution can be found in the metadata
head(leaf.dataset@meta.data)

#and you can change the active clusters in the Seurat object any time 
summary(leaf.dataset@active.ident)
leaf.dataset <- SetIdent(leaf.dataset, value = leaf.dataset@meta.data$SCT_snn_res.0.2)
summary(leaf.dataset@active.ident)

#Aim 3 Subset the data
?SeuratObject::subset

#Creating a Seurat object that consists of the cells from cluster 1 and 2
mesophyll <- subset(leaf.dataset, idents = c("1", "2"))
DimPlot(mesophyll, label = TRUE, pt.size = 1.5, label.size = 10) + NoLegend()

rest <- subset(leaf.dataset, idents = c(1, 2), invert = TRUE)

# What are the cell names of all cells  belonging to the cluster?
WhichCells(leaf.dataset, idents = "2")

#Saving the expression data for the chosen cells
c2.raw.data <- as.matrix(GetAssayData(leaf.dataset, slot = "counts")[, WhichCells(leaf.dataset, ident = "2")])



#Creating a Seurat object based on a gene expression
ATML1_pos <- subset(leaf.dataset, subset = AT4G21750 > 0.1)
DimPlot(ATML1_pos, label = TRUE, pt.size = 1.5, label.size = 10) + NoLegend()

#Creating a Seurat object based on genes expression
ATML1_PDF2 <- subset(leaf.dataset, subset = (AT4G21750 > 0.1 & AT4G04890 >0.1))
DimPlot(ATML1_PDF2, label = TRUE, pt.size = 1.5, label.size = 10) + NoLegend()

# Visualize co-expression of two features simultaneously
FeaturePlot(leaf.dataset, features = c("AT4G21750", "AT4G04890"), order = T, pt.size = 1.5, blend = TRUE)                     

#Aim 4 Integrate the data

#Step 1: Loading counts data
GSM5097889.counts <- Read10X("Data/GSM5097889", gene.column = 1)
GSM5097890.counts <- Read10X("Data/GSM5097890", gene.column = 1)
GSM5097891.counts <- Read10X("Data/GSM5097891", gene.column = 1)

#Step 2: creating Seurat objects
GSM5097889<- CreateSeuratObject(counts = GSM5097889.counts, project = "r1")
GSM5097890 <- CreateSeuratObject(counts = GSM5097890.counts, project = "r2")
GSM5097891 <- CreateSeuratObject(counts = GSM5097891.counts, project = "r3")

#Step 3: Merging Seurat objects
dim(GSM5097889)
dim(GSM5097890)
dim(GSM5097891)
leaf.dataset.merged<- merge(GSM5097889,c(GSM5097890,GSM5097891))
dim(leaf.dataset.merged)

#Step 4a: Standard pipeline: quality control
leaf.dataset.merged[["percent.mt"]] <- PercentageFeatureSet(leaf.dataset.merged, pattern = "^ATM")
leaf.dataset.merged[["percent.ct"]] <- PercentageFeatureSet(leaf.dataset.merged, pattern = "^ATC")
VlnPlot(object = leaf.dataset.merged, features = c("nFeature_RNA", "nCount_RNA","percent.ct","percent.mt"), ncol = 5)
leaf.dataset <- subset(leaf.dataset.merged, subset = percent.mt <= 20 & percent.ct <= 20 & nCount_RNA >=500)
dim(leaf.dataset.merged)

leaf.dataset.merged <- subset(leaf.dataset.merged, subset = percent.mt <= 20 & percent.ct <= 20 & nCount_RNA >=500)
VlnPlot(object = leaf.dataset.merged, features = c("nFeature_RNA", "nCount_RNA","percent.ct","percent.mt"), ncol = 5)
dim(leaf.dataset.merged)


#Step 4b: Standard processing procedures: normalization, clustering, saving

leaf.dataset.merged <- SCTransform(leaf.dataset.merged)
leaf.dataset.merged <- RunPCA(leaf.dataset.merged,verbose = FALSE)
leaf.dataset.merged <- RunUMAP(leaf.dataset.merged, dims = 1:50, verbose = FALSE)
leaf.dataset.merged <- FindNeighbors(leaf.dataset.merged, dims = 1:50, verbose = FALSE)
leaf.dataset.merged <- FindClusters(leaf.dataset.merged, resolution = 1,verbose = FALSE)
saveRDS(leaf.dataset.merged, 'Data/leaf.dataset.merged.rds')

#Step 5: Standard processing procedures: vizualization 
DimPlot(leaf.dataset.merged, group.by = "orig.ident")
DimPlot(leaf.dataset.merged, label = TRUE)

