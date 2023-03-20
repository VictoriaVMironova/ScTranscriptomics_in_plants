#The script for computer practice 1:
#Seurat: From counts to clusters

#Step 0. Preparatory. Loading libraries.
#Use _packages.R to install all the necessary packages prior loading the libraries

library(Seurat)
library(tidyverse)

#To make your data calculation reproducible between independent runs, set the seed for the random Number generator
set.seed(42)

#Step 1. Getting data into R
#GSM509788 data is stored in [Google Drive](https://drive.google.com/drive/folders/11hHUQ3nabTj230urqKnTXPCwoqWhZ8Zc?usp=sharing).
#GSM509789 is already present in the repository `Data/GSM509789/`.
#This is the version shown on presentation:
#path <- 'Data/GSM5097888/' 
#This is a light version:
path <- 'Data/GSM5097889/' 

leaf.counts<-Read10X(path, gene.column = 1)
dim(leaf.counts)

#Step 2. Making Seurat object
leaf.dataset <- CreateSeuratObject(counts = leaf.counts, project = "leaf")

head(leaf.dataset)
dim(leaf.dataset)

#Step 2.1. Testing the dataset quality
FeatureScatter(leaf.dataset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

leaf.dataset[["percent.mt"]] <- PercentageFeatureSet(leaf.dataset, pattern = "^ATM")
leaf.dataset[["percent.ct"]] <- PercentageFeatureSet(leaf.dataset, pattern = "^ATC")

VlnPlot(object = leaf.dataset, features = c("nFeature_RNA", "nCount_RNA","percent.ct","percent.mt"), ncol = 5)

leaf.dataset <- subset(leaf.dataset, subset = percent.mt <= 20 & percent.ct <= 20 & nCount_RNA >=500)
dim(leaf.dataset)

#Step 3. Data normalization 
leaf.dataset <- SCTransform(leaf.dataset)

#Step 4. Data clustering and visualization
leaf.dataset <- RunPCA(leaf.dataset,verbose = FALSE)
leaf.dataset <- RunUMAP(leaf.dataset, dims = 1:50, verbose = FALSE)

leaf.dataset <- FindNeighbors(leaf.dataset, dims = 1:50, verbose = FALSE)
leaf.dataset <- FindClusters(leaf.dataset, resolution = 1,verbose = FALSE)

DimPlot(leaf.dataset, label = TRUE)
DimPlot(leaf.dataset, label = TRUE, pt.size = 1.5, label.size = 10) + NoLegend()

#Step 5. Saving/reading Seurat object

saveRDS(leaf.dataset,'Data/leaf.dataset.rds')

leaf.dataset <- readRDS('Data/leaf.dataset.rds')

#Step 6. Gene expression visualization

FeaturePlot(leaf.dataset,features = 'AT3G24140')

FeaturePlot(leaf.dataset,features = 'AT3G24140',order = T, label = T, pt.size = 3) +
  scale_color_gradientn(colors = c('lightgray','yellow','red','darkred'))+
  ggtitle("FAMA")

FeaturePlot(leaf.dataset,features = c('AT3G24140', 'AT1G11850', 'AT2G05100', 'AT5G38410'), order = T, label = F,pt.size = 3)

marker_genes<- c('AT1G11850', 'AT3G24140', 'AT2G05100', 'AT5G01530', 'AT4G10340', 'AT1G67090', 'AT5G38410', 'AT3G54890', 'AT1G29920', 'AT3G08940', 'AT3G27690', 'AT5G38430', 'AT1G76570', 'AT4G21750', 'AT2G26330', 'AT3G45640', 'AT2G42840', 'AT5G53210', 'AT5G60880', 'AT2G26330', 'AT1G80080', 'AT1G34245', 'AT3G26744', 'AT1G12860', 'AT2G42840', 'AT5G53210', 'AT5G60880', 'AT2G26330', 'AT5G62230', 'AT5G07180', 'AT1G80080', 'AT1G34245', 'AT1G04110', 'AT4G31805', 'AT3G26744', 'AT1G12860', 'AT3G45640', 'AT3G06120', 'AT5G53210', 'AT2G26330', 'AT5G62230', 'AT5G07180', 'AT1G80080', 'AT1G34245', 'AT2G20875', 'AT1G04110', 'AT4G31805', 'AT3G26744', 'AT1G12860', 'AT3G45640', 'AT3G06120', 'AT2G26330', 'AT5G62230', 'AT5G07180', 'AT2G20875', 'AT1G04110', 'AT4G31805', 'AT3G26744', 'AT1G12860', 'AT1G14350', 'AT3G45640', 'AT3G24140', 'AT3G06120', 'AT2G26330', 'AT5G62230', 'AT5G07180', 'AT2G20875', 'AT1G04110', 'AT3G26744', 'AT1G12860', 'AT1G14350', 'AT3G45640', 'AT3G24140', 'AT2G26330', 'AT5G62230', 'AT5G07180', 'AT1G04110', 'AT1G14350', 'AT5G46240', 'AT3G45640', 'AT3G17300', 'AT5G46240', 'AT2G26330', 'AT3G45640', 'AT4G12970', 'AT2G17950', 'AT2G27250', 'AT2G45190', 'AT2G27250', 'AT4G17710', 'AT1G55580', 'AT2G22850', 'AT3G28730', 'AT4G32880', 'AT5G16560')

DoHeatmap(leaf.dataset, features = marker_genes) + NoLegend() 

DoHeatmap(leaf.dataset, features = marker_genes) + NoLegend() + scale_y_discrete(labels = c("AT3G24140" = "FAMA"))

