#The script for computer practice 2:
#Seurat: From clusters to marker genes

#Step 0. Preparatory. Loading libraries.
#Use _packages.R to install all the necessary packages prior loading the libraries

library(Seurat)
library(tidyverse)

#To make your data calculation reproducible between independent runs, set the seed for the random Number generator
set.seed(42)

#Lets reuse the Seurat object we have created in computer practice 2.
leaf.dataset <- readRDS('Data/leaf.dataset.rds')

#Or you can create the object again using the summary script below
path <- 'Data/GSM5097889/'
leaf.dataset<-Read10X(path, gene.column = 1)
leaf.dataset <- CreateSeuratObject(counts = leaf.dataset, project = "leaf")
leaf.dataset[["percent.mt"]] <- PercentageFeatureSet(leaf.dataset, pattern = "^ATM")
leaf.dataset[["percent.ct"]] <- PercentageFeatureSet(leaf.dataset, pattern = "^ATC")
leaf.dataset <- subset(leaf.dataset, subset = percent.mt <= 20 & percent.ct <= 20 & nCount_RNA >=500)
leaf.dataset <- SCTransform(leaf.dataset)
leaf.dataset <- RunPCA(leaf.dataset,verbose = FALSE)
leaf.dataset <- RunUMAP(leaf.dataset, dims = 1:50, verbose = FALSE)
leaf.dataset <- FindNeighbors(leaf.dataset, dims = 1:50, verbose = FALSE)
leaf.dataset <- FindClusters(leaf.dataset, resolution = 1, verbose = FALSE)

#The cell populations:
DimPlot(leaf.dataset, label = TRUE, pt.size = 1.5, label.size = 10) + NoLegend()

#Aim 1: Finding marker genes
leaf.markers <- FindAllMarkers(leaf.dataset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

write.csv(leaf.markers, file = "Data/leaf_markers.csv")
#leaf.markers<- read.csv(file = "Data/leaf_markers.csv", header = TRUE)  

#visualizing marker genes
top2<-leaf.markers %>%
        group_by(cluster) %>%
        slice_max(n = 2, order_by = avg_log2FC)

top2

DoHeatmap(leaf.dataset, features = top2$gene) + NoLegend()

top5<- leaf.markers %>% 
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)%>%
  filter(cluster == 0)
VlnPlot(leaf.dataset, features = top5$gene, pt.size = 0.2, ncol = 5)

DotPlot(leaf.dataset, features = unique(top2$gene), cols = c("blue", "red"), dot.scale = 8)+ coord_flip()

#Aim 2: Differentially expressed genes between two clusters

cluster2.markers <- FindMarkers(leaf.dataset, ident.1 = 2, ident.2 = 3, min.pct = 0.25)
head(cluster2.markers, n = 8)

top10 <- cluster2.markers %>% 
            slice_max(n = 10, order_by = avg_log2FC) %>% 
            row.names()
VlnPlot(leaf.dataset, features = top10, pt.size = 0.2, ncol = 5)

#FindMarkers(leaf.dataset, ident.1 = 2, ident.2 = c(1,6), min.pct = 0.25)
#FindMarkers(leaf.dataset, ident.1 = 2, ident.2 = c(1,6), min.pct = 0.25, test.use = "roc")


#Aim 3: Annotate clusters.
marker_genes<- c('AT1G11850', 'AT3G24140', 'AT2G05100', 'AT5G01530', 'AT4G10340', 'AT1G67090', 'AT5G38410', 'AT3G54890', 'AT1G29920', 'AT3G08940', 'AT3G27690', 'AT5G38430', 'AT1G76570', 'AT4G21750', 'AT2G26330', 'AT3G45640', 'AT2G42840', 'AT5G53210', 'AT5G60880', 'AT2G26330', 'AT1G80080', 'AT1G34245', 'AT3G26744', 'AT1G12860', 'AT2G42840', 'AT5G53210', 'AT5G60880', 'AT2G26330', 'AT5G62230', 'AT5G07180', 'AT1G80080', 'AT1G34245', 'AT1G04110', 'AT4G31805', 'AT3G26744', 'AT1G12860', 'AT3G45640', 'AT3G06120', 'AT5G53210', 'AT2G26330', 'AT5G62230', 'AT5G07180', 'AT1G80080', 'AT1G34245', 'AT2G20875', 'AT1G04110', 'AT4G31805', 'AT3G26744', 'AT1G12860', 'AT3G45640', 'AT3G06120', 'AT2G26330', 'AT5G62230', 'AT5G07180', 'AT2G20875', 'AT1G04110', 'AT4G31805', 'AT3G26744', 'AT1G12860', 'AT1G14350', 'AT3G45640', 'AT3G24140', 'AT3G06120', 'AT2G26330', 'AT5G62230', 'AT5G07180', 'AT2G20875', 'AT1G04110', 'AT3G26744', 'AT1G12860', 'AT1G14350', 'AT3G45640', 'AT3G24140', 'AT2G26330', 'AT5G62230', 'AT5G07180', 'AT1G04110', 'AT1G14350', 'AT5G46240', 'AT3G45640', 'AT3G17300', 'AT5G46240', 'AT2G26330', 'AT3G45640', 'AT4G12970', 'AT2G17950', 'AT2G27250', 'AT2G45190', 'AT2G27250', 'AT4G17710', 'AT1G55580', 'AT2G22850', 'AT3G28730', 'AT4G32880', 'AT5G16560')

DoHeatmap(leaf.dataset, features = marker_genes) + NoLegend() 

leaf.dataset <- RenameIdents(leaf.dataset,
                             `1`="Cluster A",
                             `2`='cluster B')
leaf.dataset$CellTypes<-Idents(leaf.dataset)

DoHeatmap(leaf.dataset, features = marker_genes) + NoLegend() 

#Your turn...

