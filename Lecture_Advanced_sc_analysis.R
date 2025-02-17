#Lecture on advanced scNA-Seq analysis


set.seed(100)
#BiocManager::install("slingshot")
#BiocManager::install("SingleCellExperiment",force = TRUE)
#BiocManager::install("scater")
library(scater)
library(slingshot)
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)

leaf.dataset <- readRDS('Data/leaf.dataset.rds')
DimPlot(leaf.dataset, label = TRUE, pt.size = 1.5, label.size = 10) + NoLegend()


#slingshot

# one curve


sce <- as.SingleCellExperiment(leaf.dataset)
sce.sling <- slingshot(sce, reducedDim='PCA')
head(sce.sling$slingPseudotime_1)


embedded <- embedCurves(sce.sling, "UMAP")
embedded <- slingCurves(embedded)[[1]]
embedded <- data.frame(embedded$s[embedded$ord,])
plotReducedDim(sce.sling, dimred = "UMAP", colour_by = "slingPseudotime_1") +
  geom_path(data = embedded, aes(x = umap_1, y = umap_2), size = 1.2, color = "black")

#with starting point

colData(sce)

sce.sling2 <- slingshot(sce, start.clus= "3", cluster=sce$seurat_clusters, reducedDim='PCA')
pseudo.paths <- slingPseudotime(sce.sling2)
head(pseudo.paths)

# Taking the rowMeans just gives us a single pseudo-time for all cells. Cells
# in segments that are shared across paths have similar pseudo-time values in 
# all paths anyway, so taking the rowMeans is not particularly controversial.
shared.pseudo <- rowMeans(pseudo.paths, na.rm=TRUE)

# Need to loop over the paths and add each one separately.
gg <- plotUMAP(sce.sling2, colour_by=I(shared.pseudo))
ggsave(filename = "Figures/UMAP_pseudotime_slingshot.png", plot = gg, width = 6, height = 5, dpi = 150)


embedded <- embedCurves(sce.sling2, "UMAP")
embedded <- slingCurves(embedded)
for (path in embedded) {
  embedded <- data.frame(path$s[path$ord,])
  gg <- gg + geom_path(data=embedded, aes(x=umap_1, y=umap_2), size=1.2)
}

gg
ggsave(filename = "Figures/UMAP_pseudotime_slingshot_traj.png", plot = gg, width = 6, height = 5, dpi = 150)

