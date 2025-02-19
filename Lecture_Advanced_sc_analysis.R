#Lecture on advanced scNA-Seq analysis

#good thing to read
#https://broadinstitute.github.io/2019_scWorkshop/functional-pseudotime-analysis.html


set.seed(100)
#for slingshot:
#BiocManager::install("slingshot")
#BiocManager::install("SingleCellExperiment",force = TRUE)
#BiocManager::install("scater")

#to install monocle3
#BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
#                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
#                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
#                       'terra', 'ggrastr'))
#devtools::install_github('cole-trapnell-lab/monocle3')

library(tradeSeq)
library(scater)
library(slingshot)
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(monocle3)
library(SeuratWrappers)
library(CytoTRACE)

leaf.dataset <- readRDS('Data/leaf.dataset.rds')
DimPlot(leaf.dataset, label = TRUE, pt.size = 1.5, label.size = 10) + NoLegend()

# Monocle3 ----

#Monocle3 offers clustering and analysis procedures as well as Seurat. You can right about it more here:
#https://cole-trapnell-lab.github.io/monocle3/docs/getting_started/
#but you can also apply it on seurat object http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/monocle3.html

cds <- as.cell_data_set(leaf.dataset)

cds <- cluster_cells(cds)
plot_cells(cds, show_trajectory_graph = TRUE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = TRUE)
wrap_plots(p1, p2)


cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)

ggsave("Figures/monocle3_traj.png", plot = last_plot(), width = 2.5, height = 2.5, dpi = 150)
#ggsave("Figures/monocle3_traj_2024.png", plot = last_plot(), width = 2.5, height = 2.5, dpi = 150)


cds <- order_cells(cds)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE)
ggsave("Figures/monocle3_traj_pseudo.png", plot = last_plot(), width = 4, height = 2.5, dpi = 150)

# To explore DEG
# Fit a GAM model to test for DE genes along pseudotime
de_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=4)

#Select significant genes (adjust p-value cutoff as needed)
sig_genes <- rownames(subset(de_genes, q_value < 0.0001))
head(sig_genes)

#Plot gene expression of top significant genes over pseudotime

rowData(cds)$gene_short_name <- rownames(cds)

plot_cells(cds, genes=c("AT1G01050", "AT1G01120", "AT1G04110", "AT1G08380"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)
ggsave("Figures/monocle3_sig_genes.png", plot = last_plot(), width = 4, height = 4, dpi = 150)
#DEG_genes <- c('AT1G12480', 'AT1G11850', 'AT2G05100', 'AT5G38410')
DEG_genes <- c("AT1G01050", "AT1G01120", "AT1G04110", "AT1G08380")
DEG_lineage_cds <- cds[rowData(cds)$gene_short_name %in% DEG_genes]
DEG_lineage_cds <- order_cells(DEG_lineage_cds)

plot_genes_in_pseudotime(DEG_lineage_cds, min_expr=0.00001)

ggsave("Figures/monocle3_sig_genes_profiles.png", plot = last_plot(), width = 4, height = 4, dpi = 150)


# # Ensure your CDS object is properly preprocessed
# cds <- preprocess_cds(cds, num_dim = 50)  # Perform PCA
# cds <- reduce_dimension(cds)  # Reduce dimensions (UMAP / PCA / TSNE)
# cds <- cluster_cells(cds, resolution=1e-5)
# plot_cells(cds, color_cells_by = "cluster")


gene_module_df<-find_gene_modules(cds[sig_genes,], resolution=c(10^seq(-6,-1)))

cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group = clusters(cds))
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2", filename = "Figures/monocle3_modules.png")

#slingshot ----

# one curve
sce <- as.SingleCellExperiment(leaf.dataset)
sce.sling <- slingshot(sce, reducedDim='PCA')
head(sce.sling$slingPseudotime_1)


embedded <- embedCurves(sce.sling, "UMAP")
embedded <- slingCurves(embedded)[[1]]
embedded <- data.frame(embedded$s[embedded$ord,])
plotReducedDim(sce.sling, dimred = "UMAP", colour_by = "slingPseudotime_1") +
  geom_path(data = embedded, aes(x = umap_1, y = umap_2), size = 1.2, color = "black")
ggsave(filename = "Figures/UMAP_pseudotime_slingshot_traj.png", plot = last_plot(), width = 6, height = 5, dpi = 150)

#with starting point

colData(sce)

sce.sling2 <- slingshot(sce, start.clus= "9", cluster=sce$seurat_clusters, reducedDim='PCA')
pseudo.paths <- slingPseudotime(sce.sling2)
head(pseudo.paths)

# Taking the rowMeans just gives us a single pseudo-time for all cells. Cells
# in segments that are shared across paths have similar pseudo-time values in 
# all paths anyway, so taking the rowMeans is not particularly controversial.
shared.pseudo <- rowMeans(pseudo.paths, na.rm=TRUE)

# Need to loop over the paths and add each one separately.
gg <- plotUMAP(sce.sling2, colour_by=I(shared.pseudo))
gg
ggsave(filename = "Figures/UMAP_pseudotime_slingshot_9.png", plot = gg, width = 6, height = 5, dpi = 150)


embedded <- embedCurves(sce.sling2, "UMAP")
embedded <- slingCurves(embedded)
for (path in embedded) {
  embedded <- data.frame(path$s[path$ord,])
  gg <- gg + geom_path(data=embedded, aes(x=umap_1, y=umap_2), size=0.5)
}

gg
ggsave(filename = "Figures/UMAP_pseudotime_slingshot_traj_9.png", plot = gg, width = 6, height = 5, dpi = 150)

#To explore the DEG
# Extract pseudotime values
pseudotime <- slingPseudotime(sce)

# Step 3: Fit Generalized Additive Models (GAM) using tradeSeq
# Prepare for differential expression analysis
# Define lineage structure from Slingshot
lineages <- slingLineages(sce)

# Fit negative binomial GAM for each gene across pseudotime
counts <- counts(sce)  # Extract raw counts for modeling
gam_model <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = slingCurveWeights(sce))

# Step 4: Test for differential expression
de_res <- associationTest(gam_model)  # Tests for significant changes along pseudotime
de_genes <- rownames(subset(de_res, pval < 0.05))  # Select significant genes

# Step 5: Visualize gene expression trends
plotSmoothers(gam_model, counts, gene = de_genes[1])  # Example for one gene


# Cytotrace ----

mat <- GetAssayData(object = leaf.dataset, assay = "RNA", slot = "counts")

as_matrix <- function(mat){
  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x
  
  for (i in seq_along(val)){
    tmp[row_pos[i],col_pos[i]] <- val[i]
  }
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}


mat2 <- as_matrix(mat)


cyto_ra <- CytoTRACE(mat2,ncores = 1,enableFast = T)

cyto_ra$CytoTRACE <- 1-cyto_ra$CytoTRACE

leaf.dataset <- AddMetaData(
  object = leaf.dataset,
  metadata = cyto_ra$CytoTRACE,
  col.name = "differentiation_level"
)

FeaturePlot(leaf.dataset, features = c('differentiation_level'))+
  scale_color_gradientn(colors = c('#018571','#dfc27d','#a6611a'))+theme(text=element_text(size=10))+NoAxes()

ggsave(filename = "Figures/Cytotrace_Diff_2024.png", plot = last_plot(), width = 4, height = 3.5, dpi = 150)



