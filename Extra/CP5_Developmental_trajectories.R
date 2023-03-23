#The script for computer practice 4:
#Developmental trajectories, developmental states

#Step 0. Preparatory. Loading libraries and the dataset.
#Use _packages.R to install all the necessary packages prior loading the libraries

library(Seurat)
library(tidyverse)
library(patchwork)
library(monocle3)
library(SeuratWrappers)
library(SeuratData)
set.seed(42)

leaf.dataset<-readRDS('Data/leaf.dataset.rds')

#Monocle3 offers clustering and analysis procedures as well as Seurat. You can right about it more here:
#https://cole-trapnell-lab.github.io/monocle3/docs/getting_started/
#but you can also apply it on seurat object http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/monocle3.html

cds <- as.cell_data_set(leaf.dataset)
cds <- cluster_cells(cds)
plot_cells(cds, show_trajectory_graph = TRUE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = TRUE)
wrap_plots(p1, p2)

#this does not work
#cds.sub <- subset(as.Seurat(cds, assay = "NULL"), monocle3_partitions == 1)
# cds.sub <- learn_graph(cds.sub)
# plot_cells(cds.sub, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)

cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)

cds <- order_cells(cds)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE)

#This also does not work
# cds_subset <- choose_cells(cds)
# subset_pr_test_res <- graph_test(cds_subset, neighbor_graph="principal_graph")
# pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.05))
