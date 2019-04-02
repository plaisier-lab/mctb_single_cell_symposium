
##########################################################
## MCTB Single Cell Symposium: 6_Seurat_visualization.R ##
##  __  __  ___ _____ ___                               ##
## |  \/  |/ __|_   _| _ )                              ##
## | |\/| | (__  | | | _ \                              ##
## |_|  |_|\___| |_| |___/                              ##
##                                                      ##
## @Course:  MCTB Single Cell Symposium                 ##
## @Developed by: Plaisier Lab                          ##
##   (https://plaisierlab.engineering.asu.edu/)         ##
## @Developed by: Banovich Lab                          ##
##   (https://www.banovichlab.org/)                     ##
## @Authors:  Chris Plaisier and Nick Banovich          ##
## @License:  GNU GPLv3                                 ##
##                                                      ##
## If this program is used in your analysis please      ##
## mention who built it. Thanks. :-)                    ##
##########################################################

# Find all the differntially expressed genes between celltypes - using max.cells.per.ident to speed up for the workshop
lung_celltype_markers = FindAllMarkers(lung1, max.cells.per.ident = 100)

# This uses dplyr to pull the top 20 DE genes from each cell-type based on absolute value of the avg_logFC
top20 = lung_celltype_markers %>% group_by(cluster) %>% top_n(20, abs(avg_logFC))

# Make a heatmap of these genes 
DoHeatmap(lung1, genes.use = top20$gene)

# View the heatmap options
?DoHeatmap

# Remake heatmap with better parameters
DoHeatmap(lung1, genes.use = top20$gene, slim.col.label = T, title = "Heatmap of DE genes between cell types", group.label.rot = 2)

# How does the heatmap look if we don't use the scaled data
DoHeatmap(lung1, genes.use = top20$gene, slim.col.label = T, title = "Heatmap of DE genes between cell types", group.label.rot = 2, use.scaled = F)

# Make a violn plot of gene expression for SRGN

VlnPlot(lung1, "SRGN")

# Make a new object with the top 2 DE genes from each celltype
top2 = lung_celltype_markers %>% group_by(cluster) %>% top_n(2, abs(avg_logFC))

# Make a violin plot of the top2 genes
VlnPlot(lung1, top2$gene)

# Make FeaturePlots of the top2 genes
FeaturePlot(lung1, top2$gene)
