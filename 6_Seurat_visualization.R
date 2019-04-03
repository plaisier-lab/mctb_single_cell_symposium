
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

# Set final cell types as ident
lungs <- SetAllIdent(lungs, id = "celltype_final")

# Find all the differntially expressed genes between celltypes - using max.cells.per.ident to speed up for the workshop
lung_celltype_markers = FindAllMarkers(lungs, max.cells.per.ident = 100, min.cells.group = 50)

# This uses dplyr to pull the top 20 DE genes from each cell-type based on absolute value of the avg_logFC
top20 = lung_celltype_markers %>% group_by(cluster) %>% top_n(20, abs(avg_logFC))

# Make a heatmap of these genes 
DoHeatmap(lungs, genes.use = top20$gene)

# View the heatmap options
?DoHeatmap

# Remake heatmap with better parameters
DoHeatmap(lungs, genes.use = top20$gene, slim.col.label = T, title = "Heatmap of DE genes between cell types", group.label.rot = 2, cells.use = row.names(lungs@meta.data[lungs@meta.data$celltype_final %in% lung_celltype_markers$cluster,]))

# How does the heatmap look if we don't use the scaled data
DoHeatmap(lungs, genes.use = top20$gene, slim.col.label = T, title = "Heatmap of DE genes between cell types", group.label.rot = 2, cells.use = row.names(lungs@meta.data[lungs@meta.data$celltype_final %in% lung_celltype_markers$cluster,]), use.scaled = F)

# Make a violn plot of gene expression for LYZ

VlnPlot(lungs, "LYZ", x.lab.rot = 2)
VlnPlot(lungs, "LYZ", x.lab.rot = 2, use.scaled = T)
VlnPlot(lungs, "LYZ", x.lab.rot = 2, use.raw = T)


# Identify DE genes between lung1 (diseased) and lung3 (control)
lungs = SetAllIdent(lungs, id = "orig.ident")
lung_disease_markers = FindMarkers(lungs, ident.1 = "lung1", ident.2 = "lung3")

# VlnPlot of top most DE genes based on absolute logFC

VlnPlot(lungs, tail(rownames(lung_disease_markers[order(abs(lung_disease_markers$avg_logFC)),]), 10))

# Identify cell type specific DE genes for secretory, cilliated, and fibroblast cells.

clubcell_lungs = SubsetData(lungs, cells.use = row.names(lungs@meta.data[lungs@meta.data$celltype_final == "Club cells",]))
endothelial_lungs = SubsetData(lungs, cells.use = row.names(lungs@meta.data[lungs@meta.data$celltype_final == "Endothelial",]))
fibroblast_lungs = SubsetData(lungs, cells.use = row.names(lungs@meta.data[lungs@meta.data$celltype_final == "Fibroblast",]))

clubcell_disease_markers = FindMarkers(clubcell_lungs, ident.1 = "lung1", ident.2 = "lung3")
endothelial_disease_markers = FindMarkers(endothelial_lungs, ident.1 = "lung1", ident.2 = "lung3")
fibroblast_disease_markers = FindMarkers(fibroblast_lungs, ident.1 = "lung1", ident.2 = "lung3")

# View top 5 disease most UPREGULATED genes in diseased lungs for each cell type with split dot plot
top5_disease = c(head(rownames(clubcell_disease_markers[order(clubcell_disease_markers$avg_logFC),]), 5), head(rownames(endothelial_disease_markers[order(endothelial_disease_markers$avg_logFC),]), 5), head(rownames(fibroblast_disease_markers[order(fibroblast_disease_markers$avg_logFC),]), 5))

# Look for duplicated values
table(top5_disease)

# Remove the duplicates
top5_disease = unique(top5_disease)

# Visualization of DE genes across all cells with SplitDotPlot
SplitDotPlotGG(lungs, genes.plot = top5_disease, grouping.var = "orig.ident", group.by = "celltype_final", x.lab.rot = 2, plot.legend = T)

# Visualization of DE genes in fibroblasts with SplitDotPlot
SplitDotPlotGG(endothelial_lungs, genes.plot = top5_disease, grouping.var = "orig.ident", group.by = "celltype_final", x.lab.rot = 2, plot.legend = T)

# Look at featureplot of these genes
FeaturePlot(lungs, top5_disease)
TSNEPlot(lungs)

# Look at a featureplot of these genes just in secretory cells
FeaturePlot(secretory_lungs, top5_disease)
TSNEPlot(secretory_lungs)


