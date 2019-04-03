##########################################################
## MCTB Single Cell Symposium: 4_Seurat_clustering.R    ##
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

# Run TSNE using PCs 1:6
lung1 = RunTSNE(object = lung1, dims.use = 1:6, do.fast = TRUE, perplexity = 30)
TSNEPlot(lung1)

# Run TSNE using PCs 1:10
lung1 = RunTSNE(object = lung1, dims.use = 1:10, do.fast = TRUE, perplexity = 30)
TSNEPlot(lung1)

# Run TSNE using PCs 1:10 with perplexity of 10
lung1 = RunTSNE(object = lung1, dims.use = 1:10, do.fast = TRUE, perplexity = 10)
TSNEPlot(lung1)

# Run TSNE using PCs 1:10 with perplexity of 100
lung1 = RunTSNE(object = lung1, dims.use = 1:10, do.fast = TRUE, perplexity = 100)
TSNEPlot(lung1)


# Find clusters
lung1 = FindClusters(object = lung1, reduction.type = "pca", dims.use = 1:10, resolution = 0.3, print.output = 0, save.SNN = TRUE)
TSNEPlot(lung1)

# Find clusters at varying resolutions
lung1 = FindClusters(object = lung1, reduction.type = "pca", resolution = c(0.02, 0.1, .5, 1, 3), print.output = 0, reuse.SNN  = TRUE)
TSNEPlot(lung1)

# View meta data
head(lung1@meta.data)

# Plot TSNE for differnt resolutions
TSNEPlot(lung1, group.by = "res.0.02")
TSNEPlot(lung1, group.by = "res.0.1")
TSNEPlot(lung1, group.by = "res.0.3")
TSNEPlot(lung1, group.by = "res.0.5")
TSNEPlot(lung1, group.by = "res.1")

#Set chosen resolution as ident

lung1 = SetIdent(lung1, ident.use = lung1@meta.data$res.0.1)

# Labeling clusters using known marker genes
FeaturePlot(lung1, c("EPCAM", "PTPRC", "PECAM1", "PDGFRB")) # Epithelial, Immune, Endothelial, Fibroblast

# Immune cell types
FeaturePlot(lung1, c("CD14", "CD3E", "MS4A1")) # Myeloid, T-cell, B-cell

# Epithelial cell types
FeaturePlot(lung1, c("MUC5B", "FOXJ1")) # Secretory cells, cilliated cells

# Dot plot to see marker gene expression by cluser

DotPlot(lung1, c("PECAM1", "PDGFRB", "MUC5B", "FOXJ1", "CD14"))

# Create new meta data column with cell type annotations using plyr
current.cluster.ids = c(0, 1, 2, 3, 4, 5, 6)
new.cluster.ids = c("Fibroblast", "Myeloid", "Myeloid", "Myeloid",
                     "Secretory", "Cilliated", "Endothelial")
lung1@meta.data$celltype = plyr::mapvalues(x = lung1@ident, from = current.cluster.ids, to = new.cluster.ids)
lung1 = SetIdent(lung1, ident.use = lung1@meta.data$celltype)
TSNEPlot(lung1)
DotPlot(lung1, c("PECAM1", "PDGFRB", "MUC5B", "FOXJ1", "CD14"))

#Find genes with differntial expression between clusters
lung_celltype_markers = FindMarkers(lung1, ident.1 = "Cilliated", ident.2 = "Secretory")

head(lung_celltype_markers)

# Use negbinom test rather than wilcox test
lung_celltype_markers_negbinom = FindMarkers(lung1, ident.1 = "Cilliated", ident.2 = "Secretory", test.use = "negbinom")

head(lung_celltype_markers_negbinom)

# Look at the difference between two DE methods
dim(lung_celltype_markers)
dim(lung_celltype_markers_negbinom)

# What genes are in the default method but not the negbinom
dim(lung_celltype_markers[!(row.names(lung_celltype_markers) %in% row.names(lung_celltype_markers_negbinom)),])
head(lung_celltype_markers[!(row.names(lung_celltype_markers) %in% row.names(lung_celltype_markers_negbinom)),])

# What genes are in the negbinom but not the default
dim(lung_celltype_markers_negbinom[!(row.names(lung_celltype_markers_negbinom) %in% row.names(lung_celltype_markers)),])
