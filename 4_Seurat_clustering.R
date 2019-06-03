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

# Run UMAP using PCs 1:6
lung1 = RunUMAP(lung1, dims = 1:6)
DimPlot(lung1, reduction = 'umap')

# Run UMAP using PCs 1:10
lung1 = RunUMAP(lung1, dims = 1:10)
DimPlot(lung1, reduction = 'umap')

# Run TSNE using PCs 1:6
lung1 = RunTSNE(lung1, dims = 1:6, do.fast = TRUE, perplexity = 30)
DimPlot(lung1, reduction = 'tsne')

# Run TSNE using PCs 1:10
lung1 = RunTSNE(lung1, dims = 1:10, do.fast = TRUE, perplexity = 30)
DimPlot(lung1, reduction = 'tsne')

# Run TSNE using PCs 1:10 with perplexity of 10
lung1 = RunTSNE(lung1, dims = 1:10, do.fast = TRUE, perplexity = 10)
DimPlot(lung1, reduction = 'tsne')

# Run TSNE using PCs 1:10 with perplexity of 100
lung1 = RunTSNE(lung1, dims = 1:10, do.fast = TRUE, perplexity = 100)
DimPlot(lung1, reduction = 'tsne')


# Find clusters
lung1 = FindNeighbors(lung1, dims = 1:10)
lung1 = FindClusters(lung1, resolution = 0.3)
DimPlot(lung1, reduction='umap')
DimPlot(lung1, reduction='tsne')

# Find clusters at varying resolutions
lung1 = FindClusters(lung1, resolution = c(0.02, 0.1 , 0.5, 1, 3))
DimPlot(lung1, reduction='umap')
DimPlot(lung1, reduction='tsne')

# View meta data
head(lung1@meta.data)

# Plot UMAP for differnt resolutions
DimPlot(lung1, group.by = "RNA_snn_res.0.02", reduction = 'umap')
DimPlot(lung1, group.by = "RNA_snn_res.0.1", reduction = 'umap')
DimPlot(lung1, group.by = "RNA_snn_res.0.3", reduction = 'umap')
DimPlot(lung1, group.by = "RNA_snn_res.0.5", reduction = 'umap')
DimPlot(lung1, group.by = "RNA_snn_res.1", reduction = 'umap')

# Plot tSNE for differnt resolutions
DimPlot(lung1, group.by = "RNA_snn_res.0.02", reduction = 'tsne')
DimPlot(lung1, group.by = "RNA_snn_res.0.1", reduction = 'tsne')
DimPlot(lung1, group.by = "RNA_snn_res.0.3", reduction = 'tsne')
DimPlot(lung1, group.by = "RNA_snn_res.0.5", reduction = 'tsne')
DimPlot(lung1, group.by = "RNA_snn_res.1", reduction = 'tsne')

#Set chosen resolution as ident
colnames(lung1[[]])
levels(lung1)
Idents(lung1) = 'RNA_snn_res.0.1'
levels(lung1)

# Labeling clusters using known marker genes
FeaturePlot(lung1, features = c("EPCAM", "PTPRC", "PECAM1", "PDGFRB")) # Epithelial, Immune, Endothelial, Fibroblast

# Immune cell types
FeaturePlot(lung1, features = c("CD14", "CD3E", "MS4A1")) # Myeloid, T-cell, B-cell

# Epithelial cell types
FeaturePlot(lung1, features = c("MUC5B", "FOXJ1")) # Secretory cells, cilliated cells

# Dot plot to see marker gene expression by cluser
DotPlot(lung1, features = c("PECAM1", "PDGFRB", "MUC5B", "FOXJ1", "CD14"))

# Create new meta data column with cell type annotations using plyr
new.cluster.ids = c("Fibroblast", "Myeloid", "Myeloid", "Secretory", "Cilliated", "Endothelial", "Fibroblast")
names(new.cluster.ids) = levels(lung1)
lung1 = RenameIdents(lung1, new.cluster.ids)
DimPlot(lung1, reduction = 'umap', label = TRUE, pt.size = 0.5) + NoLegend()

# Dot plot to see marker gene expression by cluser
DotPlot(lung1, features = c("PECAM1", "PDGFRB", "MUC5B", "FOXJ1", "CD14"))

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

