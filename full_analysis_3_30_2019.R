##########################################################
## MCTB Single Cell Symposium: 3_Seurat_basics.R        ##
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

# Import Seurat single cell analysis software
library(Seurat)

# Set your working directory
setwd('C:/Users/cplaisie/Dropbox (ASU)/mctb_single_cell_symposium_data')

# Load up data
lung1_counts = Read10X(data.dir = "Lung1/outs/filtered_gene_bc_matrices/GRCh38/")

# Take a look at the column and row labels
dim(lung1_counts)
colnames(lung1_counts)[1:10]
rownames(lung1_counts)[1:10]

# Create Seurat object
lung1 = CreateSeuratObject(raw.data = lung1_counts)

# Structure of a Seurat object
str(lung1)

# Raw data storage
lung1@raw.data[1:20,1:20]

# Data storage
lung1@data[1:20,1:20] #Currently set as raw.data

# Select mitochondrial genes - mitochondria (MT) have a double membrane, their own genome, and transcriptases.
# Ratio of MT genes to general cell genes is a great way to assess cell integrity as MT genes will have a harder
# time escaping the cell if the membrane is compromised.
mito.genes = grep(pattern = "^MT-", x = rownames(x = lung1@data), value = TRUE)
length(mito.genes)
mito.genes
percent.mito = colSums(lung1@data[mito.genes, ]) / colSums(lung1@data)
lung1 = AddMetaData(object = lung1, metadata = percent.mito, col.name = "percent.mito")

# Violin plot of three parameters
VlnPlot(object = lung1, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

# Plot of Percent.mito vs. nGene or nUMI
# Best guess data was filtered with nUMI >= 2000
par(mfrow = c(1,2))
plot(lung1@meta.data$percent.mito, lung1@meta.data$nUMI, col = rgb(0,0,1,0.5), ylim = c(0,38000), xlim = c(0,0.35))
abline(v = 0.1, col = 'red', lwd = 2)
abline(h = 2000, col = 'red', lwd = 2)
plot(lung1@meta.data$percent.mito, lung1@meta.data$nGene, col=rgb(0,0,1,0.5), ylim = c(0,6000), xlim = c(0,0.35))
abline(v = 0.1, col = 'red', lwd = 2)

# Filtering cells
#  - nUMI >= 2000 UMI
#  - percent.mito <= 0.1
lung1 = FilterCells(object = lung1, subset.names = c("nUMI", "percent.mito"), low.thresholds = c(2000, -Inf), high.thresholds = c(Inf, 0.1))
dim(lung1@data)

# Normalize data
lung1 = NormalizeData(object = lung1)
lung1@data[1:20,1:20] # Now normalized data

# Find the most variable genes in the dataset
lung1 = FindVariableGenes(object = lung1)
length(lung1@var.genes)

# Scale the data and regress out percent of mitochondrial and nUMI effects
lung1 = ScaleData(object = lung1, vars.to.regress = c('percent.mito','nUMI'), genes.use = lung1@var.genes, model.use = 'negbinom')

# Run PCA analysis
lung1 = RunPCA(object = lung1, pc.genes = lung1@var.genes, pcs.compute = 40, pcs.print = 1:30, maxit = 500, weight.by.var = FALSE)

# Plot PCAs
PCHeatmap(object = lung1, pc.use = 1:12, cells.use = 300, do.balanced = TRUE, label.columns = FALSE)

# Run TSNE
lung1 = RunTSNE(object = lung1, dims.use = 1:6, do.fast = TRUE, perplexity = 30)
TSNEPlot(lung1)

# Find clusters
lung1 = FindClusters(object = lung1, reduction.type = "pca", dims.use = 1:6, resolution = 0.3, print.output = 0, save.SNN = TRUE)
TSNEPlot(lung1)

# Find all marker genes
lung1.markers = FindAllMarkers(object = lung1, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
lung1.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
write.csv(lung1.markers,'lung1.markers.csv')

# Build classifier for clusters - 
lung1.rf_mg = BuildRFClassifier(object = lung1, training.genes = (lung1.markers %>% group_by(cluster) %>% top_n(100, avg_logFC))$gene, training.classes = lung1@ident)
lung1.rf_mg

lung1.rf = BuildRFClassifier(object = lung1, training.genes = lung1@var.genes, training.classes = lung1@ident)
lung1.rf


########################
## Load up lung2 data ##
########################
lung2_counts = Read10X(data.dir = "Lung2/outs/filtered_gene_bc_matrices/GRCh38/")
lung2 = CreateSeuratObject(raw.data = lung2_counts)
mito.genes = grep(pattern = "^MT-", x = rownames(x = lung2@data), value = TRUE)
percent.mito = colSums(lung2@data[mito.genes, ]) / colSums(lung2@data)
lung2 = AddMetaData(object = lung2, metadata = percent.mito, col.name = "percent.mito")

# Violin plot of three parameters
VlnPlot(object = lung2, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

# Plot of Percent.mito vs. nGene or nUMI
# Best guess data was filtered with nUMI >= 200
par(mfrow = c(1,2))
plot(lung2@meta.data$percent.mito, lung2@meta.data$nUMI, col = rgb(0,0,1,0.5), ylim = c(0,38000), xlim = c(0,0.35))
abline(v = 0.2, col = 'red', lwd = 2)
abline(h = 200, col = 'red', lwd = 2)
plot(lung2@meta.data$percent.mito, lung2@meta.data$nGene, col=rgb(0,0,1,0.5), ylim = c(0,6000), xlim = c(0,0.35))
abline(v = 0.2, col = 'red', lwd = 2)
dim(lung2@data)
lung2 = FilterCells(object = lung2, subset.names = c("nUMI", "percent.mito"), low.thresholds = c(200, -Inf), high.thresholds = c(Inf, 0.2))
dim(lung2@data)
lung2 = NormalizeData(object = lung2)
lung2 = FindVariableGenes(object = lung2)
lung2 = ScaleData(object = lung2, vars.to.regress = c('percent.mito','nUMI'), genes.use = lung2@var.genes, model.use = 'negbinom')
lung2 = RunPCA(object = lung2, pc.genes = lung2@var.genes, pcs.compute = 40, pcs.print = 1:30, maxit = 500, weight.by.var = FALSE)

PCHeatmap(object = lung2, pc.use = 1:12, cells.use = 300, do.balanced = TRUE, label.columns = FALSE)

# Apply classifier
#lung1 = FindClusters(object = lung1, reduction.type = "pca", dims.use = 1:6, resolution = 0.3, print.output = 0, save.SNN = TRUE)
clusts.lung2 = ClassifyCells(object = lung2, classifier = lung1.rf,new.data = lung2@data)
names(clusts.lung2) = colnames(lung2@data)
lung2 = AddMetaData(object = lung2, metadata = clusts.lung2, col.name = "clusts_lung1")
table(lung2@meta.data$clusts_lung1)
lung2 = SetAllIdent(object = lung2, id='clusts_lung1')

# Plot with lung1 clusters overlaid
lung2 = RunTSNE(object = lung2, dims.use = 1:8, do.fast = TRUE, perplexity = 30)
TSNEPlot(lung2)

