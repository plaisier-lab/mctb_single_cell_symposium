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

# Import Seurat single cell analysis software and dplyr
library(Seurat)
library(dplyr)

# Set your working directory
setwd('C:/Users/cplaisie/Dropbox (ASU)/mctb_single_cell_symposium_data')

########################
## Load up lung1 data ##
########################

# Load up data
lung1_counts = Read10X(data.dir = "Lung1/outs/filtered_gene_bc_matrices/GRCh38/")

# Take a look at the column and row labels
dim(lung1_counts)
colnames(lung1_counts)[1:10]
rownames(lung1_counts)[1:10]

# Create Seurat object
lung1 = CreateSeuratObject(raw.data = lung1_counts, project = "lung1")

# Structure of a Seurat object
str(lung1)

# Raw data storage
lung1@raw.data[1:20,1:20]

# Data storage
lung1@data[1:20,1:20] #Currently set as raw.data

# Select mitochondrial genes - mitochondria (MT) have a double membrane, their own genome, and transcriptases.
# Ratio of MT genes to general cell genes is a great way to assess cell integrity as MT genes will have a harder
# time escaping the cell if the membrane is compromised.
mito_genes = grep(pattern = "^MT-", x = rownames(x = lung1@data), value = TRUE)
length(mito_genes)
mito_genes

# Calculate percentage of mitochondria
percent_mito = colSums(lung1@data[mito_genes, ]) / colSums(lung1@data)

# Add to lung1 Seurat object meta data
lung1 = AddMetaData(object = lung1, metadata = percent_mito, col.name = "percent_mito")

# Violin plot of three parameters
VlnPlot(object = lung1, features.plot = c("nGene", "nUMI", "percent_mito"), nCol = 3)

# Plot of percent_mito vs. nGene or nUMI
# Best guess data was filtered with nUMI >= 2000
par(mfrow = c(1,2))
plot(lung1@meta.data$percent_mito, lung1@meta.data$nUMI, col = rgb(0,0,1,0.5), ylim = c(0,38000), xlim = c(0,0.35))
abline(v = 0.1, col = 'red', lwd = 2)
abline(h = 2000, col = 'red', lwd = 2)
plot(lung1@meta.data$percent_mito, lung1@meta.data$nGene, col=rgb(0,0,1,0.5), ylim = c(0,6000), xlim = c(0,0.35))
abline(v = 0.1, col = 'red', lwd = 2)

# Filtering cells
#  - nUMI >= 2000 UMI
#  - percent_mito <= 0.1
lung1 = FilterCells(object = lung1, subset.names = c("nUMI", "percent_mito"), low.thresholds = c(2000, -Inf), high.thresholds = c(Inf, 0.1))
dim(lung1@data)

# Normalize data
lung1 = NormalizeData(object = lung1)
lung1@data[1:20,1:20] # Now normalized data

# Find the most variable genes in the dataset
lung1 = FindVariableGenes(object = lung1)
length(lung1@var.genes)

# Scale the data and regress out percent of mitochondrial and nUMI effects
lung1 = ScaleData(object = lung1, vars.to.regress = c('percent_mito','nUMI'), genes.use = lung1@var.genes)

# Run PCA analysis
lung1 = RunPCA(object = lung1, pc.genes = lung1@var.genes, pcs.compute = 40, pcs.print = 1:30, maxit = 500, weight.by.var = FALSE)

# Plot PCAs
PCHeatmap(object = lung1, pc.use = 1:12, cells.use = 300, do.balanced = TRUE, label.columns = FALSE)


########################
## Load up lung2 data ##
########################
lung2_counts = Read10X(data.dir = "Lung2/outs/filtered_gene_bc_matrices/GRCh38/")
lung2 = CreateSeuratObject(raw.data = lung2_counts, project = "lung2")

# Select mitochondrial genes - mitochondria (MT) have a double membrane, their own genome, and transcriptases.
# Ratio of MT genes to general cell genes is a great way to assess cell integrity as MT genes will have a harder
# time escaping the cell if the membrane is compromised.
mito_genes = grep(pattern = "^MT-", x = rownames(x = lung2@data), value = TRUE)
percent_mito = colSums(lung2@data[mito_genes, ]) / colSums(lung2@data)
lung2 = AddMetaData(object = lung2, metadata = percent_mito, col.name = "percent_mito")

# Violin plot of three parameters
VlnPlot(object = lung2, features.plot = c("nGene", "nUMI", "percent_mito"), nCol = 3)

# Plot of percent_mito vs. nGene or nUMI
# Best guess data was filtered with nUMI >= 200
par(mfrow = c(1,2))
plot(lung2@meta.data$percent_mito, lung2@meta.data$nUMI, col = rgb(0,0,1,0.5), ylim = c(0,38000), xlim = c(0,0.35))
abline(v = 0.2, col = 'red', lwd = 2)
abline(h = 200, col = 'red', lwd = 2)
plot(lung2@meta.data$percent_mito, lung2@meta.data$nGene, col=rgb(0,0,1,0.5), ylim = c(0,6000), xlim = c(0,0.35))
abline(v = 0.2, col = 'red', lwd = 2)
dim(lung2@data)

# Filtering cells
#  - nUMI >= 2000 UMI
#  - percent_mito <= 0.2
lung2 = FilterCells(object = lung2, subset.names = c("nUMI", "percent_mito"), low.thresholds = c(200, -Inf), high.thresholds = c(Inf, 0.2))
dim(lung2@data)

# Normalize lung2
lung2 = NormalizeData(object = lung2)

# Find the most variable genes in the dataset
lung2 = FindVariableGenes(object = lung2)

# Scale the data and regress out percent of mitochondrial and nUMI effects
lung2 = ScaleData(object = lung2, vars.to.regress = c('percent_mito','nUMI'), genes.use = lung2@var.genes)

# Run PCA analysis
lung2 = RunPCA(object = lung2, pc.genes = lung2@var.genes, pcs.compute = 40, pcs.print = 1:30, maxit = 500, weight.by.var = FALSE)

# Plot PCAs
PCHeatmap(object = lung2, pc.use = 1:12, cells.use = 300, do.balanced = TRUE, label.columns = FALSE)
