##########################################################
## MCTB Single Cell Symposium: 3_Seurat_basics.R        ##
## Seurat V3.0                                          ##
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
setwd('/files/Dropbox (ASU)/mctb_single_cell_symposium_data')

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
lung1 = CreateSeuratObject(counts = lung1_counts, project = "lung1")

# Structure of a Seurat object
str(lung1)

# Raw data storage
lung1@assays$RNA@counts[1:20,1:20]

# Data storage
lung1@assays$RNA@data[1:20,1:20] #Currently set as raw.data

# Add mitochondrial gene percentatge. Use mitochondrial genes - mitochondria (MT) have a double membrane,
# their own genome, and transcriptases. Ratio of MT genes to general cell genes is a great way to assess
# cell integrity as MT genes will have a harder time escaping the cell if the membrane is compromised.
lung1[['percent.mt']] = PercentageFeatureSet(lung1, pattern='^MT-')

# Violin plot of three parameters
VlnPlot(lung1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Plot of percent_mito vs. nGene or nUMI
# Best guess data was filtered with nUMI >= 2000
plot1 = FeatureScatter(lung1, feature1 = 'nCount_RNA', feature2 = 'percent.mt')
plot2 = FeatureScatter(lung1, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')
CombinePlots(plots = list(plot1, plot2))

# Filtering cells
#  - nUMI >= 2000 UMI
#  - percent_mito <= 10%
dim(lung1@assays$RNA@data)
lung1 = subset(lung1, subset = nCount_RNA > 2000 & percent.mt < 10)
dim(lung1@assays$RNA@data)

# Normalize data
lung1 = NormalizeData(lung1)
lung1@assays$RNA@data[1:20,1:20] # Now normalized data

# Find the most variable genes in the dataset
lung1 = FindVariableFeatures(lung1, selection.method = 'vst', nfeatures=2000)
length(VariableFeatures(lung1))

# Scale the data and regress out percent of mitochondrial and nUMI effects
lung1 = ScaleData(lung1, vars.to.regress = c('percent.mt','nCounts_RNA'))

# Run PCA analysis
lung1 = RunPCA(lung1, features = VariableFeatures(lung1), npcs = 30)

# Plot PCAs
DimHeatmap(lung1, dims = 1:12, cells = 500, balanced = TRUE)


########################
## Load up lung2 data ##
########################
# Load up data
lung2_counts = Read10X(data.dir = "Lung2/outs/filtered_gene_bc_matrices/GRCh38/")

# Take a look at the column and row labels
dim(lung2_counts)
colnames(lung2_counts)[1:10]
rownames(lung2_counts)[1:10]

# Create Seurat object
lung2 = CreateSeuratObject(counts = lung2_counts, project = "lung2")

# Structure of a Seurat object
str(lung2)

# Raw data storage
lung2@assays$RNA@counts[1:20,1:20]

# Data storage
lung2@assays$RNA@data[1:20,1:20] #Currently set as raw.data

# Add mitochondrial gene percentatge. Use mitochondrial genes - mitochondria (MT) have a double membrane,
# their own genome, and transcriptases. Ratio of MT genes to general cell genes is a great way to assess
# cell integrity as MT genes will have a harder time escaping the cell if the membrane is compromised.
lung2[['percent.mt']] = PercentageFeatureSet(lung2, pattern='^MT-')

# Violin plot of three parameters
VlnPlot(lung2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Plot of percent_mito vs. nGene or nUMI
# Best guess data was filtered with nUMI >= 2000
plot1 = FeatureScatter(lung2, feature1 = 'nCount_RNA', feature2 = 'percent.mt')
plot2 = FeatureScatter(lung2, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')
CombinePlots(plots = list(plot1, plot2))

# Filtering cells
#  - nUMI >= 2000 UMI
#  - percent_mito <= 10%
dim(lung2@assays$RNA@data)
lung2 = subset(lung2, subset = nCount_RNA > 200 & percent.mt < 20)
dim(lung2@assays$RNA@data)

# Normalize data
lung2 = NormalizeData(lung2)
lung2@assays$RNA@data[1:20,1:20] # Now normalized data

# Find the most variable genes in the dataset
lung2 = FindVariableFeatures(lung2, selection.method = 'vst', nfeatures=2000)
length(VariableFeatures(lung2))

# Scale the data and regress out percent of mitochondrial and nUMI effects
lung2 = ScaleData(lung2, vars.to.regress = c('percent.mt','nCounts_RNA'))

# Run PCA analysis
lung2 = RunPCA(lung2, features = VariableFeatures(lung2), npcs = 30)

# Plot PCAs
DimHeatmap(lung2, dims = 1:12, cells = 500, balanced = TRUE)

