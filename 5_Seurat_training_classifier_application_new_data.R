##########################################################
## MCTB Single Cell Symposium:                          ##
## 5_Seurat_training_classifier_application_new_data.R  ##
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

# Build classifier for cell types found in lung1 sample
lung1_rf = BuildRFClassifier(object = lung1, training.genes = lung1@var.genes, training.classes = lung1@ident)
lung1_rf

# Confusion matrix for classifier
lung1_rf$confusion.matrix

# Calculate predictive value of classifier
confusion1 = tb.rf$confusion[c(2,1),c(2,1)]

# Run classifier on lung2 dataset
celltype_lung2 = ClassifyCells(object = lung2, classifier = lung1_rf,new.data = lung2@data)
table(celltype_lung2)

# Add names of cells to the calls
names(celltype_lung2) = colnames(lung2@data)

# Put lung1 cell type names into lung2 Seurat object
lung2 = AddMetaData(object = lung2, metadata = celltype_lung2, col.name = "celltype")

# Make sure numbers are the same
table(lung2@meta.data$celltype)

# Set this as the main identifier for labeling and clustering
lung2 = SetAllIdent(object = lung2, id='celltype')

# Plot with lung1 clusters overlaid
lung2 = RunTSNE(object = lung2, dims.use = 1:8, do.fast = TRUE, perplexity = 100)
TSNEPlot(lung2, do.label = TRUE)

###################################
## Classify with PBMC identities ##
###################################

# Build classifier from Seurat PBMC dataset
# Load the PBMC dataset
pbmc.data = Read10X(data.dir = "pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/")
pbmc = CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 200, project = "10X_PBMC")
mito_genes = grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent_mito = Matrix::colSums(pbmc@raw.data[mito_genes, ])/Matrix::colSums(pbmc@raw.data)
pbmc = AddMetaData(object = pbmc, metadata = percent_mito, col.name = "percent_mito")
pbmc = FilterCells(object = pbmc, subset.names = c("nGene", "percent_mito"), low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))
pbmc = NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc = FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
pbmc = ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent_mito"))
pbmc = RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
pbmc = FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = TRUE)
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7)
new.cluster.ids <- c("CD4 T cells", "CD14+ Monocytes", "B cells", "CD8 T cells", "FCGR3A+ Monocytes", "NK cells", "Dendritic cells", "Megakaryocytes")
pbmc@ident <- plyr::mapvalues(x = pbmc@ident, from = current.cluster.ids, to = new.cluster.ids)
pbmc = RunTSNE(object = pbmc, dims.use = 1:10, do.fast = TRUE)
TSNEPlot(object = pbmc, do.label = TRUE, pt.size = 0.5)
pbmc_rf = BuildRFClassifier(object = pbmc, training.genes = pbmc@var.genes, training.classes = pbmc@ident)
pbmc_rf
save(pbmc_rf, file = 'pbmc_rf.RData')

# WE COULD JUST LOAD THE pbmc_rf AS A PRECOMPUTED CLASSIFIER WITH THE ABOVE CODE IN A SEPARATE SUPPLEMENTARY FILE

## lung1 PBMC classification
# Subset data so we only classify myeloid cells
lung1_myeloid = SubsetData(object = lung1, cells.use = rownames(lung1@meta.data)[which(lung1@meta.data$celltype=='Myeloid')])

# Number of cell in subset
dim(lung1@data)
dim(lung1_myeloid@data)

# Classify based on PBMC myeloid cell types
celltype_lung1_myeloid = ClassifyCells(object = lung1_myeloid, classifier = pbmc_rf, new.data = lung1_myeloid@data)

# What is the breakdown of PBMC cell types found in lung2
table(celltype_lung1_myeloid)

# Label the cell type classifications
names(celltype_lung1_myeloid) = colnames(lung1_myeloid@data)

# Update main lung2
celltype_lung1_w_myeloid = as.character(lung1@meta.data$celltype)
names(celltype_lung1_w_myeloid) = rownames(lung1@meta.data)

# Change myeloid cells to PBMC cell types
celltype_lung1_w_myeloid[names(celltype_lung1_myeloid)] = as.character(celltype_lung1_myeloid)

# Add meta data to lung2 (celltype_w_myeloid)
lung1 = AddMetaData(object = lung1, metadata = celltype_lung1_w_myeloid, col.name = "celltype_w_myeloid")

# Check number match up
table(lung1@meta.data$celltype_w_myeloid)

# Set ident for plotting
lung1 = SetAllIdent(object = lung1, id='celltype_w_myeloid')

# Plot TSNE
TSNEPlot(lung1, do.label=T)

## lung2 PBMC classification
# Subset data so we only classify myeloid cells
lung2_myeloid = SubsetData(object = lung2, cells.use = rownames(lung2@meta.data)[which(lung2@meta.data$celltype=='Myeloid')])

# Classify based on PBMC myeloid cell types
celltype_lung2_myeloid = ClassifyCells(object = lung2_myeloid, classifier = pbmc_rf, new.data = lung2_myeloid@data)

# What is the breakdown of PBMC cell types found in lung2
table(celltype_lung2_myeloid)

# Update main lung2
celltype_lung2_w_myeloid = as.character(celltype_lung2)
names(celltype_lung2_w_myeloid) = names(celltype_lung2)

# Change myeloid cells to PBMC cell types
celltype_lung2_w_myeloid[names(celltype_lung2_myeloid)] = as.character(celltype_lung2_myeloid)

# Add meta data to lung2 (celltype_w_myeloid)
lung2 = AddMetaData(object = lung2, metadata = celltype_lung2_w_myeloid, col.name = "celltype_w_myeloid")

# Check number match up
table(lung2@meta.data$celltype_w_myeloid)

# Set ident for plotting
lung2 = SetAllIdent(object = lung2, id='celltype_w_myeloid')

# Plot TSNE
TSNEPlot(lung2, do.label=T)
