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

lung3 = readRDS("Lung3.rds")

# Mito genes
mito_genes = grep(pattern = "^MT-", x = rownames(x = lung3@data), value = TRUE)
percent_mito = colSums(lung3@data[mito_genes, ]) / colSums(lung3@data)

# Add to lung3 Seurat object meta data
lung3 = AddMetaData(object = lung3, metadata = percent_mito, col.name = "percent_mito")

# Filter cells
lung3 = FilterCells(object = lung3, subset.names = c("nUMI", "percent_mito"), low.thresholds = c(2000, -Inf), high.thresholds = c(Inf, 0.1))

# Normalize data
lung3 = NormalizeData(object = lung3)

# Find the most variable genes in the dataset
lung3 = FindVariableGenes(object = lung3)

# Scale the data and regress out percent of mitochondrial and nUMI effects
lung3 = ScaleData(object = lung3, vars.to.regress = c('percent_mito','nUMI'), genes.use = lung3@var.genes)

# Run PCA analysis
lung3 = RunPCA(object = lung3, pc.genes = lung3@var.genes, pcs.compute = 40, pcs.print = 1:30, maxit = 500, weight.by.var = FALSE)
# Run TSNE
lung3 = RunTSNE(object = lung3, dims.use = 1:10, do.fast = TRUE, perplexity = 100)


# Run classifier on lung3 dataset
celltype_lung3 = ClassifyCells(object = lung3, classifier = lung1_rf,new.data = lung3@data)

# Add names of cells to the calls
names(celltype_lung3) = colnames(lung3@data)

# Put lung3 cell type names into lung3 Seurat object
lung3 = AddMetaData(object = lung3, metadata = celltype_lung3, col.name = "celltype")

# Set this as the main identifier for labeling and clustering
lung3 = SetAllIdent(object = lung3, id='celltype')

# Plot with lung3 clusters overlaid
TSNEPlot(lung3, do.label = TRUE)

# Subset "myeloid cells"
lung3_myeloid = SubsetData(object = lung3, cells.use = rownames(lung3@meta.data)[which(lung3@meta.data$celltype=='Myeloid')])

# Classify based on PBMC myeloid cell types
celltype_lung3_myeloid = ClassifyCells(object = lung3_myeloid, classifier = pbmc_rf, new.data = lung3_myeloid@data)

# Label the cell type classifications
names(celltype_lung3_myeloid) = colnames(lung3_myeloid@data)

# Update main lung3
celltype_lung3_w_myeloid = as.character(lung3@meta.data$celltype)
names(celltype_lung3_w_myeloid) = rownames(lung3@meta.data)

# Change myeloid cells to PBMC cell types
celltype_lung3_w_myeloid[names(celltype_lung3_myeloid)] = as.character(celltype_lung3_myeloid)

# Add meta data to lung2 (celltype_w_myeloid)
lung3 = AddMetaData(object = lung3, metadata = celltype_lung3_w_myeloid, col.name = "celltype_w_myeloid")

# Update the orig.ident field in the meta.data
lung3@meta.data$orig.ident <- "lung3"

# Save R object
saveRDS(lung3, file = "lung3_processed.rds")
