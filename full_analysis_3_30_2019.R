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
library(dplyr)

# Set your working directory
setwd('C:/Users/cplaisie/Dropbox (ASU)/mctb_single_cell_symposium_data')
#setwd('C:/Users/plais/Dropbox (ASU)/mctb_single_cell_symposium_data')

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
mito_genes = grep(pattern = "^MT-", x = rownames(x = lung1@data), value = TRUE)
length(mito_genes)
mito_genes
percent_mito = colSums(lung1@data[mito_genes, ]) / colSums(lung1@data)
lung1 = AddMetaData(object = lung1, metadata = percent_mito, col.name = "percent_mito")

# Violin plot of three parameters
VlnPlot(object = lung1, features.plot = c("nGene", "nUMI", "percent_mito"), nCol = 3)

# Plot of Percent.mito vs. nGene or nUMI
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
lung1 = ScaleData(object = lung1, vars.to.regress = c('percent_mito','nUMI'), genes.use = lung1@var.genes, model.use = 'negbinom')

# Run PCA analysis
lung1 = RunPCA(object = lung1, pc.genes = lung1@var.genes, pcs.compute = 40, pcs.print = 1:30, maxit = 500, weight.by.var = FALSE)

# Plot PCAs
PCHeatmap(object = lung1, pc.use = 1:12, cells.use = 300, do.balanced = TRUE, label.columns = FALSE)

# Run TSNE using PCs 1:10 with perplexity of 100
lung1 = RunTSNE(object = lung1, dims.use = 1:10, do.fast = TRUE, perplexity = 100)
TSNEPlot(lung1)

# Find clusters
lung1 = FindClusters(object = lung1, dims.use = 1:10, reduction.type = "pca", resolution = 0.1, print.output = 0, save.SNN  = TRUE)
TSNEPlot(lung1)

# Create new meta data column with cell type annotations using plyr
current.cluster.ids = c(0, 1, 2, 3, 4, 5, 6)
new.cluster.ids = c("Fibroblast", "Myeloid", "Myeloid", "Myeloid",
                     "Secretory", "Cilliated", "Endothelial")
lung1@meta.data$celltype = plyr::mapvalues(x = lung1@ident, from = current.cluster.ids, to = new.cluster.ids)
lung1 = SetIdent(lung1, ident.use = lung1@meta.data$celltype)
TSNEPlot(lung1)

# Find all marker genes
lung1.markers = FindAllMarkers(object = lung1, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
lung1.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
write.csv(lung1.markers,'lung1.markers.csv')

# Build classifier for clusters -
#lung1_rf_mg = BuildRFClassifier(object = lung1, training.genes = (lung1.markers %>% #group_by(cluster) %>% top_n(100, avg_logFC))$gene, training.classes = lung1@ident)
#lung1_rf_mg

lung1_rf = BuildRFClassifier(object = lung1, training.genes = lung1@var.genes, training.classes = lung1@ident)
lung1_rf


########################
## Load up lung2 data ##
########################
lung2_counts = Read10X(data.dir = "Lung2/outs/filtered_gene_bc_matrices/GRCh38/")
lung2 = CreateSeuratObject(raw.data = lung2_counts)
mito_genes = grep(pattern = "^MT-", x = rownames(x = lung2@data), value = TRUE)
percent_mito = colSums(lung2@data[mito_genes, ]) / colSums(lung2@data)
lung2 = AddMetaData(object = lung2, metadata = percent_mito, col.name = "percent_mito")

# Violin plot of three parameters
VlnPlot(object = lung2, features.plot = c("nGene", "nUMI", "percent_mito"), nCol = 3)

# Plot of Percent.mito vs. nGene or nUMI
# Best guess data was filtered with nUMI >= 200
par(mfrow = c(1,2))
plot(lung2@meta.data$percent_mito, lung2@meta.data$nUMI, col = rgb(0,0,1,0.5), ylim = c(0,38000), xlim = c(0,0.35))
abline(v = 0.2, col = 'red', lwd = 2)
abline(h = 200, col = 'red', lwd = 2)
plot(lung2@meta.data$percent_mito, lung2@meta.data$nGene, col=rgb(0,0,1,0.5), ylim = c(0,6000), xlim = c(0,0.35))
abline(v = 0.2, col = 'red', lwd = 2)
dim(lung2@data)
lung2 = FilterCells(object = lung2, subset.names = c("nUMI", "percent_mito"), low.thresholds = c(200, -Inf), high.thresholds = c(Inf, 0.2))
dim(lung2@data)
lung2 = NormalizeData(object = lung2)
lung2 = FindVariableGenes(object = lung2)
lung2 = ScaleData(object = lung2, vars.to.regress = c('percent_mito','nUMI'), genes.use = lung2@var.genes, model.use = 'negbinom')
lung2 = RunPCA(object = lung2, pc.genes = lung2@var.genes, pcs.compute = 40, pcs.print = 1:30, maxit = 500, weight.by.var = FALSE)

PCHeatmap(object = lung2, pc.use = 1:12, cells.use = 300, do.balanced = TRUE, label.columns = FALSE)

# Apply classifier
#lung1 = FindClusters(object = lung1, reduction.type = "pca", dims.use = 1:6, resolution = 0.3, print.output = 0, save.SNN = TRUE)
clusts_lung2 = ClassifyCells(object = lung2, classifier = lung1_rf,new.data = lung2@data)
names(clusts_lung2) = colnames(lung2@data)

lung2 = AddMetaData(object = lung2, metadata = clusts_lung2, col.name = "clusts_lung1")
table(lung2@meta.data$clusts_lung1)
lung2 = SetAllIdent(object = lung2, id='clusts_lung1')
# Plot with lung1 clusters overlaid
lung2 = RunTSNE(object = lung2, dims.use = 1:8, do.fast = TRUE, perplexity = 100)

TSNEPlot(lung2)

# Subset data based on myeloid cell types
lung2_myeloid = SubsetData(object = lung2, cells.use=rownames(lung2@meta.data)[which(lung2@meta.data$clusts_lung1=='Myeloid')])
clusts_lung2_myeloid = ClassifyCells(object = lung2_myeloid, classifier = pbmc_rf, new.data = lung2_myeloid@data)
names(clusts_lung2_myeloid) = colnames(lung2_myeloid@data)
lung2_myeloid = AddMetaData(object = lung2_myeloid, metadata = clusts_lung2_myeloid, col.name = "clusts_pbmc")
table(lung2_myeloid@meta.data$clusts_pbmc)
lung2_myeloid = SetAllIdent(object = lung2, id='clusts_pbmc')

# Update main lung2
clusts_lung2_w_myeloid = as.character(clusts_lung2)
names(clusts_lung2_w_myeloid) = names(clusts_lung2)
clusts_lung2_w_myeloid[names(clusts_lung2_myeloid)] = as.character(clusts_lung2_myeloid)

lung2 = AddMetaData(object = lung2, metadata = clusts_lung2_w_myeloid, col.name = "clusts_lung1_w_myeloid")
table(lung2@meta.data$clusts_lung1_w_myeloid)
lung2 = SetAllIdent(object = lung2, id='clusts_lung1_w_myeloid')
TSNEPlot(lung2)



##########################################################
## Load up mouse data and convert IDs so we can compare ##
## PMID = 30318149                                      ##
##########################################################

mouse_data = read.csv('mouse_lung/mouse_lung_cell_type_PMID30318149.csv', header = TRUE, row.names = 1)

mouse_lung = CreateSeuratObject(raw.data = mouse_data[,1:125])
mouse_lung = NormalizeData(object = mouse_lung)
mouse_lung = FindVariableGenes(object = mouse_lung)
mouse_lung = ScaleData(object = mouse_lung, genes.use = mouse_lung@var.genes) #, model.use = 'negbinom')
mouse_lung = AddMetaData(object = mouse_lung, metadata = sapply(colnames(mouse_lung@data), function(x) { strsplit(x,"_")[[1]][1] }), col.name = "cell_type")

mouse_lung.markers = FindAllMarkers(object = mouse_lung, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

# Build classifier for clusters -

mouse_lung.rf = BuildRFClassifier(object = mouse_lung, training.genes = mouse_lung@var.genes, training.classes = mouse_lung@meta.data$cell_type)
mouse_lung.rf
# Convert from mouse gene symbols to human gene symbols
library(biomaRt)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genes = mouse_lung.rf$forest$independent.variable.names
genes = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = genes ,mart = mouse, attributesL = c("hgnc_symbol","chromosome_name", "start_position"), martL = human, uniqueRows=T)

#
tmp = genes[,'MGI.symbol']
names(tmp) = genes[,'HGNC.symbol']
common_names = intersect(genes[,'HGNC.symbol'], rownames(lung1@data))
tmp1 = lung1@scale.data[common_names,]
rownames(tmp1) = tmp[rownames(tmp1)]

cell_types_lung1 = ClassifyCells(object = lung1, classifier = mouse_lung.rf, new.data = tmp1)
names(cell_types_lung1) = colnames(lung1@data)
lung1 = AddMetaData(object = lung1, metadata = cell_types_lung1, col.name = "cell_types_lung1")
table(lung1@meta.data$cell_types_lung1)
lung1 = SetAllIdent(object = lung1, id='cell_types_lung1')
TSNEPlot(lung1)

################################
## Load PBMC data 3k from 10X ##
################################

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
