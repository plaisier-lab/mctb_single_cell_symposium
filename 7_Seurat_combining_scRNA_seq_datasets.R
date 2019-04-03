##########################################################
## MCTB Single Cell Symposium:                          ##
##     7_Seurat_combining_scRNA_seq_datasets.R          ##
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

lung3 = readRDS("lung3_processed.rds")

lungs = MergeSeurat(lung1, lung3, add.cell.id1 = "lung1", add.cell.id2 = "lung3")

lungs = NormalizeData(object = lungs)

# Find the most variable genes in the dataset
lungs = FindVariableGenes(object = lungs)
length(lungs@var.genes)

# Scale the data and regress out percent of mitochondrial and nUMI effects
lungs = ScaleData(object = lungs, vars.to.regress = c('percent_mito','nUMI'), genes.use = lungs@var.genes)

# Run PCA analysis
lungs = RunPCA(object = lungs, pc.genes = lungs@var.genes, pcs.compute = 40, pcs.print = 1:30, maxit = 500, weight.by.var = FALSE)

# Run TSne
lungs = RunTSNE(object = lungs, dims.use = 1:10, do.fast = TRUE, perplexity = 100)
TSNEPlot(lungs)

# Set celltype with classifier annotations as ident
lungs = SetAllIdent(lungs, id = "celltype_w_myeloid")
TSNEPlot(lungs)

# Identify variable genes for canonical correlation analysis to harmonize data set
g.1 = head(lung1@var.genes, 2000)
g.2 = head(lung3@var.genes, 2000)

# Remove duplicates from the list and find genes that are in both data sets
genes.use = unique(c(g.1, g.2))
genes.use = intersect(genes.use, rownames(lung1@scale.data))
genes.use = intersect(genes.use, rownames(lung3@scale.data))

# Run the CCA
lungs_CCA = RunCCA(lung1, lung3, genes.use = genes.use, num.cc = 30, add.cell.id1  = "lung1", add.cell.id2 = "lung3")

# Compare PCA plot to CCA plot
PCAPlot(lungs, group.by = "orig.ident")

DimPlot(object = lungs_CCA, reduction.use = "cca", group.by = "orig.ident", pt.size = 0.5)

# Find out how many dimensions we should use this takes a long time so I've commented it out
#Metageneplot <- MetageneBicorPlot(lungs_CCA, grouping.var = "orig.ident", dims.eval = 1:30, display.progress = FALSE)

# Align the data using the CCA
lungs_CCA <- AlignSubspace(lungs_CCA, reduction.type = "cca", grouping.var = "orig.ident", dims.align = 1:20)

# Run TSNE on the new data
lungs_CCA <- RunTSNE(lungs_CCA, reduction.use = "cca.aligned", dims.use = 1:20, do.fast = T, perplexity = 80)

# Compare merged vs CCA alighend data
TSNEPlot(lungs_CCA, group.by = "orig.ident")
TSNEPlot(lungs, group.by = "orig.ident")

# Set the ident to celltype
lungs_CCA <- SetAllIdent(lungs_CCA, id = "celltype_w_myeloid")
lungs <- SetAllIdent(lungs, id = "celltype_w_myeloid")

# Demonstrate that merging using cca vs non-cca doesn't affect the underlying data
# Table of celltype idents
table(lungs_CCA@ident)
table(lungs@ident)

# View data matix
lungs@data[1:10,1:30]

lungs_CCA@data[1:10,1:30]
