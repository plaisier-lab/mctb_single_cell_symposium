setwd('C:/Users/cplaisie/Dropbox (ASU)/mctb_single_cell_symposium_data')

# Run 2
#pdf('2_results.R')
#source('C:/Users/cplaisie/Dropbox (ASU)/mctb_single_cell_symposium/2_GettingAcquaintedWithR.R')
#dev.off()

# Run 3
pdf('3_results.pdf')
source('C:/Users/cplaisie/Dropbox (ASU)/mctb_single_cell_symposium/3_Seurat_basics.R')
dev.off()

# Run 4
pdf('4_results.pdf')
source('C:/Users/cplaisie/Dropbox (ASU)/mctb_single_cell_symposium/4_Seurat_clustering.R')
dev.off()

# Run 5
pdf('5_results.pdf')
source('C:/Users/cplaisie/Dropbox (ASU)/mctb_single_cell_symposium/5_Seurat_training_classifier_application_new_data.R')
dev.off()

# Run 6
pdf('6_results.pdf')
source('C:/Users/cplaisie/Dropbox (ASU)/mctb_single_cell_symposium/6_Seurat_visualization.R')
dev.off()

# Run 7
pdf('7_results.pdf')
source('C:/Users/cplaisie/Dropbox (ASU)/mctb_single_cell_symposium/7_Seurat_combining_scRNA_seq_datasets.R')
dev.off()
