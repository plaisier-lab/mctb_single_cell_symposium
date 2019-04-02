
##########################################################
## MCTB Single Cell Symposium: 6_Seurat_visualization.R ##
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

# Find all the differntially expressed genes between celltypes - using max.cells.per.ident to speed up for the workshop
lung_celltype_markers = FindAllMarkers(lung1, max.cells.per.ident = 100)


top10 <- lung_celltype_markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

VlnPlot(lung1, "KRT18")
DoHeatmap(lung1, genes.use = head(row.names(lung_celltype_markers)))