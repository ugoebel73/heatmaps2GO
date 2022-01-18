# heatmaps2go
This repository contains R code to extract gene ID <-> GO term associations from different databases and to display the set of genes annotated to a given set of GO terms as a heatmap.
Two non-trivial aspects of this simple endeavour are: 
1. to break the input gene list down into non-overlapping gene sets, such that each set is primarily composed of one genes annotated to a single one of the input GO terms. 
2. to display the gene sets defined by (1.) as as stacked sub-matrices (and in general to experiment with features of ComplexHeatmap)
