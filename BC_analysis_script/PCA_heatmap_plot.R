#!/usr/bin/env Rscript

#################################################################################
#     Author: Olena Maiakovska
#     Date: May, 2022
#     Availability: 
#     Usage:
#     PCA_heatmap_plot.R <relativeconcetration_file.xls>
################################################################################

args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.png"
}

#Install missing packages if needed 
list.of.packages <- c("readxl", "ggfortify", "ggplot2", "pheatmap", "RColorBrewer")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load required libraries
library("xlsx")
library("ggfortify")
library(pheatmap)
library(RColorBrewer)

# Load xls table from provided 
RC_values = read.xlsx(args[1], 1, header=T)
rownames(RC_values) = RC_values[,1]
RC_values_matrix = data.matrix(RC_values[,-1])

#To zoom in into cluster and find differences between its members, the values of several variants have to be excluded from the matix and visualization:
#delete them by following command (uncomment the next line and re-run the entire script)
#Deleting variants AAV3_L1 and AAV5_WT from the dataset provided in the manuscript Rapti et al, 2022.
#RC_values_matrix = RC_values_matrix[-c(13,37),] 
# Perform a principal components analysis on the given data matrix 
pca_res = prcomp(na.omit(RC_values_matrix), scale. = TRUE)

#If your dataset comprise o more than 50 capsid variants its recommended to resize the png plot, by modifying width and height arguments in autoplot function. 
# Plot PCA result, defining the clusters of similar/shared properties.
#To change PCs edit x and y to other PC number, in total there are ncol(pca_res$x) number of PCs. 
png(filename = 'PCA_all_capsid_variants.png', width = 768, height = 768, units = 'px')
autoplot(pca_res,x = 1, y=2,  label = TRUE, data=na.omit(RC_values_matrix), colour = "darkblue")
dev.off()
## 

#creating color palette 
myCol <- c( "azure2", brewer.pal(9, "YlOrRd"))

# Formatting RC matrix by adding the absolute constant value to eliminate zeros before log transformation
png(filename = 'Heatmap_all_capsids', width = 768, height = 768, units = 'px')
minimum_value = min(apply(RC_values_matrix[,], 1, function(x) min(x[x>0])))
RC_values_matrix_const <- RC_values_matrix + minimum_value
# if the gaps between the branches are needed add "cutree_rows = 5" or "cutree_rows = 5" arguments with any other number needed
# if axes have to be changes transform the matrix by adding "t" to the matrix, like: "t(RC_values_matrix_const)"
pheatmap(log2(RC_values_matrix_const), col=myCol, cluster_rows = T, cluster_cols = T, border_color = "white", main = "Hierarchical clustering of log2-transformed RC values")
dev.off()
