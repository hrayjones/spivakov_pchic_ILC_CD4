library("ComplexHeatmap")
library(dataRetrieval)
library(dplyr) # for `rename` & `select`
library(tidyr) # for `gather`
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(ggrepel)
library(ggExtra)
library(gridExtra)
library(grid)
library(ggpubr)
library(dtwclust)
library(viridis)
library(circlize)
library(dplyr)
library(reshape2)

#---------------------------------------------------#
# Paths
#---------------------------------------------------#
cor_path = "/Users/caz3so/scratch/20220428_RELI_spivakov_rerun/output/corr/Spearman_Corr_bigwigScores_skipZeroes_removeOutliers.tab"
meta = "/Users/caz3so/scratch/20220428_RELI_spivakov_rerun/meta/spivakov_detailed_sample_meta.tsv"

df = read.csv(cor_path, sep="\t", row.names = "srx")

order = colnames(df)


# Import the meta data
meta_df = read.csv(meta, sep="\t", header=1, row.names="srx")
meta_df <- meta_df[order,]

top_anno = HeatmapAnnotation("layout" = meta_df$layout,
                             "condition" = meta_df$condition,
                             col = list("layout" = c("SINGLE" = "#d1495b", "PAIRED" = "#66a182"),
                                        "condition" = c("K4me3" = "#8dd3c7", "K27ac" = "#ffffb3", "Input" = "#bebada", "ATAC" = "#fb8072")))

left_anno = rowAnnotation("layout" = meta_df$layout,
                             "condition" = meta_df$condition,
                             col = list("layout" = c("SINGLE" = "#d1495b", "PAIRED" = "#66a182"),
                                        "condition" = c("K4me3" = "#8dd3c7", "K27ac" = "#ffffb3", "Input" = "#bebada", "ATAC" = "#fb8072")))


mat = as.matrix(df)

#---------------------------------------------------#
# Set up Heatmaps
#---------------------------------------------------#

ht = Heatmap(mat,
             name="Pearson Correlation", 
             cluster_rows = TRUE, 
             cluster_columns = TRUE, 
             col=heat_col, 
             column_names_gp = grid::gpar(fontsize = 12), 
             row_names_gp = grid::gpar(fontsize = 12),
             top_annotation=top_anno,
             left_annotation = left_anno,
             border = TRUE
             )

ht
