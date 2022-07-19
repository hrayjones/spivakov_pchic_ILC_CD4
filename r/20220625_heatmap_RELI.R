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

#---------------------------------------------------#
# Dat Imports
#---------------------------------------------------#

path = "/Users/caz3so/workspaces/tacazares/pchic/data/RELI/20220625_RELI_98phenotypes_5samples.tsv"
  
# Import the file as a dataframe
data <- read.csv(path,
                 sep = "\t",
                 row.names = "disease")

#---------------------------------------------------#
# Plotting Figure 2
#---------------------------------------------------#
Heatmap(t(data))


        