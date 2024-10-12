rm(list=ls())
library(spacexr)
library(CARD)
library(rhdf5)
library(SingleCellExperiment)
library(pbmcapply)


load("D:/Tools/Data/MOB_ST/CARD/MOB.dge.sceset.RData")
### load the spatial data 
load("D:/Tools/Data/MOB_ST/CARD/Rep12_MOB_count_matrix.RData")
mob_sc_meta <- h5read("D:/Tools/Data/MOB_ST/CARD/mob_sc_celltype.h5ad", name='obs')
mapped <- as.character(mob_sc_meta$cellType$categories[match(mob_sc_meta$cellType$codes, 1:length(mob_sc_meta$cellType$categories)-1)])
mapped <- data.frame(cellType = mapped, cellID = sce@colData@rownames)
matrix <- assay(sce)
multipliedMatrix <- matrix * 100
multipliedMatrix@x <- round(multipliedMatrix@x, 0)
multipliedMatrix@x <- as(multipliedMatrix@x, "double")

cell_types <- mapped$cellType
cell_types <- factor(cell_types)
names(cell_types) <- mapped$cellID
# Error in check_cell_types(cell_types) : 
#   Reference: levels(cell_types) contains a cell type with name containing prohibited character /. Please rename this cell type.
levels(cell_types) <- c("EPL-IN", "GC", "MTC", "OSNs", "PGC")

reference <- Reference(multipliedMatrix, cell_types)

mob_spatial <- read.table("D:/Tools/Data/MOB_ST/CARD/mob_cor.txt", header = TRUE)

rownames(mob_spatial) <- colnames(MOB_raw)

spaceRNA <- SpatialRNA(mob_spatial, MOB_raw)

myRCTD <- create.RCTD(spaceRNA, reference)

myRCTD <- run.RCTD(myRCTD)

normalized_weights <- normalize_weights(myRCTD@results$weights)

rctd_out <- as.matrix(normalized_weights)

write.csv(rctd_out, file = "D:/Projects/CARD-master/data/画饼图/MOB/RCTD_output_mob.csv")

p1 <- CARD.visualize.pie(proportion = rctd_out, spatial_location = myRCTD@spatialRNA@coords, colors = colors, radius = 0.49) 

colors = c("#cfa6b2","#fcd4c7","#6d889f","#f47c50","#a8af7f","#7FC97F","#BEAED4",
           "#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77",
           "#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D")

print(p1)
