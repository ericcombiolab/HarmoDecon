rm(list=ls())
library(SingleCellExperiment)
library(pbmcapply)
library(CARD)
#> Warning: replacing previous import 'RcppML::nmf' by 'NMF::nmf' when loading
#> 'CARD'
### load the single cell RNA-seq data used as the reference for deconvolution
load("D:/Tools/Data/MOB_ST/CARD/MOB.dge.sceset.RData")
### load the spatial data 
load("D:/Tools/Data/MOB_ST/CARD/Rep12_MOB_count_matrix.RData")
### set paramaters
ct.varname = "cellType"
sample.varname = "sampleInfo"
ct.select = as.character(unique(colData(sce)$cellType))
location <- cbind.data.frame(x=as.numeric(sapply(strsplit(colnames(MOB_raw),split="x"),"[",1)),y=as.numeric(sapply(strsplit(colnames(MOB_raw),split="x"),"[",2)))
rownames(location) = colnames(MOB_raw)
##### create CARD object
CARD_obj = createCARDObject(
  sc_count = assays(sce)$counts,
  sc_meta = colData(sce),
  spatial_count = MOB_raw,
  spatial_location = location,
  ct.varname = ct.varname,
  ct.select = ct.select,
  sample.varname = sample.varname,
  minCountGene = 100,
  minCountSpot =5) 

CARD_obj = CARD_deconvolution(CARD_obj)

print(apply(CARD_obj@Proportion_CARD,2,summary))

CARD_obj = CARD.imputation(CARD_obj,NumGrids = 2000,ineibor = 10,exclude = NULL)## print output

print(apply(CARD_obj@refined_prop,2,summary))

write.csv(CARD_obj@Proportion_CARD, file = "D:/Projects/CARD-master/data/画饼图/MOB/CARD_output_mob.csv", row.names = FALSE)

colors = c("#cfa6b2","#fcd4c7","#6d889f","#f47c50","#a8af7f","#7FC97F","#BEAED4",
           "#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77","#D95F02",
           "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D")

CARD_obj@spatial_location = round(CARD_obj@spatial_location, 0)

p1 <- CARD.visualize.pie(proportion = CARD_obj@Proportion_CARD, 
                         spatial_location = CARD_obj@spatial_location, colors = colors, radius = 0.49) 



print(p1)