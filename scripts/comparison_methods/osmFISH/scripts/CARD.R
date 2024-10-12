rm(list=ls())
library(SingleCellExperiment)
library(pbmcapply)
library(CARD)

library(rhdf5)

osm_st_var <- h5read("D:/Data/osmFISH/st_osm.h5ad", name='var')

osm_st_count <- h5read("D:/Data/osmFISH/st_osm.h5ad", name='X')

osm_st_count <- as(osm_st_count, "dgCMatrix")

osm_st_spatial <- read.table("D:/Data/osmFISH/osm_cor.txt", header = TRUE)

osm_sc_meta <- read.table("D:/Data/osmFISH/sc_card.txt", header = TRUE, sep = "\t")

rownames(osm_st_count) <- as.character(osm_st_var$`_index`)

# colnames(osm_st_count) <- colomns

colnames(osm_st_count) <- rownames(osm_st_spatial)

osm_sc_count <- h5read("D:/Data/osmFISH/sc_osm.h5ad", name='X')

osm_sc_var <- h5read("D:/Data/osmFISH/sc_osm.h5ad", name='var')

osm_sc_count <- as(osm_sc_count, "dgCMatrix")

rownames(osm_sc_count) <- as.character(osm_sc_var$`_index`)

col_range <- 1:ncol(osm_sc_count)

colnames(osm_sc_count) <- col_range

osm_sc_meta$sampleInfo = "sample1"

osm_sc_meta$cellID = as.character(1:ncol(osm_sc_count))

CARD_obj = createCARDObject(
  sc_count = osm_sc_count,
  sc_meta = osm_sc_meta,
  spatial_count = osm_st_count,
  spatial_location = osm_st_spatial,
  ct.varname = "cellType",
  ct.select = unique(osm_sc_meta$cellType),
  sample.varname = "sampleInfo",
  minCountGene = 100,
  minCountSpot = 5)

CARD_deconv = CARD_deconvolution(CARD_object = CARD_obj)

colors = c("#FFD92F","#4DAF4A","#FCCDE5","#D9D9D9","#377EB8","#7FC97F","#BEAED4",
           "#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77",
           "#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D")


# CARD_deconv@spatial_location$y = CARD_deconv@spatial_location$x * 1.75

p1 <- CARD.visualize.pie(
  proportion = CARD_deconv@Proportion_CARD,
  spatial_location = CARD_deconv@spatial_location, 
  colors = colors, 
  radius = 0.49)

print(p1)

write.csv(CARD_deconv@Proportion_CARD, file = "D:/Projects/CARD-master/data/CARD_output.csv", row.names = FALSE)
