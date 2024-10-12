rm(list=ls())
library(SingleCellExperiment)
library(pbmcapply)
library(CARD)

library(rhdf5)

mela_st_var <- h5read("D:/Data/Melanoma/st_mela.h5ad", name='var')

mela_st_count <- h5read("D:/Data/Melanoma/st_mela.h5ad", name='X')

mela_st_count <- as(mela_st_count, "dgCMatrix")

mela_st_spatial <- read.table("D:/Data/Melanoma/Mel_cor.txt", header = TRUE)

mela_sc_meta <- read.table("D:/Data/Melanoma/sc_card.txt", header = TRUE, sep = "\t")

rownames(mela_st_count) <- as.character(mela_st_var$gene)

# colnames(mela_st_count) <- colomns

colnames(mela_st_count) <- rownames(mela_st_spatial)

mela_sc_count <- h5read("D:/Data/Melanoma/sc_mela.h5ad", name='X')

mela_sc_var <- h5read("D:/Data/Melanoma/sc_mela.h5ad", name='var')

mela_sc_count <- as(mela_sc_count, "dgCMatrix")

rownames(mela_sc_count) <- as.character(mela_sc_var$Cell)

col_range <- 1:ncol(mela_sc_count)

colnames(mela_sc_count) <- col_range

mela_sc_meta$sampleInfo = "sample1"

mela_sc_meta$cellID = as.character(1:ncol(mela_sc_count))

CARD_obj = createCARDObject(
  sc_count = mela_sc_count,
  sc_meta = mela_sc_meta,
  spatial_count = mela_st_count,
  spatial_location = mela_st_spatial,
  ct.varname = "cellType",
  ct.select = unique(mela_sc_meta$cellType),
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

write.csv(CARD_deconv@Proportion_CARD, file = "D:/Projects/CARD-master/data/画饼图/Mela/CARD_output.csv", row.names = FALSE)
