rm(list=ls())

library(Redeconve)

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

colnames(mela_sc_count) <- as.character(mela_sc_meta$X)

# col_range <- 1:ncol(mela_sc_count)
# 
# colnames(mela_sc_count) <- col_range

mela_sc_meta$sampleInfo = "sample1"

mela_sc_meta$cellID = as.character(1:ncol(mela_sc_count))

ref = get.ref(mela_sc_count,mela_sc_meta,dopar = F)

res.ct = deconvoluting(ref,mela_st_count,genemode="def",hpmode="auto",dopar=T,ncores=8)

res.prop  =  to.proportion (res.ct)

write.csv(t(res.prop), file = "D:/Projects/CARD-master/data/画饼图/Mela/Redeconve.csv")