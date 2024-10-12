rm(list=ls())
library(ggplot2)
library(SPOTlight)
library(SingleCellExperiment)
library(SpatialExperiment)
library(scater)
library(scran)
# library(CARD)
library(rhdf5)

# load("D:/Tools/Data/MOB_ST/CARD/MOB.dge.sceset.RData")
# ### load the spatial data 
# load("D:/Tools/Data/MOB_ST/CARD/Rep12_MOB_count_matrix.RData")

mela_st_var <- h5read("D:/Data/Breast_1/gm_input_st_r.h5ad", name='var')

mela_st_count <- h5read("D:/Data/Breast_1/gm_input_st_r.h5ad", name='X')

mela_st_count <- as(mela_st_count, "dgCMatrix")

mela_st_spatial <- read.table("D:/Data/Breast_1/bc_cor.txt", header = TRUE)

mela_sc_meta <- read.table("D:/Data/Breast_1/sc_card.txt", header = TRUE, sep = "\t")

rownames(mela_st_count) <- as.character(mela_st_var$`_index`)

# colnames(mela_st_count) <- colomns

colnames(mela_st_count) <- rownames(mela_st_spatial)

mela_sc_count <- h5read("D:/Data/Breast_1/gm_input_sc.h5ad", name='X')

mela_sc_var <- h5read("D:/Data/Breast_1/gm_input_sc.h5ad", name='var')

mela_sc_count <- as(mela_sc_count, "dgCMatrix")

rownames(mela_sc_count) <- as.character(mela_sc_var$genes)

col_range <- 1:ncol(mela_sc_count)

colnames(mela_sc_count) <- col_range

mela_sc_meta$cellID = as.character(1:ncol(mela_sc_count))

cell_types <- mela_sc_meta$cellType
cell_types <- factor(cell_types)
names(cell_types) <- mela_sc_meta$cellID

sce <- SingleCellExperiment(assays = list(counts = as(mela_sc_count, "sparseMatrix")))

spe <- SpatialExperiment(mela_st_count)

sce <- logNormCounts(sce)

dec <- modelGeneVar(sce)
hvg <- getTopHVGs(dec, n = 3000)

colLabels(sce) <- mela_sc_meta$cellType
genes <- !grepl(pattern = "^Rp[l|s]|Mt", x = rownames(sce))
mgs <- scoreMarkers(sce, subset.row = genes)

mgs_fil <- lapply(names(mgs), function(i) {
  x <- mgs[[i]]
  # Filter and keep relevant marker genes, those with AUC > 0.8
  x <- x[x$mean.AUC > 0, ]
  # Sort the genes from highest to lowest weight
  x <- x[order(x$mean.AUC, decreasing = TRUE), ]
  # Add gene and cluster id to the dataframe
  x$gene <- rownames(x)
  x$cluster <- i
  data.frame(x)
})

mgs_df <- do.call(rbind, mgs_fil)

assayNames(spe) <- "counts"

# split cell indices by identity
idx <- split(seq(ncol(sce)), colLabels(sce))
# downsample to at most 20 per identity & subset
n_cells <- 200
cs_keep <- lapply(idx, function(i) {
  n <- length(i)
  if (n < n_cells)
    n_cells <- n
  sample(i, n_cells)
})
sce <- sce[, unlist(cs_keep)]

res <- SPOTlight(
  x = sce,
  y = spe,
  groups = colData(sce)@listData$cellType,
  mgs = mgs_df,
  hvg = hvg,
  weight_id = "mean.AUC",
  group_id = "cluster",
  gene_id = "gene",
  slot_sp = "counts",
  model="ns",
  )

write.csv(res$mat, file = "D:/Projects/CARD-master/data/画饼图/breast_1/SPOTlight_output.csv", row.names = TRUE, col.names = TRUE)