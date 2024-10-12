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

osm_sc_meta$cellID = as.character(1:ncol(osm_sc_count))

cell_types <- osm_sc_meta$cellType
cell_types <- factor(cell_types)
names(cell_types) <- osm_sc_meta$cellID

sce <- SingleCellExperiment(assays = list(counts = as(osm_sc_count, "sparseMatrix")))

spe <- SpatialExperiment(osm_st_count)

sce <- logNormCounts(sce)

dec <- modelGeneVar(sce)
hvg <- getTopHVGs(dec, n = 3000)

colLabels(sce) <- osm_sc_meta$cellType
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

write.csv(res$mat, file = "D:/Projects/CARD-master/data/画饼图/osmFISH/SPOTlight_output.csv", row.names = TRUE, col.names = TRUE)