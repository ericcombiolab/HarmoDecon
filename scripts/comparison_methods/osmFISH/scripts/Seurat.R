rm(list=ls())
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(rhdf5)
library(SingleCellExperiment)
library(CARD)
library(SeuratData)

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
# Error in check_cell_types(cell_types) :
#   Reference: levels(cell_types) contains a cell type with name containing prohibited character /. Please rename this cell type.
# levels(cell_types) <- c("Astro", "Endo", "Excitatory L23", "Excitatory L4", "Excitatory L5",
#                         "Excitatory L6", "Inhibitory Other", "Inhibitory Pvalb", "Inhibitory Sst",
#                         "Inhibitory Vip", "Micro", "Neuron Other", "Olig", "Other", "Smc")

osm_sc <- CreateSeuratObject(counts = osm_sc_count, meta.data = osm_sc_meta)
osm_st <- CreateSeuratObject(counts = osm_st_count, assay="Spatial")

# osm_sc@meta.data$cellType <- osm_sc@colData@listData$cellType

osm_sc <- SCTransform(osm_sc)

# osm_sc <- FindVariableFeatures(osm_sc)

# osm_sc <- NormalizeData(osm_sc, normalization.method = "LogNormalize", scale.factor = 10000)

# osm_sc <- ScaleData(osm_sc, assay = 'RNA')

# osm_sc <- RunPCA(osm_sc, verbose = FALSE)

osm_st <- SCTransform(osm_st, assay = 'Spatial')

# osm_st <- NormalizeData(osm_st, normalization.method = "LogNormalize", scale.factor = 10000)

# osm_st <- ScaleData(osm_st, assay = 'Spatial')

anchors <- FindTransferAnchors(
  reference = osm_sc,
  query = osm_st,
  normalization.method = "SCT",
  dims = 1:30,
)

osm_st <- TransferData(
  anchorset = anchors, 
  reference = osm_sc,
  query = osm_st,
  refdata = list(
    cellType = "cellType"),
  k.weight = 10,
  dims = 1:30 
)

# osm_st <- IntegrateEmbeddings(
#   anchorset = anchors,
#   reference = osm_sc,
#   query = osm_st, 
#   new.reduction.name = "ref.pca"
# )

Seurat_output <- as.matrix(osm_st@assays$prediction.score.cellType@data)

Seurat_output <- t(Seurat_output)

write.csv(Seurat_output, file = "D:/Projects/CARD-master/data/画饼图/osmFISH/Seurat_output.csv", row.names = TRUE, col.names = TRUE)

colors = c("#FCCDE5","#D9D9D9","#FFD92F","#4DAF4A","#377EB8","#7FC97F","#BEAED4",
           "#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77",
           "#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D")

osm_st_spatial <- read.table("D:/Data/osmFISH/osm_cor.txt", header = TRUE)

row.names(Seurat_output) <- row.names(osm_st_spatial)

p1 <- CARD.visualize.pie(proportion = Seurat_output, spatial_location = osm_st_spatial, colors = colors, radius = 0.49)

print(p1)


