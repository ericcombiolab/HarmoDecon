rm(list=ls())
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(rhdf5)
library(SingleCellExperiment)
library(CARD)
library(SeuratData)

load("D:/Tools/Data/MOB_ST/CARD/MOB.dge.sceset.RData")
### load the spatial data 
load("D:/Tools/Data/MOB_ST/CARD/Rep12_MOB_count_matrix.RData")
mob_sc_meta <- h5read("D:/Tools/Data/MOB_ST/CARD/mob_sc_celltype.h5ad", name='obs')
mapped <- as.character(mob_sc_meta$cellType$categories[match(mob_sc_meta$cellType$codes, 1:length(mob_sc_meta$cellType$categories)-1)])
mapped <- data.frame(cellType = mapped, cellID = sce@colData@rownames)
matrix <- assay(sce)
mob_sc <- CreateSeuratObject(counts = assays(sce)$counts, meta.data = colData(sce))
MOB_st = CreateSeuratObject(counts = MOB_raw, assay="Spatial")

mob_sc@meta.data$cellType <- sce@colData@listData$cellType

mob_sc <- SCTransform(mob_sc, verbose = FALSE)

mob_sc <- FindVariableFeatures(mob_sc)

mob_sc <- NormalizeData(mob_sc)

mob_sc <- ScaleData(mob_sc, assay = 'RNA')

mob_sc <- RunPCA(mob_sc, verbose = FALSE)

MOB_st <- SCTransform(MOB_st, verbose = FALSE, assay = 'Spatial')

mob_st <- NormalizeData(mob_st)

MOB_st <- ScaleData(MOB_st, assay = 'Spatial')

anchors <- FindTransferAnchors(
  reference = mob_sc,
  query = MOB_st,
  normalization.method = "SCT",
  reference.reduction = "pca",
  dims = 1:30
  )

MOB_st <- TransferData(
  anchorset = anchors, 
  reference = mob_sc,
  query = MOB_st,
  refdata = list(
    cellType = "cellType")
)

MOB_st <- IntegrateEmbeddings(
  anchorset = anchors,
  reference = mob_sc,
  query = MOB_st, 
  new.reduction.name = "ref.pca"
)

Seurat_output <- as.matrix(MOB_st@assays$prediction.score.cellType@data)

Seurat_output <- t(Seurat_output)

write.csv(Seurat_output, file = "D:/Projects/CARD-master/data/画饼图/MOB/mob_seurat.csv", row.names = TRUE, col.names = TRUE)

colors = c("#cfa6b2","#fcd4c7","#6d889f","#f47c50","#a8af7f","#7FC97F","#BEAED4",
           "#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77","#D95F02",
           "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D")

mob_spatial <- read.table("D:/Tools/Data/MOB_ST/CARD/mob_cor.txt", header = TRUE)

row.names(Seurat_output) <- row.names(mob_spatial)

p1 <- CARD.visualize.pie(proportion = Seurat_output, spatial_location = mob_spatial, colors = colors, radius = 0.49)

print(p1)

