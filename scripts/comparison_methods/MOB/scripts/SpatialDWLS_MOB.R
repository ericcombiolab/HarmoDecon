rm(list=ls())
library(reticulate)
library(rhdf5)
library(SingleCellExperiment)
library(Giotto)
library(CARD)

load("D:/Tools/Data/MOB_ST/CARD/MOB.dge.sceset.RData")
### load the spatial data 
load("D:/Tools/Data/MOB_ST/CARD/Rep12_MOB_count_matrix.RData")
mob_sc_meta <- h5read("D:/Tools/Data/MOB_ST/CARD/mob_sc_celltype.h5ad", name='obs')
mapped <- as.character(mob_sc_meta$cellType$categories[match(mob_sc_meta$cellType$codes, 1:length(mob_sc_meta$cellType$categories)-1)])
mapped <- data.frame(cellType = mapped, cellID = sce@colData@rownames)
matrix <- assay(sce)
mob_spatial <- read.table("D:/Tools/Data/MOB_ST/CARD/mob_cor.txt", header = TRUE)
multipliedMatrix <- matrix * 100
multipliedMatrix@x <- round(multipliedMatrix@x, 0)
multipliedMatrix@x <- as(multipliedMatrix@x, "double")

my_instructions = createGiottoInstructions(python_path = 'D:/Tools/Softwares/Anaconda/envs/dl_transformer')
my_giotto_object = createGiottoObject(raw_exprs = MOB_raw, spatial_locs = mob_spatial, instructions = my_instructions)
my_giotto_object <- normalizeGiotto(gobject = my_giotto_object)

my_giotto_object <- calculateHVG(gobject = my_giotto_object)


cell_annotations = as.vector(mapped$cellType)

rank_matrix = makeSignMatrixRank(sc_matrix = multipliedMatrix, sc_cluster_ids = cell_annotations)

sign_matrix = makeSignMatrixDWLSfromMatrix(matrix = multipliedMatrix, cell_type_vector = cell_annotations, sign_gene = row.names(multipliedMatrix))

my_giotto_object = doHclust(my_giotto_object, k = 4, name = 'hier_clus')

my_giotto_object = runDWLSDeconv(gobject = my_giotto_object, n_cell=20, sign_matrix = sign_matrix, cluster_column = "hier_clus")

DWLS_output = as.matrix(my_giotto_object@spatial_enrichment$DWLS)

write.csv(DWLS_output, file = "D:/Projects/CARD-master/data/画饼图/MOB/mob_dwls.csv", row.names = TRUE, col.names = TRUE)

colors = c("#cfa6b2","#fcd4c7","#6d889f","#f47c50","#a8af7f","#7FC97F","#BEAED4",
           "#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77","#D95F02",
           "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D")

row.names(DWLS_output) <- row.names(mob_spatial)

new_DWLS_output <- DWLS_output[,c(2,3,4,5,6)]

new_DWLS_output <- apply(new_DWLS_output, 2, as.numeric)

p1 <- CARD.visualize.pie(proportion = DWLS_output, 
                         spatial_location = mob_spatial, colors = colors, radius = 0.49)
