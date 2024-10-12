rm(list=ls())

library(spacexr)
library(CARD)

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
levels(cell_types) <- c("L2-3 IT CTX-1", "L4-5 IT CTX", "L5 IT CTX", "L5 NP CTX", "L5 PT CTX",
                        "L6 CT CTX", "L6 IT CTX", "L6b CTX", "Lamp5",
                        "Sst", "Vip")

reference_osm <- Reference(osm_sc_count, cell_types)

spaceRNA_osm <- SpatialRNA(osm_st_spatial, osm_st_count)

myRCTD <- create.RCTD(spaceRNA_osm, reference_osm)

myRCTD <- run.RCTD(myRCTD)

normalized_weights <- normalize_weights(myRCTD@results$weights)

rctd_out <- as.matrix(normalized_weights)

write.csv(rctd_out, file = "D:/Projects/CARD-master/data/画饼图/osmFISH/RCTD_output.csv")



# 
# rm(list=ls())
# 
# library(spacexr)
# library(Matrix)
# library(data.table)
# library(Seurat)
# library(SeuratDisk)
# library(rhdf5)
# 
# args<-commandArgs(T)
# snrna_path = "D:/Projects/CARD-master/data/starmap_sc_rna.h5ad"
# spatial_path = "D:/Projects/CARD-master/data/starmap_spatial.h5ad"
# celltype_final = "celltype"
# output_path = "D:/Projects/CARD-master/data/画饼图/starmap"
# 
# sc_obj <- LoadH5Seurat(snrna_path)
# diff_list <- list()
# c = 0
# for (i in seq_along(unique(sc_obj@meta.data[,celltype_final]))){
#   if(sum(sc_obj@meta.data[,celltype_final] == unique(sc_obj@meta.data[,celltype_final])[i]) < 25){
#     c = c+1
#     diff_list[[c]] <- unique(sc_obj@meta.data[,celltype_final])[i]
#     print(unique(sc_obj@meta.data[,celltype_final])[i])
#   }
# }
# for (i in diff_list){
#   sc_obj = sc_obj[,sc_obj@meta.data[,celltype_final]!=i]
# }
# sc_obj@meta.data[,celltype_final] <- as.factor(as.character(sc_obj@meta.data[,celltype_final]))
# ### Load in/preprocess your data, this might vary based on your file type
# print('prepare data')
# counts <- data.frame(sc_obj@assays$RNA@counts)
# colnames(counts) <- colnames(sc_obj)
# meta_data <- data.frame(sc_obj@meta.data[,celltype_final])
# cell_types <- meta_data[,1]
# names(cell_types) <- rownames(sc_obj@meta.data)
# cell_types <- as.factor(cell_types)
# nUMI_df <- data.frame(colSums(sc_obj@assays$RNA@counts))
# nUMI <- nUMI_df$colSums.sc_obj.assays.RNA
# names(nUMI) <- rownames(nUMI_df)
# 
# ### Create the Reference object
# reference <- Reference(counts, cell_types, nUMI)
# #> Warning in Reference(counts, cell_types, nUMI): Reference: nUMI does not match colSums of counts. If this is unintended, please correct this discrepancy. If this
# #>             is intended, there is no problem.
# spatial_obj <- LoadH5Seurat(spatial_path)
# coords <- data.frame(colnames(spatial_obj))
# colnames(coords) <- 'barcodes'
# coords$xcoord <- seq_along(colnames(spatial_obj))
# coords$ycoord <- seq_along(colnames(spatial_obj))
# counts <- data.frame(spatial_obj@assays$RNA@counts) # load in counts matrix
# colnames(counts) <- colnames(spatial_obj)
# coords <- data.frame(colnames(spatial_obj))
# colnames(coords) <- 'barcodes'
# coords$xcoord <- seq_along(colnames(spatial_obj))
# coords$ycoord <- seq_along(colnames(spatial_obj))
# rownames(coords) <- coords$barcodes; coords$barcodes <- NULL # Move barcodes to rownames
# nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI
# 
# ### Create SpatialRNA object
# puck <- SpatialRNA(coords, counts, nUMI)
# myRCTD <- create.RCTD(puck, reference, max_cores = 8)
# myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
# results <- myRCTD@results
# # normalize the cell type proportions to sum to 1.
# norm_weights = sweep(results$weights, 1, rowSums(results$weights), '/') 
# cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
# spatialRNA <- myRCTD@spatialRNA
# write.csv(norm_weights, paste0(output_path, '/RCTD_result.csv'))



