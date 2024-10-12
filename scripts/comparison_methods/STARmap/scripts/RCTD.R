rm(list=ls())

library(spacexr)
library(CARD)

star_st_var <- h5read("D:/Projects/CARD-master/data/starmap_spatial.h5ad", name='var')

star_st_count <- h5read("D:/Projects/CARD-master/data/starmap_spatial.h5ad", name='X')

star_st_count <- as(star_st_count, "dgCMatrix")

star_st_spatial <- read.table("D:/Projects/CARD-master/data/combined_Locations.txt", header = TRUE)

colomns <- c('0x0', '0x1', '0x2', '0x3', '0x4', '0x5', '0x6', '0x7', '0x8',
             '0x9', '0x10', '0x11', '0x12', '0x13', '0x14', '0x15', '0x16',
             '0x17', '0x18', '1x0', '1x1', '1x2', '1x3', '1x4', '1x5', '1x6',
             '1x7', '1x8', '1x9', '1x10', '1x11', '1x12', '1x13', '1x14', '1x15',
             '1x16', '1x17', '1x18', '2x0', '2x1', '2x2', '2x3', '2x4', '2x5',
             '2x6', '2x7', '2x8', '2x9', '2x10', '2x11', '2x12', '2x13', '2x14',
             '2x15', '2x16', '2x17', '2x18', '3x0', '3x1', '3x2', '3x3', '3x4',
             '3x5', '3x6', '3x7', '3x8', '3x9', '3x10', '3x11', '3x12', '3x13',
             '3x14', '3x15', '3x16', '3x17', '3x18', '4x0', '4x1', '4x2', '4x3',
             '4x4', '4x5', '4x6', '4x7', '4x8', '4x9', '4x10', '4x11', '4x12',
             '4x13', '4x14', '4x15', '4x16', '4x17', '4x18', '5x0', '5x1',
             '5x2', '5x3', '5x4', '5x5', '5x6', '5x7', '5x8', '5x9', '5x10',
             '5x11', '5x12', '5x13', '5x14', '5x15', '5x16', '5x17', '6x0',
             '6x1', '6x2', '6x3', '6x4', '6x5', '6x6', '6x7', '6x8', '6x9',
             '6x10', '6x11', '6x12', '6x13', '6x14', '6x15', '6x16', '6x17',
             '6x18', '7x0', '7x1', '7x2', '7x3', '7x4', '7x5', '7x6', '7x7',
             '7x8', '7x9', '7x10', '7x11', '7x12', '7x13', '7x14', '7x15',
             '7x16', '7x17', '7x18', '8x0', '8x1', '8x2', '8x3', '8x4', '8x5',
             '8x6', '8x7', '8x8', '8x9', '8x10', '8x11', '8x12', '8x13', '8x14',
             '8x15', '8x16', '8x17', '8x18', '9x0', '9x1', '9x2', '9x3', '9x4',
             '9x5', '9x6', '9x7', '9x8', '9x9', '9x10', '9x11', '9x12', '9x13',
             '9x14', '9x15', '9x16', '9x17', '9x18')

rownames(star_st_spatial) <- colomns

star_sc_meta <- read.table("D:/Projects/CARD-master/data/forcard.txt", header = TRUE, sep = "\t")

rownames(star_st_count) <- as.character(star_st_var$`_index`)

colnames(star_st_count) <- colomns

star_sc_count <- h5read("D:/Projects/CARD-master/data/starmap_sc_rna.h5ad", name='X')

star_sc_var <- h5read("D:/Projects/CARD-master/data/starmap_sc_rna.h5ad", name='var')

star_sc_count <- as(star_sc_count, "dgCMatrix")

rownames(star_sc_count) <- as.character(star_sc_var$`_index`)

col_range <- 1:ncol(star_sc_count)

colnames(star_sc_count) <- col_range

# star_sc_meta$sampleInfo = "sample1"

star_sc_meta$cellID = as.character(1:ncol(star_sc_count))

cell_types <- star_sc_meta$cellType
cell_types <- factor(cell_types)
names(cell_types) <- star_sc_meta$cellID
# Error in check_cell_types(cell_types) :
#   Reference: levels(cell_types) contains a cell type with name containing prohibited character /. Please rename this cell type.
levels(cell_types) <- c("Astro", "Endo", "Excitatory L23", "Excitatory L4", "Excitatory L5",
                        "Excitatory L6", "Inhibitory Other", "Inhibitory Pvalb", "Inhibitory Sst",
                        "Inhibitory Vip", "Micro", "Neuron Other", "Olig", "Other", "Smc")

reference_star <- Reference(star_sc_count, cell_types)

spaceRNA_star <- SpatialRNA(star_st_spatial, star_st_count)

myRCTD <- create.RCTD(spaceRNA_star, reference_star)

myRCTD <- run.RCTD(myRCTD)

normalized_weights <- normalize_weights(myRCTD@results$weights)

rctd_out <- as.matrix(normalized_weights)

write.csv(rctd_out, file = "D:/Projects/CARD-master/data/画饼图/starmap/RCTD_output.csv")



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



