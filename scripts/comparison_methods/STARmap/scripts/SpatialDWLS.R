rm(list=ls())
library(reticulate)
library(rhdf5)
library(SingleCellExperiment)
library(Giotto)
library(CARD)

# load("D:/Tools/Data/MOB_ST/CARD/MOB.dge.sceset.RData")
# ### load the spatial data 
# load("D:/Tools/Data/MOB_ST/CARD/Rep12_MOB_count_matrix.RData")
# mob_sc_meta <- h5read("D:/Tools/Data/MOB_ST/CARD/mob_sc_celltype.h5ad", name='obs')
# mapped <- as.character(mob_sc_meta$cellType$categories[match(mob_sc_meta$cellType$codes, 1:length(mob_sc_meta$cellType$categories)-1)])
# mapped <- data.frame(cellType = mapped, cellID = sce@colData@rownames)

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

# matrix <- assay(sce)
# mob_spatial <- read.table("D:/Tools/Data/MOB_ST/CARD/mob_cor.txt", header = TRUE)
# multipliedMatrix <- matrix * 100
# multipliedMatrix@x <- round(multipliedMatrix@x, 0)
# multipliedMatrix@x <- as(multipliedMatrix@x, "double")

my_instructions = createGiottoInstructions(python_path = 'D:/Tools/Softwares/Anaconda/envs/dl_transformer')

st_data <- createGiottoObject(
  raw_exprs = star_st_count,
  instructions = my_instructions
)

sc_data <- createGiottoObject(
  raw_exprs = star_sc_count,
  instructions = my_instructions
)

st_data <- normalizeGiotto(gobject = st_data)
st_data <- calculateHVG(gobject = st_data)
gene_metadata = fDataDT(st_data)
featgenes = gene_metadata[hvg == 'yes']$gene_ID
st_data <- runPCA(gobject = st_data, genes_to_use = featgenes, scale_unit = F)
signPCA(st_data, genes_to_use = featgenes, scale_unit = F)
st_data <- runUMAP(st_data, dimensions_to_use = 1:10)
st_data <- createNearestNetwork(gobject = st_data, dimensions_to_use = 1:10, k = 15)
st_data <- doLeidenCluster(gobject = st_data, resolution = 0.4, n_iterations = 1000)

sc_data <- normalizeGiotto(gobject = sc_data)
sc_data <- calculateHVG(gobject = sc_data)
gene_metadata = fDataDT(sc_data)
featgenes = gene_metadata[hvg == 'yes']$gene_ID
sc_data <- runPCA(gobject = sc_data, genes_to_use = featgenes, scale_unit = F)
signPCA(sc_data, genes_to_use = featgenes, scale_unit = F)
sc_data@cell_metadata$leiden_clus <- as.character(star_sc_meta$cellType)
scran_markers_subclusters = findMarkers_one_vs_all(gobject = sc_data,
                                                   method = 'scran',
                                                   expression_values = 'normalized',
                                                   cluster_column = 'leiden_clus')
Sig_scran <- unique(scran_markers_subclusters$genes[which(scran_markers_subclusters$ranking <= 100)])
norm_exp<-2^(sc_data@norm_expr)-1
id<-sc_data@cell_metadata$leiden_clus
ExprSubset<-norm_exp[Sig_scran,]
Sig_exp<-NULL
for (i in unique(id)){
  Sig_exp<-cbind(Sig_exp,(apply(ExprSubset,1,function(y) mean(y[which(id==i)]))))
}
colnames(Sig_exp)<-unique(id)
st_data <- runDWLSDeconv(st_data,sign_matrix = Sig_exp, n_cell = 20)

write.csv(st_data@spatial_enrichment$DWLS, paste0('D:/Projects/CARD-master/data/画饼图/starmap/SpatialDWLS_result.csv'))

# my_giotto_object = createGiottoObject(raw_exprs = star_st_count, spatial_locs = star_st_spatial, instructions = my_instructions)
# my_giotto_object <- normalizeGiotto(gobject = my_giotto_object)
# 
# my_giotto_object <- calculateHVG(gobject = my_giotto_object)
# 
# gene_metadata = fDataDT(my_giotto_object)
# 
# featgenes = gene_metadata[hvg == 'yes']$gene_ID
# 
# cell_annotations = as.vector(star_sc_meta$cellType)
# 
# rank_matrix = makeSignMatrixRank(sc_matrix = star_sc_count, sc_cluster_ids = cell_annotations)
# 
# sign_matrix = makeSignMatrixDWLSfromMatrix(matrix = star_sc_count, cell_type_vector = cell_annotations, sign_gene = row.names(star_sc_count))
# 
# my_giotto_object = doHclust(my_giotto_object, k = 8, name = 'hier_clus')
# 
# my_giotto_object = runDWLSDeconv(gobject = my_giotto_object, n_cell=200, sign_matrix = sign_matrix, cluster_column = "hier_clus")
# 
# DWLS_output = as.matrix(my_giotto_object@spatial_enrichment$DWLS)

# write.csv(DWLS_output, file = "D:/Projects/CARD-master/data/画饼图/starmap/DWLS_output.csv", row.names = TRUE, col.names = TRUE)

# colors = c("#cfa6b2","#fcd4c7","#6d889f","#f47c50","#a8af7f","#7FC97F","#BEAED4",
#            "#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77","#D95F02",
#            "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D")
# 
# row.names(DWLS_output) <- row.names(mob_spatial)
# 
# new_DWLS_output <- DWLS_output[,c(2,3,4,5,6)]
# 
# new_DWLS_output <- apply(new_DWLS_output, 2, as.numeric)
# 
# p1 <- CARD.visualize.pie(proportion = DWLS_output, 
#                          spatial_location = mob_spatial, colors = colors, radius = 0.49)