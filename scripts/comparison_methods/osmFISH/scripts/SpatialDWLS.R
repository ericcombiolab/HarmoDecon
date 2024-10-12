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

# matrix <- assay(sce)
# mob_spatial <- read.table("D:/Tools/Data/MOB_ST/CARD/mob_cor.txt", header = TRUE)
# multipliedMatrix <- matrix * 100
# multipliedMatrix@x <- round(multipliedMatrix@x, 0)
# multipliedMatrix@x <- as(multipliedMatrix@x, "double")

my_instructions = createGiottoInstructions(python_path = 'D:/Tools/Softwares/Anaconda/envs/dl_transformer')

st_data <- createGiottoObject(
  raw_exprs = osm_st_count,
  instructions = my_instructions
)

sc_data <- createGiottoObject(
  raw_exprs = osm_sc_count,
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
# featgenes = gene_metadata[hvg == 'yes']$gene_ID
featgenes = gene_metadata$gene_ID
sc_data <- runPCA(gobject = sc_data, genes_to_use = featgenes, scale_unit = F)
signPCA(sc_data, genes_to_use = featgenes, scale_unit = F)
sc_data@cell_metadata$leiden_clus <- as.character(osm_sc_meta$cellType)
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
st_data <- runDWLSDeconv(st_data,sign_matrix = Sig_exp, n_cell = 200)

write.csv(st_data@spatial_enrichment$DWLS, paste0('D:/Projects/CARD-master/data/画饼图/osmFISH/SpatialDWLS_result.csv'))

# my_giotto_object = createGiottoObject(raw_exprs = osm_st_count, spatial_locs = star_st_spatial, instructions = my_instructions)
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