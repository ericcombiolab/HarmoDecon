rm(list=ls())
library(reticulate)
library(rhdf5)
library(SingleCellExperiment)
library(Giotto)
library(CARD)

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

my_instructions = createGiottoInstructions(python_path = 'D:/Tools/Softwares/Anaconda/envs/dl_transformer')

st_data <- createGiottoObject(
  raw_exprs = mela_st_count,
  instructions = my_instructions
)

sc_data <- createGiottoObject(
  raw_exprs = mela_sc_count,
  instructions = my_instructions
)

st_data <- normalizeGiotto(gobject = st_data)
st_data <- calculateHVG(gobject = st_data)
gene_metadata = fDataDT(st_data)
featgenes = gene_metadata[hvg == 'yes']$gene_ID
st_data <- Giotto::runPCA(gobject = st_data, genes_to_use = featgenes, scale_unit = F)
signPCA(st_data, genes_to_use = featgenes, scale_unit = F)
st_data <- Giotto::runUMAP(st_data, dimensions_to_use = 1:10)
st_data <- createNearestNetwork(gobject = st_data, dimensions_to_use = 1:10, k = 15)
st_data <- doLeidenCluster(gobject = st_data, resolution = 0.4, n_iterations = 1000)

sc_data <- normalizeGiotto(gobject = sc_data)
sc_data <- calculateHVG(gobject = sc_data)
gene_metadata = fDataDT(sc_data)
# featgenes = gene_metadata[hvg == 'yes']$gene_ID
featgenes = gene_metadata$gene_ID
sc_data <- Giotto::runPCA(gobject = sc_data, genes_to_use = featgenes, scale_unit = F)
signPCA(sc_data, genes_to_use = featgenes, scale_unit = F)
sc_data@cell_metadata$leiden_clus <- as.character(mela_sc_meta$cellType)
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

write.csv(st_data@spatial_enrichment$DWLS, paste0('D:/Projects/CARD-master/data/画饼图/breast_1/SpatialDWLS_result.csv'))