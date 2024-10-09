import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from tqdm import tqdm
import random

import multiprocessing

import warnings

warnings.filterwarnings("ignore")

def generate_a_spot_poisson(adata_scrna, lam, max_cell_types_in_spot, library):
    cell_num = np.random.poisson(lam=lam) + 1
    cell_type_num = random.randint(1, max_cell_types_in_spot)
    cell_type_list_selected = np.random.choice(adata_scrna.obs['celltype'].value_counts().keys(), size=cell_type_num,
                                               replace=False)

    picked_cell_type = np.unique(np.random.choice(cell_type_list_selected, size=cell_num), return_counts=True)
    picked_cells = [np.random.choice(library[picked_cell_type[0][i]], picked_cell_type[1][i], replace=False) for i in
                    range(picked_cell_type[0].size)]
    picked_cells = [x for xs in picked_cells for x in xs]
    return adata_scrna[picked_cells]

def pseudo_spot_generation(sc_exp,
                           spot_num,
                           lam,
                           max_cell_types_in_spot,
                           generation_method,
                           n_jobs=-1):
    if n_jobs > 0:
        pool = multiprocessing.Pool(processes=n_jobs)
    cell_types = sc_exp.obs['celltype'].unique()
    word_to_idx_celltype = {word: i for i, word in enumerate(cell_types)}

    cell_type_num = len(cell_types)

    generated_spots = []
    library = {i: sc_exp[sc_exp.obs['celltype'] == i].obs_names for i in sc_exp.obs['celltype'].unique()}

    if n_jobs > 0:
        args = [(sc_exp, lam, max_cell_types_in_spot, library) for i in range(spot_num)]
        generated_spots = pool.starmap(generate_a_spot_passion, tqdm(args, desc='Generating pseudo-spots'))
    else:
        # generated_spots = []
        # for f in tqdm(args, desc='Generating pseudo-spots'):
        #     generated_spots.append(generate_a_spot_passion(*f))
        for _ in tqdm(range(spot_num), desc='Generating pseudo-spots'):
            generated_spots.append(generate_a_spot_poisson(sc_exp, lam, max_cell_types_in_spot, library))

    pse_srt_table = np.zeros((spot_num, sc_exp.shape[1]), dtype=float)
    pse_fraction_table = np.zeros((spot_num, cell_type_num), dtype=float)
    for i in tqdm(range(spot_num), desc='Obtain labels'):
        one_spot = generated_spots[i]
        pse_srt_table[i] = one_spot.X.sum(axis=0)
        # for j in one_spot.obs.index:
        #     cell_type = one_spot.obs.loc[j, 'celltype']
        #     type_idx = word_to_idx_celltype[cell_type]
        #     pse_fraction_table[i, type_idx] += 1
        for j in one_spot.obs.index:
            cell_type = one_spot.obs.loc[j, 'celltype']
            if type(cell_type) == str:
                type_idx = word_to_idx_celltype[cell_type]
            else:
                type_idx = cell_type.map(word_to_idx_celltype)
            pse_fraction_table[i, type_idx] += 1
    pse_srt_table = pd.DataFrame(pse_srt_table, columns=sc_exp.var.index.values)
    adata_pse_srt = sc.AnnData(X=pse_srt_table.values)
    adata_pse_srt.obs.index = pse_srt_table.index
    adata_pse_srt.var.index = pse_srt_table.columns
    pse_fraction_table = pd.DataFrame(pse_fraction_table, columns=cell_types)
    pse_fraction_table['cell_num'] = pse_fraction_table.sum(axis=1)
    for i in pse_fraction_table.columns[:-1]:
        pse_fraction_table[i] = pse_fraction_table[i] / pse_fraction_table['cell_num']
    adata_pse_srt.obs = adata_pse_srt.obs.join(pse_fraction_table)
    return adata_pse_srt

# def generate_a_spot_passion(sc_exp,
#                             lam,
#                             max_cell_types_in_spot,
#                             generation_method,
#                             ):
#     if generation_method == 'cell':
#         cell_num = np.random.poisson(lam=lam) + 1
#         cell_list = list(sc_exp.obs.index.values)
#         picked_cells = random.choices(cell_list, k=cell_num)
#         return sc_exp[picked_cells]
#     elif generation_method == 'celltype':
#         cell_num = np.random.poisson(lam=lam) + 1
#         cell_type_list = list(sc_exp.obs['celltype'].unique())
#         cell_type_num = random.randint(1, max_cell_types_in_spot)
#
#         while (True):
#             cell_type_list_selected = random.choices(sc_exp.obs['celltype'].value_counts().keys(), k=cell_type_num)
#             if len(set(cell_type_list_selected)) == cell_type_num:
#                 break
#         sc_exp_filter = sc_exp[sc_exp.obs['celltype'].isin(cell_type_list_selected)]
#
#         picked_cell_type = random.choices(cell_type_list_selected, k=cell_num)
#         picked_cells = []
#         for i in picked_cell_type:
#             data = sc_exp[sc_exp.obs['celltype'] == i]
#             cell_list = list(data.obs.index.values)
#             picked_cells.append(random.sample(cell_list, 1)[0])
#         return sc_exp_filter[picked_cells]
#     else:
#         print('generation_method should be "cell" or "celltype" ')
#
#
# def pseudo_spot_generation(sc_exp,
#                            spot_num,
#                            lam,
#                            max_cell_types_in_spot,
#                            generation_method,
#                            n_jobs=-1
#                            ):
#     if n_jobs > 0:
#         pool = multiprocessing.Pool(processes=n_jobs)
#     cell_types = sc_exp.obs['celltype'].unique()
#
#     idx_to_word_celltype = {i: word for i, word in enumerate(cell_types)}
#
#     word_to_idx_celltype = {word: i for i, word in enumerate(cell_types)}
#
#     print(idx_to_word_celltype)
#
#     print(word_to_idx_celltype)
#
#     cell_type_num = len(sc_exp.obs['celltype'].unique())
#
#
#     args = [(sc_exp, lam, max_cell_types_in_spot, generation_method) for i in range(spot_num)]
#
#     if n_jobs > 0:
#         generated_spots = pool.starmap(generate_a_spot_passion, tqdm(args, desc='Generating pseudo-spots'))
#     else:
#         generated_spots = []
#         for f in tqdm(args, desc='Generating pseudo-spots'):
#             generated_spots.append(generate_a_spot_passion(*f))
#
#
#     # pseudo_spots = []
#     pseudo_spots_table = np.zeros((spot_num, sc_exp.shape[1]), dtype=float)
#     pseudo_fraction_table = np.zeros((spot_num, cell_type_num), dtype=float)
#     for i in tqdm(range(spot_num), desc='Obtain labels'):
#         one_spot = generated_spots[i]
#         # pseudo_spots.append(one_spot)
#         pseudo_spots_table[i] = one_spot.X.sum(axis=0)
#         for j in one_spot.obs.index:
#             cell_type = one_spot.obs.loc[j, 'celltype']
#             if type(cell_type) == str:
#                 type_idx = word_to_idx_celltype[cell_type]
#             else:
#                 type_idx = cell_type.map(word_to_idx_celltype)
#             pseudo_fraction_table[i, type_idx] += 1
#     pseudo_spots_table = pd.DataFrame(pseudo_spots_table, columns=sc_exp.var.index.values)
#     pseudo_spots = sc.AnnData(X=pseudo_spots_table.iloc[:, :].values)
#     pseudo_spots.obs.index = pseudo_spots_table.index[:]
#     pseudo_spots.var.index = pseudo_spots_table.columns[:]
#     type_list = [idx_to_word_celltype[i] for i in range(cell_type_num)]
#     pseudo_fraction_table = pd.DataFrame(pseudo_fraction_table, columns=type_list)
#     pseudo_fraction_table['cell_num'] = pseudo_fraction_table.sum(axis=1)
#     for i in pseudo_fraction_table.columns[:-1]:
#         pseudo_fraction_table[i] = pseudo_fraction_table[i] / pseudo_fraction_table['cell_num']
#     pseudo_spots.obs = pseudo_spots.obs.join(pseudo_fraction_table)
#
#
#     return pseudo_spots

if __name__=="__main__":
    # pseudo_path = "/home/comp/cszrwang/project/Geneformer_test/STARmap/toy/starmap_sc_rna.h5ad"
    # pseudo_path = "./data/starmap_sc_rna.h5ad"
    pseudo_path = "D:/Data/BGI/YANG_Chao/ESCC_Public_anno.h5ad"
    pseudo_adata = sc.read_h5ad(pseudo_path)
    # ST_path = "/home/comp/cszrwang/project/Geneformer_test/STARmap/toy/starmap_spatial_ST.h5ad"
    # ST_path = "./data/starmap_spatial_ST.h5ad"
    ST_path = "D:/Data/BGI/YANG_Chao/C03031E2_bin500.h5ad"
    ST_adata = sc.read_h5ad(ST_path)

    ST_genes = ST_adata.var.index.values
    pseudo_genes = pseudo_adata.var.index.values
    common_genes = set(ST_genes).intersection(set(pseudo_genes))
    ST_adata_filter = ST_adata[:,list(common_genes)]

    pseudo_adata_filter = pseudo_adata[:,list(common_genes)]

    spots_toy = pseudo_spot_generation(sc_exp=pseudo_adata_filter, spot_num=50000, lam=5,generation_method='celltype',
                           max_cell_types_in_spot=4, n_jobs = -1)
    # spots_toy.write('/home/comp/cszrwang/project/Geneformer_test/STARmap/toy/pseudo_50000.h5ad')
    spots_toy.write('./data/pseudo.h5ad')