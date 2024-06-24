import torch
import numpy as np
from torch.utils.data import Dataset, DataLoader
import scanpy as sc
from utils.makeGraph import ST_preprocess, intra_exp_adj
from tqdm import tqdm
import scipy.spatial.distance as sp_distance
import pandas as pd
import os
import csv

# class ActiveDataset(Dataset):
#     def __init__(self, data_path, node_num, scale, prev_data, seed=None, k=6, neighbor_method="MNN"):
#         super(ActiveDataset, self).__init__()
#         self.data = sc.read_h5ad(data_path)
#         self.node_num = node_num
#         self.seed = seed
#         self.scale = scale
#         self.k = k
#         self.neighbor_method = neighbor_method
#         self.prev_data = prev_data
#         if self.seed:
#             np.random.seed(self.seed)
#             shuffled_indices = np.random.permutation(self.data.n_obs)
#             self.data = self.data[shuffled_indices, :]
#         self.x, self.y, self.ex_adjs = self.build_graph()
#
#     def __len__(self):
#         return len(self.x)
#
#     def __getitem__(self, idx):
#         return self.x[idx], self.y[idx], self.ex_adjs[idx]
#
#     def build_graph(self):
#         node_x_ls = []
#         node_y_ls = []
#         adj_ls = []
#         ori_data = self.data.copy()
#         pre_data = ST_preprocess(ori_data, scale=self.scale)
#         self.prev_data.var = pre_data.var
#
#         pre_data = sc.concat([pre_data, self.prev_data])
#         pre_data.obs_names = range(len(pre_data.obs_names))
#         num_graphs = int(len(pre_data) / self.node_num)
#         #         x = self.data.X
#         #         y = np.array(self.data.obs)[:,:-1]
#         for i in tqdm(range(num_graphs), desc='Generating pseudo-graphs'):
#             #         for i in range(num_graphs):
#             node = pre_data[i * 200: (i + 1) * 200]
#             node_x, node_y = node.X, np.array(node.obs)[:, :-1]
#             ex_adj = intra_exp_adj(node, dist_method="cosine", corr_dist_neighbors=self.k,
#                                    PCA_dimensionality_reduction=False,
#                                    find_neighbor_method=self.neighbor_method)
#             ex_adj = np.array(ex_adj)
#             node_x, node_y, ex_adj = torch.tensor(node_x).float(), torch.tensor(node_y).float(), torch.tensor(
#                 ex_adj).float()
#             node_x_ls.append(node_x)
#             node_y_ls.append(node_y)
#             adj_ls.append(ex_adj)
#         return node_x_ls, node_y_ls, adj_ls

class PseudoDataset(Dataset):
    def __init__(self, data_path, node_num, scale, seed=None, k=6, neighbor_method="MNN"):
        super(PseudoDataset, self).__init__()
        self.data = sc.read_h5ad(data_path)
        self.node_num = node_num
        self.seed = seed
        self.scale = scale
        self.k = k
        self.neighbor_method = neighbor_method
        if self.seed:
            np.random.seed(self.seed)
            shuffled_indices = np.random.permutation(self.data.n_obs)
            self.data = self.data[shuffled_indices, :]
        self.x, self.y, self.ex_adjs = self.build_graph()

    def __len__(self):
        return len(self.x)

    def __getitem__(self, idx):
        return self.x[idx], self.y[idx], self.ex_adjs[idx]

    def build_graph(self):
        node_x_ls = []
        node_y_ls = []
        adj_ls = []
        ori_data = self.data.copy()
        pre_data = ST_preprocess(ori_data, scale=self.scale)
        num_graphs = int(len(self.data.X) / self.node_num)
        #         x = self.data.X
        #         y = np.array(self.data.obs)[:,:-1]
        for i in tqdm(range(num_graphs), desc='Generating pseudo-graphs'):
            #         for i in range(num_graphs):
            node = pre_data[i * self.node_num: (i + 1) * self.node_num]
            node_x, node_y = node.X, np.array(node.obs)[:, :-1]
            ex_adj = intra_exp_adj(node, dist_method="cosine", corr_dist_neighbors=self.k,
                                   PCA_dimensionality_reduction=False,
                                   find_neighbor_method=self.neighbor_method)
            ex_adj = np.array(ex_adj)
            node_x, node_y, ex_adj = torch.tensor(node_x).float(), torch.tensor(node_y).float(), torch.tensor(
                ex_adj).float()
            node_x_ls.append(node_x)
            node_y_ls.append(node_y)
            adj_ls.append(ex_adj)
        return node_x_ls, node_y_ls, adj_ls


class StDataset(Dataset):
    def __init__(self, data_path, location_path, hvg, scale, pseudo_st_path, spatial_dist, k=6, neighbor_method="MNN", marker_path=None):
        super(StDataset, self).__init__()
        self.data = sc.read_h5ad(data_path)
        self.pseudo_st_data = sc.read_h5ad(pseudo_st_path)
        self.k = k
        self.spatial_dist = spatial_dist
        self.scale = scale
        self.hvg = hvg
        self.neighbor_method = neighbor_method
        self.marker_path = marker_path
        self.loc_mat = pd.read_csv(location_path, sep='\t')
        self.sp_adjs = self.build_dist_adj()
        self.x, self.ex_adjs = self.build_graph()
        
    def __len__(self):
        return len(self.x)

    def __getitem__(self, idx):
        return self.x[idx], self.ex_adjs[idx], self.sp_adjs[idx]

    def build_graph(self):
        node_x_ls = []
        node_y_ls = []
        adj_ls = []
        ori_data = self.data.copy()
        pseudo_data = self.pseudo_st_data.copy()
        # pre_data = ST_preprocess(ori_data, highly_variable_genes=self.hvg, scale=self.scale)
        if self.marker_path:
            # Define an empty list to store the cell marker values
            cell_markers = []

            # Specify the directory path
            folder_dir = self.marker_path

            # Iterate over each file in the directory
            for filename in os.listdir(folder_dir):
                if filename.endswith('.csv'):
                    # Construct the full file path
                    file_path = os.path.join(folder_dir, filename)

                    # Read the CSV file
                    with open(file_path, 'r') as file:
                        reader = csv.DictReader(file)

                        # Iterate over each row in the CSV file
                        for row in reader:
                            # Extract the cell marker value from the 'cell marker' column
                            cell_marker = row['Cell marker']

                            # Add the cell marker value to the list if it's not already present
                            if cell_marker not in cell_markers:
                                cell_markers.append(cell_marker)
            pre_data = ST_preprocess(ori_data, highly_variable_genes=False, scale=self.scale)
            st_genes = pre_data.var_names
            common_genes = set(st_genes).intersection(set(cell_markers))
            pre_data = pre_data[:, list(common_genes)]
            print(f"Select {len(common_genes)} marker genes")
        else:
            pre_data = ST_preprocess(ori_data, highly_variable_genes=self.hvg, scale=self.scale)
        st_genes = pre_data.var_names
        pseudo_spots_genes = pseudo_data.var_names
        if not all(gene in pseudo_spots_genes for gene in st_genes):
            print("Not all the genes in ST recorded in pseudo spots")
            common_genes = set(st_genes).intersection(set(pseudo_spots_genes))
            pseudo_data = pseudo_data[:, list(common_genes)]
            pre_data = pre_data[:, list(common_genes)]
        else:
            pseudo_data = pseudo_data[:, pre_data.var_names]
        path = "./pseudo_graph_tmp.h5ad"
        pseudo_data.write(path)
        print(f"Select {len(pre_data.var_names)} HVGs")
        #         x = self.data.X
        #         y = np.array(self.data.obs)[:,:-1]
        #         for i in tqdm(range(num_graphs), desc='Generating pseudo-graphs'):
        #         for i in range(num_graphs):
        node = pre_data
        node_x = node.X
        ex_adj = intra_exp_adj(node, dist_method="cosine", corr_dist_neighbors=self.k, PCA_dimensionality_reduction=False,
                               find_neighbor_method=self.neighbor_method)
        ex_adj = np.array(ex_adj)
        node_x, ex_adj = torch.tensor(node_x).float(), torch.tensor(ex_adj).float()
        node_x_ls.append(node_x)
        adj_ls.append(ex_adj)
        return node_x_ls, adj_ls

    def build_dist_adj(self):
        # Define the distance threshold
        distance_threshold = self.spatial_dist
        gd_loc = self.loc_mat.copy()
        gd_loc = np.array(gd_loc)

        # Calculate pairwise distances
        distances = sp_distance.squareform(sp_distance.pdist(gd_loc))

        # Create the adjacency matrix
        sp_adj = np.where((distances <= distance_threshold) & (distances > 0), 1, 0)
        sp_adj = torch.tensor(sp_adj).float()
        sp_adjs = [sp_adj]

        return sp_adjs