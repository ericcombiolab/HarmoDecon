import torch
from torch.utils.data import Dataset, DataLoader
from stdatasets.datasets import PseudoDataset, StDataset
from models.gmgat_model import GMGATModel
import scanpy as sc

import json

import random
import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt
from tqdm import tqdm
import matplotlib.patches as mpatches
import os

import warnings
warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser(description='Process JSON input')
parser.add_argument('--input', '-i', type=str, help='Path to the JSON file')
parser.add_argument('--epoch', '-e', type=str, help='Which epoch the state is')
parser.add_argument('--seed', '-s', type=str, help='Which seed the state is')
args = parser.parse_args()

with open(args.input) as file:
    config = json.load(file)

# gpu_id = config['gpu_id']
# proj_name = config['proj_name']
sub_name = config['sub_name']
real_st_path = config['real_st_path']
real_location_path = config['real_location_path']
try:
    pseudo_st_path = config['pseudo_st_path']
except:
    pseudo_st_path = None

try:
    hvg = config['hvg']
except:
    hvg = True

try:
    scale = config['scale']
except:
    scale = True

try:
    domain_classes = config['domain_classes']
except:
    domain_classes = 16

try:
    spatial_dist = config['spatial_dist']
except:
    spatial_dist = 1.5

try:
    gpu_id = config['gpu_id']
except:
    gpu_id = None

epoch = args.epoch
seed = int(args.seed)
model_path = f"./checkpoints/{sub_name}/{sub_name}_{epoch}_seed{seed}.model"


def plot_heatmap(loc, result, sub_name, save_dir):
    loc = loc.astype(int)
    width = loc['x'].max() - loc['x'].min()
    height = loc['y'].max() - loc['y'].min()
    cell_types = result.columns.to_list()

    for celltype in cell_types:
        celltype_ind = np.where(celltype == result.columns)[0][0]
        #         print(celltype_ind)
        c_ind = celltype_ind
        fig, axes = plt.subplots(ncols=1, nrows=1, figsize=(int(8 * width / height), 8))
        s = 20
        axes.scatter(x=loc['x'], y=loc['y'], c=result.iloc[:, c_ind].to_list(), s=s, marker='s')
        axes.scatter(x=15, y=60, c=0, s=s, marker='s', )
        axes.set_xlim(loc['x'].min() - 1, loc['x'].max() + 1)
        axes.set_ylim(loc['y'].min() - 1, loc['y'].max() + 1)
        axes.axis('off')
        title = f'{sub_name}_{celltype}'
        axes.set_title(title)

        title = title.replace('/', '_')

        plt.savefig(os.path.join(save_dir, title), bbox_inches='tight',dpi=300)
        plt.show()


def draw_pie(dist, xpos, ypos, size, colors, ax):
    cumsum = np.cumsum(dist)
    cumsum = cumsum / cumsum[-1]
    pie = [0] + cumsum.tolist()
    i = 0
    for r1, r2 in zip(pie[:-1], pie[1:]):
        angles = np.linspace(2 * np.pi * r1, 2 * np.pi * r2, num=100)
        x = [0] + np.cos(angles).tolist()
        y = [0] + np.sin(angles).tolist()

        xy = np.column_stack([x, y])
        ax.scatter([xpos], [ypos], marker=xy, s=size, c=colors[i], edgecolors='none')
        i += 1

    return ax


def plot_frac_results(predict, cell_type_list, coordinates, file_name=None, point_size=1000, size_coefficient=0.0009,
                      if_show=True, color_dict=None):
    coordinates.columns = ['coor_X', 'coor_Y']
    labels = cell_type_list
    if color_dict != None:
        colors = []
        for i in cell_type_list:
            colors.append(color_dict[i])
    else:
        if len(labels) <= 10:
            colors = plt.rcParams["axes.prop_cycle"].by_key()['color'][:len(labels)]
        else:
            import matplotlib
            color = plt.get_cmap('rainbow', len(labels))
            colors = []
            for x in color([range(len(labels))][0]):
                colors.append(matplotlib.colors.to_hex(x, keep_alpha=False))

    str_len = 0
    for item in cell_type_list:
        str_len = max(str_len, len(item))
    extend_region = str_len / 15 + 1

    fig, ax = plt.subplots(figsize=(len(coordinates['coor_X'].unique()) * point_size * size_coefficient + extend_region,
                                    len(coordinates['coor_Y'].unique()) * point_size * size_coefficient))

    for i in tqdm(range(predict.shape[0]), desc="Plotting pie plots:"):
        ax = draw_pie(predict[i], coordinates['coor_X'].values[i], coordinates['coor_Y'].values[i],
                      size=point_size, ax=ax, colors=colors)

    for spine in ax.spines.values():
        spine.set_edgecolor('lightgrey')  # Set edge line color of all spines to grey

    patches = [mpatches.Patch(color=colors[i], label="{:s}".format(labels[i])) for i in range(len(colors))]
    fontsize = max(predict.shape[0] / 100, 10)
    fontsize = min(fontsize, 30)
    ax.legend(handles=patches, fontsize=fontsize, bbox_to_anchor=(1, 1), loc="upper left")
    plt.axis("equal")
    plt.xticks([])
    plt.yticks([])
    plt.tight_layout()
    if file_name != None:
        plt.savefig(file_name,dpi=300,)
    if if_show == True:
        plt.show()
    plt.close('all')


def plot_scatter_by_type(predict, cell_type_list, coordinates, point_size=400, size_coefficient=0.0009, file_path=None,
                         if_show=True):
    coordinates.columns = ['coor_X', 'coor_Y']

    for i in tqdm(range(len(cell_type_list)), desc="Plotting cell type scatter plot:"):

        fig, ax = plt.subplots(figsize=(len(coordinates['coor_X'].unique()) * point_size * size_coefficient + 1,
                                        len(coordinates['coor_Y'].unique()) * point_size * size_coefficient))
        cm = plt.cm.get_cmap('Reds')
        ax = plt.scatter(coordinates['coor_X'], coordinates['coor_Y'], s=point_size, vmin=0, vmax=1, c=predict[:, i],
                         cmap=cm)

        cbar = plt.colorbar(ax, fraction=0.05)
        labelsize = max(predict.shape[0] / 100, 10)
        labelsize = min(labelsize, 30)
        cbar.ax.tick_params(labelsize=labelsize)
        plt.axis("equal")
        plt.xticks([])
        plt.yticks([])
        plt.xlim(coordinates['coor_X'].min() - 0.5, coordinates['coor_X'].max() + 0.5)
        plt.ylim(coordinates['coor_Y'].min() - 0.5, coordinates['coor_Y'].max() + 0.5)
        plt.tight_layout()
        if file_path != None:
            name = cell_type_list[i].replace('/', '_')
            plt.savefig(file_path + '/{}.jpg'.format(name), dpi=300, bbox_inches='tight')
        if if_show == True:
            plt.show()
        plt.close('all')



if __name__ == "__main__":
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)

    if not pseudo_st_path:
        pseudo_st_path = './pseudo_spots_tmp.h5ad'
    print(f"Obtain cell types label referring to {pseudo_st_path}")
    pseudo_spots_adata = sc.read_h5ad(pseudo_st_path)

    num_cell_types = len(pseudo_spots_adata.obs.columns) - 1

    real_data = StDataset(data_path=real_st_path, location_path=real_location_path, pseudo_st_path=pseudo_st_path
                          , hvg=hvg, scale=scale, spatial_dist=spatial_dist, k=int(6*spatial_dist/1.5))

    num_hvg = len(real_data[0][0][-1])

    real_loader = DataLoader(real_data, batch_size=1)

    model = GMGATModel(block_type="GCN", num_heads=2, st_encoder_in_channels=[num_hvg, 256, 256]
                       , num_classes=domain_classes, num_cell_type=num_cell_types)

    model.load_state_dict(torch.load(model_path))

    if gpu_id:
        model = model.cuda()

    for x, ex_adj, sp_adj in real_loader:
        model.eval()
        ex_adj, sp_adj = torch.squeeze(ex_adj), torch.squeeze(sp_adj)
        ex_adj = ex_adj.fill_diagonal_(1.)
        sp_adj = sp_adj.fill_diagonal_(1.)
        if gpu_id:
            x, ex_adj, sp_adj = x.cuda(), ex_adj.cuda(), sp_adj.cuda()

        feature_dict, cls, ratio = model(x, [ex_adj, sp_adj], mode="st")

        ex_adj_out = ratio[0].cpu().detach().numpy()
        ex_adj_out = np.squeeze(ex_adj_out)


        sp_adj_out = ratio[1].cpu().detach().numpy()
        sp_adj_out = np.squeeze(sp_adj_out)

        header = pseudo_spots_adata.obs.columns[:-1]

        ex_adj_out = pd.DataFrame(ex_adj_out, columns=header)  # Create a DataFrame with the tensor data and header

        ex_adj_out = ex_adj_out[sorted(ex_adj_out.columns)]

        out_path = f"./outputs/{sub_name}"
        if not os.path.exists(out_path):
            os.makedirs(out_path)
            print("Directory created:", out_path)
        else:
            print("Directory already exists:", out_path)

        ex_adj_out.to_csv(f'{out_path}/ex_adj_out_{model_path.split("/")[-1]}.csv', index=False)  # Save the DataFrame as a CSV file without the index

        sp_adj_out = pd.DataFrame(sp_adj_out, columns=header)  # Create a DataFrame with the tensor data and header

        sp_adj_out = sp_adj_out[sorted(sp_adj_out.columns)]

        sp_adj_out.to_csv(f'{out_path}/sp_adj_out_{model_path.split("/")[-1]}.csv', index=False)  # Save the DataFrame as a CSV file without the index

        mean_out = (ex_adj_out + sp_adj_out) / 2

        mean_out.to_csv(f'{out_path}/mean_out_{model_path.split("/")[-1]}.csv', index=False)

        fig_path = f"./figures/{sub_name}"
        if not os.path.exists(fig_path):
            os.makedirs(fig_path)
            print("Directory created:", fig_path)
        else:
            print("Directory already exists:", fig_path)

        loc = pd.read_csv(real_location_path, sep='\t')

        print("Plotting heatmaps...")

        plot_heatmap(loc=loc, result=mean_out, sub_name=sub_name, save_dir=fig_path)

        print("Plotting pie charts...")

        plot_frac_results(np.array(mean_out), mean_out.columns, loc, point_size=300,
                          size_coefficient=0.0009
                          , file_name=fig_path + f'/{sub_name}_pie_plot.jpg')
