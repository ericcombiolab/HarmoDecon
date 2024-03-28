import torch.nn as nn
import numpy as np
import torch
import wandb
from torch.optim import Adam
from models.modules.gmgat import VAE
from models.modules.layers import GRL
from models.dnnTokenizer import DnnTokenizer
# from utils.lossFunctions import LossFunctions
# from utils.util import (
#     summary_bin_list_from_batch,
#     get_long_contig_logits_from_batch,
#     log_ag_graph,
#     log_knn_graph,
#     evaluate,
#     log_tsne_figure,
#     refine_gmm,
# )


class GMGATModel(nn.Module):
    def __init__(
            self,
            st_encoder_in_channels=[512, 128, 256],
            st_encoder_out_channels=[128, 256, 512],
            # ex_encoder_in_channels=[512, 128, 256],
            # ex_encoder_out_channels=[128, 256, 512],
            gaussian_size=512,
            num_classes=10,
            # lr=0.0001,
            num_heads=4,
            num_blocks=3,
            # pseudo_w_rec=None,
            # domain_w_ce=None,
            # st_sp_w_rec=None,
            # st_ex_w_rec=None,
            # st_w_mse=None,
            block_type="transformer",
            log_path="",
            # k=6,
            # use_gmm=False,
            use_bias=True,
            dropout=0.1,
            input_dimension=882,
            num_cell_type=15,
            *args,
            **kwargs,
    ):
        super().__init__()
        self.GMGAT = VAE(
            encoder_in_channels=st_encoder_in_channels,
            encoder_out_channels=st_encoder_out_channels,
            num_heads=num_heads,
            y_dim=num_classes,
            latent_dim=gaussian_size,
            dropout=dropout,
            num_blocks=num_blocks,
            use_bias=use_bias,
            block_type=block_type,
        )
        self.domain_model = nn.Sequential(
            GRL(),
            nn.Linear(512, 40),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(40, 2),
            nn.LogSoftmax()
        )
        self.devonv_model = nn.Sequential(
            GRL(),
            nn.Linear(512, 256),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(256, num_cell_type),
            nn.LogSoftmax()
        )
        self.tokenizer = DnnTokenizer(input_dimension=input_dimension)
        self.log_path = log_path

    def forward(self, node, adj):
        self.tokenizer
        pass

    def gmvae_loss(self, data, out_net, w_cat, w_gauss, w_rec):
        z, data_recon = out_net["gaussian"], out_net["reconstruct_graph"]
        logits, prob_cat = out_net["logits"], out_net["prob_cat"]
        y_mu, y_var = out_net["y_mean"], out_net["y_var"]
        mu, var = out_net["mean"], out_net["var"]

        loss_rec = self.losses.reconstruction_graph_loss(data, data_recon)
        loss_gauss = self.losses.gaussian_loss(z, mu, var, y_mu, y_var)
        loss_cat = -self.losses.entropy(logits, prob_cat) - np.log(0.1)
        loss_total = w_rec * loss_rec + w_gauss * loss_gauss + w_cat * loss_cat

        predicted_clusters = prob_cat.argmax(-1)
        highest_probs = prob_cat.max(-1).values
        loss_dict = {
            "total": loss_total,
            "predicted_clusters": predicted_clusters,
            "reconstruction": loss_rec * w_rec,
            "gaussian": loss_gauss,
            "categorical": loss_cat,
            "highest_prob": highest_probs,
            "logits": logits,
        }
        return loss_dict

    def training_step(self, batch, batch_idx):
        ag_opt, knn_opt = self.optimizers()
        knn_attributes = batch["knn_feature"]
        ag_pe_attributes = batch["ag_pe_feature"]
        ag_adj_matrix = batch["ag_adj_matrix"]
        knn_adj_matrix = batch["knn_adj_matrix"]
        mask_matrix = batch["ag_mask_matrix"]
        ag_output_dict = self.AGGMGAT(
            h=ag_pe_attributes,
            adj_matrix=ag_adj_matrix,
            mask_matrix=mask_matrix,
        )
        knn_output_dict = self.KNNGMGAT(
            h=knn_attributes,
            adj_matrix=knn_adj_matrix,
            mask_matrix=mask_matrix,
        )

        ##########################
        # Optimize AGGMGAT #
        ##########################
        knn_logits = knn_output_dict["logits"].detach()
        long_contig_knn_logits = get_long_contig_logits_from_batch(batch, knn_logits)
        long_contig_ag_logits = get_long_contig_logits_from_batch(batch, ag_output_dict["logits"])
        ag_loss_dict = self.gmvae_loss(ag_adj_matrix, ag_output_dict, self.ag_w_cat, self.ag_w_gauss, self.ag_w_rec)
        ag_reconstruction_loss = ag_loss_dict["reconstruction"]
        ag_gaussian_loss = ag_loss_dict["gaussian"]
        ag_categorical_loss = ag_loss_dict["categorical"]
        ag_cross_entropy_loss = nn.CrossEntropyLoss()(
            long_contig_ag_logits,
            long_contig_knn_logits,
        )
        ag_loss = ag_loss_dict["total"] + self.ag_w_ce * ag_cross_entropy_loss
        ag_opt.zero_grad()
        self.manual_backward(ag_loss)
        ag_opt.step()

        ##########################
        # Optimize KNNGMGAT #
        ##########################
        ag_output_dict = self.AGGMGAT(
            h=ag_pe_attributes,
            adj_matrix=ag_adj_matrix,
            mask_matrix=mask_matrix,
        )  # use updated AGGMGAT to recompute the logits.
        ag_logits = ag_output_dict["logits"].detach()
        long_contig_ag_logits = get_long_contig_logits_from_batch(batch, ag_logits)
        knn_loss_dict = self.gmvae_loss(knn_adj_matrix, knn_output_dict, self.knn_w_cat, self.knn_w_gauss,
                                        self.knn_w_rec)
        knn_reconstruction_loss = knn_loss_dict["reconstruction"]
        knn_gaussian_loss = knn_loss_dict["gaussian"]
        knn_categorical_loss = knn_loss_dict["categorical"]
        knn_cross_entropy_loss = nn.CrossEntropyLoss()(
            long_contig_ag_logits,
            long_contig_knn_logits,
        )
        knn_loss = knn_loss_dict["total"] + self.knn_w_ce * knn_cross_entropy_loss
        knn_opt.zero_grad()
        self.manual_backward(knn_loss)
        knn_opt.step()

        # logging loss curve
        loss = ag_loss + knn_loss
        self.log("train/loss", loss, on_step=False, on_epoch=True, prog_bar=False)
        self.log("train/ag_reconstruction_loss", ag_reconstruction_loss, on_step=False, on_epoch=True, prog_bar=False)
        self.log("train/ag_gaussian_loss", ag_gaussian_loss, on_step=False, on_epoch=True, prog_bar=False)
        self.log("train/ag_categorical_loss", ag_categorical_loss, on_step=False, on_epoch=True, prog_bar=False)
        self.log("train/ag_crossentropy_loss", ag_cross_entropy_loss, on_step=False, on_epoch=True, prog_bar=False)
        self.log("train/knn_reconstruction_loss", knn_reconstruction_loss, on_step=False, on_epoch=True, prog_bar=False)
        self.log("train/knn_gaussian_loss", knn_gaussian_loss, on_step=False, on_epoch=True, prog_bar=False)
        self.log("train/knn_categorical_loss", knn_categorical_loss, on_step=False, on_epoch=True, prog_bar=False)
        self.log("train/knn_crossentorpy_loss", knn_cross_entropy_loss, on_step=False, on_epoch=True, prog_bar=False)
        return {"loss": loss}

    def test_step(self):
        # TODO: update test step function here.
        pass

    def validation_step(self, batch, batch_idx):
        # TODO: add the visualization of ag and knn two lines output.
        attributes = batch["knn_feature"]
        mask_matrix = batch["ag_mask_matrix"]
        knn_adj_matrix = batch["knn_adj_matrix"]
        output_dict = self.KNNGMGAT(
            h=attributes,
            adj_matrix=knn_adj_matrix,
            mask_matrix=mask_matrix,
        )
        prob_cat = output_dict["prob_cat"]
        latent = output_dict["gaussian"]
        bin_tensor = prob_cat.argmax(-1)
        gd_bin_list, result_bin_list, non_labeled_id_list = summary_bin_list_from_batch(batch, bin_tensor)

        # Compute metrics.
        precision, recall, ARI, F1 = evaluate(
            gd_bin_list=gd_bin_list,
            result_bin_list=result_bin_list,
            non_labeled_id_list=non_labeled_id_list,
            unclassified=0,
        )

        # logging graph for visualization.
        contig_id_list = [int(id) for index, id in enumerate(torch.squeeze(batch["id"]))]
        plotting_contig_list = contig_id_list[:self.plot_graph_size]
        gd_ag_graph_path, result_ag_graph_path = log_ag_graph(
            plotting_graph_size=self.plot_graph_size,
            processed_zarr_dataset_path=self.processed_zarr_dataset_path,
            plotting_contig_list=plotting_contig_list,
            log_path=self.log_path,
            gd_bin_list=gd_bin_list,
            result_bin_list=result_bin_list,
        )
        gd_knn_graph_path, result_knn_graph_path = log_knn_graph(
            plotting_graph_size=self.plot_graph_size,
            plotting_contig_list=plotting_contig_list,
            k=self.k,
            batch=batch,
            log_path=self.log_path,
            gd_bin_list=gd_bin_list,
            result_bin_list=result_bin_list,
        )

        # Visualize latent space.
        result_tsne_figure_path = log_tsne_figure(
            batch=batch,
            latent=torch.squeeze(latent),
            log_path=self.log_path,
        )
        self.log("val/acc", attributes.shape[0], on_step=False, on_epoch=True, prog_bar=False)
        self.log("val/precision", precision, on_step=False, on_epoch=True, prog_bar=False)
        self.log("val/recall", recall, on_step=False, on_epoch=True, prog_bar=False)
        self.log("val/F1", F1, on_step=False, on_epoch=True, prog_bar=False)
        self.log("val/ARI", ARI, on_step=False, on_epoch=True, prog_bar=False)
        wandb.log({"val/ground_truth_ag_subgraph": wandb.Image(gd_ag_graph_path)})
        wandb.log({"val/result_ag_subgraph": wandb.Image(result_ag_graph_path)})
        wandb.log({"val/ground_truth_knn_subgraph": wandb.Image(gd_knn_graph_path)})
        wandb.log({"val/result_knn_subgraph": wandb.Image(result_knn_graph_path)})
        wandb.log({"val/tsne_figure": wandb.Image(result_tsne_figure_path)})

    def configure_optimizers(self):
        ag_opt = Adam(self.AGGMGAT.parameters(), lr=self.lr)
        knn_opt = Adam(self.KNNGMGAT.parameters(), lr=self.lr)
        return ag_opt, knn_opt