import torch
import torch.nn.functional as F
import torch.nn as nn
import numpy as np


class LossFunctions:
    def __init__(self):
        self.eps = 1e-8

    # def reconstruction_loss(A, Z):
    #     n, f = Z.size()
    #
    #     A = A.fill_diagonal_(1.)
    #
    #     # Calculate the loss
    #     dot_products = torch.matmul(Z, Z.t())
    #     sigmoid_dot_products = torch.sigmoid(dot_products)
    #
    #     loss = torch.sum(torch.matmul(A, torch.log(sigmoid_dot_products)) + torch.matmul((1 - A), torch.log(1 - sigmoid_dot_products)))
    #     loss = -loss / 2 / n
    #
    #     return loss

    def kl_loss(self, P, Q):
        kl_divergence = torch.sum(P * torch.log(P / Q))
        return kl_divergence

    def js_loss(self, P, Q):
        M = (P + Q) / 2
        kl_pm = self.kl_loss(P, M)
        kl_qm = self.kl_loss(Q, M)
        js_divergence = (kl_pm + kl_qm) / 2
        return js_divergence

    def self_entropy(self, tensor, dim=-1, epsilon=1e-8):
        # Compute the probability distribution
        # probabilities = F.softmax(tensor, dim=dim)

        tensor = torch.clamp(tensor, min=epsilon)

        # Calculate the entropy for each probability distribution
        entropies = -torch.sum(tensor * torch.log2(tensor), dim=dim)

        # Take the mean of the entropies
        self_entropy = torch.mean(entropies)

        return self_entropy

    # def mse_loss(A, B):
    #     mse = F.mse_loss(A, B)
    #     return mse

    def ce_loss(self, A, B):
        ce = F.cross_entropy(A, B)
        return ce

    # def mean_squared_error(self, real, predictions):
    #     """Mean Squared Error between the true and predicted outputs
    #        loss = (1/n)*Σ(real - predicted)^2
    #
    #     Args:
    #         real: (array) corresponding array containing the true labels
    #         predictions: (array) corresponding array containing the predicted labels
    #
    #     Returns:
    #         output: (array/float) depending on average parameters the result will be the mean
    #                               of all the sample losses or an array with the losses per sample
    #     """
    #     loss = (real - predictions).pow(2)
    #     return loss.sum(-1).mean()
    #
    # def reconstruction_loss(self, real, predicted):
    #     """Reconstruction loss between the true and predicted outputs
    #        mse = (1/n)*Σ(real - predicted)^2
    #
    #     Args:
    #         real: (array) corresponding array containing the true labels
    #         predictions: (array) corresponding array containing the predicted labels
    #
    #     Returns:
    #         output: (array/float) depending on average parameters the result will be the mean
    #                               of all the sample losses or an array with the losses per sample
    #     """
    #     loss = (real - predicted).pow(2)
    #     return loss.sum(-1).mean()
    #
    # def log_normal(self, x, mu, var):
    #     """Logarithm of normal distribution with mean=mu and variance=var
    #        log(x|μ, σ^2) = loss = -0.5 * Σ log(2π) + log(σ^2) + ((x - μ)/σ)^2
    #
    #     Args:
    #        x: (array) corresponding array containing the input
    #        mu: (array) corresponding array containing the mean
    #        var: (array) corresponding array containing the variance
    #
    #     Returns:
    #        output: (array/float) depending on average parameters the result will be the mean
    #                               of all the sample losses or an array with the losses per sample
    #     """
    #     if self.eps > 0.0:
    #         var = var + self.eps
    #     return -0.5 * torch.sum(
    #         np.log(2.0 * np.pi) + torch.log(var) + torch.pow(x - mu, 2) / var, dim=-1)
    #
    # def gaussian_loss(self, z, z_mu, z_var, z_mu_prior, z_var_prior):
    #     """Variational loss when using labeled data without considering reconstruction loss
    #        loss = log q(z|x,y) - log p(z) - log p(y)
    #
    #     Args:
    #        z: (array) array containing the gaussian latent variable
    #        z_mu: (array) array containing the mean of the inference model
    #        z_var: (array) array containing the variance of the inference model
    #        z_mu_prior: (array) array containing the prior mean of the generative model
    #        z_var_prior: (array) array containing the prior variance of the generative mode
    #
    #     Returns:
    #        output: (array/float) depending on average parameters the result will be the mean
    #                               of all the sample losses or an array with the losses per sample
    #     """
    #     loss = self.log_normal(z, z_mu, z_var) - self.log_normal(z, z_mu_prior, z_var_prior)
    #     return loss.mean()
    #
    # def entropy(self, logits, targets):
    #     """Entropy loss
    #         loss = (1/n) * -Σ targets*log(predicted)
    #
    #     Args:
    #         logits: (array) corresponding array containing the logits of the categorical variable
    #         real: (array) corresponding array containing the true labels
    #
    #     Returns:
    #         output: (array/float) depending on average parameters the result will be the mean
    #                               of all the sample losses or an array with the losses per sample
    #     """
    #     log_q = F.log_softmax(logits, dim=-1)
    #     return -torch.mean(torch.sum(targets * log_q, dim=-1))
    #
    def reconstruction_graph_loss(self, adj_matrix, reconstruct_graph):
        """Reconstruction graph loss
           loss = Σi Σj [(A_ij * log(reconstruct_ij)) + (1 - A_ij) * log(1 - reconstruct_ij)]

        Args:
            adj_matrix (tensor): original normalized adjacency matrix.
            reconstruct_graph (tensor): reconstruct graph by dot product.

        Returns:
            loss (tensor): loss
        """
        adj_matrix = adj_matrix.squeeze()

        # adj_matrix = adj_matrix.fill_diagonal_(1.)

        n, f = reconstruct_graph.size()

        loss = -(torch.mul(adj_matrix, torch.log(reconstruct_graph + 1e-4)) + \
                 torch.mul((1 - adj_matrix), torch.log(1 - reconstruct_graph + 1e-4))).sum()
        return loss / n / n

