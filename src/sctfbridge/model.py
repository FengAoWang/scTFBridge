import torch
import torch.nn as nn
import torch.nn.functional as F
import anndata
from array_api_compat import device

from .utils.loss_function import KL_loss, Reconstruction_loss
from .CLUB import MIEstimator
import pandas as pd
from tqdm import tqdm
import warnings
from torch.utils.data import Dataset, random_split, DataLoader
import numpy as np
from .pathexplainer import PathExplainerTorch
from pathlib import Path
import json
import os
import time

class scMultiDataset(Dataset):
    def __init__(self, anndata_list, batch_info: pd.DataFrame):
        self.rna_tensor = torch.Tensor(anndata_list[0])
        self.atac_tensor = torch.Tensor(anndata_list[1])
        self.TF_tensor = torch.Tensor(anndata_list[2])
        self.batch_info = torch.Tensor(batch_info.to_numpy())

    def __len__(self):
        return self.rna_tensor.shape[0]

    def __getitem__(self, idx):
        rna_data = self.rna_tensor[idx, :]
        atac_data = self.atac_tensor[idx, :]
        TF_data = self.TF_tensor[idx, :]
        batch_info = self.batch_info[idx, :]
        return rna_data, atac_data, TF_data, batch_info


class scOmicsDataset(Dataset):
    def __init__(self,
                 input_data,
                 batch_info: pd.DataFrame):
        self.input_tensor = torch.Tensor(input_data)
        self.batch_info = torch.Tensor(batch_info.to_numpy())

    def __len__(self):
        return self.input_tensor.shape[0]

    def __getitem__(self, idx):
        return self.input_tensor[idx, :], self.batch_info[idx, :]


def reparameterize(mean, logvar):
    std = torch.exp(logvar / 2)
    device = mean.device
    epsilon = torch.randn_like(std, device=device)
    return epsilon * std + mean


def product_of_experts(mu_set_, log_var_set_):
    tmp = 0
    for i in range(len(mu_set_)):
        tmp += torch.div(1, torch.exp(log_var_set_[i]))

    poe_var = torch.div(1., tmp)
    poe_log_var = torch.log(poe_var)

    tmp = 0.
    for i in range(len(mu_set_)):
        tmp += torch.div(1., torch.exp(log_var_set_[i])) * mu_set_[i]
    poe_mu = poe_var * tmp
    return poe_mu, poe_log_var

def gram_matrix(x, sigma=1):
    pairwise_distances = x.unsqueeze(1) - x
    return torch.exp(-pairwise_distances.norm(2, dim=2) / (2 * sigma * sigma))


class MaskedLinear(nn.Linear):
    def __init__(self, n_in, n_out, mask, latent_dim, bias=False):
        # mask 应该和转置后的权重维度相同
        # n_input x n_output_nodes
        if latent_dim != mask.shape[0] or n_out != mask.shape[1]:
            raise ValueError('Incorrect shape of the mask.')

        super().__init__(n_in, n_out, bias)

        self.register_buffer('mask', mask.t())
        self.latent_dim = latent_dim

        # 初始化权重时，只保留 mask 不为 0 的部分
        with torch.no_grad():
            self.weight[:, :self.latent_dim] *= self.mask  # 仅初始化非零的 mask 部分
            # torch.nn.init.zeros_(self.bias)

    def forward(self, x):
        # 克隆权重来避免原地操作
        masked_weight = self.weight.clone()

        # 仅保留 mask 不为 0 的地方的权重值
        masked_weight[:, :self.latent_dim] = masked_weight[:, :self.latent_dim] * self.mask + (
                1 - self.mask) * self.weight[:, :self.latent_dim].detach()

        # 计算前向传播
        x = nn.functional.linear(x, masked_weight)
        return x


class LinearLayer(nn.Module):
    def __init__(self,
                 input_dim: int,
                 output_dim: int,
                 dropout: float = 0.2,
                 batchnorm: bool = False,
                 activation=None,
                 mask=None):
        super(LinearLayer, self).__init__()
        self.linear_layer = nn.Linear(input_dim, output_dim)

        self.dropout = nn.Dropout(dropout) if dropout > 0 else None
        self.batchnorm = nn.BatchNorm1d(output_dim) if batchnorm else None

        self.activation = None
        self.mask = mask
        if activation is not None:
            if activation == 'relu':
                self.activation = F.relu
            elif activation == 'sigmoid':
                self.activation = torch.sigmoid
            elif activation == 'tanh':
                self.activation = torch.tanh
            elif activation == 'leakyrelu':
                self.activation = torch.nn.LeakyReLU()
            elif activation == 'selu':
                self.activation = torch.nn.SELU()

    def forward(self, input_x):
        # 如果 mask 存在，将它应用到线性层的权重上
        if self.mask is not None:
            masked_weight = self.linear_layer.weight * self.mask.T
            x = F.linear(input_x, masked_weight, self.linear_layer.bias)
        else:
            x = self.linear_layer(input_x)

        if self.batchnorm is not None:
            x = self.batchnorm(x)
        if self.activation is not None:
            x = self.activation(x)
        if self.dropout is not None:
            x = self.dropout(x)
        return x


class Scaler(nn.Module):
    def __init__(self, feature_dim, tau=0.5):
        super(Scaler, self).__init__()
        self.tau = tau
        # 直接在 __init__ 中注册 scale 参数
        self.scale = nn.Parameter(torch.zeros(feature_dim))

    def forward(self, inputs, mode='positive'):
        if mode == 'positive':
            scale = self.tau + (1 - self.tau) * torch.sigmoid(self.scale)
        else:
            scale = (1 - self.tau) * torch.sigmoid(-self.scale)
        device = inputs.device
        scale = scale.to(device)
        return inputs * torch.sqrt(scale + 1e-8)


class ModalVAEEncoder(nn.Module):
    def __init__(self,
                 input_dim,
                 hidden_dims,
                 latent_dim,
                 tau=0.5,
                 activation=None):
        super(ModalVAEEncoder, self).__init__()
        self.FeatureEncoder = nn.ModuleList([LinearLayer(input_dim, hidden_dims[0],
                                                         batchnorm=True, activation=activation)])
        for i in range(len(hidden_dims) - 1):
            self.FeatureEncoder.append(LinearLayer(hidden_dims[i], hidden_dims[i + 1],
                                                   batchnorm=True, activation=activation))

        #   trick: Add batchnorm layer to mu/logvar prediction layer
        self.mu_predictor = LinearLayer(hidden_dims[-1],
                                        latent_dim,
                                        batchnorm=True)
        self.logVar_predictor = LinearLayer(hidden_dims[-1],
                                            latent_dim,
                                            batchnorm=True)

        # 均值和对数方差路径
        # self.z_mean_predictor = nn.Linear(hidden_dims[-1], latent_dim)
        # self.z_mean_bn = nn.BatchNorm1d(latent_dim, affine=False, eps=1e-8)
        # self.z_mean_scaler = Scaler(tau=tau, feature_dim=latent_dim)
        #
        # self.z_logvar_predictor = nn.Linear(hidden_dims[-1], latent_dim)
        # self.z_logvar_bn = nn.BatchNorm1d(latent_dim, affine=False, eps=1e-8)
        # self.z_logvar_scaler = Scaler(tau=tau, feature_dim=latent_dim)

    def forward(self, input_x):
        for layer in self.FeatureEncoder:
            input_x = layer(input_x)

        # 均值路径
        # z_mean = self.z_mean_predictor(input_x)
        # z_mean = self.z_mean_bn(z_mean)
        # z_mean = self.z_mean_scaler(z_mean, mode='positive')
        #
        # # 对数方差路径
        # z_logvar = self.z_logvar_predictor(input_x)
        # z_logvar = self.z_logvar_bn(z_logvar)
        # z_logvar = self.z_logvar_scaler(z_logvar, mode='negative')
        #
        # # 重参数化采样
        # latent_z = reparameterize(z_mean, z_logvar)

        z_mean = self.mu_predictor(input_x)
        z_logvar = self.logVar_predictor(input_x)
        latent_z = reparameterize(z_mean, z_logvar)
        return z_mean, z_logvar, latent_z


class Encoder(nn.Module):
    def __init__(self,
                 input_dim,
                 hidden_dims,
                 activation=None,):
        super(Encoder, self).__init__()

        self.FeatureEncoder = nn.ModuleList([LinearLayer(input_dim, hidden_dims[0], batchnorm=True, activation=activation)])

        for i in range(len(hidden_dims) - 1):
            self.FeatureEncoder.append(
                LinearLayer(hidden_dims[i], hidden_dims[i + 1], batchnorm=True, activation=activation))

    def forward(self, input_x):
        for layer in self.FeatureEncoder:
            input_x = layer(input_x)
        return input_x


class ModalVAEDecoder(nn.Module):
    def __init__(self,
                 input_dim,
                 hidden_dims,
                 output_dim,
                 activation=None,
                 ):
        super(ModalVAEDecoder, self).__init__()
        self.FeatureDecoder = nn.ModuleList([LinearLayer(input_dim, hidden_dims[0],
                                                         dropout=0.2, batchnorm=True,
                                                         activation=activation)])
        for i in range(len(hidden_dims) - 1):
            self.FeatureDecoder.append(LinearLayer(hidden_dims[i], hidden_dims[i + 1],
                                                   dropout=0.2, batchnorm=True,
                                                   activation=activation))
        self.ReconsPredictor = LinearLayer(hidden_dims[-1], output_dim)

    def forward(self, input_x):
        for layer in self.FeatureDecoder:
            input_x = layer(input_x)
        DataRecons = self.ReconsPredictor(input_x)
        return DataRecons


class scTFBridge(nn.Module):
    def __init__(self,
                 input_dims,
                 encoder_hidden_dims,
                 decoder_hidden_dims,
                 latent_dim,
                 kl_weight,
                 dist,
                 batch_dims,
                 con_weight,
                 mi_weight,
                 mask,
                 batch_key=None,
                 device: torch.device = torch.device('cpu'),
                 ):
        super(scTFBridge, self).__init__()
        self.batch_key = batch_key
        self.kl_weight = kl_weight
        self.con_weight = con_weight
        self.mi_weight = mi_weight
        self.device = device

        self.temperature = 1
        self.input_dims = input_dims
        self.batch_dims = batch_dims
        self.encoder_hidden_dims = encoder_hidden_dims
        self.decoder_hidden_dims = decoder_hidden_dims

        self.RNADataDist = dist[0]
        self.ATACDataDist = dist[1]

        #   modality private Encoder
        self.RNAEncoder = Encoder(input_dims[0],
                                  encoder_hidden_dims,
                                  activation='relu')
        self.ATACEncoder = Encoder(input_dims[1],
                                   encoder_hidden_dims,
                                   activation='relu')

        self.RNAPrivateEncoder = ModalVAEEncoder(encoder_hidden_dims[-1],
                                                 [256],
                                                 latent_dim,
                                                 activation='relu')
        self.ATACPrivateEncoder = ModalVAEEncoder(encoder_hidden_dims[-1],
                                                  [256],
                                                  latent_dim,
                                                  activation='relu')

        #   modality share Encoder
        self.ModalShareEncoder = ModalVAEEncoder(encoder_hidden_dims[-1],
                                                 [256],
                                                 latent_dim,
                                                 activation='relu')

        #   modality private Decoder
        self.RNADecoder = ModalVAEDecoder(latent_dim * 2 + batch_dims,
                                          decoder_hidden_dims,
                                          input_dims[0],
                                          activation='relu')

        self.ATACDecoder = MaskedLinear(latent_dim * 2 + batch_dims,
                                        input_dims[1],
                                        mask.T,
                                        latent_dim)

        #   MI estimator
        self.MI = MIEstimator(latent_dim)

        self.latent_mode = 'mu'
        self.latent_dim = latent_dim

        self._initialize_weights()

    @classmethod
    def from_anndata(
            cls,
            rna_adata: anndata.AnnData,
            atac_adata: anndata.AnnData,
            tf_adata: anndata.AnnData,
            mask_path: str,
            batch_key: str,
            device: torch.device,
            # Model hyperparameters can be passed here with sensible defaults
            encoder_hidden_dims=None,
            decoder_hidden_dims=None,
            kl_weight: float = 1.0,
            con_weight: float = 1.0,
            mi_weight: float = 1.0,
            temperature: float = 0.1
    ):
        """
        A factory method to conveniently initialize the scMulti model from AnnData objects.

        This method automatically infers input dimensions, latent dimension, and batch
        dimensions from the provided data and file paths.

        Args:
            rna_adata (anndata.AnnData): AnnData object for RNA-seq.
            atac_adata (anndata.AnnData): AnnData object for ATAC-seq.
            tf_adata (anndata.AnnData): AnnData object for TF activities. The number of TFs
                                        is used as the latent dimension.
            mask_path (str): Path to the TF-peak binding mask file (e.g., a tab-separated file).
            batch_key (str): The column name in `.obs` that contains batch information.
            encoder_hidden_dims (List[int]): Hidden layer dimensions for the encoders.
            decoder_hidden_dims (List[int]): Hidden layer dimensions for the decoders.
            kl_weight (float): Weight for the KL divergence loss.
            con_weight (float): Weight for the contrastive loss.
            temperature (float): Temperature parameter for the contrastive loss.
            device (torch.device): Device to use.

        Returns:
            An instance of the scMulti model.
        """
        if encoder_hidden_dims is None:
            encoder_hidden_dims = [1024]
        if decoder_hidden_dims is None:
            decoder_hidden_dims = [1024]
        print("🔧 Initializing model from AnnData objects...")

        # 1. Infer input dimensions from data shapes
        input_dims = (rna_adata.shape[1], atac_adata.shape[1])
        print(f"  - RNA input features: {input_dims[0]}")
        print(f"  - ATAC input features: {input_dims[1]}")

        # 2. Infer latent dimension from the number of Transcription Factors
        latent_dim = tf_adata.shape[1]
        print(f"  - Latent dimension (from TFs): {latent_dim}")

        # 3. Infer batch dimensions from the batch_key
        if batch_key in rna_adata.obs.columns:
            batches_info = rna_adata.obs[batch_key]
            one_hot_encoded_batches = pd.get_dummies(batches_info, prefix='batch', dtype=float)
            batch_dims = one_hot_encoded_batches.shape[1]
            print(f"  - Found {batch_dims} batches from key '{batch_key}'")
        else:
            warnings.warn(
                f"Batch key '{batch_key}' not found in `rna_adata.obs`. Assuming a single batch (batch_dims=1).")
            batch_dims = 1

        # 4. Load the TF-peak mask
        # NOTE: The logic here assumes the mask file's shape is (n_peaks, n_tfs) or similar.
        # The MaskedLinear layer expects a mask of shape (output_dim, input_dim).
        # This is a critical step that requires careful dimension management.
        # For this example, we'll load it and create a correctly-sized zero-mask
        # and place the TF-peak info into it.

        print(f"  - Loading TF-peak mask from: {mask_path}")
        try:
            # This is the original TF->Peak mask, shape (n_peaks, n_tfs)
            tf_peak_mask_np = pd.read_csv(mask_path, sep='\t', header=None).values
        except FileNotFoundError:
            raise FileNotFoundError(f"Mask file not found at path: {mask_path}")

        mask_tensor = torch.tensor(tf_peak_mask_np).float()

        # 5. Call the main constructor with all inferred and provided parameters
        print("✅ Model ready.")
        return cls(
            input_dims=input_dims,
            encoder_hidden_dims=encoder_hidden_dims,
            decoder_hidden_dims=decoder_hidden_dims,
            latent_dim=latent_dim,
            batch_dims=batch_dims,
            mask=mask_tensor,
            kl_weight=kl_weight,
            mi_weight=mi_weight,
            con_weight=con_weight,
            dist=['gaussian', 'bernoulli'],
            batch_key=batch_key,
            device=device
        )

    def save(self, dir_path: str, overwrite: bool = False):
        """
        将模型配置、权重和 mask 保存到指定目录。

        Args:
            dir_path (str): 保存模型的目录路径。
            overwrite (bool): 如果目录已存在，是否覆盖。
        """
        print(f"💾 Saving model to {dir_path}...")
        p = Path(dir_path)
        if p.exists() and not overwrite:
            raise ValueError(f"Directory {dir_path} already exists. Use `overwrite=True` to overwrite.")
        p.mkdir(parents=True, exist_ok=True)

        # 1. 保存初始化所需的配置参数
        # 注意：不直接保存张量，只保存可以序列化的配置
        config = {
            "input_dims": self.input_dims,
            "encoder_hidden_dims": self.encoder_hidden_dims,  # 假设您的Encoder保存了此信息
            "decoder_hidden_dims": self.decoder_hidden_dims,  # 假设您的Decoder保存了此信息
            "latent_dim": self.latent_dim,
            "kl_weight": self.kl_weight,
            "dist": (self.RNADataDist, self.ATACDataDist),
            "batch_dims": self.batch_dims,
            "con_weight": self.con_weight,
            "mi_weight": self.mi_weight,
        }
        with open(p / "model_config.json", "w") as f:
            json.dump(config, f, indent=4)

        # 2. 保存模型权重 (state_dict)
        torch.save(self.state_dict(), p / "model_weights.pt")

        # 3. 单独保存大的张量，如 mask
        # 假设 mask 是在 __init__ 时传入的，并且保存在 self.ATACDecoder.mask
        torch.save(self.ATACDecoder.mask, p / "mask.pt")

        print("✅ Model saved successfully.")

    @classmethod
    def load(cls, dir_path: str, device: torch.device = None):
        """
        从目录中高速加载模型。

        Args:
            dir_path (str): 保存模型的目录路径。
            device (torch.device): 要将模型加载到的目标设备 (例如 torch.device('cuda:0'))。
                                   如果为 None，则自动选择 GPU 或 CPU。

        Returns:
            An instance of the scMulti model, loaded and ready for inference.
        """
        print(f"🚀 Loading model from {dir_path}...")
        p = Path(dir_path)
        if not p.exists():
            raise FileNotFoundError(f"Directory {dir_path} not found.")

        # 确定目标设备
        if device is None:
            device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        print(f"  - Using device: {device}")

        # 1. 加载配置和大的张量
        with open(p / "model_config.json", "r") as f:
            config = json.load(f)

        config["device"] = device
        # 将 mask 直接加载到目标设备
        mask = torch.load(p / "mask.pt", map_location='cpu')

        # 2. 使用配置初始化模型结构
        # 注意：模型首先在CPU上创建结构，但稍后会整体移动到device
        model = cls(mask=mask, **config)

        # 3. 高速加载模型权重
        # 这是关键！使用 map_location 直接将权重张量加载到目标设备内存
        weights_path = p / "model_weights.pt"
        state_dict = torch.load(weights_path, map_location=device)
        model.load_state_dict(state_dict)

        # 4. 将整个模型（包括在 __init__ 中创建的、但不在 state_dict 中的其他参数或缓冲区）
        #    移动到目标设备，并设置为评估模式。
        model.to(device)
        model.eval()

        print("✅ Model loaded and ready for inference.")
        return model

    def _initialize_weights(self):
        # 遍历模型的每一层并应用Kaiming初始化
        for m in self.modules():
            if isinstance(m, nn.Linear):
                nn.init.kaiming_normal_(m.weight, nonlinearity='relu')
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0)

    def get_latent_distribution(self, input_x):
        input_x_rna = input_x[0]
        input_x_atac = input_x[1]
        input_x_TF = input_x[2]

        rna_embedd = self.RNAEncoder(input_x_rna)
        atac_embedd = self.ATACEncoder(input_x_atac)

        #   share feature
        rna_share_mu, rna_share_logvar, rna_share_latent_z = self.ModalShareEncoder(rna_embedd)
        atac_share_mu, atac_share_logvar, atac_share_latent_z = self.ModalShareEncoder(atac_embedd)

        joint_mu, joint_logvar = product_of_experts([rna_share_mu, atac_share_mu],
                                                    [rna_share_logvar, atac_share_logvar])
        return {'joint_mu': joint_mu, 'joint_logvar': joint_logvar}

    def forward(self, input_x, batch_id):
        input_x_rna = input_x[0]
        input_x_atac = input_x[1]
        input_x_TF = input_x[2]

        rna_embedd = self.RNAEncoder(input_x_rna)
        atac_embedd = self.ATACEncoder(input_x_atac)

        #   share feature
        rna_share_mu, rna_share_logvar, rna_share_latent_z = self.ModalShareEncoder(rna_embedd)
        atac_share_mu, atac_share_logvar, atac_share_latent_z = self.ModalShareEncoder(atac_embedd)

        joint_mu, joint_logvar = product_of_experts([rna_share_mu, atac_share_mu],
                                                    [rna_share_logvar, atac_share_logvar])
        joint_latent_z = reparameterize(joint_mu, joint_logvar)

        #   private feature
        rna_private_mu, rna_private_logvar, rna_private_latent_z = self.RNAPrivateEncoder(rna_embedd)
        atac_private_mu, atac_private_logvar, atac_private_latent_z = self.ATACPrivateEncoder(atac_embedd)
        if self.latent_mode == 'mu':
            output = {
                'share_embedding': joint_mu,
                'RNA_private_embedding': rna_private_mu,
                'ATAC_private_embedding': atac_private_mu,
                'RNA_share_embedding': rna_share_mu,
                'ATAC_share_embedding': atac_share_mu
            }
        else:
            output = {
                'share_embedding': joint_latent_z,
                'RNA_private_embedding': rna_private_latent_z,
                'ATAC_private_embedding': atac_private_latent_z,
                'RNA_share_embedding': rna_share_latent_z,
                'ATAC_share_embedding': atac_share_latent_z
            }
        return output

    def cross_modal_generation(self, input_x, batch_id):
        # input_x_rna = torch.concat((input_x[0], batch_id), dim=1)
        # input_x_atac = torch.concat((input_x[1], batch_id), dim=1)

        input_x_rna = input_x[0]
        input_x_atac = input_x[1]
        input_x_TF = input_x[2]

        rna_embedd = self.RNAEncoder(input_x_rna)
        atac_embedd = self.ATACEncoder(input_x_atac)

        #   share feature
        rna_share_mu, rna_share_logvar, rna_share_latent_z = self.ModalShareEncoder(rna_embedd)
        atac_share_mu, atac_share_logvar, atac_share_latent_z = self.ModalShareEncoder(atac_embedd)

        prior_z_rna = torch.randn_like(atac_share_latent_z)
        prior_z_atac = torch.randn_like(rna_share_latent_z)

        device = batch_id.device
        recon_batch_id = torch.zeros_like(batch_id).to(device)

        #   cross-modal generation
        rna_recon_cross_from_atac = self.RNADecoder(
            torch.cat([atac_share_latent_z, prior_z_rna, recon_batch_id], dim=1))
        atac_recon_cross_from_rna = self.ATACDecoder(
            torch.cat([rna_share_latent_z, prior_z_atac, recon_batch_id], dim=1))

        return rna_recon_cross_from_atac, atac_recon_cross_from_rna

    def rna_generation_from_atac(self, input_x, batch_id):
        input_x_atac = input_x
        atac_embedd = self.ATACEncoder(input_x_atac)

        #   share feature
        atac_share_mu, atac_share_logvar, atac_share_latent_z = self.ModalShareEncoder(atac_embedd)

        prior_z_rna = torch.randn_like(atac_share_latent_z)
        device = batch_id.device

        recon_batch_id = torch.zeros_like(batch_id).to(device)

        #   cross-modal generation
        rna_recon_cross_from_atac = self.RNADecoder(
            torch.cat([atac_share_latent_z, prior_z_rna, recon_batch_id], dim=1))

        return rna_recon_cross_from_atac

    def rna_generation_from_rna_private(self, rna_private_latent_z, batch_id):
        latent_share = torch.randn_like(rna_private_latent_z)
        device = batch_id.device

        recon_batch_id = torch.zeros_like(batch_id).to(device)
        rna_recon_cross_from_atac = self.RNADecoder(
            torch.cat([latent_share, rna_private_latent_z, batch_id], dim=1))

        return rna_recon_cross_from_atac

    def atac_generation_from_rna(self, input_x, batch_id):
        input_x_rna = input_x
        rna_embedd = self.RNAEncoder(input_x_rna)

        #   share feature
        atac_share_mu, atac_share_logvar, atac_share_latent_z = self.ModalShareEncoder(rna_embedd)

        prior_z_atac = torch.randn_like(atac_share_latent_z)
        device = batch_id.device
        recon_batch_id = torch.zeros_like(batch_id).to(device)

        #   cross-modal generation
        atac_recon_cross_from_rna = self.ATACDecoder(
            torch.cat([atac_share_latent_z, prior_z_atac, recon_batch_id], dim=1))

        return atac_recon_cross_from_rna

    def modal_generation_from_latent(self, share_embedding, private_embedding, batch_id, modality_name):
        device = batch_id.device
        recon_batch_id = torch.zeros_like(batch_id).to(device)

        latent_z = torch.cat([share_embedding, private_embedding, recon_batch_id], dim=1)
        if modality_name == 'rna':
            modal_generation = self.RNADecoder(latent_z)
        else:
            modal_generation = self.ATACDecoder(latent_z)
        return modal_generation

    def modal_generation_from_latentTF(self, TF_embedding, batch_id, modality_name):
        #
        device = batch_id.device
        recon_batch_id = torch.zeros_like(batch_id).to(device)
        atac_share_mu, atac_share_logvar, atac_share_latent_z = self.ModalShareEncoder(TF_embedding)
        prior_z_atac = torch.randn_like(atac_share_latent_z)
        latent_z = torch.cat([atac_share_latent_z, prior_z_atac, recon_batch_id], dim=1)

        modal_generation = self.RNADecoder(latent_z)

        return modal_generation

    def modal_generation(self, input_x, batch_id):
        input_x_rna = input_x[0]
        input_x_atac = input_x[1]

        rna_embedd = self.RNAEncoder(input_x_rna)
        atac_embedd = self.ATACEncoder(input_x_atac)

        #   share feature
        rna_share_mu, rna_share_logvar, rna_share_latent_z = self.ModalShareEncoder(rna_embedd)
        atac_share_mu, atac_share_logvar, atac_share_latent_z = self.ModalShareEncoder(atac_embedd)

        joint_mu, joint_logvar = product_of_experts([rna_share_mu, atac_share_mu],
                                                    [rna_share_logvar, atac_share_logvar])

        #   private feature
        rna_private_mu, rna_private_logvar, rna_private_latent_z = self.RNAPrivateEncoder(rna_embedd)
        atac_private_mu, atac_private_logvar, atac_private_latent_z = self.ATACPrivateEncoder(atac_embedd)

        joint_latent_z = reparameterize(joint_mu, joint_logvar)

        device = batch_id.device
        recon_batch_id = torch.zeros_like(batch_id).to(device)

        #   cross-modal generation
        rna_recon = self.RNADecoder(torch.cat([joint_latent_z, rna_private_latent_z, recon_batch_id], dim=1))
        atac_recon = self.ATACDecoder(torch.cat([joint_latent_z, atac_private_latent_z, recon_batch_id], dim=1))

        return rna_recon, atac_recon

    def compute_hsic(self, x, y, sigma=1):
        m = x.shape[0]
        K = gram_matrix(x, sigma=sigma)
        L = gram_matrix(y, sigma=sigma)
        H = torch.eye(m) - 1.0 / m * torch.ones((m, m))
        device = x.device

        H = H.float().to(device)
        HSIC = torch.trace(torch.mm(L, torch.mm(H, torch.mm(K, H)))) / ((m - 1) ** 2)
        return HSIC

    def compute_loss(self, input_x, batch_id):
        # input_x_rna = torch.concat((input_x[0], batch_id), dim=1)
        # input_x_atac = torch.concat((input_x[1], batch_id), dim=1)

        input_x_rna = input_x[0]
        input_x_atac = input_x[1]
        input_x_TF = input_x[2]

        rna_embedd = self.RNAEncoder(input_x_rna)
        atac_embedd = self.ATACEncoder(input_x_atac)

        # atac_embedd = atac_embedd

        #   share feature
        rna_share_mu, rna_share_logvar, rna_share_latent_z = self.ModalShareEncoder(rna_embedd)
        atac_share_mu, atac_share_logvar, atac_share_latent_z = self.ModalShareEncoder(atac_embedd)

        joint_mu, joint_logvar = product_of_experts([rna_share_mu, atac_share_mu],
                                                    [rna_share_logvar, atac_share_logvar])
        joint_latent_z = reparameterize(joint_mu, joint_logvar)

        #   private feature
        rna_private_mu, rna_private_logvar, rna_private_latent_z = self.RNAPrivateEncoder(rna_embedd)
        atac_private_mu, atac_private_logvar, atac_private_latent_z = self.ATACPrivateEncoder(atac_embedd)

        #   multimodal latent recon
        rna_latent_z = torch.concat((joint_latent_z, rna_private_latent_z, batch_id), dim=1)
        rna_recon = self.RNADecoder(rna_latent_z)

        atac_latent_z = torch.concat((joint_latent_z, atac_private_latent_z, batch_id), dim=1)
        atac_recon = self.ATACDecoder(atac_latent_z)

        #   recon loss
        rna_recon_loss = Reconstruction_loss(rna_recon, input_x[0], 1, self.RNADataDist)
        atac_recon_loss = Reconstruction_loss(atac_recon, input_x[1], 1, self.ATACDataDist)

        #  self-modal latent recon
        rna_self_latent_z = torch.concat((rna_share_latent_z, rna_private_latent_z, batch_id), dim=1)
        rna_self_recon = self.RNADecoder(rna_self_latent_z)
        atac_self_latent_z = torch.concat((atac_share_latent_z, atac_private_latent_z, batch_id), dim=1)
        atac_self_recon = self.ATACDecoder(atac_self_latent_z)

        #   self recon loss
        rna_self_recon_loss = Reconstruction_loss(rna_self_recon, input_x[0], 1, self.RNADataDist)
        atac_self_recon_loss = Reconstruction_loss(atac_self_recon, input_x[1], 1, self.ATACDataDist)

        #   kl loss
        rna_kl_loss = KL_loss(rna_private_mu, rna_private_logvar, 1)
        atac_kl_loss = KL_loss(atac_private_mu, atac_private_logvar, 1)

        rna_prior_z = torch.randn_like(joint_latent_z)
        atac_prior_z = torch.randn_like(joint_latent_z)

        #   cross-modal generation
        rna_recon_cross_from_atac = self.RNADecoder(torch.cat([atac_share_latent_z, rna_prior_z, batch_id], dim=1))
        atac_recon_cross_from_rna = self.ATACDecoder(torch.cat([rna_share_latent_z, atac_prior_z, batch_id], dim=1))

        rna_cross_recon_loss = Reconstruction_loss(rna_recon_cross_from_atac, input_x[0], 1, self.RNADataDist)
        atac_cross_recon_loss = Reconstruction_loss(atac_recon_cross_from_rna, input_x[1], 1, self.ATACDataDist)

        # Contrastive loss with similarity matrix
        logits = (rna_share_latent_z @ atac_share_latent_z.T) / self.temperature
        rna_similarity = rna_share_latent_z @ rna_share_latent_z.T
        atac_similarity = atac_share_latent_z @ atac_share_latent_z.T

        targets = F.softmax(
            (rna_similarity + atac_similarity) / 2 * self.temperature, dim=-1
        )
        rna_loss = F.cross_entropy(logits, targets, reduction='mean')
        atac_loss = F.cross_entropy(logits.T, targets.T, reduction='mean')
        multimodal_contrastive_loss = (rna_loss + atac_loss)

        #   contrastive loss
        rna_contrastive_loss = self.contrastive_loss(rna_share_latent_z, atac_share_latent_z, rna_private_latent_z)
        atac_contrastive_loss = self.contrastive_loss(atac_share_latent_z, rna_share_latent_z, atac_private_latent_z)
        contrastive_loss = rna_contrastive_loss + atac_contrastive_loss

        share_kl_loss = (KL_loss(joint_mu, joint_logvar, 1) +
                         KL_loss(rna_share_mu, rna_share_logvar, 1) +
                         KL_loss(atac_share_mu, atac_share_logvar, 1))
        #   MI loss
        PIB_mimin = self.MI.learning_loss(rna_private_latent_z, atac_private_latent_z,
                                          joint_latent_z, rna_share_latent_z, atac_share_latent_z)
        PIB_mi = self.MI(rna_private_latent_z, atac_private_latent_z, joint_latent_z, rna_share_latent_z,
                         atac_share_latent_z)

        # hsic_loss = self.compute_hsic(joint_latent_z, batch_id)

        loss = (
                (rna_self_recon_loss + atac_self_recon_loss) +
                (rna_recon_loss + atac_recon_loss) +
                (rna_cross_recon_loss + atac_cross_recon_loss) +
                self.kl_weight * (share_kl_loss + rna_kl_loss + atac_kl_loss) +
                self.con_weight * (contrastive_loss + multimodal_contrastive_loss) +
                self.con_weight * (PIB_mi + PIB_mimin)
        )
        return loss

    def fit(self,
            adata_list: list,
            device=None,
            batch_size=128,
            lr=1e-3,
            epochs=150,
            use_early_stopping=True,
            early_stopping_patience=10,
            weight_decay=1e-4,
            val_ratio=0.2
            ):
        """
        训练模型

        Args:
            adata_list (list): 包含 [rna_data, atac_data, TF_data] 的 anndata 对象列表。
            device: training dev
            batch_size (int): 训练批次大小。
            lr (float): 优化器的学习率。
            epochs (int): 最大的训练周期数。
            use_early_stopping (bool): 是否启用早停机制。
            early_stopping_patience (int): 早停的耐心值。
            weight_decay (float): Adam 优化器的权重衰减 (L2 正则化)。
            val_ratio (float): 用于验证集的样本比例。
        """
        print(f"🔥 Starting model training for up to {epochs} epochs...")
        # 1. Device setup
        if device is None:
            device = self.device
        self.to(device)
        print(f"  - Using device: {device}")

        rna_data = adata_list[0]
        atac_data = adata_list[1]
        TF_data = adata_list[2]

        # 2. Data Preparation
        print("  - Preparing datasets and data loaders...")
        if self.batch_key in rna_data.obs.columns:
            batches_info = rna_data.obs[self.batch_key].values.tolist()
        else:
            warnings.warn("batch_key not found or not defined. Using a single batch for all cells.")
            batches_info = np.zeros(rna_data.shape[0])

        # print(batches_info)

        one_hot_encoded_batches = pd.get_dummies(batches_info, prefix='batch', dtype=float)

        rna_data_np = rna_data.X.toarray()
        atac_data_np = atac_data.X.toarray()
        TF_adata_np = TF_data.X.toarray()

        # 3. Create DataLoaders with optional validation split
        dataset = scMultiDataset([rna_data_np, atac_data_np, TF_adata_np], one_hot_encoded_batches)

        # 3. 根据 use_early_stopping 的值，条件性地创建 DataLoader
        val_loader = None
        if use_early_stopping:
            if val_ratio <= 0 or val_ratio >= 1:
                raise ValueError("val_ratio must be between 0 and 1 for early stopping.")

            val_size = int(len(dataset) * val_ratio)
            train_size = len(dataset) - val_size

            if val_size == 0:
                warnings.warn("Validation set size is 0 based on val_ratio. Disabling early stopping.")
                use_early_stopping = False
                train_loader = DataLoader(dataset, batch_size=batch_size, shuffle=True, num_workers=8, pin_memory=True)
            else:
                print(f"  - Splitting data: {train_size} training samples, {val_size} validation samples.")
                train_dataset, val_dataset = random_split(dataset, [train_size, val_size])
                train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True, drop_last=True)
                val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)
        else:
            print("  - Using all data for training. Early stopping is disabled.")
            train_loader = DataLoader(dataset, batch_size=batch_size, shuffle=True, drop_last=True)

        optimizer = torch.optim.Adam(self.parameters(), lr=lr)

        if use_early_stopping:
            best_val_loss = float('inf')
            patience_counter = 0

            # 5. 训练循环 (使用单个 epoch 级别的进度条)
            with tqdm(range(epochs), unit="epoch") as tepoch:
                for epoch in tepoch:
                    tepoch.set_description(f"Epoch {epoch + 1}/{epochs}")

                    # --- 训练阶段 ---
                    self.train()
                    total_train_loss = 0.0
                    for batch, data in enumerate(train_loader):
                        rna_data, atac_data, TF_data, batch_id = data
                        rna_data = rna_data.to(device)
                        atac_data = atac_data.to(device)
                        TF_data = TF_data.to(device)
                        batch_id = batch_id.to(device)

                        optimizer.zero_grad()
                        loss = self.compute_loss([rna_data, atac_data, TF_data], batch_id)
                        loss.backward()
                        optimizer.step()
                        total_train_loss += loss.item()

                    avg_train_loss = total_train_loss / len(train_loader)
                    loss_postfix = {'train_loss': f'{avg_train_loss:.4f}'}

                    # --- 验证与早停阶段 (仅在启用时执行) ---
                    if use_early_stopping:
                        self.eval()
                        total_val_loss = 0.0
                        with torch.no_grad():
                            for batch, data in enumerate(val_loader):
                                rna_data, atac_data, TF_data, batch_id = data
                                rna_data = rna_data.to(device)
                                atac_data = atac_data.to(device)
                                TF_data = TF_data.to(device)
                                batch_id = batch_id.to(device)
                                loss = self.compute_loss([rna_data, atac_data, TF_data], batch_id)
                                total_val_loss += loss.item()

                        avg_val_loss = total_val_loss / len(val_loader)
                        loss_postfix['val_loss'] = f'{avg_val_loss:.4f}'

                        # 早停判断
                        if avg_val_loss < best_val_loss:
                            best_val_loss = avg_val_loss
                            patience_counter = 0
                        else:
                            patience_counter += 1
                            if patience_counter >= early_stopping_patience:
                                # 使用 tepoch.write 代替 print, 避免破坏进度条格式
                                tepoch.write(f"\n🛑 Early stopping triggered at epoch {epoch + 1}.")
                                break  # 提前终止训练循环

                    # 更新进度条后缀信息
                    tepoch.set_postfix(loss_postfix)

            print("✅ Training finished successfully.")

    def get_embeddings(self,
                       adata_list: list[anndata.AnnData],
                       device: torch.device = None,
                       batch_size: int = 128,
                       latent_mode: str = 'latent_z'
                       ):
        """
        Extracts latent embeddings from the model for all modalities.

        Args:
            adata_list (list): A list of anndata objects [rna_data, atac_data, TF_data].
            device: The device to use for inference.
            batch_size (int): The batch size for processing data.
            latent_mode (str): The type of latent representation to return ('mu' or 'latent_z').

        Returns:
            A dictionary containing numpy arrays for joint and private embeddings.
        """

        print("✨ Extracting latent embeddings from the model...")

        if device is None:
            device = self.device
        self.to(device)
        self.eval()
        print(f"  - Using device: {device}")
        rna_data = adata_list[0]
        atac_data = adata_list[1]
        TF_data = adata_list[2]

        if self.batch_key in rna_data.obs.columns:
            batches_info = rna_data.obs[self.batch_key].values.tolist()
        else:
            warnings.warn("batch_key not found or not defined. Using a single batch for all cells.")
            batches_info = np.zeros(rna_data.shape[0])

        one_hot_encoded_batches = pd.get_dummies(batches_info, prefix='batch', dtype=float)

        rna_data_np = rna_data.X.toarray()
        atac_data_np = atac_data.X.toarray()
        TF_adata_np = TF_data.X.toarray()

        # 3. 创建 Dataset 和 DataLoader (使用全部数据)
        dataset = scMultiDataset([rna_data_np, atac_data_np, TF_adata_np], one_hot_encoded_batches)
        data_loader = DataLoader(dataset, batch_size=batch_size, shuffle=False, num_workers=4, pin_memory=True)

        joint_share_representations = []
        rna_private_representations = []
        atac_private_representations = []
        rna_share_representations = []
        atac_share_representations = []
        self.latent_mode = 'latent_z'

        with torch.no_grad():
            with tqdm(data_loader, unit='batch') as tepoch:
                for batch, data in enumerate(tepoch):
                    rna_data, atac_data, TF_data, batch_id = data
                    rna_data = rna_data.to(device)
                    atac_data = atac_data.to(device)
                    batch_id = batch_id.to(device)
                    TF_data = TF_data.to(device)

                    output = self([rna_data, atac_data, TF_data], batch_id)
                    share_embedding = output['share_embedding'].cpu()
                    rna_embedding = output['RNA_private_embedding'].cpu()
                    atac_embedding = output['ATAC_private_embedding'].cpu()
                    rna_share = output['RNA_share_embedding'].cpu()
                    atac_share = output['ATAC_share_embedding'].cpu()

                    joint_share_representations.append(share_embedding)
                    rna_private_representations.append(rna_embedding)
                    atac_private_representations.append(atac_embedding)
                    rna_share_representations.append(rna_share)
                    atac_share_representations.append(atac_share)
        output = {
            'joint_share_representations': torch.cat(joint_share_representations, dim=0).detach().cpu().numpy(),
            'RNA_private_representations': torch.cat(rna_private_representations, dim=0).detach().cpu().numpy(),
            'ATAC_private_representations': torch.cat(atac_private_representations, dim=0).detach().cpu().numpy(),
            'RNA_share_representations': torch.cat(rna_share_representations, dim=0).detach().cpu().numpy(),
            'ATAC_share_representations': torch.cat(atac_share_representations, dim=0).detach().cpu().numpy(),
        }
        print("✅ Embeddings extracted successfully. Returning a dictionary of representations.")
        return output

    def predict_cross_omics(self,
                           input_adata: anndata.AnnData,
                           omicsGenerated: str,
                           device: torch.device = None,
                           batch_size: int = 128, ):
        """
        Performs cross-modal prediction, generating one omics data from another.

        Args:
            input_adata (anndata.AnnData): The AnnData object of the source modality.
            omicsGenerated (str): The target modality to generate ('RNA' or 'ATAC').
            device: The device to use for inference.
            batch_size (int): The batch size for processing data.

        Returns:
            A numpy array of the generated omics data.
        """

        #   check omicsGenerated type
        valid_omics = {'RNA', 'ATAC'}
        if omicsGenerated not in valid_omics:
            raise ValueError(f"omicsGenerated must be one of {valid_omics}, but got '{omicsGenerated}'")
        input_omics = 'ATAC' if omicsGenerated == 'RNA' else 'RNA'

        print(f"🧬 Starting cross-omics prediction: {input_omics} -> {omicsGenerated}...")

        if device is None:
            device = self.device
        print(f"  - Using device: {device}")
        self.to(device)
        self.eval()

        # 2. Data Preparation
        print("  - Preparing data loader for prediction...")
        if self.batch_key in input_adata.obs.columns:
            batches_info = input_adata.obs[self.batch_key].values.tolist()
        else:
            warnings.warn("batch_key not found or not defined. Using a single batch for all cells.")
            batches_info = np.zeros(input_adata.shape[0])

        one_hot_encoded_batches = pd.get_dummies(batches_info, prefix='batch', dtype=float)
        input_data_np = input_adata.X.toarray()
        dataset = scOmicsDataset(input_data_np, one_hot_encoded_batches)
        data_loader = DataLoader(dataset, batch_size=batch_size, shuffle=False)

        all_omics_reconstructed = []
        with torch.no_grad():
            tepoch = tqdm(data_loader, desc="Predicting", unit="batch")
            for input_data_batch, batch_id_batch in tepoch:
                # 将数据移动到指定设备
                input_data_batch = input_data_batch.to(device)
                batch_id_batch = batch_id_batch.to(device)

                # 根据 omicsGenerated 调用相应的生成函数
                if omicsGenerated == 'RNA':
                    # 假设此方法存在于您的模型类中：输入ATAC，生成RNA
                    omics_reconstructed = self.rna_generation_from_atac(input_data_batch, batch_id_batch)
                elif omicsGenerated == 'ATAC':
                    # 假设此方法存在于您的模型类中：输入ATAC，重构ATAC
                    omics_reconstructed = self.atac_generation_from_rna(input_data_batch, batch_id_batch)

                # 将结果移回CPU并存入列表
                all_omics_reconstructed.append(omics_reconstructed.cpu())

        # 6. 合并结果并返回
        # 在循环外进行一次拼接，效率更高
        final_reconstruction = torch.cat(all_omics_reconstructed, dim=0).numpy()
        print(f"✅ Cross-omics prediction complete. Returning generated {omicsGenerated} matrix.")
        return final_reconstruction

    @staticmethod
    def contrastive_loss(anchor, positive, negative, margin=0.5):
        pos_dist = torch.norm(anchor - positive, p=2, dim=1)
        neg_dist = torch.norm(anchor - negative, p=2, dim=1)
        loss = torch.mean(torch.relu(pos_dist - neg_dist + margin))
        return loss

    @staticmethod
    def contrastive_loss_with_similarity(anchor, positive, negative, margin=1.0):
        pos_sim = F.cosine_similarity(anchor, positive, dim=1)
        neg_sim = F.cosine_similarity(anchor, negative, dim=1)
        loss = torch.exp(pos_sim) / (torch.exp(pos_sim) + torch.exp(neg_sim))
        loss = - torch.log(loss).mean()
        # loss = torch.mean(torch.relu(neg_sim - pos_sim + margin))
        return loss


class explainModel(nn.Module):
    def __init__(self,
                 scModel: scTFBridge,
                 mode,
                 latent_dim,
                 gene_ids,
                 cluster_num=None,
                 ):
        super(explainModel, self).__init__()
        self.scDM = scModel
        self.mode = mode
        self.gene_ids = gene_ids
        if cluster_num is not None:
            self.classifier = LinearLayer(latent_dim, cluster_num)

    def forward(self, x):
        rna_input, atac_input, batch_id = torch.split(x, [self.scDM.input_dims[0], self.scDM.input_dims[1],
                                                          self.scDM.batch_dims], dim=1)

        #   cross modal translation
        if self.mode == 'RNA2ATAC':
            atac_cross_recon = self.scDM.atac_generation_from_rna(rna_input, batch_id)
            return atac_cross_recon

        if self.mode == 'ATAC2RNA':
            rna_cross_recon = self.scDM.rna_generation_from_atac(atac_input, batch_id)
            return rna_cross_recon[:, self.gene_ids]

        #   regulatory generation
        if self.mode == 'RNARegulatory':
            rna_recon, _ = self.scDM.modal_generation([rna_input, atac_input], batch_id)
            return rna_recon[:, self.gene_ids]

        #   cluster mark gene
        if self.mode == 'joint_share':
            output = self.scDM([rna_input, atac_input], batch_id)['share_embedding']
            # output = output.detach().cpu().numpy()
            # prediction_output = self.classifier(output)
            return output[:, 4]

        if self.mode == 'rna_private':
            output = self.scDM([rna_input, atac_input], batch_id)['rna_private_embedding']
            # output = output.detach().cpu().numpy()
            # prediction_output = self.classifier(output)
            return output

        if self.mode == 'atac_private':
            output = self.scDM([rna_input, atac_input], batch_id)['atac_private_embedding']
            output = output.detach().cpu().numpy()
            prediction_output = self.classifier(output)
            return prediction_output


class explainModelLatentZ(nn.Module):
    def __init__(self,
                 scModel: scTFBridge,
                 mode,
                 latent_dim,
                 dimension_num=None
                 ):
        super(explainModelLatentZ, self).__init__()
        self.scDM = scModel
        self.mode = mode
        self.latent_dim = latent_dim
        self.dimension_num = dimension_num

    def forward(self, x):
        share_embedding, private_embedding, batch_id = torch.split(x, [self.latent_dim, self.latent_dim,
                                                                       self.scDM.batch_dims], dim=1)
        output = self.scDM.modal_generation_from_latent(share_embedding, private_embedding, batch_id, self.mode)
        if self.dimension_num is not None:
            output = output[:, [self.dimension_num]]
        return output


class explainModelLatentTF(nn.Module):
    def __init__(self,
                 scModel: scTFBridge,
                 mode,
                 latent_dim,
                 dimension_num=None
                 ):
        super(explainModelLatentTF, self).__init__()
        self.scDM = scModel
        self.mode = mode
        self.latent_dim = latent_dim
        self.dimension_num = dimension_num

    def forward(self, x):
        TF_embedding, batch_id = torch.split(x, [514, self.scDM.batch_dims], dim=1)
        output = self.scDM.modal_generation_from_latentTF(TF_embedding, batch_id, self.mode)

        if self.dimension_num is not None:
            output = output[:, [self.dimension_num]]
        return output


def get_latent_z(model: scTFBridge,
                 dataloader,
                 modal,
                 num_samples=4):
    device = next(model.parameters()).device
    share_embedding, private_embedding, sample_batch = [], [], []
    input = []
    with tqdm(dataloader, unit='batch') as tepoch:
        with torch.no_grad():
            for batch, data in enumerate(tepoch):
                rna_data, atac_data, TF_data, batch_id = data
                rna_data = rna_data.to(device)
                atac_data = atac_data.to(device)
                TF_data = TF_data.to(device)
                batch_id = batch_id.to(device)

                output = model([rna_data, atac_data, TF_data], batch_id)

                share_embedding.append(output[f'share_embedding'])
                private_embedding.append(output[f'{modal}_private_embedding'])
                sample_batch.append(batch_id)
                if modal == 'RNA':
                    input.append(rna_data)
                else:
                    input.append(atac_data)

                if len(share_embedding) * atac_data.shape[0] >= num_samples:
                    break

            share_data = torch.cat(share_embedding, dim=0)
            private_data = torch.cat(private_embedding, dim=0)
            batch_id = torch.cat(sample_batch, dim=0)
            all_data = torch.cat((share_data, private_data, batch_id), dim=1)
            rna_input = torch.cat(input, dim=0)
            return all_data, rna_input


def get_sample_data(dataloader, num_samples=4, device=torch.device('cpu')):
    samples_rna, samples_atac, sample_batch = [], [], []
    with tqdm(dataloader, unit='batch') as tepoch:
        for batch, data in enumerate(tepoch):
            rna_data, atac_data, TF_data, batch_id = data
            samples_rna.append(rna_data)
            samples_atac.append(atac_data)
            sample_batch.append(batch_id)
            if len(samples_rna) * rna_data.shape[0] >= num_samples:
                break
        rna_data = torch.cat(samples_rna, dim=0)
        atac_data = torch.cat(samples_atac, dim=0)
        batch_id = torch.cat(sample_batch, dim=0)
        all_data = torch.cat((rna_data, atac_data, batch_id), dim=1).to(device)
        rna_data = rna_data.to(device)
        return all_data, rna_data


def explain_TF2TG(scModel: scTFBridge,
                  adata_list: list[anndata.AnnData],
                  cell_type: str,
                  cell_key: str,
                  batch_key: str,
                  device=torch.device('cpu'),
                  using_sample=100
                  ):
    """
    Calculates the regulatory influence of Transcription Factors (TFs) on Target Genes (TGs)
    using attribution methods.
    """
    print(f"🧬 Starting TF-to-Target Gene (TF2TG) explanation for cell type: '{cell_type}'...")

    TF_length = scModel.latent_dim
    scModel.to(device)
    print(f"  - Using device: {device}")

    # --- Data Preparation ---
    rna_data = adata_list[0]
    atac_data = adata_list[1]
    TF_data = adata_list[2]

    scModel.latent_mode = 'latent_z'

    if cell_type != 'cellular':
        print(f"  - Filtering data for cell type: '{cell_type}'")
        rna_data = rna_data[rna_data.obs[cell_key] == cell_type]
        atac_data = atac_data[atac_data.obs[cell_key] == cell_type]
        TF_data = TF_data[TF_data.obs[cell_key] == cell_type]

    if batch_key in rna_data.obs.columns:
        batches_info = rna_data.obs[batch_key].values.tolist()
    else:
        batches_info = [0 for _ in range(rna_data.shape[0])]

    print("  - Preparing data loaders...")
    one_hot_encoded_batches = pd.get_dummies(batches_info, prefix='batch', dtype=float)

    rna_data_np = rna_data.X.toarray()
    atac_data_np = atac_data.X.toarray()
    TF_adata_np = TF_data.X.toarray()

    dim1 = rna_data.X.shape[1]

    dataset = scMultiDataset([rna_data_np, atac_data_np, TF_adata_np], one_hot_encoded_batches)
    data_loader = DataLoader(dataset, batch_size=5, shuffle=True)

    # --- Model and Explainer Initialization ---
    print("  - Initializing explanation model and background samples...")
    explain_model = explainModelLatentZ(scModel, 'rna', scModel.latent_dim, 0)
    explain_model.eval()
    explain_model.to(device)
    all_attributions = []

    test_latent_z, test_rna = get_latent_z(scModel, data_loader, num_samples=using_sample, modal='RNA')

    # --- Attribution Calculation ---
    print(f"  - Calculating attributions for {dim1} target genes. This may take a while...")
    for dim in range(dim1):
        explain_model.dimension_num = dim

        def model_loss_wrapper(z):
            rna_recon = explain_model(z)
            return F.mse_loss(rna_recon, test_rna[:, [dim]], reduction='none').mean(1).view(-1, 1)

        explainer = PathExplainerTorch(model_loss_wrapper)

        # define a baseline, in this case the zeros vector
        baseline_data = test_latent_z[0, :]
        baseline_data = torch.zeros_like(baseline_data)
        baseline_data.requires_grad = True
        attributions = explainer.attributions(test_latent_z,
                                              baseline=baseline_data,
                                              num_samples=using_sample,
                                              use_expectation=False)
        # 将attributions合并到all_attributions
        attributions = attributions.cpu().detach()
        if len(all_attributions) == 0:
            all_attributions = attributions.unsqueeze(0)  # 添加第一个维度
        else:
            all_attributions = torch.cat((all_attributions, attributions.unsqueeze(0)), dim=0)
        torch.cuda.empty_cache()


    # --- Final Processing ---
    print("  - Aggregating results...")
    all_attributions = all_attributions.numpy()
    shap_values = np.array(all_attributions)
    shap_values = shap_values[:, :, :TF_length]
    mean_shap_values = np.mean(np.abs(shap_values), axis=1)

    print(f"✅ TF2TG explanation complete for '{cell_type}'. Returning mean absolute SHAP values (Genes x TFs).")
    return mean_shap_values


def explain_RE2TG(scModel: scTFBridge,
                  adata_list: list[anndata.AnnData],
                  dataset_name: str,
                  use_gene_list: list[str],
                  cell_type: str,
                  cell_key: str,
                  batch_key: str,
                  tf_binding_path:str,
                  device=torch.device('cpu'), 
                  ):
    """
    Infers cis-regulatory links between Regulatory Elements (REs) and Target Genes (TGs)
    by integrating attribution scores, genomic distance, and molecular data.
    """
    import pybedtools
    print(f"🔗 Starting RE-to-Target Gene (RE2TG) explanation for cell type: '{cell_type}'...")

    scModel.to(device)
    print(f"  - Using device: {device}")

    # --- Data Preparation ---
    rna_data = adata_list[0]
    atac_data = adata_list[1]
    TF_data = adata_list[2]

    scModel.latent_mode = 'latent_z'

    if cell_type != 'cellular':
        print(f"  - Filtering data for cell type: '{cell_type}'")
        rna_data = rna_data[rna_data.obs[cell_key] == cell_type].copy()
        atac_data = atac_data[atac_data.obs[cell_key] == cell_type].copy()
        TF_data = TF_data[TF_data.obs[cell_key] == cell_type].copy()

    if batch_key in rna_data.obs.columns:
        batches_info = rna_data.obs[batch_key].values.tolist()
    else:
        batches_info = [0 for _ in range(rna_data.shape[0])]

    print("  - Preparing data loaders...")
    one_hot_encoded_batches = pd.get_dummies(batches_info, prefix='batch', dtype=float)
    rna_data_np = rna_data.X.toarray()
    atac_data_np = atac_data.X.toarray()
    TF_adata_np = TF_data.X.toarray()

    dim1 = rna_data.X.shape[1]
    dim2 = atac_data.X.shape[1]
    dataset = scMultiDataset([rna_data_np, atac_data_np, TF_adata_np], one_hot_encoded_batches)
    data_loader = DataLoader(dataset, batch_size=5, shuffle=True)

    # --- Model and Explainer Initialization ---
    print("  - Initializing explanation model (ATAC2RNA) and background samples...")
    explain_model = explainModel(scModel, 'ATAC2RNA', 128, 0)
    explain_model.eval()
    test_all_data, test_rna_data = get_sample_data(data_loader, num_samples=50, device=device)


    # --- Attribution Calculation ---
    gene_ids = []
    common_geneSymbol = []
    gex_data_var = rna_data.var.copy()
    for gene in use_gene_list:
        # TSS = int(gene_eqtl[gene_eqtl['GeneSymbol'] == gene]['GenePos'].tolist()[0])
        if gene in rna_data.var.index:
            gene_id = gex_data_var.index.get_loc(gene)
            gene_ids.append(gene_id)
            common_geneSymbol.append(gene)
    gene_ids = sorted(gene_ids)

    print(f"  - Calculating attributions for {len(gene_ids)} specified target genes...")
    all_attributions = []
    for gene in gene_ids:
        explain_model.gene_ids = [gene]
        def model_loss_wrapper(z):
            rna_recon = explain_model(z)
            return rna_recon
            # return F.mse_loss(rna_recon, test_rna_data[:, [gene]], reduction='none').mean(1).view(-1, 1)

        explainer = PathExplainerTorch(model_loss_wrapper)
        baseline_data = test_all_data[0, :]  # define a baseline, in this case the zeros vector
        baseline_data = torch.zeros_like(baseline_data)
        baseline_data.requires_grad = True
        attributions = explainer.attributions(test_all_data,
                                              baseline=baseline_data,
                                              num_samples=50,
                                              use_expectation=False)
        # 将attributions合并到all_attributions
        attributions = attributions.cpu().detach()
        if len(all_attributions) == 0:
            all_attributions = attributions.unsqueeze(0)  # 添加第一个维度
        else:
            all_attributions = torch.cat((all_attributions, attributions.unsqueeze(0)), dim=0)
        torch.cuda.empty_cache()  # Clear cache to reduce memory consumption

    # --- Score Integration ---
    print("  - Aggregating SHAP values (Genes x Peaks)...")
    all_attributions = all_attributions.numpy()
    shap_values = np.array(all_attributions)
    shap_values = shap_values[:, :, dim1:dim1 + dim2]
    shap_values = np.abs(shap_values)
    shap_values = shap_values.mean(axis=1)

    print("  - Loading genomic locations and calculating RE-TSS distances...")
    RE_region = pybedtools.example_bedtool(os.path.abspath(f'{tf_binding_path}/Region.bed')).sort()
    gene_infor = pd.read_csv('gene_genome_info.csv')

    gene_feature = rna_data.var.index.values
    atac_peak_feature = atac_data.var.index.values.tolist()
    peak_mean_access = atac_data.X.mean(axis=0).tolist()
    gene_mean_exp = rna_data.X.mean(axis=0).tolist()

    # print(peak_mean_access)
    peak_info = atac_data.var
    peak_info['mean_access'] = peak_mean_access[0]
    peak_info['peak'] = peak_info.index

    gene_info = rna_data.var
    gene_info['mean_exp'] = gene_mean_exp[0]

    # print(peak_mean_access)
    gene_eqtl_feature = gene_infor['GeneSymbol'].unique()
    common_features = list(set(common_geneSymbol) & set(gene_eqtl_feature))

    gene_infos = gene_infor[gene_infor['GeneSymbol'].isin(common_features)]
    gene_infos['GenePos+'] = gene_infos['GenePos'] + 1
    gene_infos_bed = gene_infos[['GeneChr', 'GenePos', 'GenePos+', 'GeneSymbol']]
    gene_infos_bed['GeneChr'] = 'chr' + gene_infos_bed['GeneChr'].astype(str)
    gene_infos_bed.sort_values(['GeneChr'], inplace=True)
    gene_infos_bed.to_csv(f'cis_regulatory/{cell_type}_gene_eqtl_bed.bed', index=False, header=False, sep='\t')
    # print(gene_eqtl_bed)

    gene_eqtl_bed = pybedtools.example_bedtool(os.path.abspath(f'cis_regulatory/{cell_type}_gene_eqtl_bed.bed')).sort()

    print('computing tss distance')
    closest = RE_region.closest(gene_eqtl_bed, d=True, k=len(common_features))
    closest.saveas(f'cis_regulatory/{cell_type}_gene_RE_distance.bed')

    shap_all = pd.DataFrame(shap_values.T, columns=common_geneSymbol, index=atac_peak_feature)

    TSS_distance = pd.read_csv(f'cis_regulatory/{cell_type}_gene_RE_distance.bed', sep='\t')

    TSS_distance.columns = ['RE_chr', 'RE_start', 'RE_end', 'TSS_chr', 'TSS_start', 'TSS_end', 'gene_symbol',
                            'TSS_distance']
    TSS_distance['peak'] = TSS_distance['RE_chr'] + ':' + TSS_distance['RE_start'].astype(str) + '-' + TSS_distance[
        'RE_end'].astype(str)
    TSS_distance.replace('.', np.nan, inplace=True)
    TSS_distance.dropna(inplace=True)

    print("  - Merging SHAP scores with distance and expression data...")
    # print("Reshaping SHAP values for efficient merge...")
    shap_long = shap_all.stack().reset_index()
    shap_long.columns = ['peak', 'gene_symbol', 'score']

    # 2. 使用pd.merge将分数合并到TSS_distance中
    # print("Merging SHAP scores...")
    TSS_distance = pd.merge(TSS_distance, shap_long, on=['peak', 'gene_symbol'], how='inner')
    # 使用 'inner' join可以只保留在两个DataFrame中都存在的peak-gene对，更高效

    # 3. 使用.map()高效添加mean_access和mean_gene_exp
    # print("Mapping mean accessibility and gene expression...")

    # 创建映射关系 (Series)
    peak_access_mapping = peak_info.set_index('peak')['mean_access']
    gene_exp_mapping = gene_info['mean_exp']  # gene_info的索引已经是基因名

    # 执行映射
    TSS_distance['mean_access'] = TSS_distance['peak'].map(peak_access_mapping)
    TSS_distance['mean_gene_exp'] = TSS_distance['gene_symbol'].map(gene_exp_mapping)

    TSS_distance['score'] = (TSS_distance['score'] * np.exp(-TSS_distance['TSS_distance'] / 30000)) * TSS_distance[
        'mean_access'] * TSS_distance['mean_gene_exp']
    print(f"✅ RE2TG explanation complete for '{cell_type}'. Returning DataFrame with integrated scores.")
    return TSS_distance


def explain_DisLatent(scModel: scTFBridge,
                      adata_list: list[anndata.AnnData],
                      omics: str,
                      cell_type: str,
                      cell_key: str,
                      batch_key: str,
                      device=torch.device('cpu'), ):
    #   check omicsGenerated type
    valid_omics = {'RNA', 'ATAC'}
    if omics not in valid_omics:
        raise ValueError(f"omicsGenerated must be one of {valid_omics}, but got '{omics}'")

    print(f"🔍 Starting latent space contribution analysis for '{omics}' omics and cell type: '{cell_type}'...")
    latent_dim = scModel.latent_dim
    scModel.to(device)
    print(f"  - Using device: {device}")
    rna_data = adata_list[0]
    atac_data = adata_list[1]
    TF_data = adata_list[2]

    scModel.latent_mode = 'latent_z'

    if cell_type != 'cellular':
        print(f"Calculating {omics} private-share contribution value for: {cell_type}")
        rna_data = rna_data[rna_data.obs[cell_key] == cell_type]
        atac_data = atac_data[atac_data.obs[cell_key] == cell_type]
        TF_data = TF_data[TF_data.obs[cell_key] == cell_type]

    if batch_key in rna_data.obs.columns:
        batches_info = rna_data.obs[batch_key].values.tolist()
    else:
        batches_info = [0 for _ in range(rna_data.shape[0])]

    print("  - Preparing data loaders...")
    one_hot_encoded_batches = pd.get_dummies(batches_info, prefix='batch', dtype=float)
    rna_data_np = rna_data.X.toarray()
    atac_data_np = atac_data.X.toarray()
    TF_adata_np = TF_data.X.toarray()

    dim1 = rna_data.X.shape[1]
    dim2 = atac_data.X.shape[1]

    gene_name = rna_data.var.index.tolist()
    peak_name = atac_data.var.index.tolist()

    # 3. 创建 Dataset 和 DataLoader (使用全部数据)
    dataset = scMultiDataset([rna_data_np, atac_data_np, TF_adata_np], one_hot_encoded_batches)
    data_loader = DataLoader(dataset, batch_size=5, shuffle=True)

    # --- Model and Explainer Initialization ---
    print("  - Initializing explanation model and background samples...")
    train_latent_z, train_omics = get_latent_z(scModel, data_loader, omics, num_samples=10)
    test_latent_z, test_omics = get_latent_z(scModel, data_loader, omics, num_samples=50)

    if omics == 'RNA':
        explain_model = explainModelLatentZ(scModel, 'rna', latent_dim, 0)
        omics_dim = dim1
        feature_list = gene_name

    else:
        explain_model = explainModelLatentZ(scModel, 'atac', latent_dim, 0)
        omics_dim = dim2
        feature_list = peak_name

    print(f"  - Calculating private vs. shared contributions for {omics_dim} features. This may take a while...")
    explain_model.to(device)
    private_contribute = []
    share_contribute = []

    for i in range(omics_dim):
        explain_model.dimension_num = i

        def model_rna_loss_wrapper(z):
            rna_recon = explain_model(z)
            return F.mse_loss(rna_recon, test_omics[:, [i]], reduction='none').mean(1).view(-1, 1)

        def model_atac_loss_wrapper(z):
            rna_recon = explain_model(z)
            return F.binary_cross_entropy_with_logits(rna_recon, test_omics[:, [i]], reduction='none').mean(1).view(-1,
                                                                                                                    1)

        baseline_data = train_latent_z[0, :]  # define a baseline, in this case the zeros vector
        baseline_data.requires_grad = True

        if omics == 'RNA':
            explainer = PathExplainerTorch(model_rna_loss_wrapper)
        else:
            explainer = PathExplainerTorch(model_atac_loss_wrapper)

        attributions = explainer.attributions(test_latent_z,
                                              baseline=baseline_data,
                                              num_samples=50,
                                              # number of samples to use when calculating the path integral
                                              use_expectation=False)

        shap_values = attributions.detach().cpu().numpy()

        shap_share = shap_values[:, :latent_dim]
        shap_private = shap_values[:, latent_dim:latent_dim * 2]
        # TF_contribution = np.abs(shap_share).mean(0)

        shap_share_sum = np.sum(np.abs(shap_share).mean(0))
        shap_private_sum = np.sum(np.abs(shap_private).mean(0))
        share_contribute.append(shap_share_sum / (shap_share_sum + shap_private_sum))
        private_contribute.append(shap_private_sum / (shap_private_sum + shap_share_sum))

        # print(f'-----{i}-----')
        # print('share', shap_share_sum, shap_share_sum / (shap_share_sum + shap_private_sum))
        # print('private', shap_private_sum, shap_private_sum / (shap_private_sum + shap_share_sum))
        # print('input', shap_input_sum)

    df = {
        'feature_name': feature_list,
        'private_embedding_contribute': private_contribute,
        'share_embedding_contribute': share_contribute,
    }
    df = pd.DataFrame(df)
    print(f"✅ Latent space contribution analysis complete for '{omics}' omics.")
    return df


def calculate_r_squared_torch(y_true, y_pred):
    # 计算总平方和 (TSS)
    tss = torch.sum((y_true - torch.mean(y_true, axis=0)) ** 2, axis=0)
    # 计算残差平方和 (RSS)
    rss = torch.sum((y_true - y_pred) ** 2, axis=0)
    # 计算R平方值
    r_squared = 1 - (rss / tss)
    return r_squared


def calculate_pcc_torch(x, y):
    """
    Calculates the Pearson Correlation Coefficient between two PyTorch tensors.

    Args:
        x (torch.Tensor): The first tensor.
        y (torch.Tensor): The second tensor.

    Returns:
        torch.Tensor: The Pearson Correlation Coefficient.
    """
    if x.ndim > 1:
        x = x.squeeze()
    if y.ndim > 1:
        y = y.squeeze()

    if x.shape != y.shape:
        raise ValueError("Input tensors must have the same shape")

    vx = x - torch.mean(x)
    vy = y - torch.mean(y)

    numerator = torch.sum(vx * vy)
    denominator = torch.sqrt(torch.sum(vx ** 2)) * torch.sqrt(torch.sum(vy ** 2))

    if denominator == 0:
        # Return 0 if the standard deviation of either variable is 0
        # (e.g., if all values in a tensor are the same)
        # This avoids division by zero.
        # Other conventions might return NaN or raise an error.
        return torch.tensor(0.0)

    return numerator / denominator
