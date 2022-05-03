
import os
import argparse
import anndata as ad
from collections import defaultdict
import scanpy as sc
import pandas as pd
import STAGATE_pyG
from sklearn.metrics import adjusted_rand_score
import numpy as np
import torch
from memory_profiler import memory_usage, profile

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--gpu',  type=int, default=0)
parser.add_argument('--n_samples', type=int)

args = parser.parse_args()

SAMPLES = {'151507': 7, 
           '151508': 7, 
           '151509': 7, 
           '151510': 7, 
           '151669': 5, 
           '151670': 5, 
           '151671': 5, 
           '151672': 5, 
           '151673': 7, 
           '151674': 7, 
           '151675': 7, 
           '151676': 7
}

hvgs=2000
hidden_dim=512
n_latent=30

adata_list = []
for sample, n_clusters in SAMPLES.items():
    input_dir = os.path.join('/work/FAC/FBM/DBC/gciriell/spacegene/Data/jhpce_human_pilot_10x', sample)
    adata = ad.read_h5ad(os.path.join(input_dir, f'{sample}.h5ad'))
    adata.obs['sample'] = [sample]*adata.shape[0]
    adata.obs['layer_guess'] = pd.read_csv(f'{input_dir}/metadata.tsv', sep='\t').loc[adata.obs_names, 'layer_guess']
    STAGATE_pyG.Cal_Spatial_Net(adata, rad_cutoff=150)
    adata_list.append(adata)



def run(adata_list, n_samples):
    adata = ad.concat(adata_list[:n_samples], pairwise=True)
    adata.uns['Spatial_Net'] = pd.concat([a.uns['Spatial_Net'] for a in adata_list]) 

    sc.pp.filter_genes(adata, min_counts=3)
    sc.pp.filter_cells(adata, min_counts=3)

    adata.layers["counts"] = adata.X.copy()

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=hvgs,
        subset=True,
        layer="counts",
        flavor="seurat_v3",

    )

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    print('N. samples', len(adata.obs['sample'].unique()))

    # adata = ad.concat(adata_list[:n_samples], pairwise=True)

    # adata.obsm['spatial'][:, 0][adata.obs['sample'] == "151509"] = 15000 + adata.obsm['spatial'][:, 0][adata.obs['sample'] == "151509"]

    # adata.obsm['spatial'][:, 1][adata.obs['sample'] == "151510"] = 15000 + adata.obsm['spatial'][:, 1][adata.obs['sample'] == "151510"]

    # adata.obsm['spatial'][:, 0][adata.obs['sample'] == "151669"] = 15000 + adata.obsm['spatial'][:, 0][adata.obs['sample'] == "151669"]
    # adata.obsm['spatial'][:, 1][adata.obs['sample'] == "151669"] = 15000 + adata.obsm['spatial'][:, 1][adata.obs['sample'] == "151669"]

    # adata.obsm['spatial'][:, 0][adata.obs['sample'] == "151670"] = 30000 + adata.obsm['spatial'][:, 0][adata.obs['sample'] == "151670"]

    # adata.obsm['spatial'][:, 1][adata.obs['sample'] == "151671"] = 30000 + adata.obsm['spatial'][:, 1][adata.obs['sample'] == "151671"]

    # adata.obsm['spatial'][:, 0][adata.obs['sample'] == "151674"] = 30000 + adata.obsm['spatial'][:, 0][adata.obs['sample'] == "151674"]
    # adata.obsm['spatial'][:, 1][adata.obs['sample'] == "151674"] = 15000 + adata.obsm['spatial'][:, 1][adata.obs['sample'] == "151674"]


    # adata.obsm['spatial'][:, 0][adata.obs['sample'] == "151675"] = 15000 + adata.obsm['spatial'][:, 0][adata.obs['sample'] == "151675"]
    # adata.obsm['spatial'][:, 1][adata.obs['sample'] == "151675"] = 30000 + adata.obsm['spatial'][:, 1][adata.obs['sample'] == "151675"]

    # adata.obsm['spatial'][:, 0][adata.obs['sample'] == "151676"] = 30000 + adata.obsm['spatial'][:, 0][adata.obs['sample'] == "151676"]
    # adata.obsm['spatial'][:, 1][adata.obs['sample'] == "151676"] = 30000 + adata.obsm['spatial'][:, 1][adata.obs['sample'] == "151676"]

    # adata.obsm['spatial'][:, 0][adata.obs['sample'] == "151507"] = 45000 + adata.obsm['spatial'][:, 0][adata.obs['sample'] == "151507"]

    # adata.obsm['spatial'][:, 1][adata.obs['sample'] == "151672"] = 45000 + adata.obsm['spatial'][:, 1][adata.obs['sample'] == "151672"]

    # adata.obsm['spatial'][:, 0][adata.obs['sample'] == "151673"] = 45000 + adata.obsm['spatial'][:, 0][adata.obs['sample'] == "151673"]
    # adata.obsm['spatial'][:, 1][adata.obs['sample'] == "151673"] = 15000 + adata.obsm['spatial'][:, 1][adata.obs['sample'] == "151673"]

    # sc.pp.filter_genes(adata, min_counts=3)
    # sc.pp.filter_cells(adata, min_counts=3)

    # adata.layers["counts"] = adata.X.copy()

    # sc.pp.highly_variable_genes(
    #     adata,
    #     n_top_genes=hvgs,
    #     subset=True,
    #     layer="counts",
    #     flavor="seurat_v3",

    # )

    # sc.pp.normalize_total(adata, target_sum=1e4)
    # sc.pp.log1p(adata)

    adata = STAGATE_pyG.train_STAGATE(adata, hidden_dims=[hidden_dim, n_latent])
    adata = STAGATE_pyG.mclust_R(adata, used_obsm='STAGATE', num_cluster=n_clusters)


mem = max(memory_usage(proc=(run, [adata_list, args.n_samples])))
mems_path = f"/work/FAC/FBM/DBC/gciriell/spacegene/Packages/cellcharter_analyses/results/dlpfc/STAGATE/memory/memory_hvgs{hvgs}_hidden_dim{hidden_dim}_nlatent{n_latent}_{'gpu' if args.gpu else 'cpu'}.csv"
row = pd.DataFrame([mem], index=[args.n_samples], )
if os.path.exists(mems_path):
    mems_df = pd.read_csv(mems_path, index_col=0)
    mems_df.columns = mems_df.columns.astype(int)
    mems_df = pd.concat((mems_df, row))
else:
    mems_df = row

mems_df.to_csv(mems_path)


