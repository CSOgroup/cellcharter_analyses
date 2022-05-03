import anndata as ad
from collections import defaultdict
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import cellchart as cc
from time import time
from sklearn.metrics import adjusted_rand_score
from sklearn.decomposition import TruncatedSVD
import argparse
import squidpy as sq
from memory_profiler import memory_usage, profile
import gc
import os

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

hvgs = 5000
n_latent = 5
nhood_layers = 4

adata_list = []
for sample, n_clusters in SAMPLES.items():
    adata = ad.read_h5ad(f'/work/FAC/FBM/DBC/gciriell/spacegene/Data/jhpce_human_pilot_10x/{sample}/{sample}.h5ad')
    adata.obs['sample'] = [sample]*adata.shape[0]
    adata_list.append(adata)

def run(n_samples):
    print('N. samples', n_samples)

    adata = ad.concat(adata_list[:n_samples], pairwise=True)    

    sc.pp.filter_genes(adata, min_counts=3)
    sc.pp.filter_cells(adata, min_counts=3)

    adata.layers["counts"] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=hvgs,
        subset=True,
        layer="counts",
        flavor="seurat_v3",

    )

    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key='sample')
    model = scvi.model.SCVI(adata, n_latent=n_latent)
    model.train(early_stopping=True)
    adata.obsm['X_scVI'] = model.get_latent_representation(adata).astype(np.float32)

    cc.tl.SpatialCluster.aggregate_neighbors(adata, nhood_layers, X_key='X_scVI', out_key='X_cellcharter')

    cls = cc.tl.SpatialCluster(7, random_state=42, gpus=args.gpu)
    cls.fit(adata, X_key='X_cellcharter')
    adata.obs[f'cluster_cellcharter'] = cls.predict(adata, X_key='X_cellcharter')


mem = max(memory_usage(proc=(run, [args.n_samples])))
mems_path = f"/work/FAC/FBM/DBC/gciriell/spacegene/Packages/cellcharter_analyses/results/dlpfc/CellCharter/memory/memory_hvgs{hvgs}_nlatent{n_latent}_nhoodlayers{nhood_layers}_{'gpu' if args.gpu else 'cpu'}.csv"
row = pd.DataFrame([mem], index=[args.n_samples], )
if os.path.exists(mems_path):
    mems_df = pd.read_csv(mems_path, index_col=0)
    mems_df.columns = mems_df.columns.astype(int)
    mems_df = pd.concat((mems_df, row))
else:
    mems_df = row

mems_df.to_csv(mems_path)


