import anndata as ad
from collections import defaultdict
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import cellcharter as cc
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
parser.add_argument('--n_cpus', type=int)

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
times = defaultdict(list)

times_path = f"/work/FAC/FBM/DBC/gciriell/spacegene/Packages/cellcharter_analyses/results/dlpfc/CellCharter/time/time_hvgs{hvgs}_nlatent{n_latent}_nhoodlayers{nhood_layers}_{'gpu' if args.gpu else 'cpu'}_ncpus{args.n_cpus}.csv"

adata_list = []
for sample, n_clusters in SAMPLES.items():
    adata = ad.read_h5ad(f'/work/FAC/FBM/DBC/gciriell/spacegene/Data/jhpce_human_pilot_10x/{sample}/{sample}.h5ad')
    adata.obs['sample'] = [sample]*adata.shape[0]
    adata_list.append(adata)


adata = ad.concat(adata_list[:args.n_samples], pairwise=True)    

sc.pp.filter_genes(adata, min_counts=3)
sc.pp.filter_cells(adata, min_counts=3)

adata.layers["counts"] = adata.X.copy()

sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)

start_hvg = time()
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=hvgs,
    subset=True,
    layer="counts",
    flavor="seurat_v3",

)
time_hvg = time() - start_hvg

rng = np.random.default_rng(12345)
seeds = rng.integers(low=0, high=32768, size=10)
for i in seeds:
    scvi.settings.seed = i 
    start_scvi = time()
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key='sample')
    model = scvi.model.SCVI(adata, n_latent=n_latent)
    model.train(early_stopping=True)
    adata.obsm['X_scVI'] = model.get_latent_representation(adata).astype(np.float32)
    time_scvi = time() - start_scvi

    start_neigh = time()
    cc.tl.SpatialCluster.aggregate_neighbors(adata, nhood_layers, X_key='X_scVI')
    time_neigh = time() - start_neigh

    start_cls = time()
    cls = cc.tl.SpatialCluster(7, random_state=i, gpus=args.gpu)
    cls.fit(adata, X_key='X_cellcharter')
    adata.obs[f'cluster_cellcharter'] = cls.predict(adata, X_key='X_cellcharter')
    time_cls = time() - start_cls

    row = pd.DataFrame([[args.n_samples, time_hvg, time_scvi, time_neigh, time_cls]], columns=['n_samples', 'HVGs', 'scVI', 'nhood_aggregation', 'clustering'])
    if os.path.exists(times_path):
        times_df = pd.read_csv(times_path, index_col=0)
        times_df = pd.concat((times_df, row))
    else:
        times_df = row

    times_df.to_csv(times_path)

