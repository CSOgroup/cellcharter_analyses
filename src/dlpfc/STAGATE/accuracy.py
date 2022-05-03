
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

os.environ['R_HOME'] = '/dcsrsoft/spack/hetre/v1.2/spack/opt/spack/linux-rhel8-zen2/gcc-9.3.0/r-4.0.5-ipkqr7imy7lxmyo2gfyhpccb3ooinj25/rlib/R'
os.environ['R_USER'] = '/work/FAC/FBM/DBC/gciriell/spacegene/envs/stagate/lib/python3.9/site-packages/rpy2'

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--gpu',  type=int, default=0)

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


hvgs=5000
hidden_dim=1024
n_latent=30

aris = defaultdict(list)

for sample, n_clusters in SAMPLES.items():
    input_dir = os.path.join('/work/FAC/FBM/DBC/gciriell/spacegene/Data/jhpce_human_pilot_10x', sample)
    adata = ad.read_h5ad(f'/work/FAC/FBM/DBC/gciriell/spacegene/Data/jhpce_human_pilot_10x/{sample}.h5ad')

    adata = sc.read_visium(path=input_dir)
    adata.var_names_make_unique()

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

    metadata = pd.read_csv(f'{input_dir}/metadata.tsv', sep='\t')

    adata.obs['layer_guess'] = metadata.loc[adata.obs_names, 'layer_guess']
    

    STAGATE_pyG.Cal_Spatial_Net(adata, rad_cutoff=150)
    STAGATE_pyG.Stats_Spatial_Net(adata)

    rng = np.random.default_rng(12345)
    seeds = rng.integers(low=0, high=32768, size=10)
    for i in seeds:
        adata = STAGATE_pyG.train_STAGATE(adata, hidden_dims=[hidden_dim, n_latent], random_seed=i)
        adata = STAGATE_pyG.mclust_R(adata, used_obsm='STAGATE', num_cluster=n_clusters)

        obs_df = adata.obs.dropna()
        ari = adjusted_rand_score(obs_df['mclust'], obs_df['layer_guess'])
        
        aris[sample].append(ari)
        aris_df = pd.DataFrame.from_dict(aris, orient='index')
        aris_df.to_csv(f"/work/FAC/FBM/DBC/gciriell/spacegene/Packages/cellcharter_analyses/results/dlpfc/STAGATE/accuracy/ARI_hvgs{hvgs}_hidden_dim{hidden_dim}_n_latent{n_latent}_{'gpu' if args.gpu else 'cpu'}.csv")

    
