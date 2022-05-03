
import os
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



SAMPLES = {
    '151507': 7, 
    '151672': 5, 
    '151673': 7, 
}

hvgs = [1000, 2000, 5000]
hidden_dims = [256, 512, 1024]
n_latents = [5, 10, 30, 50]

aris = defaultdict(list)
times = defaultdict(list)

print('Loading data...')
    

for sample, n_clusters in SAMPLES.items():
    for hvg in hvgs:
        input_dir = os.path.join('/work/FAC/FBM/DBC/gciriell/spacegene/Data/jhpce_human_pilot_10x', sample)
        adata = ad.read_h5ad(f'/work/FAC/FBM/DBC/gciriell/spacegene/Data/jhpce_human_pilot_10x/{sample}.h5ad')

        adata = sc.read_visium(path=input_dir)
        adata.var_names_make_unique()

        adata.layers["counts"] = adata.X.copy()

        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=hvg,
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

        for hidden_dim in hidden_dims:
            for n_latent in n_latents:
                rng = np.random.default_rng(12345)
                seeds = rng.integers(low=0, high=32768, size=10)
                for i in seeds:
                    adata = STAGATE_pyG.train_STAGATE(adata, hidden_dims=[hidden_dim, n_latent], random_seed=i)
                    adata = STAGATE_pyG.mclust_R(adata, used_obsm='STAGATE', num_cluster=n_clusters)

                    obs_df = adata.obs.dropna()
                    ARI = adjusted_rand_score(obs_df['mclust'], obs_df['layer_guess'])

                    aris[f'{sample}_{hvg}_{hidden_dim}_{n_latent}'].append(ARI)
                

                    aris_df = pd.DataFrame.from_dict(aris, orient='index')
                    aris_df.to_csv(f"/work/FAC/FBM/DBC/gciriell/spacegene/Packages/cellcharter_analyses/results/dlpfc/STAGATE/ARI_tuning_{sample}.csv")