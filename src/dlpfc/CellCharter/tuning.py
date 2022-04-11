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

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--gpu',  type=int, default=0)

args = parser.parse_args()

SAMPLES = {
    '151507': 7, 
    '151672': 5, 
    '151673': 7, 
}

hvgs = [500, 1000, 2000, 5000]
n_latents = [5, 10, 15]
nhood_dists = [2, 3, 4]


aris = defaultdict(list)
times = defaultdict(list)

print('Loading data...')
    

for sample, n_clusters in SAMPLES.items():
    for hvg in hvgs:
        adata = ad.read_h5ad(f'/work/FAC/FBM/DBC/gciriell/spacegene/Data/jhpce_human_pilot_10x/{sample}/{sample}.h5ad')

        

        sc.pp.filter_genes(adata, min_counts=3)
        sc.pp.filter_cells(adata, min_counts=3)

        adata.layers["counts"] = adata.X.copy()

        sc.pp.normalize_total(adata, target_sum=1e6)
        sc.pp.log1p(adata)

        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=hvg,
            subset=True,
            layer="counts",
            flavor="seurat_v3",

        )

        #sc.pp.filter_genes(adata, min_counts=100)
        #sc.pp.calculate_qc_metrics(adata)

        #genes = adata.var_names[adata.var.highly_variable]
        #sq.gr.sepal(adata, max_neighs=6, genes=genes, n_jobs=1)
        #adata = adata[:, adata.uns["sepal_score"].index].copy()

        for n_latent in n_latents:
            for nhood_dist in nhood_dists:
                rng = np.random.default_rng(12345)
                seeds = rng.integers(low=0, high=32768, size=10)
                for i in seeds:
                    scvi.settings.seed = i 
                    scvi.model.SCVI.setup_anndata(adata, layer="counts")
                    model = scvi.model.SCVI(adata, n_latent=n_latent)
                    model.train(early_stopping=True, early_stopping_patience=15, plan_kwargs=dict(lr=0.01, optimizier='AdamW'))
                    adata.obsm['X_scVI'] = model.get_latent_representation(adata).astype(np.float32)

                    cc.tl.SpatialCluster.aggregate_neighbors(adata, nhood_dist, X_key='X_scVI')
                    
                    cls = cc.tl.SpatialCluster(int(n_clusters), random_state=i, gpus=args.gpu)
                    cls.fit(adata, X_key='X_cellcharter')
                    adata.obs['cluster_cellcharter'] = cls.predict(adata, X_key='X_cellcharter')
                    
                    
                    ari = adjusted_rand_score(adata[~adata.obs['sce.layer_guess'].isna()].obs['sce.layer_guess'], adata[~adata.obs['sce.layer_guess'].isna()].obs['cluster_cellcharter'])
                    print('ARI:', ari)
                    aris[f'{sample}_{hvg}_{n_latent}_{nhood_dist}'].append(ari)
                

                    aris_df = pd.DataFrame.from_dict(aris, orient='index')
                    aris_df.to_csv(f"/work/FAC/FBM/DBC/gciriell/spacegene/Packages/cellcharter_analyses/results/dlpfc/CellCharter/ARI_tuning_{sample}_lr.csv")