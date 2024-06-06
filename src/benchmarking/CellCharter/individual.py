import anndata as ad
from collections import defaultdict
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import cellcharter as cc
from sklearn.metrics import adjusted_rand_score
import argparse

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

hvgs = 5000
n_latent = 5
nhood_layers = 4

aris = defaultdict(list)

for sample, n_clusters in SAMPLES.items():
    adata = ad.read_h5ad(f'../../../data/visium_human_dlpfc/preprocessed_h5ad/{sample}.h5ad')

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

    rng = np.random.default_rng(12345)
    seeds = rng.integers(low=0, high=32768, size=10)
    for i in seeds:
        scvi.settings.seed = i 
        scvi.model.SCVI.setup_anndata(adata, layer="counts")
        model = scvi.model.SCVI(adata, n_latent=n_latent)
        model.train(early_stopping=True, plan_kwargs=dict(optimizer='AdamW', reduce_lr_on_plateau=True))
        adata.obsm['X_scVI'] = model.get_latent_representation(adata).astype(np.float32)

        cc.tl.SpatialCluster.aggregate_neighbors(adata, nhood_layers, X_key='X_scVI', out_key='X_cellcharter')
        
        cls = cc.tl.SpatialCluster(int(n_clusters), random_state=i, gpus=args.gpu)
        cls.fit(adata, X_key='X_cellcharter')
        adata.obs['cluster_cellcharter'] = cls.predict(adata, X_key='X_cellcharter')
        
        ari = adjusted_rand_score(adata[~adata.obs['sce.layer_guess'].isna()].obs['sce.layer_guess'], adata[~adata.obs['sce.layer_guess'].isna()].obs['cluster_cellcharter'])
        aris[sample].append(ari)
    

        aris_df = pd.DataFrame.from_dict(aris, orient='index')
        aris_df.to_csv(f"../../../results/benchmarking/individual/ARI_CellCharter_hvgs{hvgs}_nlatent{n_latent}_nhoodlayers{nhood_layers}_{'gpu' if args.gpu else 'cpu'}_individual.csv")