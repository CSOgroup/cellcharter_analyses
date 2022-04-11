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

SAMPLES = {#'151507': 7, 
           '151508': 7, 
           '151509': 7, 
           '151510': 7, 
           '151669': 5, 
           '151670': 5, 
           '151671': 5, 
           #'151672': 5, 
           #'151673': 7, 
           '151674': 7, 
           '151675': 7, 
           '151676': 7
}

GROUPS = {
    '151507': 0, 
    '151508': 0, 
    '151509': 0, 
    '151510': 0, 
    '151669': 1, 
    '151670': 1, 
    '151671': 1, 
    '151672': 1, 
    '151673': 2, 
    '151674': 2, 
    '151675': 2, 
    '151676': 2

}

hvgs = 5000
n_latent = 5
nhood_layers = 4


aris = defaultdict(list)
times = defaultdict(list)
adata_list = []
for sample, n_clusters in SAMPLES.items():
    adata = ad.read_h5ad(f'/work/FAC/FBM/DBC/gciriell/spacegene/Data/jhpce_human_pilot_10x/{sample}/{sample}.h5ad')
    adata.obs['sample'] = [sample]*adata.shape[0]
    adata.obs['group'] = [GROUPS[sample]]*adata.shape[0]
    adata_list.append(adata)

adata = ad.concat(adata_list, pairwise=True)   

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
for i, seed in enumerate(seeds):
    scvi.settings.seed = seed
    start_scvi = time()
    scvi.model.SCVI.setup_anndata(adata, layer="counts", categorical_covariate_keys=["group", "sample"])
    model = scvi.model.SCVI(adata, n_latent=n_latent)
    model.train(early_stopping=True)
    adata.obsm['X_scVI'] = model.get_latent_representation(adata).astype(np.float32)
    time_scvi = time() - start_scvi

    start_neigh = time()
    cc.tl.SpatialCluster.aggregate_neighbors(adata, nhood_layers, X_key='X_scVI')
    time_neigh = time() - start_neigh

    start_cls = time()
    cls = cc.tl.SpatialCluster(7, random_state=seed, gpus=args.gpu)
    cls.fit(adata, X_key='X_cellcharter')
    adata.obs[f'cluster_cellcharter_{i}'] = cls.predict(adata, X_key='X_cellcharter')
    time_cls = time() - start_cls

    times[f'{sample} (hvg)'].append(time_hvg)
    times[f'{sample} (scvi)'].append(time_scvi)
    times[f'{sample} (aggregation)'].append(time_neigh)
    times[f'{sample} (clustering)'].append(time_cls)

    times_df = pd.DataFrame.from_dict(times, orient='index')
    times_df.to_csv(f"/work/FAC/FBM/DBC/gciriell/spacegene/Packages/cellcharter_analyses/results/dlpfc/CellCharter/time/time_hvgs{hvgs}_nlatent{n_latent}_nhoodlayers{nhood_layers}_{'gpu' if args.gpu else 'cpu'}_combined_groupsamplebatch.csv")
    
    adata.write(f"/work/FAC/FBM/DBC/gciriell/spacegene/Packages/cellcharter_analyses/results/dlpfc/CellCharter/labels/ARI_hvgs{hvgs}_nlatent{n_latent}_nhoodlayers{nhood_layers}_{'gpu' if args.gpu else 'cpu'}_combined_groupsamplebatch.h5ad")

    for sample, n_clusters in SAMPLES.items():
        adata_sample = adata[adata.obs['sample'] == sample]

        ari = adjusted_rand_score(
            adata_sample[~adata_sample.obs['sce.layer_guess'].isna()].obs['sce.layer_guess'], 
            adata_sample[~adata_sample.obs['sce.layer_guess'].isna()].obs[f'cluster_cellcharter_{i}'])
        print('ARI:', ari)
        aris[f'{sample}'].append(ari)

        aris_df = pd.DataFrame.from_dict(aris, orient='index')
        aris_df.to_csv(f"/work/FAC/FBM/DBC/gciriell/spacegene/Packages/cellcharter_analyses/results/dlpfc/CellCharter/accuracy/ARI_hvgs{hvgs}_nlatent{n_latent}_nhoodlayers{nhood_layers}_{'gpu' if args.gpu else 'cpu'}_combined_groupsamplebatch.csv")

