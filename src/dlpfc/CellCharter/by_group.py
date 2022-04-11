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
from memory_profiler import memory_usage
import gc


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

N_CLUSTERS = [7,5,7]

hvgs = 5000
n_latent = 5
nhood_layers = 4

adata_list = []
for sample, n_clusters in SAMPLES.items():
    adata = ad.read_h5ad(f'/work/FAC/FBM/DBC/gciriell/spacegene/Data/jhpce_human_pilot_10x/{sample}/{sample}.h5ad')
    adata.obs['sample'] = [sample]*adata.shape[0]
    adata_list.append(adata)

rng = np.random.default_rng(12345)
seeds = rng.integers(low=0, high=32768, size=10)

aris = defaultdict(list)
for j in range(3):
    adata = ad.concat(adata_list[j*3:j*3+3], pairwise=True)    

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
    for i in seeds:
        scvi.settings.seed = i 
        scvi.model.SCVI.setup_anndata(adata, layer="counts")
        model = scvi.model.SCVI(adata, n_latent=n_latent)
        model.train(early_stopping=True)
        adata.obsm['X_scVI'] = model.get_latent_representation(adata).astype(np.float32)

        cc.tl.SpatialCluster.aggregate_neighbors(adata, nhood_layers, X_key='X_scVI')

        cls = cc.tl.SpatialCluster(N_CLUSTERS[j], random_state=42, gpus=args.gpu)
        cls.fit(adata, X_key='X_cellcharter')
        adata.obs[f'cluster_cellcharter_{i}'] = cls.predict(adata, X_key='X_cellcharter')

        adata.write(f"/work/FAC/FBM/DBC/gciriell/spacegene/Packages/cellcharter_analyses/results/dlpfc/CellCharter/labels/ARI_hvgs{hvgs}_nlatent{n_latent}_nhoodlayers{nhood_layers}_{'gpu' if args.gpu else 'cpu'}_group{j}_nobatch.h5ad")


        for sample in list(SAMPLES.keys())[j*3:j*3+3]:
            
            ari = adjusted_rand_score(
                adata[(~adata.obs['sce.layer_guess'].isna()) & (adata.obs['sample'] == sample)].obs['sce.layer_guess'],
                adata[(~adata.obs['sce.layer_guess'].isna()) & (adata.obs['sample'] == sample)].obs[f'cluster_cellcharter_{i}'])
            aris[sample].append(ari)


        aris_df = pd.DataFrame.from_dict(aris, orient='index')
        aris_df.to_csv(f"/work/FAC/FBM/DBC/gciriell/spacegene/Packages/cellcharter_analyses/results/dlpfc/CellCharter/accuracy/ARI_hvgs{hvgs}_nlatent{n_latent}_nhoodlayers{nhood_layers}_{'gpu' if args.gpu else 'cpu'}_group_no-batch.csv")

    
