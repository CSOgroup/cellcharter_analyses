import anndata as ad
from collections import defaultdict
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
from time import time
import argparse
import os
from sotip import MED, get_ground_distance, cal_ME_heterogeneity, MED_phEMD_mp, merge_cls_paga, adjusted_rand_score


def get_emd_distmat(adata,ME_knn):
    knn = ME_knn
    spatial_var='spatial'
    cls_key='leiden'
    ME_var_names_np_unique = np.array(adata.obs[cls_key].cat.categories) 

    MED(adata,use_cls=cls_key,nn=knn,copy=False,ME_var_names_np_unique=ME_var_names_np_unique,spatial_var=spatial_var) 
    sc.tl.paga(adata,groups=cls_key)
    sc.pl.paga_compare(adata,basis='X_umap',show=True)    
    gd_method = 'paga_guided_umap'
    gd = get_ground_distance(adata,method=gd_method,cls_key=cls_key,embed_key=None,connect_threshold=0.5)  
    heter_key = 'ME_heter_{0}_{1}'.format(cls_key,gd_method)
    cal_ME_heterogeneity(adata,copy=False,key_added=heter_key) 
    adata_phEMD = MED_phEMD_mp(
        adata.copy(),
        GD_method=gd_method,
        MED_knn=knn,
        CT_obs=cls_key,
        ifspatialplot=False,
        OT_method='pyemd',
        ME_precompyted=True,
        GD_precomputed=True,
        mp=200
    )
    return adata_phEMD.obsm['X_ME_EMD_mat']

def get_sotip_ari(adata,n_neighbors,LIBD_cls_num=7):
    knn_indices, knn_dists, forest = sc.neighbors.compute_neighbors_umap( adata.obsp['ME_EMD_mat'], n_neighbors=n_neighbors, metric='precomputed' )
    adata.obsp['distances'], adata.obsp['connectivities'] = sc.neighbors._compute_connectivities_umap(
        knn_indices,
        knn_dists,
        adata.shape[0],
        n_neighbors, # change to neighbors you plan to use
    )
    adata.uns['neighbors_EMD'] = adata.uns['neighbors'].copy()

    sc.tl.leiden(adata,neighbors_key='neighbors_EMD',key_added='leiden_EMD',resolution=1)
    sc.tl.paga(adata,groups='leiden_EMD',neighbors_key='neighbors_EMD')
    merge_cls_paga(adata,thresh=0,min_cls=LIBD_cls_num,paga_plot=False)

    adata_valid = adata[np.logical_not(adata.obs['sce.layer_guess'].isna())]
    cur_ari = adjusted_rand_score(adata_valid.obs['sce.layer_guess'],adata_valid.obs['leiden_EMD_merge'])
    return cur_ari

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--gpu',  type=int, default=0)
parser.add_argument('--n_samples', type=int)
parser.add_argument('--run-index', type=int, default=None)
parser.add_argument('--n_cpus', type=int)

args = parser.parse_args()

SAMPLES = {'151507': 7, 
           '151508': 7, 
           '151509': 7, 
           '151510': 7, 
           '151673': 7, 
           '151674': 7, 
           '151675': 7, 
           '151676': 7,
           '151669': 5, 
           '151670': 5, 
           '151671': 5, 
           '151672': 5, 
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

hvgs = 2000
n_latent = 10
knn=30
ME_knn=50
EMD_neighbors=500
times = defaultdict(list)

times_path = f"../../../results/benchmarking/time/time_SOTIP_hvgs{hvgs}_nlatent{n_latent}_knn{knn}_MEknn{ME_knn}_EMDneighbors{EMD_neighbors}_ncpus{args.n_cpus}.csv"

adata_list = []
for sample, n_clusters in SAMPLES.items():
    adata = ad.read_h5ad(f'../../../data/Visium_DLPFC/preprocessed_h5ad/{sample}.h5ad')
    adata.obs['sample'] = [sample]*adata.shape[0]
    adata.obs['group'] = [GROUPS[sample]]*adata.shape[0]
    adata_list.append(adata)

rng = np.random.default_rng(12345)
seeds = rng.integers(low=0, high=32768, size=5)

if args.run_index is not None:
    seeds = [seeds[args.run_index]]

print(seeds)

for i in seeds:
    print('Setting seed', i)
    scvi.settings.seed = i 

    adata = ad.concat(adata_list[:args.n_samples], pairwise=True)
    adata.obs_names_make_unique()
    print(adata)
    print(adata.X)

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
    start_scvi = time()
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key='sample')
    model = scvi.model.SCVI(adata, n_latent=n_latent)
    model.train(early_stopping=True)
    adata.obsm['X_scVI'] = model.get_latent_representation(adata).astype(np.float32)
    time_scvi = time() - start_scvi

    start_cell_typing = time()
    sc.pp.neighbors(adata, use_rep='X_scVI', n_neighbors=knn)
    sc.tl.umap(adata, random_state=i)
    sc.tl.leiden(adata,resolution=1, random_state=i)
    time_cell_typing = time() - start_cell_typing

    start_emd_distance = time()
    emd_distmat = get_emd_distmat(adata, ME_knn)
    adata.obsp['ME_EMD_mat'] = emd_distmat
    time_emd_distance = time() - start_emd_distance

    start_clustering = time()
    sotip_ari = get_sotip_ari(adata, EMD_neighbors, LIBD_cls_num=7)
    time_clustering = time() - start_clustering

    row = pd.DataFrame([[args.n_samples, time_hvg, time_scvi, time_cell_typing, time_emd_distance, time_clustering]], columns=['n_samples', 'HVGs', 'scVI', 'cell_typing', 'emd_distance', 'clustering'])
    if os.path.exists(times_path):
        times_df = pd.read_csv(times_path, index_col=0)
        times_df = pd.concat((times_df, row))
    else:
        times_df = row

    times_df.to_csv(times_path)
