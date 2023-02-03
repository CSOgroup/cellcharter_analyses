
import scanpy as sc
import numpy as np
from sotip import MED, get_ground_distance, cal_ME_heterogeneity, MED_phEMD_mp, merge_cls_paga, adjusted_rand_score
import argparse
import anndata as ad
from collections import defaultdict
import pandas as pd
import torch

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
parser.add_argument('--sample',  type=str, default=None)

args = parser.parse_args()

SAMPLES = {
    '151507': 7, 
    '151672': 5, 
    '151673': 7, 
}

hvgs = [1000, 2000, 5000]
pcs_list = [50,100]
knn_list = [15,30,50]
ME_knn_list = [10,30,50]
EMD_neighbors_list = np.arange(100,600,200)

aris = defaultdict(list)

if args.sample is not None:
    SAMPLES = {args.sample: SAMPLES[args.sample]}

for sample, n_clusters in SAMPLES.items():
    rng = np.random.default_rng(12345)
    seeds = rng.integers(low=0, high=32768, size=5)
    for i in seeds[2:]:
        np.random.seed(i)
        torch.manual_seed(i)
        for hvg in hvgs:
            for pcs in pcs_list:
                for knn in knn_list:
                    for ME_knn in ME_knn_list:
                        for EMD_neighbors in EMD_neighbors_list:
                            adata = ad.read_h5ad(f'../../../data/Visium_DLPFC/preprocessed_h5ad/{sample}.h5ad')

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

                            sc.pp.scale(adata)
                            sc.pp.pca(adata,n_comps=pcs)
                            sc.pp.neighbors(adata,n_pcs=pcs,n_neighbors=knn)
                            sc.tl.umap(adata, random_state=i)
                            sc.tl.leiden(adata,resolution=1, random_state=i)
                            emd_distmat = get_emd_distmat(adata, ME_knn)
                            adata.obsp['ME_EMD_mat'] = emd_distmat

                            pd_dict = {
                                'PCA':[],
                                'Cell_cls_knn':[],
                                'Cell_cls_res':[],
                                'ME_knn':[],
                                'EMD_neighbors':[],
                                'leiden_ari':[],
                                'sotip_ari':[]
                            }

                            sotip_ari = get_sotip_ari(adata, EMD_neighbors, LIBD_cls_num=n_clusters)
                            aris[f'{sample}_{hvg}_{pcs}_{knn}_{ME_knn}_{EMD_neighbors}'].append(sotip_ari)

                            aris_df = pd.DataFrame.from_dict(aris, orient='index')
                            aris_df.to_csv(f"../../../results/benchmarking/hyperparm_tuning/ARI_SOTIP_tuning_{sample}.csv")
