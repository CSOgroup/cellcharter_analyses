
import scanpy as sc
import pandas as pd
from pathlib import Path
from scanpy.readwrite import read_visium
from scanpy._utils  import check_presence_download
import networkx as nx
import scipy.sparse as sp
from scipy.spatial import distance
from sklearn import metrics
from collections import defaultdict

import numpy as np
import torch
import os
import anndata as ad
from SEDR.src.SEDR_train import SEDR_Train
from sedr_utils import load_ST_file, adata_preprocess, graph_construction, res_search_fixed_clus, Params

k = 10
knn_distanceType = 'euclidean'
epochs = 200
using_dec=True
using_mask=False
feat_w = 10
gcn_w = 0.1
dec_kl_w = 10
gcn_lr = 0.01
gcn_decay = 0.01
dec_cluster_n = 10
dec_interval = 20
dec_tol = 0.0
p_drop=0.1
gcn_hidden1=32
feat_hidden1=100

eval_resolution = 1
eval_graph_n = 20

pcs = [100, 200, 300]
gcn_hiddens2 = [4, 8, 16]
feat_hiddens2 = [10, 20, 30]


SAMPLES = {
    '151507': 7, 
    '151672': 5, 
    '151673': 7
}

device = 'cuda:0' if torch.cuda.is_available() else 'cpu'
print('===== Using device: ' + device)

aris = defaultdict(list)
times = defaultdict(list)

for sample, n_clusters in SAMPLES.items():
    input_dir = os.path.join('../../../data/Visium_DLPFC/raw', sample)
    adata = load_ST_file(file_fold=input_dir, load_images=False)
    adata.var_names_make_unique()

    for feat_hidden2 in feat_hiddens2:
        for gcn_hidden2 in gcn_hiddens2:
            for pc in pcs:
                params = Params(dict(
                    k=k,
                    knn_distanceType=knn_distanceType,
                    epochs=epochs,
                    using_dec=using_dec,
                    using_mask=using_mask,
                    feat_w=feat_w,
                    gcn_w=gcn_w,
                    dec_kl_w=dec_kl_w,
                    gcn_lr=gcn_lr,
                    gcn_decay=gcn_decay,
                    dec_cluster_n=dec_cluster_n,
                    dec_interval=dec_interval,
                    dec_tol=dec_tol,
                    eval_resolution=eval_resolution,
                    eval_graph_n=eval_graph_n,
                    cell_feat_dim=pc,
                    gcn_hidden2=gcn_hidden2,
                    gcn_hidden1=gcn_hidden1,
                    feat_hidden2=feat_hidden2,
                    feat_hidden1=feat_hidden1,
                    p_drop=p_drop
                ))

                

                for i in range(5):
                    np.random.seed(i)
                    torch.manual_seed(i)
                    torch.cuda.manual_seed(i)

                    adata_X = adata_preprocess(adata, min_cells=5, pca_n_comps=params.cell_feat_dim)
                    graph_dict = graph_construction(adata.obsm['spatial'], adata.shape[0], params)
                    params.cell_num = adata.shape[0]
                    params.device = device

                    sedr_net = SEDR_Train(adata_X, graph_dict, params)
                    if params.using_dec:
                        sedr_net.train_with_dec()
                    else:
                        sedr_net.train_without_dec()
                    sedr_feat, _, _, _ = sedr_net.process()

                    adata_sedr = ad.AnnData(sedr_feat)
                    adata_sedr.uns['spatial'] = adata.uns['spatial']
                    adata_sedr.obsm['spatial'] = adata.obsm['spatial']

                    sc.pp.neighbors(adata_sedr, n_neighbors=params.eval_graph_n)

                    eval_resolution = res_search_fixed_clus(adata_sedr, n_clusters)

                    sc.tl.leiden(adata_sedr, key_added="SEDR_leiden", resolution=eval_resolution)

                    df_meta = pd.read_csv(f'{input_dir}/metadata.tsv', sep='\t')
                    df_meta['SEDR'] = adata_sedr.obs['SEDR_leiden'].tolist()
                    #df_meta = df_meta[~pd.isnull(df_meta['layer_guess'])]
                    ari = metrics.adjusted_rand_score(
                        df_meta[~pd.isnull(df_meta['layer_guess'])]['layer_guess'], 
                        df_meta[~pd.isnull(df_meta['layer_guess'])]['SEDR']
                    )
                    print(ari)

                    aris[f'{sample}_{p_drop}_{feat_hidden1}_{feat_hidden2}_{gcn_hidden1}_{gcn_hidden2}_{pc}'].append(ari)

                    aris_df = pd.DataFrame.from_dict(aris, orient='index')
                    aris_df.to_csv(f"../../../results/benchmarking/hyperparm_tuning/ARI_SEDR_tuning_{sample}.csv")
        
