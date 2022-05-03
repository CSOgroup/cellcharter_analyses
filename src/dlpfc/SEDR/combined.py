import anndata as ad
from collections import defaultdict
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import adjusted_rand_score
import argparse
import os
import torch

from SEDR.src.SEDR_train import SEDR_Train
from sedr_utils import load_ST_file, params_sedr, adata_preprocess, graph_construction, res_search_fixed_clus, Params

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


params_sedr['feat_hidden2'] = 10
params_sedr['gcn_hidden2'] = 16
params_sedr['cell_feat_dim'] = 100 #pc

params = Params(params_sedr)
device = 'cuda:0' if args.gpu else 'cpu'

aris = defaultdict(list)

adata_list = []
for sample, n_clusters in SAMPLES.items():
    input_dir = os.path.join('/work/FAC/FBM/DBC/gciriell/spacegene/Data/jhpce_human_pilot_10x', sample)
    adata = load_ST_file(file_fold=input_dir, load_images=False)
    adata.obs['sample'] = [sample]*adata.shape[0]
    adata.obs['layer_guess'] = pd.read_csv(f'{input_dir}/metadata.tsv', sep='\t').loc[adata.obs_names, 'layer_guess']
    adata_list.append(adata)

adata = ad.concat(adata_list, pairwise=True)
adata.var_names_make_unique()
adata_X = adata_preprocess(adata, min_cells=5, pca_n_comps=params.cell_feat_dim)
graph_dict = graph_construction(adata.obsm['spatial'], adata.shape[0], params)
params.cell_num = adata.shape[0]
params.device = device

rng = np.random.default_rng(12345)
seeds = rng.integers(low=0, high=32768, size=10)
for i, seed in enumerate(seeds):
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)

    sedr_net = SEDR_Train(adata_X, graph_dict, params)
    if params.using_dec:
        sedr_net.train_with_dec()
    else:
        sedr_net.train_without_dec()
    sedr_feat, _, _, _ = sedr_net.process()

    adata_sedr = ad.AnnData(sedr_feat)
    # adata_sedr.uns['spatial'] = adata.uns['spatial']
    adata_sedr.obsm = adata.obsm
    adata_sedr.obs = adata.obs

    sc.pp.neighbors(adata_sedr, n_neighbors=params.eval_graph_n)

    eval_resolution = res_search_fixed_clus(adata_sedr, n_clusters)

    sc.tl.leiden(adata_sedr, key_added=f"SEDR_leiden_{i}", resolution=eval_resolution)

    adata.write(f"/work/FAC/FBM/DBC/gciriell/spacegene/Packages/cellcharter_analyses/results/dlpfc/SEDR/labels/ARI_feat_hidden2_{params.feat_hidden2}_gcn_hidden2_{params.gcn_hidden2}_pc{params.cell_feat_dim}_{'gpu' if args.gpu else 'cpu'}_combined.h5ad")

    for sample, n_clusters in SAMPLES.items():
        adata_sample = adata_sedr[adata_sedr.obs['sample'] == sample]

        obs_df = adata_sample.obs.dropna()
        ari = adjusted_rand_score(obs_df[f'SEDR_leiden_{i}'], obs_df['layer_guess'])
        aris[f'{sample}'].append(ari)

        aris_df = pd.DataFrame.from_dict(aris, orient='index')
        aris_df.to_csv(f"/work/FAC/FBM/DBC/gciriell/spacegene/Packages/cellcharter_analyses/results/dlpfc/SEDR/accuracy/ARI_feat_hidden2_{params.feat_hidden2}_gcn_hidden2_{params.gcn_hidden2}_pc{params.cell_feat_dim}_{'gpu' if args.gpu else 'cpu'}_combined.csv")