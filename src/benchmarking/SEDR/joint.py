import anndata as ad
from collections import defaultdict
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import adjusted_rand_score
import argparse
import os
import torch
import scvi

from SEDR.src.SEDR_train import SEDR_Train
from sedr_utils import load_ST_file, params_sedr, adata_preprocess, graph_construction, res_search_fixed_clus, Params

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--gpu',  type=int, default=0)
parser.add_argument('--scvi', default=True, action='store_false')

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


params_sedr['feat_hidden2'] = 10
params_sedr['gcn_hidden2'] = 16

# Number of PCs
params_sedr['cell_feat_dim'] = 100 if args.scvi is None else 10

hvgs=5000


params = Params(params_sedr)
device = 'cuda:0' if args.gpu else 'cpu'

aris = defaultdict(list)

adata_list = []
for sample, n_clusters in SAMPLES.items():
    input_dir = os.path.join('../../../data/Visium_DLPFC/raw', sample)
    adata = load_ST_file(file_fold=input_dir, load_images=False)
    adata.obs['sample'] = [sample]*adata.shape[0]
    adata.obs['group'] = [GROUPS[sample]]*adata.shape[0]
    adata.obs['layer_guess'] = pd.read_csv(f'{input_dir}/metadata.tsv', sep='\t').loc[adata.obs_names, 'layer_guess']
    adata_list.append(adata)

adata = ad.concat(adata_list, pairwise=True)
adata.var_names_make_unique()

if args.scvi:
    sc.pp.filter_genes(adata, min_counts=3)
    sc.pp.filter_cells(adata, min_counts=3)

    adata.layers["counts"] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=hvgs,
        subset=True,
        layer="counts",
        flavor="seurat_v3",

    )

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
else:
    adata_X = adata_preprocess(adata, min_cells=5, pca_n_comps=params.cell_feat_dim)
graph_dict = graph_construction(adata.obsm['spatial'], adata.shape[0], params)
params.cell_num = adata.shape[0]
params.device = 'cpu'

rng = np.random.default_rng(12345)
seeds = rng.integers(low=0, high=32768, size=10)
for i, seed in enumerate(seeds):
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)

    if args.scvi:
        scvi.settings.seed = seed
        scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="group")
        model = scvi.model.SCVI(adata, n_latent=params_sedr['cell_feat_dim'])
        model.train(early_stopping=True)
        adata_X = model.get_latent_representation(adata).astype(np.float32)

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

    adata.write(f"../../../results/benchmarking/joint/labels_SEDR_feat_hidden2_{params.feat_hidden2}_gcn_hidden2_{params.gcn_hidden2}_pc{params.cell_feat_dim}_{'gpu' if args.gpu else 'cpu'}_joint.h5ad")

    for sample, n_clusters in SAMPLES.items():
        adata_sample = adata_sedr[adata_sedr.obs['sample'] == sample]

        obs_df = adata_sample.obs.dropna()
        ari = adjusted_rand_score(obs_df[f'SEDR_leiden_{i}'], obs_df['layer_guess'])
        aris[f'{sample}'].append(ari)

        aris_df = pd.DataFrame.from_dict(aris, orient='index')
        aris_df.to_csv(f"../../../results/benchmarking/joint/ARI_SEDR_feat_hidden2_{params.feat_hidden2}_gcn_hidden2_{params.gcn_hidden2}_pc{params.cell_feat_dim}_{'gpu' if args.gpu else 'cpu'}_joint.csv")