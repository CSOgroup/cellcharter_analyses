import anndata as ad
from collections import defaultdict
import numpy as np
import pandas as pd
import scanpy as sc
from time import time
import argparse
import os
import torch
from SEDR.src.SEDR_train import SEDR_Train
from sedr_utils import load_ST_file, params_sedr, adata_preprocess, graph_construction, res_search_fixed_clus, Params
import scvi

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--gpu',  type=int, default=0)
parser.add_argument('--n_samples', type=int)
parser.add_argument('--n_cpus', type=int)
parser.add_argument('--run_index', type=int)
parser.add_argument('--scvi', default=True, action='store_false')

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

params_sedr['feat_hidden2'] = 10
params_sedr['gcn_hidden2'] = 16
params_sedr['cell_feat_dim'] = 100 if args.scvi is None else 10

scvi_hvgs=5000

params = Params(params_sedr)
device = 'cuda:0' if args.gpu else 'cpu'

times = defaultdict(list)

times_path = f"../../../results/benchmarking/time/time_SEDR_feat_hidden2_{params.feat_hidden2}_gcn_hidden2_{params.gcn_hidden2}_pc{params.cell_feat_dim}_{'gpu' if args.gpu else 'cpu'}_ncpus{args.n_cpus}.csv"

adata_list = []
for sample, n_clusters in SAMPLES.items():
    input_dir = os.path.join('/work/FAC/FBM/DBC/gciriell/spacegene/Data/jhpce_human_pilot_10x', sample)
    adata = load_ST_file(file_fold=input_dir, load_images=False)
    adata.obs['sample'] = [sample]*adata.shape[0]
    adata.obs['group'] = [GROUPS[sample]]*adata.shape[0]
    adata.obs['layer_guess'] = pd.read_csv(f'{input_dir}/metadata.tsv', sep='\t').loc[adata.obs_names, 'layer_guess']
    adata_list.append(adata)


adata = ad.concat(adata_list[:args.n_samples], pairwise=True)    
adata.var_names_make_unique()

start_preprocess = time()
if args.scvi:
    sc.pp.filter_genes(adata, min_counts=3)
    sc.pp.filter_cells(adata, min_counts=3)

    adata.layers["counts"] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=scvi_hvgs,
        subset=True,
        layer="counts",
        flavor="seurat_v3",
    )
else:
    adata_X = adata_preprocess(adata, min_cells=5, pca_n_comps=params.cell_feat_dim)

graph_dict = graph_construction(adata.obsm['spatial'], adata.shape[0], params)
params.cell_num = adata.shape[0]
params.device = device

time_preprocess = time() - start_preprocess

rng = np.random.default_rng(12345)
seeds = rng.integers(low=0, high=32768, size=5)
i = seeds[args.run_index]
np.random.seed(i)
torch.manual_seed(i)
torch.cuda.manual_seed(i)

start_scvi = time()
if args.scvi:
    scvi.settings.seed = i
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="group")
    model = scvi.model.SCVI(adata)
    model.train(early_stopping=True)
    adata_X = model.get_latent_representation(adata).astype(np.float32)
    
time_scvi = time() - start_scvi

start_sedr = time()

sedr_net = SEDR_Train(adata_X, graph_dict, params)
if params.using_dec:
    sedr_net.train_with_dec()
else:
    sedr_net.train_without_dec()
sedr_feat, _, _, _ = sedr_net.process()

time_sedr = time() - start_sedr

start_cls = time()
adata_sedr = ad.AnnData(sedr_feat)
# adata_sedr.uns['spatial'] = adata.uns['spatial']
adata_sedr.obsm = adata.obsm
adata_sedr.obs = adata.obs

sc.pp.neighbors(adata_sedr, n_neighbors=params.eval_graph_n)

eval_resolution = res_search_fixed_clus(adata_sedr, n_clusters)

sc.tl.leiden(adata_sedr, key_added="SEDR_leiden", resolution=eval_resolution)
time_cls = time() - start_cls

row = pd.DataFrame([[args.n_samples, time_preprocess, time_scvi, time_sedr, time_cls]], columns=['n_samples', 'preprocess', 'scVI', 'SEDR', 'clustering'])
if os.path.exists(times_path):
    times_df = pd.read_csv(times_path, index_col=0)
    times_df = pd.concat((times_df, row))
else:
    times_df = row

times_df.to_csv(times_path)