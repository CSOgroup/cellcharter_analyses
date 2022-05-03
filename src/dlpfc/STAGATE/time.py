import anndata as ad
from collections import defaultdict
import numpy as np
import pandas as pd
import scanpy as sc
from time import time
import argparse
import os
import torch
import STAGATE_pyG
from torch_geometric.loader import DataLoader
from tqdm.auto import tqdm
import torch.nn.functional as F

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--gpu',  type=int, default=0)
parser.add_argument('--n_samples', type=int)
parser.add_argument('--n_cpus', type=int)
parser.add_argument('--run_index', type=int)

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


hvgs=5000
hidden_dim=1024
n_latent=30
lr=0.001
weight_decay=1e-4
num_epoch=1000

device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')


times = defaultdict(list)

times_path = f"/work/FAC/FBM/DBC/gciriell/spacegene/Packages/cellcharter_analyses/results/dlpfc/STAGATE/time/time__hvgs{hvgs}_hidden_dim{hidden_dim}_nlatent{n_latent}_{'gpu' if args.gpu else 'cpu'}_ncpus{args.n_cpus}.csv"

adata_list = []
for sample, n_clusters in SAMPLES.items():
    input_dir = os.path.join('/work/FAC/FBM/DBC/gciriell/spacegene/Data/jhpce_human_pilot_10x', sample)
    adata = ad.read_h5ad(os.path.join(input_dir, f'{sample}.h5ad'))
    adata.obs['sample'] = [sample]*adata.shape[0]
    adata.obs['layer_guess'] = pd.read_csv(f'{input_dir}/metadata.tsv', sep='\t').loc[adata.obs_names, 'layer_guess']
    adata_list.append(adata)


adata = ad.concat(adata_list[:args.n_samples], pairwise=True)    

start_preprocess = time()

sc.pp.filter_genes(adata, min_counts=3)
sc.pp.filter_cells(adata, min_counts=3)

adata.layers["counts"] = adata.X.copy()

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=hvgs,
    subset=True,
    layer="counts",
    flavor="seurat_v3",

)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

adata_list = [adata[adata.obs['sample'] == sample] for sample in adata.obs['sample'].unique()]

for a in adata_list:
    STAGATE_pyG.Cal_Spatial_Net(a, rad_cutoff=150)


adata = ad.concat(adata_list, pairwise=True)
adata.uns['Spatial_Net'] = pd.concat([a.uns['Spatial_Net'] for a in adata_list]) 

data_list = [STAGATE_pyG.Transfer_pytorch_Data(adata) for adata in adata_list]
for temp in data_list:
    temp.to(device)

data = STAGATE_pyG.Transfer_pytorch_Data(adata)

# batch_size=1 or 2
loader = DataLoader(data_list, batch_size=1, shuffle=True)

time_preprocess = time() - start_preprocess

rng = np.random.default_rng(12345)
seeds = rng.integers(low=0, high=32768, size=5)
i = seeds[args.run_index]
np.random.seed(i)
torch.manual_seed(i)
torch.cuda.manual_seed(i)

start_stagate = time()

model = STAGATE_pyG.STAGATE(hidden_dims = [data_list[0].x.shape[1]]+[hidden_dim, n_latent]).to(device)
optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)
for epoch in tqdm(range(1, num_epoch+1)):
    for batch in loader:
        model.train()
        optimizer.zero_grad()
        z, out = model(batch.x, batch.edge_index)
        loss = F.mse_loss(batch.x, out) #F.nll_loss(out[data.train_mask], data.y[data.train_mask])
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 5.)
        optimizer.step()

data.to(device)
model.eval()
z, out = model(data.x, data.edge_index)

STAGATE_rep = z.to('cpu').detach().numpy()
adata.obsm['STAGATE'] = STAGATE_rep

time_stagate = time() - start_stagate

start_cls = time()

adata = STAGATE_pyG.mclust_R(adata, used_obsm='STAGATE', num_cluster=n_clusters)

time_cls = time() - start_cls

row = pd.DataFrame([[args.n_samples, time_preprocess, time_stagate, time_cls]], columns=['n_samples', 'preprocess', 'STAGATE', 'clustering'])
if os.path.exists(times_path):
    times_df = pd.read_csv(times_path, index_col=0)
    times_df = pd.concat((times_df, row))
else:
    times_df = row

times_df.to_csv(times_path)

