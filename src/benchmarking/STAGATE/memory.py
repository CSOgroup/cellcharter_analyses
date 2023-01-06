
import os
import argparse
import anndata as ad
from collections import defaultdict
import scanpy as sc
import pandas as pd
import STAGATE_pyG
from sklearn.metrics import adjusted_rand_score
import numpy as np
import torch
from memory_profiler import memory_usage, profile
from torch_geometric.loader import DataLoader
from tqdm.auto import tqdm
import torch.nn.functional as F
import scvi

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--gpu',  type=int, default=0)
parser.add_argument('--n_samples', type=int)
parser.add_argument('--scvi', default=False, action='store_true')

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

hvgs=2000
hidden_dim=512
n_latent=30
lr=0.001
weight_decay=1e-4
num_epoch=1000

adata_list = []
for sample, n_clusters in SAMPLES.items():
    input_dir = '../../../data/Visium_DLPFC/preprocessed_h5ad/'
    adata = ad.read_h5ad(os.path.join(input_dir, f'{sample}.h5ad'))
    adata.obs['sample'] = [sample]*adata.shape[0]
    adata.obs['layer_guess'] = pd.read_csv(f'{input_dir}/metadata.tsv', sep='\t').loc[adata.obs_names, 'layer_guess']
    STAGATE_pyG.Cal_Spatial_Net(adata, rad_cutoff=150)
    adata_list.append(adata)



def run(adata_list, n_samples):
    adata = ad.concat(adata_list, pairwise=True)  

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

    if args.scvi:
        scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="group")
        model = scvi.model.SCVI(adata, n_latent=n_latent)
        model.train(early_stopping=True)

        adata= ad.AnnData(
            X=model.get_latent_representation(adata).astype(np.float32), 
            obs=adata.obs, 
            obsm=adata.obsm, 
            obsp=adata.obsp, 
            uns=adata.uns)

    adata_list = [adata[adata.obs['sample'] == sample] for sample in SAMPLES.keys()]

    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

    # Consturcting network for each batch
    for a in adata_list:
        STAGATE_pyG.Cal_Spatial_Net(a, rad_cutoff=150)
        #STAGATE_pyG.Stats_Spatial_Net(temp_adata)
    #STAGATE_pyG.Cal_Spatial_Net(adata, rad_cutoff=150)
    adata = ad.concat(adata_list, pairwise=True)
    adata.uns['Spatial_Net'] = pd.concat([a.uns['Spatial_Net'] for a in adata_list]) 

    data_list = [STAGATE_pyG.Transfer_pytorch_Data(adata) for adata in adata_list]
    for temp in data_list:
        temp.to(device)



    #STAGATE_pyG.Cal_Spatial_Net(adata, rad_cutoff=150)
    data = STAGATE_pyG.Transfer_pytorch_Data(adata)

    # batch_size=1 or 2
    loader = DataLoader(data_list, batch_size=1, shuffle=True)

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

    adata = STAGATE_pyG.mclust_R(adata, used_obsm='STAGATE', num_cluster=n_clusters)


mem = max(memory_usage(proc=(run, [adata_list, args.n_samples])))
mems_path = f"../../../results/benchmarking/memory/memory_STAGATE_hvgs{hvgs}_hidden_dim{hidden_dim}_nlatent{n_latent}_{'gpu' if args.gpu else 'cpu'}.csv"
row = pd.DataFrame([mem], index=[args.n_samples], )
if os.path.exists(mems_path):
    mems_df = pd.read_csv(mems_path, index_col=0)
    mems_df.columns = mems_df.columns.astype(int)
    mems_df = pd.concat((mems_df, row))
else:
    mems_df = row

mems_df.to_csv(mems_path)


