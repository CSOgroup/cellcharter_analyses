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

hvgs=5000
hidden_dim=1024
n_latent=30
lr=0.001
weight_decay=1e-4
num_epoch=1000


aris = defaultdict(list)
adata_list = []
for sample, n_clusters in SAMPLES.items():
    input_dir = '../../../data/Visium_DLPFC/preprocessed_h5ad/'
    adata = ad.read_h5ad(os.path.join(input_dir, f'{sample}.h5ad'))
    adata.obs['sample'] = [sample]*adata.shape[0]
    adata.obs['group'] = [GROUPS[sample]]*adata.shape[0]
    adata.obs['layer_guess'] = pd.read_csv(f'{input_dir}/metadata.tsv', sep='\t').loc[adata.obs_names, 'layer_guess']
    adata_list.append(adata)


adata = ad.concat(adata_list, pairwise=True)  

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


rng = np.random.default_rng(12345)
seeds = rng.integers(low=0, high=32768, size=10)
for i, seed in enumerate(seeds):
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="group")
    model = scvi.model.SCVI(adata, n_latent=n_latent)
    model.train(early_stopping=True)

    adata_scvi = ad.AnnData(
        X=model.get_latent_representation(adata).astype(np.float32), 
        obs=adata.obs, 
        obsm=adata.obsm, 
        obsp=adata.obsp, 
        uns=adata.uns)

    adata_list = [adata_scvi[adata_scvi.obs['sample'] == sample] for sample in SAMPLES.keys()]

    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

    # Consturcting network for each batch
    for a in adata_list:
        STAGATE_pyG.Cal_Spatial_Net(a, rad_cutoff=150)
        #STAGATE_pyG.Stats_Spatial_Net(temp_adata)
    #STAGATE_pyG.Cal_Spatial_Net(adata, rad_cutoff=150)
    adata_scvi = ad.concat(adata_list, pairwise=True)
    adata_scvi.uns['Spatial_Net'] = pd.concat([a.uns['Spatial_Net'] for a in adata_list]) 

    data_list = [STAGATE_pyG.Transfer_pytorch_Data(a) for a in adata_list]
    for temp in data_list:
        temp.to(device)



    #STAGATE_pyG.Cal_Spatial_Net(adata, rad_cutoff=150)
    data = STAGATE_pyG.Transfer_pytorch_Data(adata_scvi)



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
    adata_scvi.obsm[f'STAGATE_{i}'] = STAGATE_rep

    adata.obsm[f'STAGATE_{i}'] = STAGATE_rep
    adata.obsm[f'scVI_{i}'] = adata_scvi.X
    adata_stagate = STAGATE_pyG.mclust_R(adata_scvi, used_obsm=f'STAGATE_{i}', num_cluster=7)
    adata.obs[f'mclust_{i}'] = adata_stagate.obs['mclust']

    adata.write(f"../../../results/benchmarking/joint/labels_STAGATE_hvgs{hvgs}_hidden_dim{hidden_dim}_nlatent{n_latent}_{'gpu' if args.gpu else 'cpu'}_joint.h5ad")

    for sample, n_clusters in SAMPLES.items():
        adata_sample = adata[adata.obs['sample'] == sample]

        obs_df = adata_sample.obs.dropna()
        ari = adjusted_rand_score(obs_df[f'mclust_{i}'], obs_df['layer_guess'])
        aris[f'{sample}'].append(ari)

        aris_df = pd.DataFrame.from_dict(aris, orient='index')
        aris_df.to_csv(f"../../../results/benchmarking/joint/ARI_STAGATE_hvgs{hvgs}_hidden_dim{hidden_dim}_nlatent{n_latent}_{'gpu' if args.gpu else 'cpu'}_joint.csv")

