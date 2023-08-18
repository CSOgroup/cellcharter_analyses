
import os
import anndata as ad
from collections import defaultdict
import pandas as pd
import STAGATE_pyG
import torch
from torch_geometric.loader import DataLoader
from tqdm.auto import tqdm
import torch.nn.functional as F
from time import time

hidden_dim=512
n_latent=10
lr=0.001
weight_decay=1e-4
num_epoch=1000
PROCESSED_PATH = '/scratch/mvarrone/goltsev_2018/Processed'
FILE_NAME = 'BALBc-MRL_scarches_noreducelr_patience5_hidden128'

adata = ad.read_h5ad(os.path.join(PROCESSED_PATH, f'{FILE_NAME}.h5ad'))
adata= ad.AnnData(
            X=adata.obsm[f'X_trVAE'], 
            obs=adata.obs, 
            obsm=adata.obsm, 
            obsp=adata.obsp, 
            uns=adata.uns)
adata.obs['X'] = adata.obsm['spatial'][:,0]
adata.obs['Y'] = adata.obsm['spatial'][:,1]

times = dict()

start_processing = time()
adata_list = [adata[adata.obs['sample'] == sample] for sample in adata.obs['sample'].unique()]

Batch_list = list()

for adata_sample in adata_list:
    Batch_list.extend(STAGATE_pyG.Batch_Data(adata_sample, num_batch_x=4, num_batch_y=6,
                                spatial_key=['X', 'Y'], plot_Stats=True))

for adata_batch in Batch_list:
    STAGATE_pyG.Cal_Spatial_Net(adata_batch, rad_cutoff=9000)

for adata_sample in adata_list:
    STAGATE_pyG.Cal_Spatial_Net(adata_sample, rad_cutoff=9000)

times['processing'] = [time() - start_processing]

start_spatial_aggregation = time()
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

data_list = [STAGATE_pyG.Transfer_pytorch_Data(a) for a in Batch_list]
for temp in data_list:
    temp.to(device)


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

for adata_sample in adata_list:
    data = STAGATE_pyG.Transfer_pytorch_Data(adata_sample)
    data.to(device)
    model.eval()
    z, out = model(data.x, data.edge_index)

    STAGATE_rep = z.to('cpu').detach().numpy()
    adata_sample.obsm['STAGATE'] = STAGATE_rep

adata = ad.concat(adata_list, pairwise=True)

times['spatial_aggregation'] = [time() - start_spatial_aggregation]


start_clustering = time()
adata = STAGATE_pyG.mclust_R(adata, used_obsm='STAGATE', num_cluster=11)
times['clustering'] = [time() - start_clustering]

times_df = pd.DataFrame.from_dict(times, orient='columns')
times_df.to_csv(f"../../results/codex_mouse_spleen/time_codex/time_STAGATE_trVAE_{'gpu' if torch.cuda.is_available() else 'cpu'}.csv")

adata.write(f"../../results/codex_mouse_spleen/BALBc-MRL_scarches_STAGATE_trVAE.h5ad")


