# CellCharter analyses

Code to reproduce the analyses of [CellCharter](https://doi.org/10.1038/s41588-023-01588-4).

## Figures
| Notebook / Script              | Figures         |
|--------------------------------|-----------------|
|  Visium Human DLPFC (benchmarking) ([GitHub](https://github.com/CSOgroup/cellcharter_analyses/blob/main/visium_human_dlpfc.ipynb), [nbviewer](https://nbviewer.org/github/CSOgroup/cellcharter_analyses/blob/main/visium_human_dlpfc.ipynb)) | Fig. 1 |
|  ATAC+RNA Mouse Brain ([GitHub](https://github.com/CSOgroup/cellcharter_analyses/blob/main/ATAC-RNA_mouse_brain.ipynb), [nbviewer](https://nbviewer.org/github/CSOgroup/cellcharter_analyses/blob/main/ATAC-RNA_mouse_brain.ipynb)) | Fig. 2 |
|  CODEX Mouse Spleen ([GitHub](https://github.com/CSOgroup/cellcharter_analyses/blob/main/codex_mouse_spleen.ipynb), [nbviewer](https://nbviewer.org/github/CSOgroup/cellcharter_analyses/blob/main/codex_mouse_spleen.ipynb)) | Figs. 2-3 |
|  CosMx Human NSCLC ([GitHub](https://github.com/CSOgroup/cellcharter_analyses/blob/main/cosmx_human_nsclc.ipynb), [nbviewer](https://nbviewer.org/github/CSOgroup/cellcharter_analyses/blob/main/cosmx_human_nsclc.ipynb)) | Figs. 4-5 |
|  MERSCOPE Human Lung Cancer ([GitHub](https://github.com/CSOgroup/cellcharter_analyses/blob/main/merscope_human_lung_cancer.ipynb), [nbviewer](https://nbviewer.org/github/CSOgroup/cellcharter_analyses/blob/main/merscope_human_lung_cancer.ipynb)) | Fig. 5 |
|  IMC Human Lung Cancer ([GitHub](https://github.com/CSOgroup/cellcharter_analyses/blob/main/imc_human_lung_cancer.ipynb), [nbviewer](https://nbviewer.org/github/CSOgroup/cellcharter_analyses/blob/main/imc_human_lung_cancer.ipynb)) | Fig. 6|
|  RNA-seq LUAD datasets ([GitHub](https://github.com/CSOgroup/cellcharter_analyses/blob/main/hypoxia_neutrophils_bulkDatasets_analyses.R)) | Fig. 6|

## Data
The datasets used are accessible in this [Figshare collection](https://figshare.com/collections/CellCharter_analyses/6414188) and must be placed (and eventually unzipped) in the respective directories inside the `data` directory.

## Installation
1. Create a conda or pyenv environment
2. Install Python < 3.11 and [PyTorch](https://pytorch.org) < 2.0.0. If you are planning to use a GPU, make sure to download and install the correct version of PyTorch first.
3. Install the library used for dimensionality reduction and batch effect removal according the data type you are planning to analyze:
    -   [scVI](https://github.com/scverse/scvi-tools) for spatial transcriptomics and/or epigenomics data such as 10x Visium and Xenium, Nanostring CosMx, Vizgen MERSCOPE, Stereo-seq, DBiT-seq, MERFISH and seqFISH data.
    -   A modified version of [scArches](https://github.com/theislab/scarches)'s TRVAE model for spatial proteomics data such as Akoya CODEX, Lunaphore COMET, CyCIF, IMC and MIBI-TOF data.
4. Install CellCharter using pip:

```bash
pip install cellcharter==0.2.0
```

We suggest using `mamba` to install the dependencies.
Installing the latest version of the dependencies (in particular `scvi-tools` and `spatialdata`) may lead to dependency conflicts. 
However, this should not be a problem because CellCharter doesn't use any of the mismatching features.

If you want to make sure to avoid conflicts, here is the list of the dependencies' versions that do not lead to a dependency conflict:
python=3.10

- pytorch=1.12.1
- jax=0.4.14
- jaxlib=0.4.14
- flax=0.7.2
- chex=0.1.7
- squidpy=1.3.0
- scikit-learn=1.3.0
- pycave=3.2.1
- spatialdata=0.0.12
- rasterio=1.3.8
- sknw=0.14
- urllib3=1.26.16
- typing-extensions=4.5.0
- numpy=1.23.4
- markdown-it-py=2.2.0
- torchmetrics=0.11.4

For the analysis of spatial transcriptomics and epigenomics data:

- scvi-tools=0.20.3

For the analysis of spatial proteomics data:
- scarches=0.5.9

We report here an example of an installation aimed at analyzing spatial transcriptomics data (and thus installing `scvi-tools`).
This example is based on a Linux CentOS 7 system with an NVIDIA A100 GPU.

```bash
conda create -n cellcharter-env -c conda-forge python=3.10 mamba=1.5.1
conda activate cellcharter-env
mamba install pytorch==1.12.1 torchvision==0.13.1 torchaudio==0.12.1 cudatoolkit=11.6 -c pytorch -c conda-forge
pip install jax==0.4.14 jaxlib==0.4.14 chex==0.1.7 flax==0.7.2 
pip install scvi-tools==0.20.3
pip install cellcharter==0.2.0
```
It may not be necessary to install `jax` and related libraries separately.
A different system may require different commands to install PyTorch and JAX. Refer to their respective documentation for more details.
