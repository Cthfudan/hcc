# Integration of endo cells using scvi-tools

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

import scanpy as sc
import scvi

from matplotlib.pyplot import rc_context

sc.settings.verbosity = 3
sc.set_figure_params(dpi=80, frameon=False, color_map='viridis_r', dpi_save=600)
palette_d = ['#378C4F', '#6DB6FFFF', '#F5CDCD', '#D9579B', '#A59ACB', '#7464AA', '#7BBC5E',  '#006DDBFF', '#B6DBFFFF', '#E2A7CC', '#924900FF']

adata_sc = scvi.data.read_h5ad('scRNA_caf.h5ad')
adata_sn = scvi.data.read_h5ad('snRNA_caf.h5ad')

adata = adata_sc.concatenate(adata_sn)

adata.layers['counts'] = adata.X.copy() # preserve counts
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

# select HVGs
sc.pp.highly_variable_genes(
    adata,
    flavor="seurat_v3",
    n_top_genes=2000,
    layer="counts",
    batch_key="sampletype",
    subset=True
)

adata = adata[~adata.obs.patient_id.isin(["snPT7", "snPT9", "snPT8", "snPT1"])] # remove patients with too few CAFs

# setup anndata
scvi.model.SCVI.setup_anndata(
    adata,
    layer='counts',
    categorical_covariate_keys=["sampletype", "patient_id"],
    continuous_covariate_keys=["mt_ratio"]
)

vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood='nb')
vae.train()

# get latent representation 
adata.obsm['X_scvi'] = vae.get_latent_representation()
adata.layers['scvi_normalized'] = vae.get_normalized_expression(library_size=10e4)

# use the corrected cell embeddings from scvi
sc.pp.neighbors(adata, use_rep='X_scvi')
sc.tl.umap(adata, min_dist=0.3)
sc.pl.umap(
    adata,
    color=["sampletype", "patient_id", "celltype"],
    frameon=False
)

# intergration using scANVI
np.unique(adata.obs["celltype"], return_counts = True)

lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=adata,
    unlabeled_category="unknown",
    labels_key="celltype",
)

lvae.train(max_epochs=25)

# predict labels in snRNA-seq, and get the latent space
adata.obs["C_scANVI"] = lvae.predict(adata)
adata.obsm["X_scANVI"] = lvae.get_latent_representation(adata)
adata.layers['scanvi_normalized'] = lvae.get_normalized_expression(library_size=10e4)

# use the corrected cell embeddings from scanvi
sc.pp.neighbors(adata, use_rep='X_scANVI')
sc.tl.umap(adata, min_dist=0.5)

sc.tl.leiden(adata, resolution=0.8, key_added="scanvi_leiden")

sc.pl.umap(adata, color=["C_scANVI"], palette=["#00468B", "#ED0000", "#42B540", "#0099B4"], save="caf_landscape.pdf")

adata.write("adata_caf_integ.h5ad")

# find markers for each cluster
adata.uns['log1p']['base'] = None # solve the key error
sc.tl.rank_genes_groups(adata, 'scanvi_leiden', method='wilcoxon', use_raw=True)

# top 20 genes in each clusters
pd.set_option('display.max_columns', None)
pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(20)

sc.pl.umap(adata, color=["CD74", "MYH11", "CD36", "TOP2A"])

# re-annotation
adata.obs['celltype_scanvi'] = (
    adata.obs["scanvi_leiden"]
    .map(lambda x: {"0": "lpCAF", "1": "mCAF", "2": "vCAF", "3": "apCAF", "4": "vCAF", "5": "mCAF", "6": "mCAF", "7": "mCAF", "8": "vCAF"}.get(x, x))
    .astype("category")
)

markers = {"apCAF": ["CD74", "CCL5"],
                      "lpCAF": ["CD36", "CPM"],
           "mCAF": ["LUM", "VCAN"], 
                      "vCAF": ["MYH11", "MCAM"] }

sc.pl.dotplot(adata, markers, "celltype_scanvi", color_map = "Reds", use_raw=True, standard_scale = "var")

adata.write("adata_caf_integ.h5ad")

