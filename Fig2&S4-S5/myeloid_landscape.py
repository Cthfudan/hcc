# Integration of myeloid cells using scvi-tools

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

adata_sc = scvi.data.read_h5ad('scRNA_myeanno.h5ad')
adata_sn = scvi.data.read_h5ad('snRNA_myeanno.h5ad')

adata = adata_sc.concatenate(adata_sn)

adata.layers['counts'] = adata.X.copy() # preserve counts
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

# select HVGs
sc.pp.highly_variable_genes(
    adata,
    flavor="seurat_v3",
    n_top_genes=3000,
    layer="counts",
    batch_key="sampletype",
    subset=True
)

# setup anndata
scvi.model.SCVI.setup_anndata(
    adata,
    layer='counts',
    categorical_covariate_keys=["batch", "patient_id"],
    continuous_covariate_keys=["mt_ratio", "ribo_ratio"]
)

vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood='nb')
vae.train(max_epochs=500, early_stopping=True)

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

# integration with scANVI
# we here assume that cells from snRNA-seq is un-annotated
adata.obs["celltype_scanvi"] = 'Unknown'
sc_mask = adata.obs['sampletype'] == "scRNA-seq"
adata.obs["celltype_scanvi"][sc_mask] = adata.obs.celltype[sc_mask].values

np.unique(adata.obs["celltype_scanvi"], return_counts = True)

lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=adata,
    unlabeled_category="Unknown",
    labels_key="celltype_scanvi",
)

# model training
lvae.train(max_epochs=50)

# predict labels in snRNA-seq, and get the latent space
adata.obs["C_scANVI"] = lvae.predict(adata)
adata.obsm["X_scANVI"] = lvae.get_latent_representation(adata)
adata.layers['scanvi_normalized'] = lvae.get_normalized_expression(library_size=10e4)

# use the corrected cell embeddings from scanvi
# the label prediction from scanvi seem to be raletively goo# 
sc.pp.neighbors(adata, use_rep='X_scANVI')
sc.tl.umap(adata, min_dist=0.2)

# clustering on the scanvi latent space
sc.tl.leiden(adata, key_added='scanvi_leiden', resolution=0.7)

# find markers for each cluster
adata.uns['log1p']['base'] = None # solve the key error
sc.tl.rank_genes_groups(adata, 'scanvi_leiden', method='wilcoxon')

# top 20 genes in each clusters
pd.set_option('display.max_columns', None)
pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(20)

# re-annotation
adata.obs['int_clusters'] = (
    adata.obs["scanvi_leiden"]
    .map(lambda x: {"0": "TREM2+ LAM", "1": "FOLR2+ TAM", "2": "Inflam-TAM", "3": "Classical monocyte", "4": "Kupffer cell", "5": "cDC2", "6": "Cycling cell", "7": "FOLR2+ TAM", "8": "IFN-TAM", "9": "TREM2+ LAM", "10": "Nonclassical monocyte", "11": "cDC1", "12": "cDC3", "13": "Classical monocyte"}.get(x, x))
    .astype("category")
)
# this expression, from github, is very useful for reassignment of certain variables in dataframe

# downstream analysis
# marker visualization
markers = {
    'Classical monocyte': ['LYZ', 'S100A9'], 
    'Cycling cell': ['TOP2A', 'MKI67'], 
    'FOLR2+ TAM': ['FOLR2', 'MRC1', 'CD163'], # FOLR2+ CD163+ onco-fetal macrophage
    'IFN-TAM': ['IFIT1', 'IFI44L'], 
    'Inflam-TAM': ['CCL4', 'CCL3L1', 'IL1B'],  # inflammation
    'Kupffer cell': ['CD5L', 'VCAM1'],  # tissue resident
    'Nonclassical monocyte': ['FCGR3A', 'LILRB2'], 
    'TREM2+ LAM': ['TREM2', 'SPP1', 'FABP5'],  # lipid metabolism
    'cDC1': ['CLEC9A', 'CPNE3'], 
    'cDC2': ['CD1C', 'CLEC10A'], 
    'cDC3': ['LAMP3', 'CCR7'], 
}

sc.pl.dotplot(adata, markers, 'int_clusters', color_map='Reds', use_raw=True, standard_scale='var', expression_cutoff=0.5, save="marker_myeloid.pdf")


