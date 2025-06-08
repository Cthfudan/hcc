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

adata_sn = scvi.data.read_h5ad('snRNA_TNK.h5ad')
adata_sc = scvi.data.read_h5ad('scRNA_TNK.h5ad')

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
    unlabeled_category="Unknown",
    labels_key="celltype",
)

# model training
lvae.train(max_epochs=25)

# predict labels in snRNA-seq, and get the latent space
adata.obs["C_scANVI"] = lvae.predict(adata)
adata.obsm["X_scANVI"] = lvae.get_latent_representation(adata)
adata.layers['scanvi_normalized'] = lvae.get_normalized_expression(library_size=10e4)

# use the corrected cell embeddings from scanvi
# the label prediction from scanvi seem to be raletively goo# 
sc.pp.neighbors(adata, use_rep='X_scANVI')
sc.tl.umap(adata, min_dist=0.2)

sc.tl.leiden(adata, resolution=1, key_added='scanvi_leiden_1.0')
sc.pl.umap(adata, color=['scanvi_leiden_1.0', "celltype"])

# find markers for each cluster
adata.uns['log1p']['base'] = None # solve the key error
sc.tl.rank_genes_groups(adata, 'celltype', method='wilcoxon', use_raw=True)

# top 20 genes in each clusters
pd.set_option('display.max_columns', None)
pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(50)

# re-annotation
adata.obs['celltype_scanvi'] = (
    adata.obs["celltype"]
    .map(lambda x: {"CCR6+ CD8+ T": "CD8+ Trm", "CCR7+ CD4+ T": "CD4+ T-naive", "CD4+ Tcm": "CD4+ Tem", "Effector CD8+ T": "CD8+ Tem", "NKT": "CD8+ Temra", "TNFRSF9+ CD8+ T": "CD8+ Tex"}.get(x, x))
    .astype("category")
)

sc.pl.umap(adata, color = 'celltype_scanvi', save="landscape_TNK.pdf")

# marker gene visualization
markers = {
    "Canonical": ['CD3D', 'CD4', 'CD8A'], 
    'Activated CD4+ T': ['FOS', 'JUNB'],
    'Activated CD8+ T': ['JUND', 'HSP90AB1'], 
    'CD8+ Trm': ['CCL20', 'CCR6'], # Th17 like chemokines
    'CD4+ T-naive': ['CCR7', 'SARAF'], 
    'CD4+ Tem': ['CD40LG', 'TXNIP'], 
    'CD16+ NK': ['FCGR3A', 'FGFBP2'], 
    'CXCL13+ CD4+ T': ['CXCL13', 'IL6ST'], 
    'CD8+ Tem': ['GZMK', 'IFNG'], 
    'CD8+ Temra': ['PRF1', 'GZMH'],  # inflammation
    'Proliferating T': ['TOP2A', 'MKI67'], 
    'CD8+ Tex': ['TNFRSF9', 'PDCD1'],  # tissue resident
    'Treg': ['FOXP3', 'CTLA4'], 
    'XCL1+ NK': ['XCL1', 'XCL2']
}

adata_CD8T = adata[adata.obs.celltype_scanvi.isin(['Activated CD8+ T', 'CD8+ Trm', 'CD8+ Tem', 'CD8+ Temra', 'CD8+ Tex'])]

markers_cd8 = {
    'Activated CD8+ T': ['JUND', 'HSP90AB1'], 
    'CD8+ Trm': ['CCL20', 'CCR6'], # Th17 like chemokines
    'CD8+ Tem': ['GZMK', 'IFNG'], 
    'CD8+ Temra': ['PRF1', 'GZMH'],  # inflammation
    'CD8+ Tex': ['TNFRSF9', 'PDCD1', 'LAG3'],  # tissue resident
}

sc.pl.dotplot(adata_CD8T, markers_cd8, 'celltype_scanvi', color_map='Reds', use_raw=True, standard_scale='var', save="CD8_marker.pdf")

adata_CD4T = adata[adata.obs.celltype_scanvi.isin(['Activated CD4+ T', 'CD4+ T-naive', 'CD4+ Tem', 'CXCL13+ CD4+ T', 'Treg'])]

markers_cd4 = {
    "Canonical": ['CD3D', 'CD4', 'CD8A'], 
    'Activated CD4+ T': ['FOS', 'JUNB'],
    'CD4+ T-naive': ['CCR7', "SARAF"], 
    'CD4+ Tem': ['KLRB1', 'TXNIP'],
    'CXCL13+ CD4+ T': ['CXCL13', 'IL6ST'], 
    'Treg': ['FOXP3', 'CTLA4'],
}

sc.pl.dotplot(adata_CD4T, markers_cd4, 'celltype_scanvi', color_map='Reds', use_raw=True, standard_scale='var', save="cd4_marker.pdf")

adata_nk = adata[adata.obs.celltype_scanvi.isin(['XCL1+ NK', 'CD16+ NK'])]

markers_NK = {
    'CD16+ NK': ['FCGR3A', 'FGFBP2'], 
    'XCL1+ NK': ['XCL1', 'XCL2']
}

sc.pl.dotplot(adata_nk, markers_NK, 'celltype_scanvi', color_map='Reds', use_raw=True, standard_scale='var', save="markers_NK.pdf")




