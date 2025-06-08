# scRNA and snRNA integration using scvi-tools
import scvi
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context

sc.settings.verbosity = 3
sc.set_figure_params(dpi=80, frameon=True, color_map='viridis_r', dpi_save=600)
palette_d = ['#378C4F', '#7BBC5E', '#E2A7CC', '#D9579B', '#A59ACB', '#006DDBFF', '#7464AA',  '#6DB6FFFF', '#B6DBFFFF', '#F5CDCD', '#924900FF'] # user defined discrete colors

adata_sc = scvi.data.read_h5ad('scRNA_anno.h5ad')
adata_sn = scvi.data.read_h5ad('snRNA_anno.h5ad')
adata = adata_sc.concatenate(adata_sn)

adata.layers['counts'] = adata.X.copy() # preserve counts
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata # freeze the process, and store the normalized counts

# select the top 1500 variable genes for data integration
sc.pp.highly_variable_genes(adata, 
                           flavor='seurat_v3', 
                           n_top_genes=1500, 
                           batch_key='sampletype', 
                            layer='counts',
                           subset=True)

# Train a scvi model
scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    categorical_covariate_keys=["batch", "patient_id"],
    continuous_covariate_keys=["mt_ratio", "ribo_ratio"]
)

model = scvi.model.SCVI(adata)
model.train(max_epochs=600, early_stopping=True)

adata.obsm['X_scvi'] = model.get_latent_representation()
adata.layers['scvi_normalized'] = model.get_normalized_expression(library_size=10e4)
adata.obsm['X_scvi'].shape

# plot the uncorrected cell embeddings
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_neighbors=20, n_pcs=30)
sc.tl.umap(adata, min_dist=0.3)
sc.pl.umap(adata, color=['sampletype', 'patient_id'], frameon=False, save='landscape_uncorrect.pdf')

# use the corrected cell embeddings from scvi
sc.pp.neighbors(adata, use_rep='X_scvi')
sc.tl.umap(adata, min_dist=0.5)

# use integration with a label transfer assuming snRNA-seq data is unknown
adata.obs["celltype_scanvi"] = 'Unknown'
sc_mask = adata.obs['sampletype'] == "scRNA-seq"
adata.obs["celltype_scanvi"][sc_mask] = adata.obs.celltype[sc_mask].values

# Integration with scANVI
lvae = scvi.model.SCANVI.from_scvi_model(model, adata=adata, labels_key='celltype_scanvi', unlabeled_category='Unknown')
lvae.train()
adata.obsm["X_scANVI"] = lvae.get_latent_representation(adata)
adata.layers['scANVI_normalized'] = lvae.get_normalized_expression(library_size=10e4)
sc.pp.neighbors(adata, use_rep='X_scANVI')
sc.tl.umap(adata, min_dist=0.3)

# create a integration cluster
adata.obs['int_clusters'] = (
    adata.obs["new_clusters_2"]
    .map(lambda x: {"B cells": "B cells", 
                    "CAF": "CAF", 
                    "Cycling cells_H": "Cycling cells", 
                    "Cycling cells_M": "Cycling cells", 
                    "Cycling cells_T": "Cycling cells", 
                    "Endothelial cells": "Endothelial cells", 
                    "Hepatocytes": "Hepatocytes", 
                    "Myeloid cells": "Myeloid cells", 
                    "NK cells": "NK cells", 
                    "Neutrophils": "Neutrophils", 
                    "Plasma cells": "Plasma cells", 
                    "T cells": "T cells"}.get(x, x)).astype("category")
)

sc.pl.umap(
    adata,
    color=['int_clusters'],
    palette=palette_d, 
    frameon=False,
    s=1,
    title='',
    save = '_landscape.pdf'
)

sc.pl.umap(adata, color=['patient_id'], frameon=True, title='', s=1, save = '_patient_id.pdf')
sc.pl.umap(adata, color = ['sampletype'], frameon=True, title = '', s=1, save='_technique.pdf')

# Marker visualization
markers = {
    'B cells': ['CD79A', 'BANK1'], 
    'CAF': ['PDGFRB', 'TAGLN'], 
    'Cycling cells': ['MKI67', 'TOP2A'], 
    'Endothelial cells': ['CD34', 'PECAM1'], 
    'Hepatocytes': ['APOB', 'TF', 'HNF4A'], 
    'Myeloid cells': ['C1QB', 'LYZ'], 
    'NK cells': ['NKG7', 'KLRB1'], 
    'Granulocytes': ['OSM', 'HCAR3'], 
    'Plasma cells': ['JCHAIN', 'MZB1'],  
    'T cells': ['IL7R', 'CD3D']
}

sc.pl.dotplot(adata, markers, 'int_clusters', color_map='Reds', use_raw=True, standard_scale='var')