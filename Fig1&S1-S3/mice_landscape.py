import scvi
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context

sc.settings.verbosity = 3
sc.set_figure_params(dpi=150, frameon=True, color_map='viridis_r', dpi_save=600)
palette_d = ['#378C4F', '#7BBC5E', '#E2A7CC', '#D9579B', '#6DB6FFFF','#A59ACB', '#006DDBFF',  '#7464AA', '#B6DBFFFF', '#F5CDCD', '#924900FF'] # user defined discrete colors

adata = sc.read_h5ad("sc_mice_new.h5ad")
print(adata.X)
adata.layers['counts'] = adata.X.copy() # preserve counts

# log transformation of the data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

sc.pp.neighbors(adata, use_rep='X_harmony')
sc.tl.umap(adata, min_dist=0.4)
sc.tl.leiden(adata, resolution=1.0, key_added='sc_leiden')
sc.tl.rank_genes_groups(adata, 'sc_leiden', method='wilcoxon')
pd.set_option('display.max_columns', None)
pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(20)

# Annotation
adata.obs['annotation'] = (
    adata.obs["sc_leiden"]
    .map(lambda x: {"0": "Hepatocytes", 
                    "1": "Granulocytes", 
                    "2": "Mono/macro", 
                    "3": "Mono/macro", 
                    "4": "Mono/macro", 
                    "5": "Hepatocytes", 
                    "6": "T cells", 
                    "7": "Endothelial cells", 
                    "8": "NK cells", 
                    "9": "B cells", 
                    "10": "T cells", 
                    "11": "Granulocytes", 
                    "12": "Cycling cells", 
                    "13": "Hepatocytes", 
                    "14": "T cells", 
                    "15": "Granulocytes", 
                    "16": "Mono/macro", 
                    "17": "Endothelial cells", 
                    "18": "Monocyte-linage", 
                    "19": "Hepatocytes", 
                    "20": "Monocyte-linage", 
                    "21": "Monocyte-linage", 
                    "22": "Monocyte-linage", 
                    "23": "Cycling cells", 
                    "24": "Granulocytes", 
                    "25": "Hepatocytes", 
                    "26": "Hepatocytes", 
                    "27": "Plasma cells", 
                    "28": "Hepatocytes", 
                    "29": "Endothelial cells", 
                    "30": "Monocyte-linage", 
                    "31": "CAF", \
                    "32": "Granulocytes"}.get(x, x))
    .astype("category")
)

adata = adata[adata.obs['RNA_snn_res.1'] != '31']

sc.pl.umap(adata, color = 'annotation', palette=palette_d, frameon=False, s=1, save='_landscape.pdf')

from pandas.api.types import CategoricalDtype
cell_type = CategoricalDtype(categories=["B cells", 
                                         "Plasma cells",
                                         "Endothelial cells",
                                         "Fibroblasts",
                                         "Hepatocytes",
                                         "Mono/macro",
                                         "Granulocytes",
                                         "NK cells",
                                         "T cells",
                                         "Cycling cells"], ordered=True)

adata.obs['annotation_dotplot'] = (
    adata.obs["annotation"]
    .map(lambda x: {"Monocyte-linage": "Mono/macro", "CAF": "Fibroblasts"}.get(x, x))
    .astype("category")
)

adata.obs["annotation_dotplot"] = adata.obs["annotation_dotplot"].astype(cell_type)
# Marker visualization
markers = {
    'B cells': ['Cd79a', 'Ms4a1'], 
    'Plasma cells': ['Jchain', 'Mzb1'],  
    'Endothelial cells': ['Ptprb', 'Maf'], 
    'Fibroblasts': ['Ecm1', 'Plvap'], 
    'Hepatocytes': ['Apob', 'Ttr'], 
    'Mono/macro': ['C1qb', 'Lyz2'], 
    'Granulocytes': ['Itgam', 'Mcl1'], 
    'NK cells': ['Xcl1', 'Klrk1'],  
    'T cells': ['Il7r', 'Cd3d'], 
    'Cycling cells': ['Mki67', 'Top2a'], 
}
sc.pl.dotplot(adata, markers,'annotation_dotplot', color_map='Reds', standard_scale='var', save="marker_gene_dotplot.pdf")

sc.pl.umap(adata, color="sample_id", s=1, save="_sample_id.pdf")
sc.pl.umap(adata, color="condition", s=1, palette=["#FF7F0E", "#1F77B4"], save="_condition.pdf")
adata.write('sc_mice_annotation.h5ad')


