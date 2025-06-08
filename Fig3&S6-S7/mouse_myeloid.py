import scvi
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context

sc.settings.verbosity = 3
sc.set_figure_params(dpi=150, frameon=True, color_map='viridis_r', dpi_save=600)
palette_d = ['#378C4F', '#6DB6FFFF', '#B6DBFFFF', '#E2A7CC', '#924900FF','#F5CDCD', '#D9579B', '#7BBC5E',  '#7464AA', '#006DDBFF', '#A59ACB'] # user defined discrete colors

adata = sc.read_h5ad('mye_filter_new.h5ad')
adata.layers['counts'] = adata.X.copy() # preserve counts

# log transformation of the data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata # freeze the process, and store the normalized counts

sc.pp.neighbors(adata, use_rep='X_harmony')
sc.tl.umap(adata, min_dist=0.6)

sc.pl.umap(
    adata,
    color=['RNA_snn_res.0.8'],
    frameon=False
)

sc.tl.leiden(adata, resolution=0.8, key_added='sc_leiden')

adata.uns['log1p']['base'] = None # solve the key error
sc.tl.rank_genes_groups(adata, 'sc_leiden', method='wilcoxon')

pd.set_option('display.max_columns', None)
pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(20)

# Annotation
adata.obs['annotation2'] = (
    adata.obs["sc_leiden"]
    .map(lambda x: {"0": "Flor2+ TAM", "1": "Classical monocytes", "2": "Trem2+ LAM", "3": "Classical monocytes", "4": "RTM-like TAM", "5": "DC1", "6": "Trem2+ LAM", "7": "Cycling cells", "8": "Inflam-TAM", "9": "RTM-like TAM", "10": "DC2", "11": "DC3", "12": "Nonclassical monocytes", "13": "Classical monocytes"}.get(x, x))
    .astype("category")
)

sc.pl.umap(adata, color = 'annotation2', palette=palette_d, s=4, save="mye_landscape.pdf")

# Marker visualization
markers = {
    'Classical monocytes': ['F13a1', 'Fn1'], 
    'Cycling cells': ['Mki67', 'Top2a'], 
    'DC1': ['Cd209a', 'Cbfa2t3'], 
    'DC2': ['Tcf4', 'Siglech'], #pDC
    'DC3': ['Clec9a', 'Xcr1'], 
    'Folr2+ TAM': ['Folr2', 'Cd163'], 
    'Inflam-TAM': ['Ccl2', 'Ccl3'], 
    'Nonclassical monocytes': ['Ace', 'Ceacam1'], 
    'RTM-like TAM': ['Itgad', 'Vcam1'], 
    'Trem2+ LAM': ['Trem2', 'Spp1', 'Fabp5']
}

sc.pl.dotplot(adata, markers, groupby='annotation2', standard_scale='var', save="mye_marker_dotplot.pdf")



