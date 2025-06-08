# SVM predictions on mouse data
import scanpy as sc
import numpy as np
import pandas as pd
import pickle

clf = pickle.load(open("svm_model_landscape.pickle", "rb"))
adata_pred = sc.read_h5ad('/home/chutianhao/R/Projects/snRNA_scRNA_hcc/project/svm/data/sc_mouse_all.h5ad')
adata_pred.obs.annotation.unique()

adata_pred = adata_pred[adata_pred.obs.annotation.isin(["B cells", 
                                                        "CAF", 
                                                        "Endothelial cells", 
                                                        "Hepatocytes", 
                                                        "Monocyte-linage", 
                                                        "NK cells", 
                                                        "Granulocytes", 
                                                        "Plasma cells",
                                                        "T cells"])]

# make sure the mouse and human celltype match
adata_pred.obs['annotation'] = (
    adata_pred.obs["annotation"]
    .map(lambda x: {"Granulocytes": "Myeloid cells", "Monocyte-linage": "Myeloid cells", "NK cells": "T/NK", "T cells": "T/NK"}.get(x, x))
    .astype("category")
)
adata_pred.obs.annotation.unique()

adata_pred.obs[['nCount_RNA']].describe()

X = adata_pred.X.copy()
Y = adata_pred.obs.copy()

X = X.toarray()
Y = Y.loc[:, 'annotation']
Y = Y.values.to_list()

# scale the data before training
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
X = scaler.fit_transform(X)

Y_pred = clf.predict(X)
Y_proba = clf.predict_proba(X)
Y_pred = Y_pred.tolist()

# calculate the model precision, recall and F1-score
from sklearn import metrics
print(metrics.classification_report(Y, Y_pred, digits=3))

Y_df = pd.DataFrame(data=Y_proba, 
                    columns=["B cells", 
                             "CAF", 
                             "Endothelial cells", 
                             "Hepatocytes", 
                             "Myeloid cells", 
                             "Plasma cells",
                             "T/NK"])
Y_df['annotation'] = Y
Y_df['prediction'] = Y_pred
Y_df.to_csv('data/cell_prob_landscape.csv')

# re-predict the TREM2+ signature to the whole mouse dataset
adata_mouse = sc.read_h5ad('data/mye_mouse_train.h5ad')
X = adata_mouse.X.copy().toarray()
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
X = scaler.fit_transform(X)

Y_pred = clf.predict(X)
Y_proba = clf.predict_proba(X)
Y_proba = pd.DataFrame(data=Y_proba, 
                       columns=['Classical monocyte', 
                             'FOLR2+ TAM', 
                             'Inflam-TAM', 
                             'Kupffer cell', 
                             'Nonclassical monocyte', 
                             'TREM2+ LAM'])

adata = sc.read_h5ad("/home/chutianhao/python/snRNA_scRNA_hcc/scvi/integration_mouse_myeloid/adata_filtered.h5ad")
sc.settings.verbosity = 3
sc.set_figure_params(dpi=150, frameon=True, color_map='viridis_r', dpi_save=600)
palette_d = ['#378C4F', '#6DB6FFFF', '#B6DBFFFF', '#E2A7CC', '#924900FF','#F5CDCD', '#D9579B', '#7BBC5E',  '#7464AA', '#006DDBFF', '#A59ACB'] # user defined discrete colors

sc.pl.umap(adata, color = 'annotation2', s=4,)
adata.obs['Trem2_proba'] = Y_proba['TREM2+ LAM'].values.tolist()
adata.obs['Folr2_proba'] = Y_proba['FOLR2+ TAM'].values.tolist()
adata.obs['RTM_proba'] = Y_proba['Kupffer cell'].values.tolist()
adata.obs['Inflam_proba'] = Y_proba['Inflam-TAM'].values.tolist()
adata.obs['clas_mono_proba'] = Y_proba['Classical monocyte'].values.tolist()
adata.obs['nocl_mono_proba'] = Y_proba['Nonclassical monocyte'].values.tolist()

sc.pl.umap(adata, color='Trem2_proba', s=4, save='Trem2_proba_umap.pdf')
# save this new anndata obj
adata.write("adata_proba_added.h5ad")