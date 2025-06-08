# train a SVM classifier on human data
import scanpy as sc
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize

# import data 
train_file = "/home/chutianhao/R/Projects/snRNA_scRNA_hcc/project/svm/data/sc_human_train.h5ad"
adata_train = sc.read_h5ad(train_file)

# prepare training data
X = adata_train.X.copy()
Y = adata_train.obs.copy()

group = "clusters_2"
Y = Y.loc[:, group]
Y = Y.values.to_list()

# scale the data before training
from sklearn.preprocessing import StandardScaler
X = X.toarray()
scaler = StandardScaler()
X = scaler.fit_transform(X)

# search the best C and gamma values, using grid search
from sklearn.svm import SVC
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.model_selection import GridSearchCV

C_range = np.logspace(-2, 7, 10)
gamma_range = np.logspace(-9, 0, 10)
param_grid = dict(gamma=gamma_range, C=C_range)
cv = StratifiedShuffleSplit(n_splits=5, test_size=0.2, random_state=42)
grid = GridSearchCV(SVC(cache_size=400000), param_grid=param_grid, cv=cv)
grid.fit(X, Y)

print(
    "The best parameters are %s with a score of %0.2f"
    % (grid.best_params_, grid.best_score_)
)

# train the classifier with RBF kernel
from sklearn import svm 

clf = svm.SVC(kernel='rbf', C = 100, gamma=0.001, cache_size=400000, probability=True)
clf.fit(X, Y)

# test the model on test dataset
adata_test = sc.read_h5ad("/home/chutianhao/R/Projects/snRNA_scRNA_hcc/project/svm/data/sc_human_test.h5ad")
X_test = adata_test.X.copy()
Y_test = adata_test.obs.copy()

X_test = X_test.toarray()
X_test = scaler.fit_transform(X_test)

Y_test = Y_test.loc[:, group]
Y_test = Y_test.values.to_list()
Y_pred = clf.predict(X_test)
Y_proba = clf.predict_proba(X_test)
Y_pred = Y_pred.tolist()

match = [i for i, j in zip(Y_test, Y_pred) if i == j]
nomatch = [i for i, j in zip(Y_test, Y_pred) if i != j]

print(len(match), len(nomatch), len(Y_pred))

# calculate the model precision, recall and F1-score
from sklearn import metrics

print(metrics.classification_report(Y_test, Y_pred, digits=3))

# save the trained svm model
import pickle
filename = "svm_model_landscape.pickle"
pickle.dump(clf, open(filename, "wb"))