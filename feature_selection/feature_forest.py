# -*- coding: utf-8 -*-
"""Use scikit-learn and frontier to identify useful parameters for quality
control by constructing and querying a series of random forests."""
__author__ = "Sam Nicholls <sn8@sanger.ac.uk>"
__copyright__ = "Copyright (c) Sam Nicholls 2014"
__version__ = "0.0.1"
__maintainer__ = "Sam Nicholls <sam@samnicholls.net>"

import sys

import numpy as np

from sklearn.cross_validation import StratifiedKFold, StratifiedShuffleSplit
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier

from data_munging.load_data import data, target, levels
from util import plot_tree


all_parameters = data.frontier_labels

NUM_FOLDS = 10
PROP_VALIDATION = 0.1 #10%

NUM_TREES = 500
NUM_PARAMETERS = 10

# Withhold some proportion of the data set for validation later (stratify)
sss = StratifiedShuffleSplit(target, n_iter=1, test_size=PROP_VALIDATION)
for train_index, test_index in sss:
    X_train, X_test = data[train_index], data[test_index]
    y_train, y_test = target[train_index], target[test_index]

# Initialise folds on remaining data and store indexes for use later
kf = StratifiedKFold(y_train, n_folds=NUM_FOLDS)
fold_indexes = []
for train_index, test_index in kf:
    fold_indexes.append({
            "train_index": train_index,
            "test_index": test_index
    })

# Init score storing structures
importance_scores = np.zeros([NUM_FOLDS, len(all_parameters)])
importance_scores_stdev = np.zeros([NUM_FOLDS, len(all_parameters)])
tree_cv_scores = np.zeros(NUM_FOLDS)
parameter_union = np.zeros(len(all_parameters), dtype=int)

# For each fold, build a forest, skim the NUM_PARAMETERS best features as
# measured by the classifier's feature_importance property, then fit and score a
# single decision tree on that feature set.
for n_fold, indexer in enumerate(fold_indexes):
    print "\n[FRST] Constructing Forest#%d" % (n_fold + 1)
    Xf_train, Xf_test = X_train[indexer["train_index"]], X_train[indexer["test_index"]]
    yf_train, yf_test = y_train[indexer["train_index"]], y_train[indexer["test_index"]]

    clf = RandomForestClassifier(n_estimators=NUM_TREES, random_state=0)
    clf.fit(Xf_train, yf_train)

    importance_scores[n_fold] = clf.feature_importances_
    importance_scores_stdev[n_fold] = np.std([tree.feature_importances_ for tree in clf.estimators_], axis=0)

    param_mask = np.zeros(len(all_parameters), dtype=int)
    print "\n[    ] µImp.\tσp.Imp.\tName"
    for i in reversed(importance_scores[n_fold].argsort()[-NUM_PARAMETERS:]):
        param_mask[i] = 1
        parameter_union[i] += 1
        print "[PARM] %.4f\t%.4f\t%s" % (
                importance_scores[n_fold][i],
                importance_scores_stdev[n_fold][i],
                all_parameters[i]
        )

    clf_t = DecisionTreeClassifier()
    print "[TREE] Fitting Tree#%d" % (n_fold + 1)
    param_index_list = np.where(param_mask > 0)[0]
    clf_t.fit(Xf_train[:, param_index_list], yf_train)
    tree_cv_scores[n_fold] = clf_t.score(Xf_test[:, param_index_list], yf_test)
    print "[TREE] CV %.2f" % tree_cv_scores[n_fold]

print "\n[ENSM] Average Single Tree CV %.2f" % np.average(tree_cv_scores)

print "\n[    ] #Used\tµImp.\tσp.Imp.\tName"
importance_indices = (np.average(importance_scores, axis=0).argsort())[::-1]
importance_averages = np.average(importance_scores, axis=0)
importance_std_devs = np.average(importance_scores_stdev, axis=0)
for i in importance_indices:
    mask = parameter_union[i]
    if mask > 0:
        print "[PARM] %d\t%.4f\t%.4f\t%s" % (
                parameter_union[i],
                importance_averages[i],
                # TODO Is it safe to do variance pooling?
                importance_std_devs[i],
                all_parameters[i]
        )
print "[    ] Sorted on Average Importance (µImp)."

# Fit and score a single decision tree using all available training data and
# the union of all previously forest selected parameters and validate against
# the withheld data set.
clf_t = DecisionTreeClassifier(criterion="entropy", min_samples_split=25, min_samples_leaf=10)
print "\n[TREE] Fit parameter union tree"
param_index_list = np.where(parameter_union > 0)[0]
clf_t.fit(X_train[:, param_index_list], y_train)
print "[TREE] CV %.2f" % clf_t.score(X_test[:, param_index_list], y_test)

parameter_names = []
for i in param_index_list:
    parameter_names.append(all_parameters[i])
plot_tree("feature_selection/plots/final_tree.pdf", clf_t, parameter_names)
