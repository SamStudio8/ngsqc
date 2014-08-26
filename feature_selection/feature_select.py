"""Use scikit-learn and frontier to identify useful parameters for quality control."""
__author__ = "Sam Nicholls <sn8@sanger.ac.uk>"
__copyright__ = "Copyright (c) Sam Nicholls 2014"
__version__ = "0.0.1"
__maintainer__ = "Sam Nicholls <sam@samnicholls.net>"

import sys

import numpy as np

from frontier import frontier
from frontier.IO.BamcheckReader import BamcheckReader
from frontier.IO.AQCReader import AQCReader

from sklearn.cross_validation import StratifiedKFold, StratifiedShuffleSplit
from sklearn.tree import DecisionTreeClassifier

from data_munging.problem_def import CLASSES, NO_VARIANCE, RAW_COUNTS

DATA_DIR = "/store/sanger/ngsqc/bamcheck/bamcheck_2013dec25_ratios_out-1000/"
TARGET_PATH = "/store/sanger/ngsqc/bamcheck/crohns-uc-table-a.2013dec25.manual_qc_update.txt"
NUM_FOLDS = 10
PROP_VALIDATION = 0.1 # 10%

def validate_data_set(data_set):
    #raise Exception("%s is not a valid option for data_set" % data_set)
    return True

def validate_param_set(parameter_set):
    #raise Exception("%s is not a valid option for parameter_set" % parameter_set)
    return True

def setup():
    try:
        validate_data_set(None)
        validate_param_set(None)
    except Exception, e:
        print e
        sys.exit(1)

    statplexer = frontier.Statplexer(
        DATA_DIR,
        TARGET_PATH,
        CLASSES,
        BamcheckReader,
        AQCReader
    )

    return statplexer

def iterate():
    pass

# Setup structures
statplexer = setup()
all_parameters = statplexer.exclude_parameters(NO_VARIANCE + RAW_COUNTS)
data, target, levels = statplexer.get_data_by_target(all_parameters, None)

# Withhold 10% (stratify)
sss = StratifiedShuffleSplit(target, n_iter=1, test_size=PROP_VALIDATION)
for train_index, test_index in sss:
    X_train, X_test = data[train_index], data[test_index]
    y_train, y_test = target[train_index], target[test_index]

# Initialise folds on remaining data and store fold indexes for use later
kf = StratifiedKFold(y_train, n_folds=NUM_FOLDS)
fold_indexes = []
for train_index, test_index in kf:
    fold_indexes.append({
            "train_index": train_index,
            "test_index": test_index
    })

# While there are some number of parameters, conduct folding and elimination
param_mask = np.ones(len(all_parameters), dtype=int)
last_cv = 0.0
while sum(param_mask) >= 5:

    # Iterate over folds
    scores = np.zeros([len(all_parameters), NUM_FOLDS])
    for n_fold, indexer in enumerate(fold_indexes):
        Xf_train, Xf_test = X_train[indexer["train_index"]], X_train[indexer["test_index"]]
        yf_train, yf_test = y_train[indexer["train_index"]], y_train[indexer["test_index"]]

        # Ignore a parameter each time
        for n_col, mask in enumerate(param_mask):
            if mask == 0:
                continue

            # Fit and score
            print "[%d][%d] Fitting and Scoring..." % (n_fold + 1, n_col + 1)
            param_mask[n_col] = 0
            clf = DecisionTreeClassifier()
            param_index_list = np.where(param_mask > 0)[0]
            clf.fit(Xf_train[:, param_index_list], yf_train)
            scores[n_col][n_fold] = clf.score(Xf_test[:, param_index_list], yf_test)
            param_mask[n_col] = 1

        n_fold += 1

    averages = np.zeros(len(all_parameters))
    for n_col, mask in enumerate(param_mask):
        if mask == 0:
            continue
        averages[n_col] = np.average(scores[n_col])
    max_avg = np.max(averages)
    print max_avg

    for n_col, mask in enumerate(param_mask):
        if mask == 0:
            continue
        if averages[n_col] >= max_avg:
            param_mask[n_col] = 0
            print "[ ][  ] Pruned parameter %s" % all_parameters[n_col]
            print "[ ][  ] Average CV moved from %.2f to %.2f" % (last_cv, max_avg)

    last_cv = max_avg

for n_col, mask in enumerate(param_mask):
    if mask > 0:
        print "[ ][  ] Maintained parameter %s" % all_parameters[n_col]

# Now validate against withheld validation set
clf = DecisionTreeClassifier()
param_index_list = np.where(param_mask > 0)[0]
clf.fit(X_train[:, param_index_list], y_train)
print clf.score(X_test[:, param_index_list], y_test)
