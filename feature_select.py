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

CLASSES = {
    "pass": {
        "class": ["pass"],
        "names": ["pass", "passed"],
        "code": 1,
    },
    "fail": {
        "class": ["fail"],
        "names": ["fail", "failed"],
        "code": -1,
    },
    "warn": {
        "class": ["warn"],
        "names": ["warn", "warning"],
        "code": 0,
    },
}

DATA_DIR = "/store/sanger/ngsqc/bamcheck/bamcheck_2013dec25_ratios_out-100/"
TARGET_PATH = "/store/sanger/ngsqc/bamcheck/crohns-uc-table-a.2013dec25.manual_qc_update.txt"
NUM_FOLDS = 5
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
all_parameters = statplexer.list_parameters()
data, target, levels = statplexer.get_data_by_target(all_parameters, None)

# Withhold 10% (stratify)
sss = StratifiedShuffleSplit(target, n_iter=1, test_size=PROP_VALIDATION)
for train_index, test_index in sss:
    X_train, X_test = data[train_index], data[test_index]
    y_train, y_test = target[train_index], target[test_index]

# While there are some number of parameters, conduct folding and elimination
param_mask = np.ones(len(all_parameters), dtype=int)
while sum(param_mask) >= 5:
    # Initialise classifier and CV folding on remaining data
    kf = StratifiedKFold(y_train, n_folds=NUM_FOLDS)

    # Iterate over folds
    scores = np.zeros([len(all_parameters), NUM_FOLDS])
    n_fold = 0
    for train_index, test_index in kf:
        Xf_train, Xf_test = X_train[train_index], X_train[test_index]
        yf_train, yf_test = y_train[train_index], y_train[test_index]

        # Ignore a parameter each time
        for n_col, mask in enumerate(param_mask):
            if mask == 0:
                continue


            # Fit and score
            print "[%d][%d] Fitting and Scoring..." % (n_fold + 1, n_col + 1)
            param_mask[n_col] = 0
            clf = DecisionTreeClassifier()
            clf.fit(Xf_train[:, param_mask], yf_train)
            scores[n_col][n_fold] = clf.score(Xf_test[:, param_mask], yf_test)
            param_mask[n_col] = 1

        n_fold += 1

    averages = np.zeros(len(all_parameters))
    for n_col, mask in enumerate(param_mask):
        if mask == 0:
            continue
        averages[n_col] = np.average(scores[n_col])
    min_avg = np.min(averages[np.nonzero(averages)])
    print min_avg

    for n_col, mask in enumerate(param_mask):
        if mask == 0:
            continue
        if averages[n_col] <= min_avg:
            print averages[n_col]
            param_mask[n_col] = 0
            print "[ ][  ] Pruned parameter %s" % all_parameters[n_col]


# Now validate against withheld validation set
clf = DecisionTreeClassifier()
clf.fit(X_train[:, param_mask], y_train)
print clf.score(X_test, y_test)
