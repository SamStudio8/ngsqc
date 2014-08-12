"""Use scikit-learn and frontier to identify useful parameters for quality control."""
__author__ = "Sam Nicholls <sn8@sanger.ac.uk>"
__copyright__ = "Copyright (c) Sam Nicholls 2014"
__version__ = "0.0.1"
__maintainer__ = "Sam Nicholls <sam@samnicholls.net>"

import sys

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

DATA_DIR = "/store/sanger/ngsqc/bamcheck/bamcheck_2013dec25_ratios_out/"
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
all_parameters = statplexer.list_parameters()
data, target, levels = statplexer.get_data_by_target(all_parameters, None)

# Withold 10% (stratify)
sss = StratifiedShuffleSplit(target, n_iter=1, test_size=PROP_VALIDATION)
for train_index, test_index in sss:
    X_train, X_test = data[train_index], data[test_index]
    y_train, y_test = target[train_index], target[test_index]

# Initialise classifier and CV folding
clf = DecisionTreeClassifier()
kf = StratifiedKFold(y_train, n_folds=NUM_FOLDS)

# Iterate over 10 folds on remaining 90%
for train_index, test_index in kf:
    Xf_train, Xf_test = X_train[train_index], X_train[test_index]
    yf_train, yf_test = y_train[train_index], y_train[test_index]
    clf.fit(Xf_train, yf_train)

    print clf.score(Xf_test, yf_test)

print clf.score(X_test, y_test)
