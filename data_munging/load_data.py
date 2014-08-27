import numpy as np

from frontier import frontier
from frontier.IO.BamcheckReader import BamcheckReader
from frontier.IO.AQCReader import AQCReader
from problem_def import CLASSES

DATA_DIR = "/store/sanger/ngsqc/bamcheck/bamcheck_2013dec25_ratios_out/"
TARGET_PATH = "/store/sanger/ngsqc/bamcheck/crohns-uc-table-a.2013dec25.manual_qc_update.txt"
USE_TARGETS = [1,-1]

statplexer = frontier.Statplexer(
    DATA_DIR,
    TARGET_PATH,
    CLASSES,
    BamcheckReader,
    AQCReader
)

all_parameters = statplexer.list_parameters()
data, target, levels = statplexer.get_data_by_target(all_parameters, USE_TARGETS)

data = data.transform({
    "error-rate":
        lambda x,f,i: x*100,
    "duplicate-mapped-ratio":
        lambda x,f,i: (f.get("reads-duplicated", i) / f.get("reads-mapped", i)) * 100.0,
    "percent-bases-mapped":
        lambda x,f,i: (f.get("bases-mapped-(cigar)", i) / f.get("total-length", i)) * 100.0,
    "percent-reads-qc-fail":
        lambda x,f,i: (f.get("reads-MQ0", i) / f.get("raw-total-sequences", i)) * 100.0,
    "quality-dropoff-fwd-mean-runmed-decline-range":
        lambda x,f,i: (f.get("quality-dropoff-fwd-mean-runmed-decline-high-value", i) - f.get("quality-dropoff-fwd-mean-runmed-decline-low-value", i)),
    "quality-dropoff-rev-mean-runmed-decline-range":
        lambda x,f,i: (f.get("quality-dropoff-rev-mean-runmed-decline-high-value", i) - f.get("quality-dropoff-rev-mean-runmed-decline-low-value", i)),
    "max-total-mean-baseline-deviation":
        # Currently get_single is required to force transform to be applied row-by-row
        lambda x,f,i: np.max([
            f.get_single("A-percent-total-mean-baseline-deviation", i),
            f.get_single("C-percent-total-mean-baseline-deviation", i),
            f.get_single("G-percent-total-mean-baseline-deviation", i),
            f.get_single("T-percent-total-mean-baseline-deviation", i)
        ]),
    "max-max-baseline-deviation":
        # Currently get_single is required to force transform to be applied row-by-row
        lambda x,f,i: np.max([
            f.get_single("A-percent-max-baseline-deviation", i),
            f.get_single("C-percent-max-baseline-deviation", i),
            f.get_single("G-percent-max-baseline-deviation", i),
            f.get_single("T-percent-max-baseline-deviation", i)
        ]),
    "range-baseline-deviation":
        # Currently get_single is required to force transform to be applied row-by-row
        lambda x,f,i: np.max([
                f.get_single("A-percent-max-baseline-deviation", i),
                f.get_single("C-percent-max-baseline-deviation", i),
                f.get_single("G-percent-max-baseline-deviation", i),
                f.get_single("T-percent-max-baseline-deviation", i)
            ]) - np.min([
                f.get_single("A-percent-max-baseline-deviation", i),
                f.get_single("C-percent-max-baseline-deviation", i),
                f.get_single("G-percent-max-baseline-deviation", i),
                f.get_single("T-percent-max-baseline-deviation", i)
            ]),
}, add_unknown=True)

data = data.exclude([
    # Converted to range
    "quality-dropoff-fwd-mean-runmed-decline-high-value",
    "quality-dropoff-fwd-mean-runmed-decline-low-value",

    # Converted to range
    "quality-dropoff-rev-mean-runmed-decline-high-value",
    "quality-dropoff-rev-mean-runmed-decline-low-value",

    # Converted to max of max
    "A-percent-total-mean-baseline-deviation",
    "C-percent-total-mean-baseline-deviation",
    "G-percent-total-mean-baseline-deviation",
    "T-percent-total-mean-baseline-deviation",

    # Convered to max of max
    "A-percent-max-baseline-deviation",
    "C-percent-max-baseline-deviation",
    "G-percent-max-baseline-deviation",
    "T-percent-max-baseline-deviation",

    # percent-max-baseline-deviation == max(percent-max-above-baseline, percent-max-below-baseline)
    "A-percent-max-above-baseline",
    "C-percent-max-above-baseline",
    "G-percent-max-above-baseline",
    "T-percent-max-above-baseline",

    "A-percent-max-below-baseline",
    "C-percent-max-below-baseline",
    "G-percent-max-below-baseline",
    "T-percent-max-below-baseline",

    # percent-total-mean-baseline-deviation == sum(percent-mean-above-baseline and percent-mean-below-baseline)
    "A-percent-mean-above-baseline",
    "C-percent-mean-above-baseline",
    "G-percent-mean-above-baseline",
    "T-percent-mean-above-baseline",

    "A-percent-mean-below-baseline",
    "C-percent-mean-below-baseline",
    "G-percent-mean-below-baseline",
    "T-percent-mean-below-baseline",

    # Raw counts
    "pairs-with-other-orientation",
    "outward-oriented-pairs",
    "reads-QC-failed",
    "pairs-on-different-chromosomes",
    "reads-duplicated",
    "reads-unmapped"
    "reads-unpaired",
    "reads-MQ0",
    "1st-fragments",
    "last-fragments"
    "inward-oriented-pairs",
    "mismatches",
    "reads-mapped",
    "reads-paired",
    "sequences",
    "raw-total-sequences",
    "bases-duplicated",
    "bases-mapped-(cigar)",
    "bases-mapped",
    "total-length",

    # No variance
    "bases-trimmed",
    "filtered-sequences",
    "is-paired",
    "is-sorted",
    "maximum-length",
    "non-primary-alignments",
    "quality-dropoff-high-iqr-threshold",
    "quality-dropoff-ignore-edge-cycles",
    "quality-dropoff-runmed-k",

    # Potentially unhelpful, refers to a read position
    "quality-dropoff-fwd-mean-runmed-decline-start-read-cycle",
    "quality-dropoff-fwd-mean-runmed-decline-end-read-cycle",
    "quality-dropoff-rev-mean-runmed-decline-start-read-cycle",
    "quality-dropoff-rev-mean-runmed-decline-end-read-cycle",

    "quality-dropoff-fwd-high-iqr-start-read-cycle",
    "quality-dropoff-fwd-high-iqr-end-read-cycle",
    "quality-dropoff-rev-high-iqr-start-read-cycle",
    "quality-dropoff-rev-high-iqr-end-read-cycle",
])

print "[DATA] %d samples" % len(data)
print "[DATA] %s levels" % str(levels)
print "[DATA] Samples by Level %s" % sorted(statplexer.count_targets_by_class(target).items())
