import numpy as np

from frontier import frontier
from frontier.IO.BamcheckReader import BamcheckReader
from frontier.IO.AQCReader import AQCReader
from problem_def import CLASSES

DATA_DIR = "/store/sanger/ngsqc/bamcheck/bamcheck_2013dec25_ratios_out-50/"
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

# Check correct...
data['error-rate'] *= 100
data['duplicate-mapped-ratio'] = (data["reads-duplicated"] / data["reads-mapped"]) * 100
data['percent-bases-mapped'] = (data["bases-mapped-(cigar)"] / data["total-length"]) * 100
data['percent-reads-qc-fail'] = (data["reads-MQ0"] / data["raw-total-sequences"]) * 100.
data["quality-dropoff-fwd-mean-runmed-decline-range"] = data["quality-dropoff-fwd-mean-runmed-decline-high-value"] - data["quality-dropoff-fwd-mean-runmed-decline-low-value"]
data["quality-dropoff-rev-mean-runmed-decline-range"] = data["quality-dropoff-rev-mean-runmed-decline-high-value"] - data["quality-dropoff-rev-mean-runmed-decline-low-value"]
data["max-total-mean-baseline-deviation"] = data[[
            "A-percent-total-mean-baseline-deviation",
            "C-percent-total-mean-baseline-deviation",
            "G-percent-total-mean-baseline-deviation",
            "T-percent-total-mean-baseline-deviation"
        ]].max(axis=1)
data["max-max-baseline-deviation"] = data[[
            "A-percent-max-baseline-deviation",
            "C-percent-max-baseline-deviation",
            "G-percent-max-baseline-deviation",
            "T-percent-max-baseline-deviation",
        ]].max(axis=1)

data["range-baseline-deviation"] = data[[
            "A-percent-max-baseline-deviation",
            "C-percent-max-baseline-deviation",
            "G-percent-max-baseline-deviation",
            "T-percent-max-baseline-deviation"
        ]].max(axis=1) - data[[
            "A-percent-max-baseline-deviation",
            "C-percent-max-baseline-deviation",
            "G-percent-max-baseline-deviation",
            "T-percent-max-baseline-deviation"
        ]].min(axis=1)

# Drop cols...
data = data.drop([
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
        "reads-unmapped",
        "reads-unpaired",
        "reads-MQ0",
        "1st-fragments",
        "last-fragments",
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
    ], axis=1)


print "[DATA] %d samples" % len(data)
print "[DATA] %s levels" % str(levels)
print "[DATA] Samples by Level %s" % sorted(statplexer.count_targets_by_class(target).items())

