from frontier import frontier
from frontier.IO.BamcheckReader import BamcheckReader
from frontier.IO.AQCReader import AQCReader
from problem_def import CLASSES, NO_VARIANCE, RAW_COUNTS

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
        lambda x, f, i: x*100,
    "percent-bases-mapped":
        lambda x,f,i: (f.get("bases-mapped-(cigar)", i) / f.get("total-length", i)) * 100.0,
    "percent-reads-qc-fail":
        lambda x,f,i: (f.get("reads-MQ0", i) / f.get("raw-total-sequences", i)) * 100.0,
    "quality-dropoff-fwd-mean-runmed-decline-range":
        lambda x,f,i: (f.get("quality-dropoff-fwd-mean-runmed-decline-high-value", i) - f.get("quality-dropoff-fwd-mean-runmed-decline-low-value", i)),
    "quality-dropoff-rev-mean-runmed-decline-range":
        lambda x,f,i: (f.get("quality-dropoff-rev-mean-runmed-decline-high-value", i) - f.get("quality-dropoff-rev-mean-runmed-decline-low-value", i))
}, add_unknown=True)

data = data.exclude([
    "quality-dropoff-fwd-mean-runmed-decline-high-value",
    "quality-dropoff-fwd-mean-runmed-decline-low-value",

    "quality-dropoff-rev-mean-runmed-decline-high-value",
    "quality-dropoff-rev-mean-runmed-decline-low-value",
])

print "[DATA] %d samples" % len(data)
print "[DATA] %s levels" % str(levels)
print "[DATA] Samples by Level %s" % sorted(statplexer.count_targets_by_class(target).items())

