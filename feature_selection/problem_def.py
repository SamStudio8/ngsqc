
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

NO_VARIANCE = [
    "bases-trimmed",
    "filtered-sequences",
    "is-paired",
    "is-sorted",
    "maximum-length",
    "non-primary-alignments",
    "quality-dropoff-high-iqr-threshold",
    "quality-dropoff-ignore-edge-cycles",
    "quality-dropoff-runmed-k parameter",
]

RAW_COUNTS = [
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
]
