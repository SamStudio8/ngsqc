# Summary Numbers Reference

## samtools stats
(see https://github.com/samtools/samtools/blob/master/stats.c for implementation)

#### raw total sequences

    stats->nreads_filtered+stats->nreads_1st+stats->nreads_2nd

This appears to include all 1st and 2nd reads as well as all filtered reads. The source note indicates it is “not counting excluded seqs (and none of the below)” but I don’t know what that means. 

#### filtered sequences

    stats->nreads_filtered
    
Number of reads that are filtered

#### sequences

    stats->nreads_1st+stats->nreads_2nd
    
Number of 1st and 2nd reads not counting filtered 

#### <strike>is paired</strike>

#### is sorted

    stats->is_sorted ? 1 : 0
    
1 if the file is sorted, 0 otherwise

#### 1st fragments

    stats->nreads_1st
    
Number of 1st fragments (e.g. forward reads in paired end reads)

#### last fragments

    stats->nreads_2nd
    
Number of last fragments (e.g. reverse reads in paired end reads)

#### reads mapped

    stats->nreads_paired_and_mapped+stats->nreads_single_mapped
    
Total number of mapped reads

#### reads mapped and paired

    stats->nreads_paired_and_mapped
    
Number of reads with paired-end flag  set and both mates mapped

#### reads unmapped

    stats->nreads_unmapped
    
Total number of unmapped reads

#### reads properly paired

    stats->nreads_properly_paired
    
Number of reads with “proper-pair” bit set (e.g. they are mapped and paired in a way which makes sense

#### <strike>reads unpaired</strike>
    
#### reads paired

    stats->nreads_paired_tech
    
Number of reads with “paired-end” flag set

#### reads duplicated

    stats->nreads_dup
    
Number of reads with PCR or optical duplicate flag set

#### reads MQ0

    stats->nreads_mq0
    
Number of reads mapped with mapping quality zero

#### reads QC failed

    stats->nreads_QCfailed
    
Number of reads flagged as QC failed

#### non-primary alignments

    stats->nreads_secondary
    
Number of reads which represent secondary alignments

#### total length

    stats->total_len
    
Number of bases represented  by all reads (ignoring clipping)

#### bases mapped

    stats->nbases_mapped
    
Number of bases represented by all mapped reads (ignoring clipping)

#### bases mapped (cigar)

    stats->nbases_mapped_cigar
    
Number of bases represented according to CIGAR strings (Match and Inserted)

#### bases trimmed

    stats->nbases_trimmed
    
Number of bases trimmed

#### bases duplicated

    stats->total_len_dup
    
Number of bases duplicated

#### mismatches

    stats->nmismatches
    
Number of mismatches (according to NM fields)

#### error rate

    stats->nbases_mapped_cigar ? (float)stats->nmismatches/stats->nbases_mapped_cigar : 0)
    
Rate of mismatch bases per bases mapped

#### average length

    (stats->nreads_1st+stats->nreads_2nd)?stats->total_len/(stats->nreads_1st+stats->nreads_2nd) : 0
    
Average read length

#### maximum length

    stats->max_len
    
Maximum read length

#### average quality

    stats->total_len?stats->sum_qual/stats->total_len:0
    
Average base quality

#### insert size average

    avg_isize
    
Average insertion size

#### insert size standard deviation

    sd_isize
    
Standard deviation of insert size

#### inward oriented pairs

    nisize_inward
    
Number of pairs oriented inwards

#### outward oriented pairs

    nisize_outward
    
Number of pairs oriented outwards

#### pairs with other orientation

    nisize_other
    
Number of pairs with other (neither inwards nor outwards) orientation

#### pairs on different chromosomes

    stats->nreads_anomalous/2
    
Number of pairs mapped to different chromosomes



## bamcheckr

### [indel-peaks](https://github.com/wtsi-hgi/seq_autoqc/blob/master/bamcheckr/R/indel-peaks.r)

#### fwd.percent.insertions.above.baseline
#### fwd.percent.insertions.below.baseline
#### rev.percent.insertions.above.baseline
#### rev.percent.insertions.below.baseline

These metrics result from fitting a baseline to the insertions vs read cycle graph, then subtracting that baseline and counting the number of insertions above the baselline and below the baseline, and dividing each by the total number of insertions to yield the percentage above and below the baseline. This is done for forward (first fragment) and reverse (last fragment) reads separately. 

#### fwd.percent.deletions.above.baseline
#### fwd.percent.deletions.below.baseline
#### rev.percent.deletions.above.baseline
#### rev.percent.deletions.below.baseline

These metrics result from fitting a baseline to the deletions vs read cycle graph, then subtracting that baseline and counting the number of deletions above the baselline and below the baseline, and dividing each by the total number of deletions to yield the percentage above and below the baseline. This is done for forward (first fragment) and reverse (last fragment) reads separately. 

N.B. the indel-peaks code apparently doesn’t output its parameters, making it hard to interpret!  Oops. We should add this (e.g. as indel.peaks.runmed.k).


### [quality-dropoff](https://github.com/wtsi-hgi/seq_autoqc/blob/master/bamcheckr/R/quality-dropoff.r)

#### quality.dropoff.fwd.high.iqr.max.contiguous.read.cycles

These represent the maximum number of contiguous read cycles in which the quality scores had an IQR (interquartile range) higher than the threshold (quality.dropoff.high.iqr.threshold). This metric is for the forward (first fragment) reads. 

#### quality.dropoff.rev.high.iqr.max.contiguous.read.cycles

The same as above but for the reverse (last fragment) reads.


#### quality.dropoff.fwd.high.iqr.start.read.cycle
#### quality.dropoff.fwd.high.iqr.end.read.cycle

These annotate the above ‘quality.dropoff.fwd.high.iqr.max.contiguous.read.cycles’ results to indicate which read cycles the contiguous block of high IQR reads started and ended at. Any two out of three of these parameters can be used to calculate the other (e.g.   quality.dropoff.fwd.high.iqr.max.contiguous.read.cycles = quality.dropoff.fwd.high.iqr.end.read.cycle - quality.dropoff.fwd.high.iqr.start.read.cycle + 1)

#### quality.dropoff.rev.high.iqr.start.read.cycle
#### quality.dropoff.rev.high.iqr.end.read.cycle

The same as above but for the reverse (last fragment) reads. 

#### quality.dropoff.fwd.mean.runmed.decline.start.read.cycle
#### quality.dropoff.fwd.mean.runmed.decline.end.read.cycle
#### quality.dropoff.fwd.mean.runmed.decline.max.contiguous.read.cycles
#### quality.dropoff.fwd.mean.runmed.decline.high.value
#### quality.dropoff.fwd.mean.runmed.decline.low.value

Similar to the high IQR test above, this test locates the longest contiguous block of read cycles in the forward reads during which the quality declines or holds steady without ever increasing. ‘*.max.contiguous.read.cycles’ gives the length of the block of contiguous read cycles, ‘*.start.read.cycle’ gives the read cycle at which the block begins, and *.end.read.cycle’ gives the read cycle at which the block ends (the following read cycle the quality should be expected to either increase or not exist). The ‘*.high.value’ gives the highest quality value in the block (which one would expect to find at the lowest read cycle) and ‘*.low.value’ gives the lowest quality value in the block.

***Probably we should have included high-low value as ‘.quality.range’ so that there would be a single metric to filter on. Oops.



#### quality.dropoff.rev.mean.runmed.decline.start.read.cycle
#### quality.dropoff.rev.mean.runmed.decline.end.read.cycle
#### quality.dropoff.rev.mean.runmed.decline.max.contiguous.read.cycles
#### quality.dropoff.rev.mean.runmed.decline.high.value
#### quality.dropoff.rev.mean.runmed.decline.low.value

The same as above but for the reverse reads. 


#### quality.dropoff.high.iqr.threshold
#### quality.dropoff.runmed.k
#### quality.dropoff.ignore.edge.cycles
These three are the parameters which were used in the quality-dropoff calculations: 


### [base-content-deviation](https://github.com/wtsi-hgi/seq_autoqc/blob/master/bamcheckr/R/base-content-deviation.r)

All of the following result from fitting a baseline to each of the base percent vs read cycle plots.

#### A.percent.mean.above.baseline
#### C.percent.mean.above.baseline
#### G.percent.mean.above.baseline
#### T.percent.mean.above.baseline

Taking all read cycles at which the N (e.g. A, C, G, or T) base percent was above the baseline and subtracting the baseline from it, these represent the mean of those values. 

#### A.percent.mean.below.baseline
#### C.percent.mean.below.baseline
#### G.percent.mean.below.baseline
#### T.percent.mean.below.baseline

As above but for the absolute value of values below the baseline. 

#### A.percent.max.above.baseline
#### C.percent.max.above.baseline
#### G.percent.max.above.baseline
#### T.percent.max.above.baseline

The base percent value at the read cycle with the highest deviation above the baseline. 

#### A.percent.max.below.baseline
#### C.percent.max.below.baseline
#### G.percent.max.below.baseline
#### T.percent.max.below.baseline

The absolute value of the base percent value at the read cycle where it is fuathest below the baseline. 

#### A.percent.max.baseline.deviation

The higher of ‘A.percent.max.above.baseline’ or ‘A.percent.max.below.baseline’

#### C.percent.max.baseline.deviation
#### G.percent.max.baseline.deviation
#### T.percent.max.baseline.deviation

Likewise for the other bases.

#### A.percent.total.mean.baseline.deviation

The sum of ‘A.percent.mean.above.baseline’ and ‘A.percent.mean.below.baseline’

#### C.percent.total.mean.baseline.deviation
#### G.percent.total.mean.baseline.deviation
#### T.percent.total.mean.baseline.deviation

Likewise for the other bases.

* Again, the base-content-deviation detectors are missing the parameters used in fitting.

