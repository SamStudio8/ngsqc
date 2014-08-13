NGSQC Pipeline
==============

## Count total lanelets for each study

    $ ls crohns/ | grep "bamcheck\$" | wc -l
    9508
    $ ls uc/ | grep "bamcheck\$" | wc -l
    3947

## Extract list of target BAM files (currently targeting crohns)

    $ awk '$3 == "crohns" { print $2 }' crohns-uc-table-a.2013dec25.manual_qc_update.txt > crohns.samples.txt
    $ sed -i 's/$/.bam.goldilocks.bam/' crohns.samples.txt
    $ wc -l crohns.samples.txt
    9508 crohns.samples.txt # Number of lanelets

    $ sort -u crohns.samples.txt > crohns.samples.sorted.txt
    $ wc -l crohns.samples.sorted.txt
    2923 crohns.samples.sorted.txt # Number of human samples

## Extraction from iRODS

    $ bsub -G hgi -J "samirods[1-2923]%50" -q small -o /lustre/scratch113/teams/hgi/goldilocks/joblog/%J.%I.o -e /lustre/scratch113/teams/hgi/goldilocks/joblog/%J.%I.e bash -c 'file=$(awk NR==$LSB_JOBINDEX /lustre/scratch113/teams/hgi/autoqc_ml/crohns/crohns.samples.sorted.txt) && echo "samtools_irods for ${file}" && /software/solexa/bin/samtools_irods view -bh irods: ${file} 3:46000001-47000000 > $(basename ${file}).goldilocks.bam && rm $(basename ${file}).bai'

## Check list against directory of extracted data

    $ vimdiff crohns.samples.sorted.txt <(ls -1 goldilocks-3:46000001-47000000/ | grep -v "bai$")

## Merge

    $ bsub -o /lustre/scratch113/teams/hgi/goldilocks/joblog/samtools_merge.%J.o -e /lustre/scratch113/teams/hgi/goldilocks/joblog/samtools_merge.%J.e -G hgi -J "samtools_merge" -M2000 -R "select[mem>2000] rusage[mem=2000]" bash -c '/nfs/users/nfs_s/sn8/samtools/samtools merge -f -u -b /lustre/scratch113/teams/hgi/autoqc_ml/crohns/goldilocks-3:46000001-47000000.fofn /lustre/scratch113/teams/hgi/goldilocks/goldilocks.merged.bam'

## Count samples and lanelets in merged output

    $ samtools view -H goldilocks.merged.bam | grep "^RG" | wc -l
    $ samtools view -H goldilocks.merged.bam | grep "^@RG" | perl -pi 's/.*SM:(.*?)\s.*/\1/;' | sort | uniq | wc -l

