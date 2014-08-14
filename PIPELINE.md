NGSQC Pipeline
==============

## Count totals for each study

#### Lanelets

    $ ls crohns/ | grep "bamcheck\$" | wc -l
    9508
    $ ls uc/ | grep "bamcheck\$" | wc -l
    3947
    ==============
    13455 lanelets

#### Samples (as process below)

    $ wc -l crohns.samples.sorted.txt
    2932
    $ wc -l uc.samples.sorted.txt
    1992
    ============
    4924 samples

## Extract list of target BAM files (currently targeting crohns)

    $ awk '$3 == "crohns" { print $2 }' crohns-uc-table-a.2013dec25.manual_qc_update.txt > crohns.samples.txt
    $ sed -i 's/$/.bam.goldilocks.bam/' crohns.samples.txt
    $ wc -l crohns.samples.txt
    9508 crohns.samples.txt

    $ sort -u crohns.samples.txt > crohns.samples.sorted.txt
    $ wc -l crohns.samples.sorted.txt
    2923 crohns.samples.sorted.txt

## Query iRODS for chosen samples

#### Acquire kerberos token

    $ kinit

#### List iRODS contents

    $ /software/irods/icommands/bin/ils /humgen/projects/crohns/20130909 | grep -v "bai$" > crohns20130909.ils
    $ /software/irods/icommands/bin/ils /humgen/projects/crohns/20131023 | grep -v "bai$" > crohns20131023.ils

#### Remove header and characters trailing sample name, then sort

    $ sed '1d; s/\s*\(.*\)\..*$/\1/' crohns20130909.ils > crohns20130909-samples.ils
    $ sed '1d; s/\s*\(.*\)\..*$/\1/' crohns20131023.ils > crohns20131023-samples.ils
    $ sort -u crohns20131023-samples.ils > crohns20130909-samples-sorted.ils
    $ sort -u crohns20131023-samples.ils > crohns20131023-samples-sorted.ils

    $ wc -l crohns20131023-samples-sorted.ils
    2697 crohns20131023-samples-sorted.ils
    $ wc -l crohns20130909-samples-sorted.ils
    1817 crohns20130909-samples-sorted.ils
    ============
    4514 samples


#### Confirm samples are unique to one release

    $ comm -12 crohns20131023-samples-sorted.ils crohns20130909-samples-sorted.ils | wc -l
    0

#### Discover that each release happens to relate to one study after all...

    $ comm -12 crohns20130909-samples-sorted.ils crohns.samples.sorted.name-only.txt | wc -l
    0
    $ comm -12 crohns20131023-samples-sorted.ils uc.samples.sorted.name-only.txt | wc -l
    0

#### Identify crohns-only BAM files missing from iRODS storage

    $ comm -13 crohns20131023-samples-sorted.ils crohns.samples.sorted.name-only.txt | wc -l
    226

## Extraction from iRODS

    $ bsub -G hgi -J "samirods[1-2923]%50" -q small -o /lustre/scratch113/teams/hgi/goldilocks/joblog/%J.%I.o -e /lustre/scratch113/teams/hgi/goldilocks/joblog/%J.%I.e bash -c 'file=$(awk NR==$LSB_JOBINDEX /lustre/scratch113/teams/hgi/autoqc_ml/crohns/crohns.samples.sorted.irods.txt) && echo "samtools_irods for ${file}" && /software/solexa/bin/samtools_irods view -bh irods: ${file} 3:46000001-47000000 > $(basename ${file}).goldilocks.bam && rm $(basename ${file}).bai'

## Check list against directory of extracted data

    $ vimdiff crohns.samples.sorted.txt <(ls -1 goldilocks-3:46000001-47000000/ | grep -v "bai$")

## Merge

    $ bsub -o /lustre/scratch113/teams/hgi/goldilocks/joblog/samtools_merge.%J.o -e /lustre/scratch113/teams/hgi/goldilocks/joblog/samtools_merge.%J.e -G hgi -J "samtools_merge" -M2000 -R "select[mem>2000] rusage[mem=2000]" bash -c '/nfs/users/nfs_s/sn8/samtools/samtools merge -f -u -b /lustre/scratch113/teams/hgi/autoqc_ml/crohns/goldilocks-3:46000001-47000000.fofn /lustre/scratch113/teams/hgi/goldilocks/goldilocks.merged.bam'

## Count samples and lanelets in merged output

    $ samtools view -H goldilocks.merged.bam | grep "^RG" | wc -l
    $ samtools view -H goldilocks.merged.bam | grep "^@RG" | perl -pi 's/.*SM:(.*?)\s.*/\1/;' | sort | uniq | wc -l

