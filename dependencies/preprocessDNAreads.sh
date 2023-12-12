#!/bin/bash
# Brian Bushnell's read preprocessing pipeline

set -o errexit
set -o nounset

READS=( $@ )

BBMAP=/opt/bioinformatics/bbmap
READDIR=$(dirname ${READS[0]} )
IN_STR="in=${READS[0]}"
if [ ${#READS[@]} -gt 1 ] ; then IN_STR="$IN_STR in2=${READS[1]}" ; fi

# Format conversion
INTERLEAVED="${READS[0]%_R1*}.interleaved"
INTERLEAVED=$(basename $INTERLEAVED)
$BBMAP/reformat.sh $IN_STR out=${INTERLEAVED}.fastq.gz

# Adapter trimming
CLEANED="${INTERLEAVED}.cleaned"
$BBMAP/bbduk.sh in="${INTERLEAVED}.fastq.gz" out="${CLEANED}.fastq.gz" ref=$BBMAP/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo ftm=5 trimpolya=10 entropy=0.1 interleaved=t ordered=t usejni

# Contaminant filtering
DECONTAMINATED="${READS[0]%_R1*}.decontaminated"
DECONTAMINATED=$(basename $DECONTAMINATED)
CONTAMINANTS="${READS[0]%_R1*}.contaminants"
CONTAMINANTS=$(basename $CONTAMINANTS)
STATS="${READS[0]%_R1*}.contaminant.stats.txt"
STATS=$(basename $STATS)
$BBMAP/bbduk.sh in="${CLEANED}.fastq.gz" out="${DECONTAMINATED}.fastq.gz" outm="${CONTAMINANTS}.fastq.gz" ref=$BBMAP/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=$STATS interleaved=t ordered=t usejni

# Error correction
CORRECTED="${READS[0]%_R1*}.corrected"
CORRECTED=$(basename $CORRECTED)
$BBMAP/tadpole.sh in="${DECONTAMINATED}.fastq.gz" out="${CORRECTED}.fastq.gz" mode=correct k=50

## Paired-read merging
#MERGED="${READS[0]%_R1*}.merged"
#UNMERGED="${READS[0]%_R1*}.unmergeable"
#IHIST="${READS[0]%_R1*}.ihist.txt"
#OUTC="${READS[0]%_R1*}.kmer-cardinality.txt"
#$BBMAP/bbmerge-auto.sh in="${CORRECTED}.fastq.gz" out="${MERGED}.fastq.gz" outu="${UNMERGED}.fastq.gz" ihist="$IHIST" outc=$OUTC ecct extend2=20 iterations=5 usejni=t 

# Split corrected read pairs for bowtie
$BBMAP/reformat.sh in="${CORRECTED}.fastq.gz" out="${CORRECTED}_R#.fastq.gz"

# Normalize corrected reads
NORMALIZED="${READS[0]%_R1*}.normalized"
NORMALIZED=$(basename $NORMALIZED)
KHIST="${READS[0]%_R1*}.khist.txt"
KHIST=$(basename $KHIST)
$BBMAP/bbnorm.sh in="${CORRECTED}.fastq.gz" out="${NORMALIZED}.fastq.gz" target=100 min=5 khist=$KHIST

