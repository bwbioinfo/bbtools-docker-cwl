#!/bin/bash
# Brian Bushnell's read preprocessing pipeline
# 2020-04-03 Added minlength filter and checks for proper pairing to final reformat step
set -o errexit
set -o nounset

READS=( $@ )

BBMAP=/opt/bioinformatics/bbmap
READDIR=$(dirname ${READS[0]} )
IN_STR="in=${READS[0]}"
if [ ${#READS[@]} -gt 1 ] ; then IN_STR="$IN_STR in2=${READS[1]}" ; fi
STEM=${READS[0]%.f*q*}
STEM=${STEM%_R1}
STEM=${STEM%_1}
STEM=${STEM%.1}
STEM=$(basename $STEM)
# Format conversion
INTERLEAVED="${STEM}.interleaved"
$BBMAP/reformat.sh $IN_STR out=${INTERLEAVED}.fastq.gz

# Adapter trimming
CLEANED="${STEM}.cleaned"
$BBMAP/bbduk.sh in="${INTERLEAVED}.fastq.gz" out="${CLEANED}.fastq.gz" ref=$BBMAP/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo ftm=5 trimpolya=0 entropy=0.1 usejni

# Contaminant filtering
DECONTAMINATED="${STEM}.decontaminated"
CONTAMINANTS="${STEM}.contaminants"
STATS="${STEM}.contaminant.stats.txt"
RRNA_DB=/opt/bioinformatics/sortmerna/rRNA_databases
$BBMAP/bbduk.sh in="${CLEANED}.fastq.gz" out="${DECONTAMINATED}.fastq.gz" outm="${CONTAMINANTS}.fastq.gz" ref=$BBMAP/resources/phix174_ill.ref.fa.gz,$RRNA_DB/silva-euk-18s-id95.fasta,$RRNA_DB/silva-euk-28s-id98.fasta minkmerfraction=0.75 k=31 stats=$STATS usejni

# Split decontaminated read pairs for mapping

# Clean up
rm ${INTERLEAVED}.fastq.gz
rm ${CLEANED}.fastq.gz

if [ ${#READS[@]} -gt 1 ] 
	then $BBMAP/reformat.sh in="${DECONTAMINATED}.fastq.gz" out="${DECONTAMINATED}_R#.fastq.gz" verifyinterleaved=t minlength=50 requirebothbad=f outsingle="${DECONTAMINATED}.singletons.fastq.gz"
#	rm ${DECONTAMINATED}.fastq.gz
fi


# Normalize decontaminated reads
NORMALIZED="${STEM}.normalized"
KHIST="${STEM}.khist.txt"
$BBMAP/bbnorm.sh in="${DECONTAMINATED}.fastq.gz" out="${NORMALIZED}.fastq.gz" target=100 min=5 khist=$KHIST


if [ ${#READS[@]} -gt 1 ] 
	then $BBMAP/reformat.sh in="${NORMALIZED}.fastq.gz" out="${NORMALIZED}_R#.fastq.gz" verifyinterleaved=t minlength=50 requirebothbad=f outsingle="${NORMALIZED}.singletons.fastq.gz"
#	rm ${DECONTAMINATED}.fastq.gz
fi


