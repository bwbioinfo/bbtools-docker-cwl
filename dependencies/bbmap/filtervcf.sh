#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified February 5, 2018

Description:  Filters VCF files by position or other attributes.
Filtering by optional fields (such as allele frequency) require VCF files
generated by CallVariants.

Usage:  filtervcf.sh in=<file> out=<file>

I/O parameters:
in=<file>       Input VCF.
out=<file>      Output VCF.
ref=<file>      Reference fasta (optional).
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.
bgzip=f         Use bgzip for gzip compression.

Position-filtering parameters:
minpos=         Ignore variants not overlapping this range.
maxpos=         Ignore variants not overlapping this range.
contigs=        Comma-delimited list of contig names to include. These 
                should have no spaces, or underscores instead of spaces.
invert=f        Invert position filters.

Type-filtering parameters:
sub=t           Keep substitutions.
del=t           Keep deletions.
ins=t           Keep insertions.

Variant-quality filtering parameters:
minreads=0              Ignore variants seen in fewer reads.
minqualitymax=0         Ignore variants with lower max base quality.
minedistmax=0           Ignore variants with lower max distance from read ends.
minmapqmax=0            Ignore variants with lower max mapq.
minidmax=0              Ignore variants with lower max read identity.
minpairingrate=0.0      Ignore variants with lower pairing rate.
minstrandratio=0.0      Ignore variants with lower plus/minus strand ratio.
minquality=0.0          Ignore variants with lower average base quality.
minedist=0.0            Ignore variants with lower average distance from ends.
minavgmapq=0.0          Ignore variants with lower average mapq.
minallelefraction=0.0   Ignore variants with lower allele fraction.  This
                        should be adjusted for high ploidies.
minid=0                 Ignore variants with lower average read identity.
minscore=0.0            Ignore variants with lower Phred-scaled score.
clearfilters            Reset all variant filters to zero.

There are additionally max filters for score, quality, mapq, allelefraction,
and identity.

Java Parameters:
-Xmx            This will be passed to Java to set memory usage, overriding the program's automatic memory detection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory.
-eoom           This flag will cause the process to exit if an out-of-memory exception occurs.  Requires Java 8u92+.
-da             Disable assertions.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
"
}

pushd . > /dev/null
DIR="${BASH_SOURCE[0]}"
while [ -h "$DIR" ]; do
  cd "$(dirname "$DIR")"
  DIR="$(readlink "$(basename "$DIR")")"
done
cd "$(dirname "$DIR")"
DIR="$(pwd)/"
popd > /dev/null

#DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx4g"
z2="-Xms4g"
EA="-ea"
EOOM=""
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	parseXmx "$@"
	if [[ $set == 1 ]]; then
		return
	fi
	freeRam 4000m 42
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

filtervcf() {
	if [[ $SHIFTER_RUNTIME == 1 ]]; then
		#Ignore NERSC_HOST
		shifter=1
	elif [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_144_64bit
		module load samtools/1.4
		module load pigz
	elif [[ $NERSC_HOST == denovo ]]; then
		module unload java
		module load java/1.8.0_144
		module load PrgEnv-gnu/7.1
		module load samtools/1.4
		module load pigz
	elif [[ $NERSC_HOST == cori ]]; then
		module use /global/common/software/m342/nersc-builds/denovo/Modules/jgi
		module use /global/common/software/m342/nersc-builds/denovo/Modules/usg
		module unload java
		module load java/1.8.0_144
		module load PrgEnv-gnu/7.1
		module load samtools/1.4
		module load pigz
	fi
	local CMD="java $EA $EOOM $z $z2 -cp $CP var2.FilterVCF $@"
	echo $CMD >&2
	eval $CMD
}

filtervcf "$@"
