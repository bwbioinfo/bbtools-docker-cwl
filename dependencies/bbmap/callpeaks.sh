#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified May 31, 2018

Description:  Calls peaks from a 2-column (x, y) tab-delimited histogram.

Usage:        callpeaks.sh in=<histogram file> out=<output file>

Peak-calling parameters:
in=<file>           'in=stdin.fq' will pipe from standard in.
out=<file>          Write the peaks to this file.  Default is stdout.
minHeight=2         (h) Ignore peaks shorter than this.
minVolume=5         (v) Ignore peaks with less area than this.
minWidth=3          (w) Ignore peaks narrower than this.
minPeak=2           (minp) Ignore peaks with an X-value below this. 
                    Useful when low-count kmers are filtered).
maxPeak=BIG         (maxp) Ignore peaks with an X-value above this.
maxPeakCount=10     (maxpc) Print up to this many peaks (prioritizing height).
countColumn=1       (col) For multi-column input, this column, zero-based, 
                    contains the counts.
ploidy=-1           Specify ploidy; otherwise it will be autodetected.
logscale=f          Transform to log-scale prior to peak-calling.  Useful
                    for kmer-frequency histograms.

Smoothing parameters:
smoothradius=0      Integer radius of triangle filter.  Set above zero to 
                    smooth data prior to peak-calling.  Higher values are 
                    smoother.
smoothprogressive=f Set to true to widen the filter as the x-coordinate
                    increases.  Useful for kmer-frequency histograms.
maxradius=10        Maximum radius of progressive smoothing function.
progressivemult=2   Increment radius each time depth increases by this factor.

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

z="-Xmx200m"
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
}
calcXmx "$@"

stats() {
	if [[ $SHIFTER_RUNTIME == 1 ]]; then
		#Ignore NERSC_HOST
		shifter=1
	elif [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_144_64bit
		module load pigz
	elif [[ $NERSC_HOST == denovo ]]; then
		module unload java
		module load java/1.8.0_144
		module load pigz
	elif [[ $NERSC_HOST == cori ]]; then
		module use /global/common/software/m342/nersc-builds/denovo/Modules/jgi
		module use /global/common/software/m342/nersc-builds/denovo/Modules/usg
		module unload java
		module load java/1.8.0_144
		module load pigz
	fi
	local CMD="java $EA $EOOM -Xmx120m -cp $CP jgi.CallPeaks $@"
#	echo $CMD >&2
	eval $CMD
}

stats "$@"
