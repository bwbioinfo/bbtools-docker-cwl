#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified May 9, 2018

Description:  Looks at accessions to see how to compress them.

Usage:  analyzeaccession.sh *accession2taxid.gz out=<output file>

Parameters:
lines=-1        If positive, stop after this many lines.

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

z="-Xmx400m"
z2="-Xms400m"
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

a_sample_mt() {
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
	local CMD="java $EA $EOOM $z -cp $CP tax.AnalyzeAccession $@"
	echo $CMD >&2
	eval $CMD
}

a_sample_mt "$@"
