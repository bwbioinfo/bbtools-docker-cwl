#!/bin/bash
set -e

#Written by Brian Bushnell
#Last updated July 11, 2018

#Fetches and renames RefSeq.
#Be sure the taxonomy server is updated first, or run with local taxonomy data!
#To use this script outside of NERSC when using local taxonomy data,
#add "taxpath=/path/to/taxonomy_directory/" to each BBTools command


#Ensure necessary executables are in your path
#module load bbtools
module load pigz

#Fetch RefSeq
time wget -nv ftp://ftp.ncbi.nlm.nih.gov/refseq/release/complete/*genomic.fna.gz

#Rename by taxID by looking up gi numbers or accessions
time gi2taxid.sh -Xmx1g *genomic.fna.gz out=renamed.fa.gz tree=null table=null accession=null zl=6 server

#(skipped) Optionally, delete the old files
rm *genomic.fna.gz

#Sort by taxonomy.
#This makes sketching by taxa use much less memory because sketches can be written to disk as soon as they are finished.
time sortbyname.sh in=renamed.fa.gz memmult=0.33 out=sorted.fa.gz zl=8 pigz=32 taxa tree=auto gi=ignore fastawrap=1023 minlen=60 readbufferlen=2 readbuffers=1

#Make a blacklist of kmers occuring in at least 250 different species.
time sketchblacklist.sh -Xmx63g in=sorted.fa.gz prepasses=1 tree=auto taxa taxlevel=species ow out=blacklist_refseq_species_250.sketch mincount=250 k=31,24

#Generate 31 sketch files, with one sketch per species.
time bbsketch.sh -Xmx63g in=sorted.fa.gz out=taxa#.sketch mode=taxa tree=auto accession=null gi=null files=31 ow unpigz minsize=400 prefilter autosize blacklist=blacklist_refseq_species_250.sketch k=31,24 depth

#A query such as contigs.fa can now be compared to the new reference sketches like this:
#comparesketch.sh in=contigs.fa k=31,24 tree=auto taxa*.sketch blacklist=blacklist_refseq_species_250.sketch

#On NERSC systems, you can then set the default path to nt by pointing /global/projectb/sandbox/gaag/bbtools/refseq/current at the path to the new sketches.
#Then you can use the default set of nt sketches like this:
#comparesketch.sh in=contigs.fa refseq tree=auto
#That command automatically adds the default path to the sketches, the blacklist, and the correct values for K.
