#!/bin/bash
#Run this on gpweb34
#module unload oracle-jdk
#module load oracle-jdk/1.8_64bit

nohup /global/projectb/sandbox/gaag/bbtools/jgi-bbtools/taxserver.sh -da -Xmx9g port=3071 verbose tree=auto sketchonly index domain=https://nt-sketch.jgi-psf.org killcode=xxxxx oldcode=xxxxx oldaddress=https://nt-sketch.jgi-psf.org/kill/ nt k=31,24 1>>ntlog18.o 2>&1 &

#simple mode, for testing:
#/global/projectb/sandbox/gaag/bbtools/jgi-bbtools/taxserver.sh -ea -Xmx9g port=3071 verbose tree=auto sketchonly nt k=31,24 index=f
