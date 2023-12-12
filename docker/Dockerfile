FROM ubuntu:20.04

ARG INSTALL_PATH=/opt
ARG DEPPATH=./Dependencies

LABEL maintainer="csfg-helpdesk@concordia.ca"
LABEL version="1.00"
LABEL description="This is the custom Docker Image for RNA-Seq and DNA-Seq Preprocessing, for use by CSFG."

# Disable Prompt During Packages Installation
ARG DEBIAN_FRONTEND=noninteractive

# Update Ubuntu Software repository
RUN apt update
RUN apt-get install -y openjdk-11-jdk

### BBMAP ###
# copy BBMAP program
RUN  mkdir ${INSTALL_PATH}/bioinformatics
RUN  mkdir ${INSTALL_PATH}/sortmerna
COPY ${DEPPATH}/bbmap ${INSTALL_PATH}/bioinformatics/bbmap
COPY ${DEPPATH}/rRNA_databases ${INSTALL_PATH}/bioinformatics/sortmerna/rRNA_databases

### DNA-SEQ ###
# copy preprocess DNA script
COPY $DEPPATH/preprocessDNAreads.sh $INSTALL_PATH/preprocessDNAreads.sh

### RNA-SEQ ###
# copy preprocess RNA script
COPY $DEPPATH/preprocessRNAreads.sh $INSTALL_PATH/preprocessRNAreads.sh

ENV  PATH="$INSTALL_PATH/:$PATH"
ENV  PATH="$INSTALL_PATH/bioinformatics/bbmap/:$PATH"

RUN  chmod -R a+rxX $INSTALL_PATH/*