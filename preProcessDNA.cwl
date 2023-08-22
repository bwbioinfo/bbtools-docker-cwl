cwlVersion: v1.2
class: CommandLineTool
baseCommand: ["/bin/bash", "preprocessing_script.sh"]
doc: |
  This CWL workflow uses the provided preprocessing script to perform quality trimming,
  adapter removal, and data cleaning on paired-end Illumina sequencing data.

inputs:
  input_forward:
    type: File
    inputBinding:
      position: 0
    doc: "Input file containing forward reads in FASTQ format."

  input_reverse:
    type: File
    inputBinding:
      position: 1
    doc: "Input file containing reverse reads in FASTQ format."

  adapters:
    type: File
    inputBinding:
      position: 2
    doc: "File containing adapter sequences in FASTA format."

outputs:
  cleaned_forward:
    type: File
    outputBinding:
      glob: cleaned_forward.fastq
    doc: "Cleaned forward reads in FASTQ format."

  cleaned_reverse:
    type: File
    outputBinding:
      glob: cleaned_reverse.fastq
    doc: "Cleaned reverse reads in FASTQ format."

requirements:
  DockerRequirement:
    dockerPull: docker pull ghcr.io/bwbioinfo/bbtools-docker-cwl:latest

stdout: preprocessing.log
stderr: preprocessing.err
