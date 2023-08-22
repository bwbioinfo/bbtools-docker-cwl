cwlVersion: v1.2
class: CommandLineTool
baseCommand: ["/bin/bash", "preprocessing_script.sh"]

inputs:
  input_forward:
    type: File
    inputBinding:
      position: 0
  input_reverse:
    type: File
    inputBinding:
      position: 1
  adapters:
    type: File
    inputBinding:
      position: 2

outputs:
  cleaned_forward:
    type: File
    outputBinding:
      glob: cleaned_forward.fastq
  cleaned_reverse:
    type: File
    outputBinding:
      glob: cleaned_reverse.fastq

requirements:
  DockerRequirement:
    dockerPull: rockylinux/rockylinux:latest

stdout: preprocessing.log
stderr: preprocessing.err