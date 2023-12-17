cwlVersion: v1.0
class: CommandLineTool
baseCommand: tadpole.sh

requirements:
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: ghcr.io/bwbioinfo/bbtools-docker-cwl:latest

inputs:
  - id: in_reads
    type: File
    inputBinding:
      position: 1
      valueFrom: $(inputs.in_reads.path)

  - id: out_reads
    type: File
    outputBinding:
      glob: "$(inputs.in_reads.basename)_tadpoled.fastq"

  # Add more input parameters as needed

outputs:
  - id: tadpoled_reads
    type: File
    outputBinding:
      glob: "$(inputs.in_reads.basename)_tadpoled.fastq"

stdout: tadpole.log
stderr: tadpole.err
