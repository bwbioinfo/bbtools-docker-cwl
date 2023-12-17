cwlVersion: v1.0
class: CommandLineTool
baseCommand: bbnorm.sh


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
      glob: "$(inputs.in_reads.basename)_normalized.fastq"

  - id: target_depth
    type: int
    inputBinding:
      position: 2
      prefix: target=

  - id: mindepth
    type: int?
    inputBinding:
      position: 3
      prefix: mindepth=

  # Add more input parameters as needed

outputs:
  - id: normalized_reads
    type: File
    outputBinding:
      glob: "$(inputs.in_reads.basename)_normalized.fastq"

stdout: bbnorm.log
stderr: bbnorm.err
