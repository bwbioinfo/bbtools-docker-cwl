cwlVersion: v1.0
class: CommandLineTool
baseCommand: bbduk.sh

requirements:
  - class: InlineJavascriptRequirement

inputs:
  - id: in_reads
    type: File
    inputBinding:
      position: 1
      valueFrom: $(inputs.in_reads.path)

  - id: out_reads
    type: File
    outputBinding:
      glob: "$(inputs.in_reads.basename)_filtered.fastq"

  # Add more input parameters as needed

outputs:
  - id: filtered_reads
    type: File
    outputBinding:
      glob: "$(inputs.in_reads.basename)_filtered.fastq"

stdout: bbduk.log
stderr: bbduk.err
