cwlVersion: v1.0
class: CommandLineTool
baseCommand: reformat.sh

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
      glob: "$(inputs.in_reads.basename)_reformatted.fastq"

  # Add more input parameters as needed

outputs:
  - id: reformatted_reads
    type: File
    outputBinding:
      glob: "$(inputs.in_reads.basename)_reformatted.fastq"

stdout: reformat.log
stderr: reformat.err
