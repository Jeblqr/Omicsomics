cwlVersion: v1.2
class: CommandLineTool

baseCommand: [fastqc]

requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: biocontainers/fastqc:v0.11.9_cv8

inputs:
  fastq_files:
    type:
      type: array
      items: File
    inputBinding:
      position: 1
      itemSeparator: ''

outputs:
  reports:
    type:
      type: array
      items: File
    outputBinding:
      glob: $(inputs.fastq_files.map(f => f.basename.replace(/\.fastq(?:\.gz)?$/, '') + '_fastqc.zip'))
