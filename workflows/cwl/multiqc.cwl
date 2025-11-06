cwlVersion: v1.2
class: CommandLineTool

baseCommand: [multiqc]

requirements:
  DockerRequirement:
    dockerPull: ewels/multiqc:1.16

inputs:
  reports:
    type:
      type: array
      items: File
    inputBinding:
      position: 1
  output_dir:
    type: string
    inputBinding:
      prefix: --outdir

outputs:
  report:
    type: File
    outputBinding:
      glob: $(inputs.output_dir + '/multiqc_report.html')
