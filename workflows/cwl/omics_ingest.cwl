cwlVersion: v1.2
class: Workflow

requirements:
  InlineJavascriptRequirement: {}

inputs:
  fastq_files:
    type:
      type: array
      items: File
  output_dir:
    type: string
    default: results/qc

outputs:
  multiqc_report:
    type: File
    outputSource: multiqc/report

steps:
  fastqc:
    run: fastqc.cwl
    in:
      fastq_files: fastq_files
    out: [reports]

  multiqc:
    run: multiqc.cwl
    in:
      reports: fastqc/reports
      output_dir: output_dir
    out: [report]
