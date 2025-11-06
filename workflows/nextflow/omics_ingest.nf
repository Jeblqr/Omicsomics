nextflow.enable.dsl=2

params.reads = params.reads ?: 'data/reads/*.fastq.gz'
params.outdir = params.outdir ?: 'results/qc'

process FASTQC {
  tag "${sample_id}"
  publishDir "${params.outdir}/qc", mode: 'copy'

  input:
  tuple val(sample_id), path(reads)

  output:
  path("${sample_id}_fastqc.zip"), emit: qc_zip
  path("${sample_id}_fastqc.html"), emit: qc_report

  script:
  """
  fastqc --threads ${task.cpus} --outdir . ${reads}
  """
}

process MULTIQC {
  publishDir "${params.outdir}/multiqc", mode: 'copy'

  input:
  path qc_reports from FASTQC.out.qc_zip.collect()
  path qc_htmls from FASTQC.out.qc_report.collect()

  output:
  path 'multiqc_report.html'
  path 'multiqc_data'

  script:
  """
  multiqc . --outdir .
  """
}

workflow {
  Channel
    .fromFilePairs(params.reads, flat: true)
    .map { sample_id, reads -> tuple(sample_id.tokenize('/').last().replaceAll(/\W+/, '_'), reads) }
    | FASTQC

  MULTIQC(FASTQC.out)
}
