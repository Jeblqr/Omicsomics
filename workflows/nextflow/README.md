# Nextflow Pipelines

This directory contains prototype Nextflow pipelines used by the Omicsomics platform.

- `omics_ingest.nf`: Minimal QC pipeline that runs FastQC on paired FASTQ files and aggregates results with MultiQC. Input glob is controlled via `--reads` and results are stored in `--outdir` (default `results/qc`).

## Running locally

```bash
nextflow run omics_ingest.nf --reads "data/reads/*_R{1,2}.fastq.gz" --outdir results/qc
```

Ensure `fastqc` and `multiqc` are available in the execution environment, or update the process containers to use existing bioinformatics images.
