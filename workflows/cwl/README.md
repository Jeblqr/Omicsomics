# CWL Workflows

Prototype Common Workflow Language (CWL) definitions for the Omicsomics ingestion/QC pipelines.

- `fastqc.cwl`: Runs FastQC against a collection of FASTQ files (BioContainers image).
- `multiqc.cwl`: Aggregates FastQC outputs into a single MultiQC HTML report.
- `omics_ingest.cwl`: Example workflow wiring `fastqc` and `multiqc` together. Accepts an array of FASTQ files and produces a combined MultiQC report.

Run with `cwltool`:

```bash
cwltool omics_ingest.cwl --fastq_files data/reads/sample1.fastq.gz --fastq_files data/reads/sample2.fastq.gz
```
