"""Pipeline template service for common workflows."""

from typing import List, Dict, Any
from enum import Enum


class TemplateCategory(str, Enum):
    """Pipeline template categories."""

    TRANSCRIPTOMICS = "transcriptomics"
    GENOMICS = "genomics"
    PROTEOMICS = "proteomics"
    METABOLOMICS = "metabolomics"
    EPIGENOMICS = "epigenomics"
    SINGLECELL = "singlecell"
    MULTIOMICS = "multiomics"
    GWAS = "gwas"
    QUALITY_CONTROL = "quality_control"
    ENRICHMENT = "enrichment"


# Common pipeline templates available to all users
COMMON_PIPELINES = [
    {
        "id": "rna-seq-basic",
        "name": "RNA-seq Basic Analysis",
        "description": "Standard RNA-seq pipeline: QC, alignment, quantification, and DEG analysis",
        "category": "transcriptomics",
        "steps": [
            {
                "name": "Quality Control",
                "tool": "fastqc",
                "version": "0.11.9",
                "parameters": {"threads": 4},
            },
            {
                "name": "Adapter Trimming",
                "tool": "trim_galore",
                "version": "0.6.7",
                "parameters": {"quality": 20, "length": 20},
            },
            {
                "name": "Read Alignment",
                "tool": "star",
                "version": "2.7.10",
                "parameters": {"threads": 8, "outSAMtype": "BAM SortedByCoordinate"},
            },
            {
                "name": "Feature Counting",
                "tool": "featureCounts",
                "version": "2.0.1",
                "parameters": {"threads": 4, "isPairedEnd": True},
            },
            {
                "name": "Differential Expression",
                "tool": "deseq2",
                "version": "1.34.0",
                "parameters": {"alpha": 0.05, "lfcThreshold": 0},
            },
        ],
        "parameters": {
            "genome_index": {
                "type": "path",
                "required": True,
                "description": "STAR genome index directory",
            },
            "annotation_gtf": {
                "type": "path",
                "required": True,
                "description": "Gene annotation GTF file",
            },
            "design_formula": {
                "type": "string",
                "default": "~ condition",
                "description": "DESeq2 design formula",
            },
            "alpha": {
                "type": "float",
                "default": 0.05,
                "description": "Significance threshold",
            },
        },
        "inputs": ["fastq_r1", "fastq_r2", "genome_index", "annotation_gtf"],
        "outputs": ["counts_matrix", "deg_results", "qc_report", "normalized_counts"],
    },
    {
        "id": "variant-calling",
        "name": "Variant Calling Pipeline",
        "description": "WGS/WES variant calling: alignment, preprocessing, variant calling, and annotation",
        "category": "genomics",
        "steps": [
            "bwa_align",
            "mark_duplicates",
            "base_recalibration",
            "haplotypecaller",
            "variant_annotation",
        ],
        "inputs": ["fastq_r1", "fastq_r2", "reference_genome", "known_variants"],
        "outputs": ["vcf_file", "annotated_variants", "qc_metrics"],
    },
    {
        "id": "chip-seq",
        "name": "ChIP-seq Analysis",
        "description": "ChIP-seq pipeline: QC, alignment, peak calling, and motif analysis",
        "category": "epigenomics",
        "steps": [
            "fastqc",
            "bowtie2_align",
            "remove_duplicates",
            "macs2_peak_calling",
            "homer_motif_analysis",
        ],
        "inputs": ["fastq", "genome_index", "control_fastq"],
        "outputs": ["peaks_bed", "bigwig", "motif_results"],
    },
    {
        "id": "single-cell-rna",
        "name": "Single Cell RNA-seq",
        "description": "scRNA-seq analysis: QC, normalization, clustering, and marker identification",
        "category": "singlecell",
        "steps": [
            "cellranger_count",
            "seurat_qc",
            "normalization",
            "clustering",
            "marker_detection",
        ],
        "inputs": ["fastq_dir", "transcriptome_ref"],
        "outputs": ["cell_by_gene_matrix", "clusters", "markers", "umap_plot"],
    },
    {
        "id": "proteomics-label-free",
        "name": "Label-free Proteomics Quantification",
        "description": "LC-MS/MS proteomics: identification, quantification, and differential analysis",
        "category": "proteomics",
        "steps": [
            {
                "name": "Database Search",
                "tool": "maxquant",
                "version": "2.0.3",
                "parameters": {"fdr": 0.01, "minPeptides": 2},
            },
            {
                "name": "Peptide Filtering",
                "tool": "peptide_filter",
                "version": "1.0.0",
                "parameters": {"minLength": 7, "removeContaminants": True},
            },
            {
                "name": "Protein Inference",
                "tool": "protein_prophet",
                "version": "5.2.0",
                "parameters": {"probability": 0.95},
            },
            {
                "name": "Differential Analysis",
                "tool": "limma",
                "version": "3.50.0",
                "parameters": {"fdr": 0.05, "logFC": 1},
            },
        ],
        "parameters": {
            "fasta_database": {
                "type": "path",
                "required": True,
                "description": "Protein FASTA database",
            },
            "fdr_threshold": {
                "type": "float",
                "default": 0.01,
                "description": "False discovery rate",
            },
            "min_peptides": {
                "type": "integer",
                "default": 2,
                "description": "Minimum peptides per protein",
            },
            "imputation_method": {
                "type": "string",
                "default": "knn",
                "description": "Missing value imputation",
            },
        },
        "inputs": ["raw_files", "fasta_database", "sample_metadata"],
        "outputs": ["protein_groups", "de_proteins", "qc_report", "intensity_matrix"],
    },
    {
        "id": "proteomics-qc",
        "name": "Proteomics Quality Control",
        "description": "Comprehensive QC for proteomics data: missing values, CV, PCA, normalization",
        "category": "quality_control",
        "steps": [
            {
                "name": "Data Loading",
                "tool": "proteus",
                "version": "1.0.0",
                "parameters": {},
            },
            {
                "name": "Missing Value Analysis",
                "tool": "qc_missing",
                "version": "1.0.0",
                "parameters": {"threshold": 0.5},
            },
            {
                "name": "CV Analysis",
                "tool": "qc_cv",
                "version": "1.0.0",
                "parameters": {"max_cv": 0.3},
            },
            {
                "name": "Normalization",
                "tool": "normalize",
                "version": "1.0.0",
                "parameters": {"method": "median"},
            },
            {
                "name": "PCA",
                "tool": "pca",
                "version": "1.0.0",
                "parameters": {"n_components": 2},
            },
        ],
        "parameters": {
            "missing_threshold": {
                "type": "float",
                "default": 0.5,
                "description": "Max missing ratio",
            },
            "cv_threshold": {"type": "float", "default": 0.3, "description": "Max CV"},
            "normalization": {
                "type": "string",
                "default": "median",
                "description": "Normalization method",
            },
        },
        "inputs": ["protein_groups", "sample_metadata"],
        "outputs": ["qc_report", "filtered_data", "normalized_data", "plots"],
    },
    {
        "id": "go-enrichment",
        "name": "GO Term Enrichment Analysis",
        "description": "Gene Ontology enrichment for differentially expressed genes/proteins",
        "category": "enrichment",
        "steps": [
            {
                "name": "ID Mapping",
                "tool": "biomart",
                "version": "2.50.0",
                "parameters": {"dataset": "hsapiens_gene_ensembl"},
            },
            {
                "name": "GO Enrichment",
                "tool": "topgo",
                "version": "2.46.0",
                "parameters": {"ontology": "BP", "algorithm": "elim"},
            },
            {
                "name": "KEGG Pathways",
                "tool": "clusterprofiler",
                "version": "4.2.0",
                "parameters": {"organism": "hsa", "pvalueCutoff": 0.05},
            },
        ],
        "parameters": {
            "organism": {
                "type": "string",
                "default": "human",
                "description": "Organism",
            },
            "ontology": {
                "type": "string",
                "default": "BP",
                "description": "GO ontology (BP/MF/CC)",
            },
            "pvalue_cutoff": {
                "type": "float",
                "default": 0.05,
                "description": "P-value threshold",
            },
        },
        "inputs": ["gene_list", "background_genes"],
        "outputs": ["enrichment_results", "pathway_plots", "network"],
    },
    {
        "id": "metabolomics-untargeted",
        "name": "Untargeted Metabolomics",
        "description": "Metabolomics analysis: feature detection, annotation, and pathway analysis",
        "category": "metabolomics",
        "steps": [
            "xcms_peak_detection",
            "camera_annotation",
            "metaboanalyst",
            "pathway_enrichment",
        ],
        "inputs": ["mzml_files", "sample_metadata"],
        "outputs": ["feature_table", "identified_compounds", "pathway_results"],
    },
    {
        "id": "gwas",
        "name": "Genome-Wide Association Study",
        "description": "GWAS pipeline: QC, imputation, association testing, and visualization",
        "category": "gwas",
        "steps": [
            "plink_qc",
            "imputation",
            "association_test",
            "manhattan_plot",
        ],
        "inputs": ["genotype_data", "phenotype_data"],
        "outputs": ["association_results", "manhattan_plot", "qq_plot"],
    },
    {
        "id": "metagenomics",
        "name": "Metagenomics Taxonomic Profiling",
        "description": "Metagenomic analysis: quality control, taxonomic classification, and diversity analysis",
        "category": "multiomics",
        "steps": [
            "quality_filtering",
            "kraken2_classification",
            "metaphlan",
            "diversity_analysis",
        ],
        "inputs": ["fastq_files"],
        "outputs": ["taxonomy_table", "diversity_metrics", "abundance_plots"],
    },
]


def get_all_pipeline_templates() -> List[Dict[str, Any]]:
    """Get all available pipeline templates."""
    return COMMON_PIPELINES


def get_pipeline_template(pipeline_id: str) -> Dict[str, Any] | None:
    """Get a specific pipeline template by ID."""
    for template in COMMON_PIPELINES:
        if template["id"] == pipeline_id:
            return template
    return None


def get_pipelines_by_category(category: str) -> List[Dict[str, Any]]:
    """Get pipeline templates filtered by category."""
    return [p for p in COMMON_PIPELINES if p["category"] == category]
