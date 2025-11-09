"""Pipeline template service for common workflows."""

from typing import List, Dict, Any

# Common pipeline templates available to all users
COMMON_PIPELINES = [
    {
        "id": "rna-seq-basic",
        "name": "RNA-seq Basic Analysis",
        "description": "Standard RNA-seq pipeline: QC, alignment, quantification, and DEG analysis",
        "category": "transcriptomics",
        "steps": [
            "fastqc",
            "trim_galore",
            "star_align",
            "featurecounts",
            "deseq2",
        ],
        "inputs": ["fastq_r1", "fastq_r2", "genome_index", "annotation_gtf"],
        "outputs": ["counts_matrix", "deg_results", "qc_report"],
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
            "maxquant",
            "peptide_filtering",
            "protein_inference",
            "limma_analysis",
        ],
        "inputs": ["raw_files", "fasta_database"],
        "outputs": ["protein_groups", "de_proteins", "qc_report"],
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
