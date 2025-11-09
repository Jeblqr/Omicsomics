"""Registry of available analysis tools with their specifications."""

from typing import Dict, List, Optional
from app.schemas.tools import (
    ToolSpecification,
    ToolCategory,
    OmicsType,
    ParameterType,
    ParameterDefinition,
    InputDefinition,
    OutputDefinition,
    ToolVersion,
    ResourceRequirement,
)


class ToolsRegistry:
    """Central registry for all available analysis tools."""

    def __init__(self):
        self._tools: Dict[str, ToolSpecification] = {}
        self._initialize_tools()

    def _initialize_tools(self):
        """Initialize the registry with common tools."""
        # Genomics tools
        self.register_tool(self._create_gatk_spec())
        self.register_tool(self._create_bwa_spec())
        self.register_tool(self._create_samtools_spec())
        self.register_tool(self._create_bcftools_spec())

        # Transcriptomics tools
        self.register_tool(self._create_star_spec())
        self.register_tool(self._create_salmon_spec())
        self.register_tool(self._create_deseq2_spec())
        self.register_tool(self._create_edger_spec())

        # Single-cell tools
        self.register_tool(self._create_seurat_spec())
        self.register_tool(self._create_scanpy_spec())
        self.register_tool(self._create_cellranger_spec())

        # Epigenomics tools
        self.register_tool(self._create_macs2_spec())
        self.register_tool(self._create_homer_spec())

        # Proteomics tools
        self.register_tool(self._create_maxquant_spec())

        # Metabolomics tools
        self.register_tool(self._create_xcms_spec())

        # Quality control tools
        self.register_tool(self._create_fastqc_spec())
        self.register_tool(self._create_multiqc_spec())

    def register_tool(self, spec: ToolSpecification):
        """Register a tool in the registry."""
        self._tools[spec.id] = spec

    def get_tool(self, tool_id: str) -> Optional[ToolSpecification]:
        """Get a tool specification by ID."""
        return self._tools.get(tool_id)

    def list_tools(
        self,
        category: Optional[ToolCategory] = None,
        omics_type: Optional[OmicsType] = None,
        tags: Optional[List[str]] = None,
    ) -> List[ToolSpecification]:
        """List tools with optional filters."""
        tools = list(self._tools.values())

        if category:
            tools = [t for t in tools if t.category == category]

        if omics_type:
            tools = [t for t in tools if omics_type in t.omics_types]

        if tags:
            tools = [t for t in tools if any(tag in t.tags for tag in tags)]

        return tools

    def search_tools(self, query: str) -> List[ToolSpecification]:
        """Search tools by name, description, or tags."""
        query_lower = query.lower()
        results = []

        for tool in self._tools.values():
            if (
                query_lower in tool.name.lower()
                or query_lower in tool.description.lower()
                or any(query_lower in tag.lower() for tag in tool.tags)
            ):
                results.append(tool)

        return results

    # Tool specifications

    def _create_gatk_spec(self) -> ToolSpecification:
        """GATK - Genome Analysis Toolkit."""
        return ToolSpecification(
            id="gatk",
            name="GATK",
            description="Genome Analysis Toolkit for variant discovery in high-throughput sequencing data",
            category=ToolCategory.VARIANT_CALLING,
            omics_types=[OmicsType.GENOMICS],
            current_version="4.4.0.0",
            versions=[
                ToolVersion(
                    version="4.4.0.0", docker_image="broadinstitute/gatk:4.4.0.0"
                ),
                ToolVersion(
                    version="4.3.0.0", docker_image="broadinstitute/gatk:4.3.0.0"
                ),
            ],
            inputs=[
                InputDefinition(
                    name="input_bam",
                    label="Input BAM",
                    description="Aligned reads in BAM format",
                    type=ParameterType.FILE,
                    required=True,
                    file_extensions=[".bam"],
                    format="bam",
                ),
                InputDefinition(
                    name="reference_fasta",
                    label="Reference Genome",
                    description="Reference genome in FASTA format",
                    type=ParameterType.FILE,
                    required=True,
                    file_extensions=[".fasta", ".fa"],
                    format="fasta",
                ),
            ],
            outputs=[
                OutputDefinition(
                    name="output_vcf",
                    label="Output VCF",
                    description="Variants in VCF format",
                    type=ParameterType.FILE,
                    file_extension=".vcf.gz",
                    format="vcf",
                ),
            ],
            parameters=[
                ParameterDefinition(
                    name="command",
                    type=ParameterType.ENUM,
                    label="GATK Command",
                    description="GATK command to run",
                    required=True,
                    enum_values=[
                        "HaplotypeCaller",
                        "GenotypeGVCFs",
                        "VariantFiltration",
                        "BaseRecalibrator",
                    ],
                    default="HaplotypeCaller",
                    ui_widget="dropdown",
                ),
                ParameterDefinition(
                    name="min_base_quality",
                    type=ParameterType.INTEGER,
                    label="Minimum Base Quality",
                    description="Minimum base quality score",
                    default=20,
                    min_value=0,
                    max_value=60,
                    ui_widget="slider",
                ),
            ],
            resources=ResourceRequirement(
                min_cpu=4, max_cpu=32, min_memory_gb=8, max_memory_gb=64
            ),
            docker_image="broadinstitute/gatk:4.4.0.0",
            documentation_url="https://gatk.broadinstitute.org/",
            citation="McKenna et al. (2010) The Genome Analysis Toolkit",
            license="BSD-3-Clause",
            tags=["variant-calling", "genomics", "snp", "indel"],
        )

    def _create_star_spec(self) -> ToolSpecification:
        """STAR - RNA-seq aligner."""
        return ToolSpecification(
            id="star",
            name="STAR",
            description="Spliced Transcripts Alignment to a Reference",
            category=ToolCategory.ALIGNMENT,
            omics_types=[OmicsType.TRANSCRIPTOMICS],
            current_version="2.7.11a",
            versions=[
                ToolVersion(
                    version="2.7.11a", docker_image="quay.io/biocontainers/star:2.7.11a"
                ),
            ],
            inputs=[
                InputDefinition(
                    name="fastq_1",
                    label="FASTQ Read 1",
                    description="Forward reads in FASTQ format",
                    type=ParameterType.FILE,
                    required=True,
                    file_extensions=[".fastq", ".fq", ".fastq.gz", ".fq.gz"],
                    format="fastq",
                ),
                InputDefinition(
                    name="fastq_2",
                    label="FASTQ Read 2",
                    description="Reverse reads in FASTQ format (for paired-end)",
                    type=ParameterType.FILE,
                    required=False,
                    file_extensions=[".fastq", ".fq", ".fastq.gz", ".fq.gz"],
                    format="fastq",
                ),
                InputDefinition(
                    name="genome_index",
                    label="Genome Index",
                    description="STAR genome index directory",
                    type=ParameterType.DIRECTORY,
                    required=True,
                ),
            ],
            outputs=[
                OutputDefinition(
                    name="aligned_bam",
                    label="Aligned BAM",
                    description="Aligned reads in BAM format",
                    type=ParameterType.FILE,
                    file_extension=".bam",
                    format="bam",
                ),
                OutputDefinition(
                    name="gene_counts",
                    label="Gene Counts",
                    description="Gene expression counts",
                    type=ParameterType.FILE,
                    file_extension=".tab",
                    format="tsv",
                ),
            ],
            parameters=[
                ParameterDefinition(
                    name="threads",
                    type=ParameterType.INTEGER,
                    label="Threads",
                    description="Number of threads to use",
                    default=8,
                    min_value=1,
                    max_value=64,
                    ui_widget="slider",
                ),
                ParameterDefinition(
                    name="out_sam_type",
                    type=ParameterType.ENUM,
                    label="Output Format",
                    description="Output file format",
                    enum_values=["BAM SortedByCoordinate", "BAM Unsorted", "SAM"],
                    default="BAM SortedByCoordinate",
                    ui_widget="dropdown",
                ),
                ParameterDefinition(
                    name="quant_mode",
                    type=ParameterType.ENUM,
                    label="Quantification Mode",
                    description="Gene quantification mode",
                    enum_values=["GeneCounts", "TranscriptomeSAM", "-"],
                    default="GeneCounts",
                    ui_widget="dropdown",
                ),
            ],
            resources=ResourceRequirement(
                min_cpu=8, max_cpu=64, min_memory_gb=32, max_memory_gb=128
            ),
            docker_image="quay.io/biocontainers/star:2.7.11a",
            documentation_url="https://github.com/alexdobin/STAR",
            citation="Dobin et al. (2013) STAR: ultrafast universal RNA-seq aligner",
            license="MIT",
            tags=["alignment", "rna-seq", "transcriptomics", "splicing"],
        )

    def _create_deseq2_spec(self) -> ToolSpecification:
        """DESeq2 - Differential expression analysis."""
        return ToolSpecification(
            id="deseq2",
            name="DESeq2",
            description="Differential gene expression analysis based on the negative binomial distribution",
            category=ToolCategory.DIFFERENTIAL_EXPRESSION,
            omics_types=[OmicsType.TRANSCRIPTOMICS],
            current_version="1.40.0",
            versions=[
                ToolVersion(
                    version="1.40.0",
                    docker_image="bioconductor/bioconductor_docker:RELEASE_3_17",
                ),
            ],
            inputs=[
                InputDefinition(
                    name="count_matrix",
                    label="Count Matrix",
                    description="Gene expression count matrix",
                    type=ParameterType.FILE,
                    required=True,
                    file_extensions=[".csv", ".tsv", ".txt"],
                    format="csv",
                ),
                InputDefinition(
                    name="sample_info",
                    label="Sample Information",
                    description="Sample metadata with conditions",
                    type=ParameterType.FILE,
                    required=True,
                    file_extensions=[".csv", ".tsv"],
                    format="csv",
                ),
            ],
            outputs=[
                OutputDefinition(
                    name="results_table",
                    label="Results Table",
                    description="Differential expression results",
                    type=ParameterType.FILE,
                    file_extension=".csv",
                    format="csv",
                ),
                OutputDefinition(
                    name="normalized_counts",
                    label="Normalized Counts",
                    description="Normalized expression values",
                    type=ParameterType.FILE,
                    file_extension=".csv",
                    format="csv",
                ),
            ],
            parameters=[
                ParameterDefinition(
                    name="design_formula",
                    type=ParameterType.STRING,
                    label="Design Formula",
                    description="Design formula (e.g., '~ condition')",
                    required=True,
                    default="~ condition",
                ),
                ParameterDefinition(
                    name="alpha",
                    type=ParameterType.FLOAT,
                    label="Significance Level",
                    description="Adjusted p-value threshold",
                    default=0.05,
                    min_value=0.0,
                    max_value=1.0,
                    ui_widget="slider",
                ),
                ParameterDefinition(
                    name="lfc_threshold",
                    type=ParameterType.FLOAT,
                    label="Log2 Fold Change Threshold",
                    description="Minimum log2 fold change",
                    default=0.0,
                    min_value=0.0,
                    max_value=10.0,
                    ui_widget="slider",
                ),
            ],
            resources=ResourceRequirement(
                min_cpu=2, max_cpu=16, min_memory_gb=8, max_memory_gb=32
            ),
            docker_image="bioconductor/bioconductor_docker:RELEASE_3_17",
            documentation_url="https://bioconductor.org/packages/DESeq2/",
            citation="Love et al. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2",
            license="LGPL",
            tags=[
                "differential-expression",
                "rna-seq",
                "statistics",
                "transcriptomics",
            ],
        )

    def _create_seurat_spec(self) -> ToolSpecification:
        """Seurat - Single-cell analysis."""
        return ToolSpecification(
            id="seurat",
            name="Seurat",
            description="R toolkit for single cell genomics",
            category=ToolCategory.SINGLE_CELL,
            omics_types=[OmicsType.SINGLE_CELL],
            current_version="5.0.0",
            versions=[
                ToolVersion(version="5.0.0", docker_image="satijalab/seurat:5.0.0"),
            ],
            inputs=[
                InputDefinition(
                    name="count_matrix",
                    label="Count Matrix",
                    description="Single-cell expression matrix",
                    type=ParameterType.FILE,
                    required=True,
                    file_extensions=[".h5", ".mtx", ".csv"],
                    format="h5ad",
                ),
            ],
            outputs=[
                OutputDefinition(
                    name="seurat_object",
                    label="Seurat Object",
                    description="Processed Seurat object",
                    type=ParameterType.FILE,
                    file_extension=".rds",
                    format="rds",
                ),
                OutputDefinition(
                    name="clusters",
                    label="Cell Clusters",
                    description="Cell clustering results",
                    type=ParameterType.FILE,
                    file_extension=".csv",
                    format="csv",
                ),
            ],
            parameters=[
                ParameterDefinition(
                    name="min_cells",
                    type=ParameterType.INTEGER,
                    label="Minimum Cells",
                    description="Minimum number of cells expressing a gene",
                    default=3,
                    min_value=0,
                    max_value=100,
                ),
                ParameterDefinition(
                    name="min_features",
                    type=ParameterType.INTEGER,
                    label="Minimum Features",
                    description="Minimum number of features per cell",
                    default=200,
                    min_value=0,
                    max_value=10000,
                ),
                ParameterDefinition(
                    name="resolution",
                    type=ParameterType.FLOAT,
                    label="Clustering Resolution",
                    description="Resolution parameter for clustering",
                    default=0.8,
                    min_value=0.1,
                    max_value=2.0,
                    ui_widget="slider",
                ),
            ],
            resources=ResourceRequirement(
                min_cpu=4, max_cpu=32, min_memory_gb=16, max_memory_gb=128
            ),
            docker_image="satijalab/seurat:5.0.0",
            documentation_url="https://satijalab.org/seurat/",
            citation="Hao et al. (2021) Integrated analysis of multimodal single-cell data",
            license="MIT",
            tags=[
                "single-cell",
                "clustering",
                "dimensionality-reduction",
                "visualization",
            ],
        )

    def _create_macs2_spec(self) -> ToolSpecification:
        """MACS2 - Peak calling."""
        return ToolSpecification(
            id="macs2",
            name="MACS2",
            description="Model-based Analysis of ChIP-Seq for peak calling",
            category=ToolCategory.PEAK_CALLING,
            omics_types=[OmicsType.EPIGENOMICS],
            current_version="2.2.9",
            versions=[
                ToolVersion(
                    version="2.2.9", docker_image="quay.io/biocontainers/macs2:2.2.9"
                ),
            ],
            inputs=[
                InputDefinition(
                    name="treatment_bam",
                    label="Treatment BAM",
                    description="ChIP-seq treatment file",
                    type=ParameterType.FILE,
                    required=True,
                    file_extensions=[".bam"],
                    format="bam",
                ),
                InputDefinition(
                    name="control_bam",
                    label="Control BAM",
                    description="ChIP-seq control/input file",
                    type=ParameterType.FILE,
                    required=False,
                    file_extensions=[".bam"],
                    format="bam",
                ),
            ],
            outputs=[
                OutputDefinition(
                    name="peaks",
                    label="Peaks",
                    description="Called peaks in narrowPeak format",
                    type=ParameterType.FILE,
                    file_extension="_peaks.narrowPeak",
                    format="bed",
                ),
                OutputDefinition(
                    name="summits",
                    label="Peak Summits",
                    description="Peak summits in BED format",
                    type=ParameterType.FILE,
                    file_extension="_summits.bed",
                    format="bed",
                ),
            ],
            parameters=[
                ParameterDefinition(
                    name="genome_size",
                    type=ParameterType.ENUM,
                    label="Genome Size",
                    description="Effective genome size",
                    enum_values=["hs", "mm", "ce", "dm"],
                    default="hs",
                    ui_widget="dropdown",
                ),
                ParameterDefinition(
                    name="qvalue",
                    type=ParameterType.FLOAT,
                    label="Q-value Cutoff",
                    description="Minimum FDR cutoff",
                    default=0.05,
                    min_value=0.0,
                    max_value=1.0,
                    ui_widget="slider",
                ),
            ],
            resources=ResourceRequirement(
                min_cpu=2, max_cpu=8, min_memory_gb=4, max_memory_gb=16
            ),
            docker_image="quay.io/biocontainers/macs2:2.2.9",
            documentation_url="https://github.com/macs3-project/MACS",
            citation="Zhang et al. (2008) Model-based Analysis of ChIP-Seq (MACS)",
            license="BSD",
            tags=["chip-seq", "peak-calling", "epigenomics"],
        )

    def _create_fastqc_spec(self) -> ToolSpecification:
        """FastQC - Quality control."""
        return ToolSpecification(
            id="fastqc",
            name="FastQC",
            description="Quality control tool for high throughput sequence data",
            category=ToolCategory.QUALITY_CONTROL,
            omics_types=[OmicsType.GENOMICS, OmicsType.TRANSCRIPTOMICS],
            current_version="0.12.1",
            versions=[
                ToolVersion(
                    version="0.12.1", docker_image="biocontainers/fastqc:v0.12.1"
                ),
            ],
            inputs=[
                InputDefinition(
                    name="fastq_files",
                    label="FASTQ Files",
                    description="Raw sequencing reads",
                    type=ParameterType.FILE,
                    required=True,
                    multiple=True,
                    file_extensions=[".fastq", ".fq", ".fastq.gz", ".fq.gz"],
                    format="fastq",
                ),
            ],
            outputs=[
                OutputDefinition(
                    name="html_report",
                    label="HTML Report",
                    description="Quality control report",
                    type=ParameterType.FILE,
                    file_extension=".html",
                    format="html",
                ),
            ],
            parameters=[
                ParameterDefinition(
                    name="threads",
                    type=ParameterType.INTEGER,
                    label="Threads",
                    description="Number of threads",
                    default=4,
                    min_value=1,
                    max_value=32,
                ),
            ],
            resources=ResourceRequirement(
                min_cpu=1, max_cpu=8, min_memory_gb=2, max_memory_gb=8
            ),
            docker_image="biocontainers/fastqc:v0.12.1",
            documentation_url="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/",
            license="GPL",
            tags=["qc", "quality-control", "fastq"],
        )

    # Placeholder methods for other tools
    def _create_bwa_spec(self) -> ToolSpecification:
        """BWA - DNA sequence alignment."""
        return ToolSpecification(
            id="bwa",
            name="BWA",
            description="Burrows-Wheeler Aligner for short-read alignment",
            category=ToolCategory.ALIGNMENT,
            omics_types=[OmicsType.GENOMICS],
            current_version="0.7.17",
            docker_image="biocontainers/bwa:0.7.17",
            tags=["alignment", "genomics"],
        )

    def _create_samtools_spec(self) -> ToolSpecification:
        return ToolSpecification(
            id="samtools",
            name="SAMtools",
            description="Tools for manipulating SAM/BAM files",
            category=ToolCategory.PREPROCESSING,
            omics_types=[OmicsType.GENOMICS, OmicsType.TRANSCRIPTOMICS],
            current_version="1.18",
            docker_image="biocontainers/samtools:1.18",
            tags=["bam", "sam", "preprocessing"],
        )

    def _create_bcftools_spec(self) -> ToolSpecification:
        return ToolSpecification(
            id="bcftools",
            name="BCFtools",
            description="Tools for manipulating VCF/BCF files",
            category=ToolCategory.PREPROCESSING,
            omics_types=[OmicsType.GENOMICS],
            current_version="1.18",
            docker_image="biocontainers/bcftools:1.18",
            tags=["vcf", "variant", "preprocessing"],
        )

    def _create_salmon_spec(self) -> ToolSpecification:
        return ToolSpecification(
            id="salmon",
            name="Salmon",
            description="Transcript quantification from RNA-seq data",
            category=ToolCategory.QUANTIFICATION,
            omics_types=[OmicsType.TRANSCRIPTOMICS],
            current_version="1.10.0",
            docker_image="combinelab/salmon:1.10.0",
            tags=["quantification", "rna-seq", "transcriptomics"],
        )

    def _create_edger_spec(self) -> ToolSpecification:
        return ToolSpecification(
            id="edger",
            name="edgeR",
            description="Differential expression analysis of RNA-seq data",
            category=ToolCategory.DIFFERENTIAL_EXPRESSION,
            omics_types=[OmicsType.TRANSCRIPTOMICS],
            current_version="3.42.0",
            docker_image="bioconductor/bioconductor_docker:RELEASE_3_17",
            tags=["differential-expression", "rna-seq"],
        )

    def _create_scanpy_spec(self) -> ToolSpecification:
        return ToolSpecification(
            id="scanpy",
            name="Scanpy",
            description="Single-cell analysis in Python",
            category=ToolCategory.SINGLE_CELL,
            omics_types=[OmicsType.SINGLE_CELL],
            current_version="1.9.6",
            docker_image="python:3.11",
            tags=["single-cell", "python"],
        )

    def _create_cellranger_spec(self) -> ToolSpecification:
        return ToolSpecification(
            id="cellranger",
            name="Cell Ranger",
            description="10x Genomics single-cell data processing",
            category=ToolCategory.PREPROCESSING,
            omics_types=[OmicsType.SINGLE_CELL],
            current_version="7.2.0",
            tags=["single-cell", "10x", "preprocessing"],
        )

    def _create_homer_spec(self) -> ToolSpecification:
        return ToolSpecification(
            id="homer",
            name="HOMER",
            description="Motif discovery and ChIP-seq analysis",
            category=ToolCategory.MOTIF_ANALYSIS,
            omics_types=[OmicsType.EPIGENOMICS],
            current_version="4.11",
            tags=["motif", "chip-seq", "epigenomics"],
        )

    def _create_maxquant_spec(self) -> ToolSpecification:
        return ToolSpecification(
            id="maxquant",
            name="MaxQuant",
            description="Quantitative proteomics software",
            category=ToolCategory.PROTEOMICS,
            omics_types=[OmicsType.PROTEOMICS],
            current_version="2.4.0",
            tags=["proteomics", "mass-spec", "quantification"],
        )

    def _create_xcms_spec(self) -> ToolSpecification:
        return ToolSpecification(
            id="xcms",
            name="XCMS",
            description="LC/MS data processing and analysis",
            category=ToolCategory.METABOLOMICS,
            omics_types=[OmicsType.METABOLOMICS],
            current_version="3.22.0",
            docker_image="bioconductor/bioconductor_docker:RELEASE_3_17",
            tags=["metabolomics", "lcms"],
        )

    def _create_multiqc_spec(self) -> ToolSpecification:
        return ToolSpecification(
            id="multiqc",
            name="MultiQC",
            description="Aggregate results from bioinformatics analyses",
            category=ToolCategory.QUALITY_CONTROL,
            omics_types=[
                OmicsType.GENOMICS,
                OmicsType.TRANSCRIPTOMICS,
                OmicsType.EPIGENOMICS,
            ],
            current_version="1.19",
            docker_image="ewels/multiqc:1.19",
            tags=["qc", "report", "aggregation"],
        )


# Global registry instance
_registry = ToolsRegistry()


def get_tools_registry() -> ToolsRegistry:
    """Get the global tools registry instance."""
    return _registry
