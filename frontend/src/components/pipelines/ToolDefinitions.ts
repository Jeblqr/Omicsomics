/**
 * Three-tier tool library structure:
 * Level 1: Omics Type (Genomics, Proteomics, Metabolomics, Transcriptomics, Multi-omics)
 * Level 2: Operation Category (Input, QC, Alignment, Peak Calling, Variant Calling, etc.)
 * Level 3: Specific Tools with parameters
 */

export interface ToolParameter {
  name: string;
  type: 'string' | 'number' | 'boolean' | 'select' | 'multiselect' | 'file';
  label: string;
  description: string;
  default?: any;
  required?: boolean;
  options?: string[]; // For select/multiselect
  min?: number; // For number
  max?: number; // For number
  unit?: string; // For display (e.g., "bp", "MB", "%")
}

export interface ToolDefinition {
  id: string;
  name: string;
  type: 'input' | 'process' | 'filter' | 'transform' | 'analysis' | 'output' | 'visualization';
  tool: string;
  version: string;
  omicsType: string; // Level 1
  operationCategory: string; // Level 2
  description: string;
  parameterTemplate: ToolParameter[];
  defaultParameters?: Record<string, any>;
  documentation?: string;
}

// ========== GENOMICS TOOLS ==========

const genomicsTools: ToolDefinition[] = [
  // Input
  {
    id: 'genomics-fastq-input',
    name: 'FASTQ Input',
    type: 'input',
    tool: 'file-input',
    version: '1.0',
    omicsType: 'Genomics',
    operationCategory: 'Input',
    description: 'Read FASTQ sequencing files (single-end or paired-end)',
    parameterTemplate: [
      { name: 'file_path', type: 'file', label: 'FASTQ File', description: 'Path to FASTQ file', required: true },
      { name: 'paired_end', type: 'boolean', label: 'Paired-end', description: 'Is this paired-end data?', default: false },
      { name: 'mate_file', type: 'file', label: 'Mate File', description: 'Path to mate file (for paired-end)', required: false },
    ],
    defaultParameters: { paired_end: false },
  },
  {
    id: 'genomics-vcf-input',
    name: 'VCF Input',
    type: 'input',
    tool: 'file-input',
    version: '1.0',
    omicsType: 'Genomics',
    operationCategory: 'Input',
    description: 'Read VCF variant call files',
    parameterTemplate: [
      { name: 'file_path', type: 'file', label: 'VCF File', description: 'Path to VCF file', required: true },
      { name: 'index_file', type: 'file', label: 'Index File', description: 'Path to VCF index (.tbi)', required: false },
    ],
    defaultParameters: {},
  },

  // QC
  {
    id: 'genomics-fastqc',
    name: 'FastQC',
    type: 'process',
    tool: 'fastqc',
    version: '0.12.1',
    omicsType: 'Genomics',
    operationCategory: 'QC',
    description: 'Comprehensive quality control checks for high-throughput sequence data',
    parameterTemplate: [
      { name: 'threads', type: 'number', label: 'Threads', description: 'Number of CPU threads', default: 4, min: 1, max: 32 },
      { name: 'kmers', type: 'number', label: 'K-mer Length', description: 'K-mer length for overrepresented sequences', default: 7, min: 2, max: 10 },
      { name: 'contaminants', type: 'file', label: 'Contaminants File', description: 'Custom contaminants list', required: false },
    ],
    defaultParameters: { threads: 4, kmers: 7 },
    documentation: 'FastQC performs quality checks on raw sequence data',
  },
  {
    id: 'genomics-trim-galore',
    name: 'Trim Galore',
    type: 'filter',
    tool: 'trim_galore',
    version: '0.6.10',
    omicsType: 'Genomics',
    operationCategory: 'QC',
    description: 'Adapter and quality trimming with automatic adapter detection',
    parameterTemplate: [
      { name: 'quality', type: 'number', label: 'Quality Threshold', description: 'Phred quality score threshold', default: 20, min: 0, max: 40, unit: 'Q' },
      { name: 'length', type: 'number', label: 'Min Length', description: 'Minimum read length after trimming', default: 20, min: 10, max: 200, unit: 'bp' },
      { name: 'stringency', type: 'number', label: 'Stringency', description: 'Adapter overlap stringency', default: 1, min: 1, max: 10 },
      { name: 'error_rate', type: 'number', label: 'Error Rate', description: 'Maximum allowed error rate', default: 0.1, min: 0, max: 1 },
      { name: 'clip_r1', type: 'number', label: 'Clip R1 5\'', description: 'Remove N bp from 5\' end of R1', default: 0, min: 0, max: 50, unit: 'bp' },
      { name: 'clip_r2', type: 'number', label: 'Clip R2 5\'', description: 'Remove N bp from 5\' end of R2', default: 0, min: 0, max: 50, unit: 'bp' },
      { name: 'three_prime_clip_r1', type: 'number', label: 'Clip R1 3\'', description: 'Remove N bp from 3\' end of R1', default: 0, min: 0, max: 50, unit: 'bp' },
    ],
    defaultParameters: { quality: 20, length: 20, stringency: 1, error_rate: 0.1 },
  },
  {
    id: 'genomics-cutadapt',
    name: 'Cutadapt',
    type: 'filter',
    tool: 'cutadapt',
    version: '4.4',
    omicsType: 'Genomics',
    operationCategory: 'QC',
    description: 'Remove adapter sequences from high-throughput sequencing reads',
    parameterTemplate: [
      { name: 'adapter_3p', type: 'string', label: '3\' Adapter', description: '3\' adapter sequence', default: '' },
      { name: 'adapter_5p', type: 'string', label: '5\' Adapter', description: '5\' adapter sequence', default: '' },
      { name: 'quality_cutoff', type: 'number', label: 'Quality Cutoff', description: 'Quality cutoff value', default: 20, min: 0, max: 40, unit: 'Q' },
      { name: 'minimum_length', type: 'number', label: 'Min Length', description: 'Discard reads shorter than this', default: 20, min: 1, max: 500, unit: 'bp' },
      { name: 'max_n', type: 'number', label: 'Max N Content', description: 'Discard reads with more than N ambiguous bases', default: 10, min: 0, max: 100 },
    ],
    defaultParameters: { quality_cutoff: 20, minimum_length: 20, max_n: 10 },
  },

  // Alignment
  {
    id: 'genomics-bwa-mem',
    name: 'BWA-MEM',
    type: 'process',
    tool: 'bwa',
    version: '0.7.17',
    omicsType: 'Genomics',
    operationCategory: 'Alignment',
    description: 'Burrows-Wheeler Aligner for short read alignment (DNA-seq)',
    parameterTemplate: [
      { name: 'reference', type: 'file', label: 'Reference Genome', description: 'Path to reference genome FASTA', required: true },
      { name: 'threads', type: 'number', label: 'Threads', description: 'Number of CPU threads', default: 8, min: 1, max: 64 },
      { name: 'min_seed_length', type: 'number', label: 'Min Seed Length', description: 'Minimum seed length', default: 19, min: 5, max: 50, unit: 'bp' },
      { name: 'band_width', type: 'number', label: 'Band Width', description: 'Band width for banded alignment', default: 100, min: 10, max: 1000 },
      { name: 'mark_shorter', type: 'boolean', label: 'Mark Shorter Splits', description: 'Mark shorter split hits as secondary', default: false },
    ],
    defaultParameters: { threads: 8, min_seed_length: 19, band_width: 100 },
  },
  {
    id: 'genomics-bowtie2',
    name: 'Bowtie2',
    type: 'process',
    tool: 'bowtie2',
    version: '2.5.1',
    omicsType: 'Genomics',
    operationCategory: 'Alignment',
    description: 'Fast and memory-efficient read aligner',
    parameterTemplate: [
      { name: 'reference_index', type: 'file', label: 'Reference Index', description: 'Path to Bowtie2 index', required: true },
      { name: 'threads', type: 'number', label: 'Threads', description: 'Number of CPU threads', default: 8, min: 1, max: 64 },
      { name: 'mode', type: 'select', label: 'Alignment Mode', description: 'End-to-end or local alignment', default: 'end-to-end', options: ['end-to-end', 'local'] },
      { name: 'sensitivity', type: 'select', label: 'Sensitivity', description: 'Alignment sensitivity preset', default: 'sensitive', options: ['very-fast', 'fast', 'sensitive', 'very-sensitive'] },
    ],
    defaultParameters: { threads: 8, mode: 'end-to-end', sensitivity: 'sensitive' },
  },

  // Variant Calling
  {
    id: 'genomics-gatk-haplotypecaller',
    name: 'GATK HaplotypeCaller',
    type: 'analysis',
    tool: 'gatk_haplotypecaller',
    version: '4.4.0',
    omicsType: 'Genomics',
    operationCategory: 'Variant Calling',
    description: 'Call germline SNPs and indels via local re-assembly of haplotypes',
    parameterTemplate: [
      { name: 'reference', type: 'file', label: 'Reference Genome', description: 'Reference genome FASTA', required: true },
      { name: 'ploidy', type: 'number', label: 'Ploidy', description: 'Sample ploidy', default: 2, min: 1, max: 10 },
      { name: 'min_base_quality', type: 'number', label: 'Min Base Quality', description: 'Minimum base quality', default: 10, min: 0, max: 40, unit: 'Q' },
      { name: 'emit_ref_confidence', type: 'select', label: 'Emit Ref Confidence', description: 'Reference confidence mode', default: 'NONE', options: ['NONE', 'BP_RESOLUTION', 'GVCF'] },
    ],
    defaultParameters: { ploidy: 2, min_base_quality: 10, emit_ref_confidence: 'NONE' },
  },
  {
    id: 'genomics-bcftools-call',
    name: 'BCFtools Call',
    type: 'analysis',
    tool: 'bcftools',
    version: '1.17',
    omicsType: 'Genomics',
    operationCategory: 'Variant Calling',
    description: 'SNP and indel calling from samtools mpileup',
    parameterTemplate: [
      { name: 'ploidy', type: 'number', label: 'Ploidy', description: 'Sample ploidy', default: 2, min: 1, max: 10 },
      { name: 'caller', type: 'select', label: 'Caller', description: 'Calling algorithm', default: 'multiallelic', options: ['multiallelic', 'consensus'] },
      { name: 'min_alt_frac', type: 'number', label: 'Min Alt Fraction', description: 'Minimum alternate allele fraction', default: 0.2, min: 0, max: 1 },
    ],
    defaultParameters: { ploidy: 2, caller: 'multiallelic', min_alt_frac: 0.2 },
  },

  // Output
  {
    id: 'genomics-vcf-output',
    name: 'VCF Output',
    type: 'output',
    tool: 'file-output',
    version: '1.0',
    omicsType: 'Genomics',
    operationCategory: 'Output',
    description: 'Save variants in VCF format',
    parameterTemplate: [
      { name: 'output_path', type: 'file', label: 'Output Path', description: 'Output VCF file path', required: true },
      { name: 'compress', type: 'boolean', label: 'Compress', description: 'Compress output with bgzip', default: true },
      { name: 'index', type: 'boolean', label: 'Create Index', description: 'Create tabix index', default: true },
    ],
    defaultParameters: { compress: true, index: true },
  },
];

// ========== TRANSCRIPTOMICS TOOLS ==========

const transcriptomicsTools: ToolDefinition[] = [
  // Alignment
  {
    id: 'trans-star',
    name: 'STAR',
    type: 'process',
    tool: 'star',
    version: '2.7.11',
    omicsType: 'Transcriptomics',
    operationCategory: 'Alignment',
    description: 'Spliced Transcripts Alignment to a Reference - RNA-seq aligner',
    parameterTemplate: [
      { name: 'genome_dir', type: 'file', label: 'Genome Directory', description: 'Path to STAR genome index', required: true },
      { name: 'threads', type: 'number', label: 'Threads', description: 'Number of CPU threads', default: 8, min: 1, max: 64 },
      { name: 'out_sam_type', type: 'select', label: 'Output Format', description: 'Output file format', default: 'BAM SortedByCoordinate', options: ['SAM', 'BAM Unsorted', 'BAM SortedByCoordinate'] },
      { name: 'quant_mode', type: 'select', label: 'Quantification', description: 'Gene/transcript quantification', default: 'GeneCounts', options: ['None', 'GeneCounts', 'TranscriptomeSAM'] },
    ],
    defaultParameters: { threads: 8, out_sam_type: 'BAM SortedByCoordinate', quant_mode: 'GeneCounts' },
  },
  {
    id: 'trans-hisat2',
    name: 'HISAT2',
    type: 'process',
    tool: 'hisat2',
    version: '2.2.1',
    omicsType: 'Transcriptomics',
    operationCategory: 'Alignment',
    description: 'Fast and sensitive alignment for RNA-seq',
    parameterTemplate: [
      { name: 'index', type: 'file', label: 'Index', description: 'Path to HISAT2 index', required: true },
      { name: 'threads', type: 'number', label: 'Threads', description: 'Number of CPU threads', default: 8, min: 1, max: 64 },
      { name: 'rna_strandness', type: 'select', label: 'Strandness', description: 'RNA strand specificity', default: 'unstranded', options: ['unstranded', 'FR', 'RF'] },
    ],
    defaultParameters: { threads: 8, rna_strandness: 'unstranded' },
  },

  // Feature Extraction
  {
    id: 'trans-feature-counts',
    name: 'featureCounts',
    type: 'process',
    tool: 'featureCounts',
    version: '2.0.6',
    omicsType: 'Transcriptomics',
    operationCategory: 'Feature Extraction',
    description: 'Count reads mapping to genomic features',
    parameterTemplate: [
      { name: 'annotation', type: 'file', label: 'GTF/GFF File', description: 'Gene annotation file', required: true },
      { name: 'threads', type: 'number', label: 'Threads', description: 'Number of CPU threads', default: 4, min: 1, max: 32 },
      { name: 'is_paired_end', type: 'boolean', label: 'Paired-end', description: 'Is data paired-end?', default: true },
      { name: 'strandness', type: 'select', label: 'Strandness', description: 'Strand specificity', default: '0', options: ['0', '1', '2'] },
      { name: 'feature_type', type: 'string', label: 'Feature Type', description: 'Feature type to count', default: 'exon' },
    ],
    defaultParameters: { threads: 4, is_paired_end: true, strandness: '0', feature_type: 'exon' },
  },
  {
    id: 'trans-htseq-count',
    name: 'HTSeq Count',
    type: 'process',
    tool: 'htseq',
    version: '2.0.3',
    omicsType: 'Transcriptomics',
    operationCategory: 'Feature Extraction',
    description: 'Count aligned reads overlapping genomic features',
    parameterTemplate: [
      { name: 'annotation', type: 'file', label: 'GTF File', description: 'Gene annotation GTF', required: true },
      { name: 'mode', type: 'select', label: 'Mode', description: 'Counting mode', default: 'union', options: ['union', 'intersection-strict', 'intersection-nonempty'] },
      { name: 'strandness', type: 'select', label: 'Strandness', description: 'Strand specificity', default: 'no', options: ['yes', 'no', 'reverse'] },
      { name: 'min_quality', type: 'number', label: 'Min Quality', description: 'Min alignment quality', default: 10, min: 0, max: 60, unit: 'Q' },
    ],
    defaultParameters: { mode: 'union', strandness: 'no', min_quality: 10 },
  },

  // Analysis
  {
    id: 'trans-deseq2',
    name: 'DESeq2',
    type: 'analysis',
    tool: 'deseq2',
    version: '1.40.0',
    omicsType: 'Transcriptomics',
    operationCategory: 'Differential Expression',
    description: 'Differential gene expression analysis using negative binomial distribution',
    parameterTemplate: [
      { name: 'alpha', type: 'number', label: 'Alpha', description: 'Significance cutoff (FDR)', default: 0.05, min: 0, max: 1 },
      { name: 'lfc_threshold', type: 'number', label: 'LFC Threshold', description: 'Log2 fold change threshold', default: 0, min: 0, max: 10 },
      { name: 'fit_type', type: 'select', label: 'Fit Type', description: 'Dispersion fit type', default: 'parametric', options: ['parametric', 'local', 'mean'] },
      { name: 'test', type: 'select', label: 'Test', description: 'Statistical test', default: 'Wald', options: ['Wald', 'LRT'] },
    ],
    defaultParameters: { alpha: 0.05, lfc_threshold: 0, fit_type: 'parametric', test: 'Wald' },
  },
  {
    id: 'trans-edger',
    name: 'edgeR',
    type: 'analysis',
    tool: 'edger',
    version: '3.42.0',
    omicsType: 'Transcriptomics',
    operationCategory: 'Differential Expression',
    description: 'Differential expression analysis of RNA-seq data',
    parameterTemplate: [
      { name: 'fdr', type: 'number', label: 'FDR', description: 'False discovery rate cutoff', default: 0.05, min: 0, max: 1 },
      { name: 'lfc', type: 'number', label: 'Log FC', description: 'Log fold change cutoff', default: 1, min: 0, max: 10 },
      { name: 'method', type: 'select', label: 'Method', description: 'Test method', default: 'glmQLFit', options: ['glmQLFit', 'glmFit', 'exactTest'] },
    ],
    defaultParameters: { fdr: 0.05, lfc: 1, method: 'glmQLFit' },
  },

  // Normalization
  {
    id: 'trans-normalize',
    name: 'RNA-seq Normalization',
    type: 'transform',
    tool: 'normalize',
    version: '1.0',
    omicsType: 'Transcriptomics',
    operationCategory: 'Normalization',
    description: 'Normalize RNA-seq count data',
    parameterTemplate: [
      { name: 'method', type: 'select', label: 'Method', description: 'Normalization method', default: 'TMM', options: ['TMM', 'RLE', 'upperquartile', 'CPM', 'TPM', 'RPKM'] },
      { name: 'log_transform', type: 'boolean', label: 'Log Transform', description: 'Apply log2 transformation', default: false },
    ],
    defaultParameters: { method: 'TMM', log_transform: false },
  },
];

// ========== PROTEOMICS TOOLS ==========

const proteomicsTools: ToolDefinition[] = [
  // Input
  {
    id: 'prot-raw-input',
    name: 'Raw File Input',
    type: 'input',
    tool: 'file-input',
    version: '1.0',
    omicsType: 'Proteomics',
    operationCategory: 'Input',
    description: 'Read mass spectrometry raw data files (.raw, .mzML, .mzXML)',
    parameterTemplate: [
      { name: 'file_path', type: 'file', label: 'File Path', description: 'Path to MS data file', required: true },
      { name: 'format', type: 'select', label: 'Format', description: 'File format', default: 'mzML', options: ['mzML', 'mzXML', 'raw'] },
    ],
    defaultParameters: { format: 'mzML' },
  },

  // Peak Calling
  {
    id: 'prot-maxquant',
    name: 'MaxQuant',
    type: 'analysis',
    tool: 'maxquant',
    version: '2.4.0',
    omicsType: 'Proteomics',
    operationCategory: 'Peak Calling',
    description: 'Quantitative proteomics analysis with MaxQuant',
    parameterTemplate: [
      { name: 'fasta', type: 'file', label: 'FASTA Database', description: 'Protein database FASTA', required: true },
      { name: 'threads', type: 'number', label: 'Threads', description: 'Number of CPU threads', default: 8, min: 1, max: 64 },
      { name: 'fdr', type: 'number', label: 'FDR', description: 'False discovery rate', default: 0.01, min: 0, max: 1 },
      { name: 'min_peptide_length', type: 'number', label: 'Min Peptide Length', description: 'Minimum peptide length', default: 7, min: 5, max: 20, unit: 'AA' },
      { name: 'label_free', type: 'boolean', label: 'Label-free Quant', description: 'Enable LFQ', default: true },
    ],
    defaultParameters: { threads: 8, fdr: 0.01, min_peptide_length: 7, label_free: true },
  },
  {
    id: 'prot-openms-peakpicker',
    name: 'OpenMS PeakPicker',
    type: 'process',
    tool: 'openms',
    version: '3.0.0',
    omicsType: 'Proteomics',
    operationCategory: 'Peak Calling',
    description: 'Detect peaks in MS spectra',
    parameterTemplate: [
      { name: 'signal_to_noise', type: 'number', label: 'Signal-to-Noise', description: 'Minimum signal-to-noise ratio', default: 3.0, min: 0, max: 100 },
      { name: 'peak_width', type: 'number', label: 'Peak Width', description: 'Expected peak width', default: 0.1, min: 0.01, max: 10, unit: 'Da' },
    ],
    defaultParameters: { signal_to_noise: 3.0, peak_width: 0.1 },
  },

  // Analysis
  {
    id: 'prot-perseus',
    name: 'Perseus',
    type: 'analysis',
    tool: 'perseus',
    version: '1.6.15',
    omicsType: 'Proteomics',
    operationCategory: 'Differential Expression',
    description: 'Statistical analysis of proteomics data',
    parameterTemplate: [
      { name: 'fdr', type: 'number', label: 'FDR', description: 'False discovery rate', default: 0.05, min: 0, max: 1 },
      { name: 's0', type: 'number', label: 'S0', description: 'Artificial within-group variance', default: 0.1, min: 0, max: 10 },
    ],
    defaultParameters: { fdr: 0.05, s0: 0.1 },
  },
];

// ========== METABOLOMICS TOOLS ==========

const metabolomicsTools: ToolDefinition[] = [
  // Peak Calling
  {
    id: 'metab-xcms',
    name: 'XCMS',
    type: 'analysis',
    tool: 'xcms',
    version: '3.22.0',
    omicsType: 'Metabolomics',
    operationCategory: 'Peak Calling',
    description: 'LC-MS peak detection, alignment, and integration',
    parameterTemplate: [
      { name: 'ppm', type: 'number', label: 'PPM', description: 'Mass accuracy (ppm)', default: 25, min: 1, max: 100, unit: 'ppm' },
      { name: 'peakwidth_min', type: 'number', label: 'Min Peak Width', description: 'Minimum peak width', default: 5, min: 1, max: 100, unit: 's' },
      { name: 'peakwidth_max', type: 'number', label: 'Max Peak Width', description: 'Maximum peak width', default: 20, min: 1, max: 100, unit: 's' },
      { name: 'snthresh', type: 'number', label: 'S/N Threshold', description: 'Signal-to-noise threshold', default: 10, min: 1, max: 100 },
    ],
    defaultParameters: { ppm: 25, peakwidth_min: 5, peakwidth_max: 20, snthresh: 10 },
  },
  {
    id: 'metab-mzmine',
    name: 'MZmine',
    type: 'analysis',
    tool: 'mzmine',
    version: '3.5.0',
    omicsType: 'Metabolomics',
    operationCategory: 'Peak Calling',
    description: 'Mass spectrometry data processing for metabolomics',
    parameterTemplate: [
      { name: 'noise_level', type: 'number', label: 'Noise Level', description: 'Noise level threshold', default: 1000, min: 0, max: 1000000 },
      { name: 'mz_tolerance', type: 'number', label: 'm/z Tolerance', description: 'm/z tolerance', default: 0.005, min: 0.001, max: 1, unit: 'Da' },
    ],
    defaultParameters: { noise_level: 1000, mz_tolerance: 0.005 },
  },

  // Analysis
  {
    id: 'metab-metaboanalyst',
    name: 'MetaboAnalyst',
    type: 'analysis',
    tool: 'metaboanalyst',
    version: '5.0',
    omicsType: 'Metabolomics',
    operationCategory: 'Statistical Analysis',
    description: 'Comprehensive metabolomics data analysis',
    parameterTemplate: [
      { name: 'normalization', type: 'select', label: 'Normalization', description: 'Normalization method', default: 'none', options: ['none', 'sum', 'median', 'quantile', 'pareto'] },
      { name: 'scaling', type: 'select', label: 'Scaling', description: 'Data scaling', default: 'auto', options: ['none', 'auto', 'pareto', 'range'] },
    ],
    defaultParameters: { normalization: 'none', scaling: 'auto' },
  },
];

// ========== MULTI-OMICS TOOLS ==========

const multiomicsTools: ToolDefinition[] = [
  // Integration
  {
    id: 'multi-mixomics',
    name: 'mixOmics',
    type: 'analysis',
    tool: 'mixomics',
    version: '6.24.0',
    omicsType: 'Multi-omics',
    operationCategory: 'Integration',
    description: 'Multi-omics data integration and analysis',
    parameterTemplate: [
      { name: 'method', type: 'select', label: 'Method', description: 'Integration method', default: 'DIABLO', options: ['PLS', 'sPLS', 'DIABLO', 'MOFA'] },
      { name: 'ncomp', type: 'number', label: 'Components', description: 'Number of components', default: 2, min: 1, max: 20 },
    ],
    defaultParameters: { method: 'DIABLO', ncomp: 2 },
  },
];

// ========== VISUALIZATION TOOLS ==========

const visualizationTools: ToolDefinition[] = [
  {
    id: 'viz-heatmap',
    name: 'Heatmap',
    type: 'visualization',
    tool: 'heatmap',
    version: '1.0',
    omicsType: 'General',
    operationCategory: 'Visualization',
    description: 'Generate heatmap visualization',
    parameterTemplate: [
      { name: 'clustering', type: 'select', label: 'Clustering', description: 'Clustering method', default: 'hierarchical', options: ['none', 'hierarchical', 'kmeans'] },
      { name: 'color_scheme', type: 'select', label: 'Color Scheme', description: 'Color palette', default: 'RdBu', options: ['RdBu', 'viridis', 'plasma', 'cividis'] },
    ],
    defaultParameters: { clustering: 'hierarchical', color_scheme: 'RdBu' },
  },
  {
    id: 'viz-volcano',
    name: 'Volcano Plot',
    type: 'visualization',
    tool: 'volcano',
    version: '1.0',
    omicsType: 'General',
    operationCategory: 'Visualization',
    description: 'Generate volcano plot for differential analysis',
    parameterTemplate: [
      { name: 'pvalue_threshold', type: 'number', label: 'P-value Threshold', description: 'Significance threshold', default: 0.05, min: 0, max: 1 },
      { name: 'fc_threshold', type: 'number', label: 'Fold Change Threshold', description: 'Log2 fold change threshold', default: 1, min: 0, max: 10 },
    ],
    defaultParameters: { pvalue_threshold: 0.05, fc_threshold: 1 },
  },
  {
    id: 'viz-pca',
    name: 'PCA Plot',
    type: 'visualization',
    tool: 'pca',
    version: '1.0',
    omicsType: 'General',
    operationCategory: 'Visualization',
    description: 'Principal Component Analysis plot',
    parameterTemplate: [
      { name: 'n_components', type: 'number', label: 'Components', description: 'Number of components', default: 2, min: 2, max: 10 },
      { name: 'scale', type: 'boolean', label: 'Scale Data', description: 'Scale data before PCA', default: true },
    ],
    defaultParameters: { n_components: 2, scale: true },
  },
  {
    id: 'viz-boxplot',
    name: 'Box Plot',
    type: 'visualization',
    tool: 'boxplot',
    version: '1.0',
    omicsType: 'General',
    operationCategory: 'Visualization',
    description: 'Generate box plot for data distribution',
    parameterTemplate: [
      { name: 'show_outliers', type: 'boolean', label: 'Show Outliers', description: 'Display outlier points', default: true },
      { name: 'notch', type: 'boolean', label: 'Notched', description: 'Draw notched box plot', default: false },
    ],
    defaultParameters: { show_outliers: true, notch: false },
  },
];

// ========== COMMON TOOLS (across omics types) ==========

const commonTools: ToolDefinition[] = [
  // Data Processing
  {
    id: 'common-filter-by-variance',
    name: 'Filter by Variance',
    type: 'filter',
    tool: 'variance_filter',
    version: '1.0',
    omicsType: 'General',
    operationCategory: 'Data Processing',
    description: 'Filter features by variance threshold',
    parameterTemplate: [
      { name: 'threshold', type: 'number', label: 'Variance Threshold', description: 'Minimum variance to keep', default: 0.01, min: 0, max: 100 },
      { name: 'percentile', type: 'number', label: 'Percentile', description: 'Keep top N% by variance', default: 100, min: 0, max: 100, unit: '%' },
    ],
    defaultParameters: { threshold: 0.01, percentile: 100 },
  },
  {
    id: 'common-missing-value',
    name: 'Handle Missing Values',
    type: 'transform',
    tool: 'missing_value',
    version: '1.0',
    omicsType: 'General',
    operationCategory: 'Data Processing',
    description: 'Handle missing values in dataset',
    parameterTemplate: [
      { name: 'method', type: 'select', label: 'Method', description: 'Imputation method', default: 'remove', options: ['remove', 'mean', 'median', 'knn', 'zero'] },
      { name: 'threshold', type: 'number', label: 'Threshold', description: 'Max missing % to keep feature', default: 20, min: 0, max: 100, unit: '%' },
    ],
    defaultParameters: { method: 'remove', threshold: 20 },
  },

  // Output
  {
    id: 'common-csv-output',
    name: 'CSV Output',
    type: 'output',
    tool: 'file-output',
    version: '1.0',
    omicsType: 'General',
    operationCategory: 'Output',
    description: 'Save results in CSV format',
    parameterTemplate: [
      { name: 'output_path', type: 'file', label: 'Output Path', description: 'Output CSV file path', required: true },
      { name: 'separator', type: 'select', label: 'Separator', description: 'Column separator', default: ',', options: [',', '\t', ';', '|'] },
      { name: 'include_header', type: 'boolean', label: 'Include Header', description: 'Write column names', default: true },
    ],
    defaultParameters: { separator: ',', include_header: true },
  },
  {
    id: 'common-report-output',
    name: 'HTML Report',
    type: 'output',
    tool: 'report-generator',
    version: '1.0',
    omicsType: 'General',
    operationCategory: 'Output',
    description: 'Generate comprehensive HTML report',
    parameterTemplate: [
      { name: 'output_path', type: 'file', label: 'Output Path', description: 'Output HTML file path', required: true },
      { name: 'include_qc', type: 'boolean', label: 'Include QC', description: 'Include QC metrics', default: true },
      { name: 'include_plots', type: 'boolean', label: 'Include Plots', description: 'Include visualization plots', default: true },
    ],
    defaultParameters: { include_qc: true, include_plots: true },
  },
];

// ========== EXPORT ALL TOOLS ==========

export const allTools: ToolDefinition[] = [
  ...genomicsTools,
  ...transcriptomicsTools,
  ...proteomicsTools,
  ...metabolomicsTools,
  ...multiomicsTools,
  ...visualizationTools,
  ...commonTools,
];

// Helper functions
export const getToolsByOmicsType = (omicsType: string): ToolDefinition[] => {
  return allTools.filter((tool) => tool.omicsType === omicsType || tool.omicsType === 'General');
};

export const getToolsByCategory = (category: string): ToolDefinition[] => {
  return allTools.filter((tool) => tool.operationCategory === category);
};

export const getToolById = (id: string): ToolDefinition | undefined => {
  return allTools.find((tool) => tool.id === id);
};

export const getOmicsTypes = (): string[] => {
  return Array.from(new Set(allTools.map((tool) => tool.omicsType)));
};

export const getCategoriesByOmicsType = (omicsType: string): string[] => {
  return Array.from(
    new Set(
      allTools
        .filter((tool) => tool.omicsType === omicsType || tool.omicsType === 'General')
        .map((tool) => tool.operationCategory)
    )
  );
};
