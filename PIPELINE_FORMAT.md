# Pipeline Format Specification

## Overview

Omicsomics uses a **standardized JSON format** for pipeline definitions. This enables:

- **Pipeline sharing** across users and projects
- **Cross-platform portability**
- **Version control** for workflows
- **Reproducible research** with documented workflows

## Core Concepts

### Pipeline vs Run

- **Pipeline**: A reusable workflow template
  - Defines sequence of analysis steps
  - Contains parameters and configurations
  - Can be shared publicly or kept private
  - Portable across projects and users
  - Supports visual editing
- **Run**: An execution instance of a pipeline
  - Uses specific pipeline configuration
  - Operates on specific data files
  - Tracks execution status and results
  - Stores complete pipeline config for reproducibility
  - Belongs to a specific project

## Standard Pipeline Format

### Complete Example

```json
{
  "format_version": "1.0.0",
  "pipeline": {
    "name": "RNA-seq Analysis Pipeline",
    "description": "Complete RNA-seq analysis from FASTQ to differential expression",
    "version": "2.1.0",
    "category": "transcriptomics",
    "author": "user@example.com",
    "is_public": true,
    "definition": {
      "nodes": [
        {
          "id": "node_1",
          "type": "input",
          "label": "Raw FASTQ Files",
          "data": {
            "accepts": ["fastq", "fastq.gz"],
            "required": true,
            "parameters": {
              "paired_end": true
            }
          },
          "position": { "x": 100, "y": 100 }
        },
        {
          "id": "node_2",
          "type": "process",
          "label": "Quality Control",
          "data": {
            "tool_config": {
              "tool_id": "fastqc",
              "tool_version": "0.12.1",
              "parameters": {
                "threads": 4
              },
              "input_mappings": {
                "fastq_files": "node_1.output"
              },
              "output_mappings": {
                "html_report": "qc_report"
              },
              "enabled": true
            }
          },
          "position": { "x": 300, "y": 100 }
        },
        {
          "id": "node_3",
          "type": "filter",
          "label": "Adapter Trimming",
          "data": {
            "tool": "Trimmomatic",
            "version": "0.39",
            "parameters": {
              "leading": 3,
              "trailing": 3,
              "minlen": 36
            }
          },
          "position": { "x": 500, "y": 100 }
        },
        {
          "id": "node_4",
          "type": "process",
          "label": "Alignment",
          "data": {
            "tool_config": {
              "tool_id": "star",
              "tool_version": "2.7.11a",
              "parameters": {
                "threads": 8,
                "out_sam_type": "BAM SortedByCoordinate",
                "quant_mode": "GeneCounts"
              },
              "input_mappings": {
                "fastq_1": "node_3.output.read1",
                "fastq_2": "node_3.output.read2",
                "genome_index": "reference.genome_index"
              },
              "output_mappings": {
                "aligned_bam": "aligned_reads",
                "gene_counts": "raw_counts"
              },
              "resources": {
                "min_cpu": 8,
                "max_cpu": 16,
                "min_memory_gb": 32,
                "max_memory_gb": 64
              },
              "enabled": true
            }
          },
          "position": { "x": 700, "y": 100 }
        },
        {
          "id": "node_5",
          "type": "analysis",
          "label": "Differential Expression",
          "data": {
            "tool": "DESeq2",
            "version": "1.34.0",
            "parameters": {
              "alpha": 0.05,
              "lfcThreshold": 0
            }
          },
          "position": { "x": 900, "y": 100 }
        },
        {
          "id": "node_6",
          "type": "output",
          "label": "Results",
          "data": {
            "outputs": ["counts", "de_genes", "plots"],
            "format": "csv"
          },
          "position": { "x": 1100, "y": 100 }
        }
      ],
      "edges": [
        {
          "id": "edge_1",
          "source": "node_1",
          "target": "node_2",
          "sourceHandle": null,
          "targetHandle": null
        },
        {
          "id": "edge_2",
          "source": "node_2",
          "target": "node_3",
          "sourceHandle": null,
          "targetHandle": null
        },
        {
          "id": "edge_3",
          "source": "node_3",
          "target": "node_4",
          "sourceHandle": null,
          "targetHandle": null
        },
        {
          "id": "edge_4",
          "source": "node_4",
          "target": "node_5",
          "sourceHandle": null,
          "targetHandle": null
        },
        {
          "id": "edge_5",
          "source": "node_5",
          "target": "node_6",
          "sourceHandle": null,
          "targetHandle": null
        }
      ],
      "parameters": {
        "reference_genome": "GRCh38",
        "annotation_gtf": "gencode.v38.annotation.gtf",
        "comparison_groups": ["treatment", "control"]
      }
    },
    "metadata": {
      "original_author": "admin@omicsomics.org",
      "created_at": "2025-01-09T10:00:00Z",
      "updated_at": "2025-01-09T12:00:00Z",
      "tags": ["RNA-seq", "differential expression", "validated"],
      "citations": [
        "DOI:10.1093/bioinformatics/btp120",
        "DOI:10.1186/s13059-014-0550-8"
      ]
    }
  }
}
```

## Node Types

### 1. Input Node

Represents data input to the pipeline.

```json
{
  "id": "input_1",
  "type": "input",
  "label": "Raw Data",
  "data": {
    "accepts": ["vcf", "bam", "fastq"],
    "required": true,
    "multiple": false,
    "parameters": {}
  },
  "position": { "x": 100, "y": 100 }
}
```

### 2. Process Node

Data processing step.

```json
{
  "id": "process_1",
  "type": "process",
  "label": "Quality Control",
  "data": {
    "tool": "FastQC",
    "version": "0.11.9",
    "container": "biocontainers/fastqc:0.11.9",
    "parameters": {
      "threads": 4
    }
  },
  "position": { "x": 300, "y": 100 }
}
```

### 3. Filter Node

Data filtering/selection.

```json
{
  "id": "filter_1",
  "type": "filter",
  "label": "Quality Filter",
  "data": {
    "tool": "custom",
    "parameters": {
      "min_quality": 30,
      "remove_duplicates": true
    }
  },
  "position": { "x": 500, "y": 100 }
}
```

### 4. Transform Node

Data transformation/conversion.

```json
{
  "id": "transform_1",
  "type": "transform",
  "label": "Normalize",
  "data": {
    "tool": "DESeq2",
    "method": "rlog",
    "parameters": {
      "blind": false
    }
  },
  "position": { "x": 700, "y": 100 }
}
```

### 5. Analysis Node

Statistical analysis/computation.

```json
{
  "id": "analysis_1",
  "type": "analysis",
  "label": "Differential Expression",
  "data": {
    "tool": "limma",
    "version": "3.50.0",
    "parameters": {
      "lfc_threshold": 1,
      "adj_pvalue": 0.05
    }
  },
  "position": { "x": 900, "y": 100 }
}
```

### 6. Output Node

Results output.

```json
{
  "id": "output_1",
  "type": "output",
  "label": "Results",
  "data": {
    "outputs": ["results.csv", "plots.pdf"],
    "format": "csv",
    "compress": true
  },
  "position": { "x": 1100, "y": 100 }
}
```

## Edges

Edges define connections between nodes:

```json
{
  "id": "edge_1",
  "source": "node_1",
  "target": "node_2",
  "type": "default",
  "sourceHandle": null,
  "targetHandle": null,
  "label": "FASTQ files"
}
```

### Special Edge Types

- `default`: Standard connection
- `bridge`: Auto-generated when merging pipelines
- `conditional`: Conditional execution (planned)
- `parallel`: Parallel execution (planned)

## Global Parameters

Pipeline-wide parameters:

```json
{
  "parameters": {
    "reference_genome": "GRCh38",
    "organism": "Homo sapiens",
    "threads": 8,
    "memory_gb": 32,
    "temp_dir": "/tmp",
    "output_dir": "/results"
  }
}
```

## API Operations

### Export Pipeline

```bash
# Export a pipeline
curl -H "Authorization: Bearer $TOKEN" \
  "http://localhost:8001/api/v1/custom-pipelines/1/export" \
  > my_pipeline.json
```

### Import Pipeline

```bash
# Import a pipeline
curl -X POST "http://localhost:8001/api/v1/custom-pipelines/import" \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d @my_pipeline.json
```

Response:

```json
{
  "id": 42,
  "name": "RNA-seq Analysis Pipeline",
  "message": "Pipeline imported successfully",
  "original_author": "user@example.com"
}
```

### Create Pipeline

```bash
curl -X POST "http://localhost:8001/api/v1/custom-pipelines/" \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{
    "name": "My Pipeline",
    "description": "Custom analysis workflow",
    "category": "genomics",
    "is_public": false,
    "definition": {
      "nodes": [...],
      "edges": [...],
      "parameters": {}
    }
  }'
```

### Merge Pipelines

```bash
# Merge multiple pipelines
curl -X POST "http://localhost:8001/api/v1/custom-pipelines/merge" \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{
    "pipeline_ids": [1, 2, 3]
  }'
```

Response:

```json
{
  "merged_definition": {
    "nodes": [...],  // All nodes with remapped IDs
    "edges": [...],  // Original edges + bridge edges
    "parameters": {},
    "metadata": {
      "merged_from": [1, 2, 3],
      "merge_strategy": "sequential"
    }
  },
  "source_count": 3,
  "total_nodes": 15,
  "total_edges": 14
}
```

## Pipeline Categories

Predefined categories:

- `genomics` - WGS, WES, variant calling
- `transcriptomics` - RNA-seq, expression analysis
- `proteomics` - Protein identification, quantification
- `metabolomics` - Metabolite profiling
- `epigenomics` - ChIP-seq, ATAC-seq, methylation
- `single_cell` - scRNA-seq, scATAC-seq
- `multi_omics` - Cross-omics integration
- `gwas` - Association studies
- `custom` - User-defined workflows

## Validation Rules

### Required Fields

- `name` - Pipeline name (string)
- `definition` - Pipeline definition (object)
- `definition.nodes` - At least one node (array)
- `definition.edges` - Edges list (array, can be empty)

### Node Validation

- Unique `id` within pipeline
- Valid `type` (input, process, filter, transform, analysis, output)
- `data` object present

### Edge Validation

- `source` and `target` must reference existing nodes
- No circular dependencies (DAG required)
- At least one path from input to output

## Use Cases

### 1. Share Pipeline Publicly

```python
# Create and share
pipeline = await create_custom_pipeline(
    db, name="QC Pipeline",
    description="Standard QC workflow",
    definition=qc_definition,
    owner_id=user.id,
    is_public=True  # Share publicly
)
```

### 2. Import Community Pipeline

```bash
# Download from repository
wget https://pipelines.omicsomics.org/rnaseq-basic.json

# Import
curl -X POST .../custom-pipelines/import -d @rnaseq-basic.json
```

### 3. Version Control

```bash
# Track pipeline in git
git add pipelines/my_pipeline_v1.0.json
git commit -m "Initial pipeline version"

# Update and tag
git add pipelines/my_pipeline_v1.1.json
git tag -a v1.1 -m "Added normalization step"
```

### 4. Merge Workflows

```python
# Combine QC + alignment + analysis
qc_pipeline = get_pipeline(1)
align_pipeline = get_pipeline(2)
analysis_pipeline = get_pipeline(3)

merged = merge_pipelines([
    qc_pipeline,
    align_pipeline,
    analysis_pipeline
])
```

## Tool Configuration

### Tool Config Structure

Each process node can specify a `tool_config` object in its `data` field:

```json
{
  "tool_config": {
    "tool_id": "star", // Tool identifier from registry
    "tool_version": "2.7.11a", // Specific version (optional)
    "parameters": {
      // Tool-specific parameters
      "threads": 8,
      "out_sam_type": "BAM SortedByCoordinate"
    },
    "input_mappings": {
      // Map inputs to data sources
      "fastq_1": "previous_node.output",
      "genome_index": "reference.index"
    },
    "output_mappings": {
      // Map outputs to destinations
      "aligned_bam": "aligned_reads",
      "gene_counts": "expression_matrix"
    },
    "resources": {
      // Resource requirements
      "min_cpu": 8,
      "max_cpu": 16,
      "min_memory_gb": 32,
      "max_memory_gb": 64,
      "gpu_required": false
    },
    "environment": {
      // Environment variables
      "TMPDIR": "/tmp",
      "OMP_NUM_THREADS": "8"
    },
    "enabled": true // Toggle tool execution
  }
}
```

### Available Tools

Query available tools via API:

```bash
# List all tools
GET /api/v1/tools/

# Search tools
GET /api/v1/tools/search?q=alignment

# Get tool details
GET /api/v1/tools/star

# Get tool parameters
GET /api/v1/tools/star/parameters
```

### Tool Categories

- **alignment**: BWA, STAR, Bowtie2, HISAT2
- **variant_calling**: GATK, FreeBayes, VarScan, Strelka
- **quality_control**: FastQC, MultiQC, Qualimap
- **quantification**: Salmon, Kallisto, featureCounts, HTSeq
- **differential_expression**: DESeq2, edgeR, limma
- **single_cell**: Seurat, Scanpy, Cell Ranger, Monocle
- **peak_calling**: MACS2, HOMER, SICER
- **functional_analysis**: GSEA, GO enrichment, KEGG pathway

### Parameter Types

Tools define parameter types for validation:

- **string**: Text values
- **integer**: Whole numbers (with min/max)
- **float**: Decimal numbers (with min/max)
- **boolean**: true/false flags
- **enum**: Select from predefined values
- **file**: File input (with extension validation)
- **directory**: Directory path
- **array**: Multiple values

### Validation

Validate pipelines before execution:

```bash
# Validate pipeline definition
POST /api/v1/custom-pipelines/validate
{
  "nodes": [...],
  "edges": [...],
  "parameters": {...}
}

# Validate existing pipeline
GET /api/v1/custom-pipelines/{id}/validate
```

Returns validation results:

```json
{
  "valid": true,
  "errors": [],
  "warnings": [
    "Node star_align: Input 'fastq_2' not mapped (optional for single-end)"
  ]
}
```

### Dependencies

Check pipeline dependencies:

```bash
GET /api/v1/custom-pipelines/{id}/dependencies
```

Returns:

```json
{
  "tools": [
    {
      "name": "STAR",
      "version": "2.7.11a",
      "docker_image": "quay.io/biocontainers/star:2.7.11a",
      "category": "alignment"
    },
    {
      "name": "DESeq2",
      "version": "1.40.0",
      "docker_image": "bioconductor/bioconductor_docker:RELEASE_3_17",
      "category": "differential_expression"
    }
  ],
  "formats": ["fastq", "bam", "csv"],
  "resources": {
    "min_cpu": 10,
    "max_cpu": 32,
    "min_memory_gb": 40,
    "max_memory_gb": 96,
    "gpu_required": false
  },
  "node_count": 5
}
```

## Best Practices

### 1. Naming

- Use descriptive names: "RNA-seq DE Analysis" not "Pipeline 1"
- Include version in name if multiple versions exist
- Use consistent naming scheme across pipelines

### 2. Documentation

- Provide clear description
- Document all parameters
- List required inputs and expected outputs
- Include citations for methods

### 3. Parameters

- Use sensible defaults from tool registry
- Document parameter ranges
- Validate parameter values before execution
- Use typed parameters for better validation

### 4. Modularity

- Keep pipelines focused (5-10 nodes ideal)
- Create reusable sub-pipelines
- Use merge for complex workflows
- Avoid deep nesting

### 5. Tool Selection

- Use registry tools when available (tested and validated)
- Specify tool versions for reproducibility
- Check resource requirements match available resources
- Validate tool configurations before saving pipeline

### 6. Input/Output Mapping

- Explicitly map all required inputs
- Use descriptive mapping names
- Check data format compatibility between nodes
- Document expected input formats in pipeline description

### 5. Metadata

- Track tool versions
- Document data sources
- Include authorship info
- Add relevant tags

## Future Enhancements

- [ ] Conditional branching
- [ ] Parallel execution paths
- [ ] Subworkflows/nested pipelines
- [ ] Dynamic parameters
- [ ] Resource requirements (CPU, RAM)
- [ ] Cost estimation
- [ ] Execution time prediction
- [ ] Pipeline validation service
- [ ] Pipeline marketplace
- [ ] Automated testing framework

## See Also

- [Data Format](DATA_FORMAT.md) - Unified data format specification
- [User Guide](docs/PIPELINE_BUILDER.md) - Visual pipeline builder guide
- [API Documentation](http://localhost:8001/docs) - Complete API reference
