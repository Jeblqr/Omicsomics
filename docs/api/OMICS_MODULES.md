# 组学分析模块 API 文档

本文档详细介绍了 Omicsomics 平台中各组学分析模块的 API 端点及其使用方法。

## 目录

- [基因组学分析](#基因组学分析-wgswes)
- [转录组学分析](#转录组学分析-bulk-rna-seq)
- [单细胞分析](#单细胞分析-scrna-seq)
- [表观组学分析](#表观组学分析-chip-seqatac-seq)
- [可视化模块](#可视化模块)

---

## 基因组学分析 (WGS/WES)

**基础路径**: `/api/v1/genomics`

### 1. 质量控制

**端点**: `POST /genomics/qc`

**描述**: 对 FASTQ 文件运行 FastQC 质量控制。

**请求体**:

```json
{
  "sample_id": 1,
  "fastq_files": ["/path/to/sample_R1.fastq.gz", "/path/to/sample_R2.fastq.gz"],
  "output_dir": "/path/to/qc_results"
}
```

**响应**:

```json
{
  "workflow_id": 123,
  "status": "queued",
  "message": "FastQC analysis queued for execution"
}
```

### 2. 接头修剪

**端点**: `POST /genomics/trim`

**描述**: 使用 fastp 或 Trimmomatic 进行接头修剪。

**请求体**:

```json
{
  "sample_id": 1,
  "fastq_files": ["/path/to/sample_R1.fastq.gz", "/path/to/sample_R2.fastq.gz"],
  "output_dir": "/path/to/trimmed",
  "tool": "fastp",
  "params": {
    "qualified_quality_phred": 20,
    "unqualified_percent_limit": 40,
    "length_required": 50
  }
}
```

**fastp 参数示例**:

- `qualified_quality_phred`: 质量值阈值 (默认 20)
- `length_required`: 最小读长 (默认 15)
- `cut_front`: 切除前端低质量碱基
- `cut_tail`: 切除尾端低质量碱基

### 3. 序列比对

**端点**: `POST /genomics/align`

**描述**: 将 FASTQ 文件比对到参考基因组。

**请求体**:

```json
{
  "sample_id": 1,
  "fastq_files": [
    "/path/to/trimmed_R1.fastq.gz",
    "/path/to/trimmed_R2.fastq.gz"
  ],
  "reference_genome": "/path/to/reference/hg38.fa",
  "output_bam": "/path/to/output/aligned.bam",
  "aligner": "bwa-mem",
  "threads": 16
}
```

**支持的比对工具**:

- `bwa-mem`: 适用于短读长 (Illumina)
- `bowtie2`: 适用于 ChIP-seq, RNA-seq
- `minimap2`: 适用于长读长 (PacBio, Nanopore)

### 4. 变异检测

**端点**: `POST /genomics/variant-calling`

**描述**: 从比对的 BAM 文件中检测变异。

**请求体**:

```json
{
  "sample_id": 1,
  "input_bam": "/path/to/aligned.bam",
  "reference_genome": "/path/to/reference/hg38.fa",
  "output_vcf": "/path/to/output/variants.vcf",
  "caller": "gatk4"
}
```

**支持的变异检测工具**:

- `gatk4`: GATK4 HaplotypeCaller (适用于 germline 变异)
- `freebayes`: 贝叶斯变异检测
- `deepvariant`: 基于深度学习的变异检测

### 5. 变异注释

**端点**: `POST /genomics/annotate-variants`

**描述**: 为变异添加功能注释信息。

**请求体**:

```json
{
  "sample_id": 1,
  "input_vcf": "/path/to/variants.vcf",
  "output_vcf": "/path/to/annotated_variants.vcf",
  "annotator": "vep"
}
```

**支持的注释工具**:

- `vep`: Ensembl Variant Effect Predictor (推荐)
- `snpeff`: 快速功能注释
- `annovar`: 综合注释工具

**VEP 提供的注释信息**:

- 变异后果 (Consequence): missense, synonymous, frameshift 等
- 基因和转录本: HGNC 基因名称, Ensembl 转录本 ID
- 蛋白质改变: HGVS 蛋白质命名
- 群体频率: gnomAD, 1000 Genomes
- 致病性预测: SIFT, PolyPhen, CADD
- 临床相关性: ClinVar 注释

### 6. 完整流水线

**端点**: `POST /genomics/complete-pipeline`

**描述**: 运行完整的基因组学分析流水线（QC → Trim → Align → Call → Annotate）。

**请求体**:

```json
{
  "sample_id": 1,
  "fastq_files": ["/path/to/sample_R1.fastq.gz", "/path/to/sample_R2.fastq.gz"],
  "reference_genome": "/path/to/reference/hg38.fa",
  "output_dir": "/path/to/output",
  "skip_qc": false,
  "skip_trimming": false,
  "aligner": "bwa-mem",
  "variant_caller": "gatk4",
  "annotator": "vep",
  "threads": 16
}
```

---

## 转录组学分析 (bulk RNA-seq)

**基础路径**: `/api/v1/transcriptomics`

### 1. 转录本定量

**端点**: `POST /transcriptomics/quantify`

**描述**: 运行 RNA-seq 定量分析。

**请求体**:

```json
{
  "sample_id": 1,
  "fastq_files": ["/path/to/sample_R1.fastq.gz", "/path/to/sample_R2.fastq.gz"],
  "reference_index": "/path/to/salmon_index",
  "output_dir": "/path/to/quant_results",
  "tool": "salmon",
  "threads": 8
}
```

**支持的定量工具**:

#### Salmon (推荐 - 快速、准确)

- 无需比对，使用伪比对
- 自动检测文库类型
- 输出: `quant.sf` (TPM, counts)

#### Kallisto

- 超快速伪比对
- 适合大规模数据集
- 输出: `abundance.tsv`

#### STAR

- 传统比对 + 定量
- 适合发现新转录本
- 输出: BAM 文件 + `ReadsPerGene.out.tab`

#### HISAT2

- 快速比对工具
- 适合转录组和基因组比对

### 2. 生成 Count 矩阵

**端点**: `POST /transcriptomics/count-matrix`

**描述**: 从 BAM 文件生成基因水平的 count 矩阵。

**请求体**:

```json
{
  "sample_id": 1,
  "input_bam": "/path/to/aligned.bam",
  "gtf_file": "/path/to/annotation.gtf",
  "output_file": "/path/to/counts.txt",
  "threads": 4
}
```

**featureCounts 特性**:

- 快速准确的基因水平计数
- 支持多样本批处理
- 自动处理多映射读长
- 支持链特异性文库

### 3. 差异表达分析

**端点**: `POST /transcriptomics/differential-expression`

**描述**: 识别不同条件下的差异表达基因。

**请求体**:

```json
{
  "sample_id": 1,
  "counts_matrix": "/path/to/counts.csv",
  "sample_metadata": "/path/to/metadata.csv",
  "output_dir": "/path/to/de_results",
  "tool": "deseq2",
  "design_formula": "~ condition"
}
```

**样本元数据格式** (`metadata.csv`):

```csv
sample,condition,batch
sample1,control,1
sample2,control,1
sample3,treated,1
sample4,treated,1
```

**DESeq2 输出**:

- `deseq2_results.csv`: 差异表达结果 (log2FC, p-value, padj)
- `MA_plot.pdf`: MA 图
- `PCA_plot.pdf`: PCA 图
- `volcano_plot.pdf`: 火山图

**结果表列说明**:

- `baseMean`: 平均表达量
- `log2FoldChange`: log2 倍数变化
- `lfcSE`: log2FC 标准误
- `pvalue`: 原始 p 值
- `padj`: FDR 校正后的 p 值

### 4. 富集分析

**端点**: `GET /transcriptomics/enrichment/gsea`

**描述**: 基因集富集分析 (GSEA)。

**状态**: 规划中 - 将集成 GSEApy 或 Enrichr API

---

## 单细胞分析 (scRNA-seq)

**基础路径**: `/api/v1/singlecell`

### 1. Cell Ranger 定量

**端点**: `POST /singlecell/cellranger`

**描述**: 运行 10x Genomics Cell Ranger 流水线。

**请求体**:

```json
{
  "sample_id": 1,
  "fastq_dir": "/path/to/fastq_dir",
  "sample_name": "sample1",
  "transcriptome": "/path/to/cellranger_ref/refdata-gex-GRCh38-2020-A",
  "output_dir": "/path/to/output",
  "chemistry": "auto",
  "threads": 16
}
```

**Cell Ranger 输出**:

- `filtered_feature_bc_matrix.h5`: 过滤后的特征-条码矩阵
- `raw_feature_bc_matrix.h5`: 原始矩阵
- `web_summary.html`: QC 报告
- `metrics_summary.csv`: 指标汇总

**Chemistry 选项**:

- `auto`: 自动检测 (推荐)
- `SC3Pv3`: 3' v3 chemistry
- `SC5P-PE`: 5' paired-end
- `SC5P-R2`: 5' R2-only

### 2. Scanpy 预处理

**端点**: `POST /singlecell/preprocess`

**描述**: 运行 Scanpy 标准预处理流水线。

**请求体**:

```json
{
  "sample_id": 1,
  "input_h5ad": "/path/to/raw.h5ad",
  "output_h5ad": "/path/to/processed.h5ad",
  "min_genes": 200,
  "min_cells": 3,
  "max_pct_mt": 20.0,
  "n_top_genes": 2000,
  "n_neighbors": 10,
  "leiden_resolution": 0.5
}
```

**预处理步骤**:

1. **质控过滤**

   - `min_genes`: 每个细胞最少检测到的基因数 (默认 200)
   - `min_cells`: 每个基因最少表达的细胞数 (默认 3)
   - `max_pct_mt`: 最大线粒体基因比例 (默认 20%)

2. **标准化**

   - Total count normalization (target_sum = 10,000)
   - Log1p transformation

3. **高变基因识别**

   - `n_top_genes`: 高变基因数量 (默认 2000)

4. **降维**

   - PCA: 50 个主成分
   - UMAP: 2D 嵌入

5. **聚类**
   - Leiden 算法
   - `leiden_resolution`: 分辨率参数 (0.1-2.0)

### 3. Seurat 批次整合

**端点**: `POST /singlecell/integrate`

**描述**: 使用 Seurat 进行批次校正和数据整合。

**请求体**:

```json
{
  "sample_id": 1,
  "input_h5_files": [
    "/path/to/batch1.h5",
    "/path/to/batch2.h5",
    "/path/to/batch3.h5"
  ],
  "output_rds": "/path/to/integrated.rds",
  "batch_key": "batch"
}
```

**Seurat Integration 特性**:

- 识别锚点细胞 (anchor cells)
- 校正批次效应
- 保留生物学变异
- 适用于多实验整合

**适用场景**:

- 不同批次的相同样本类型
- 不同技术平台的数据
- 不同时间点的样本

### 4. 细胞类型注释

**端点**: `POST /singlecell/annotate`

**描述**: 基于标记基因注释细胞类型。

**请求体**:

```json
{
  "sample_id": 1,
  "input_h5ad": "/path/to/processed.h5ad",
  "output_h5ad": "/path/to/annotated.h5ad",
  "marker_genes": {
    "T_cells": ["CD3D", "CD3E", "CD3G"],
    "B_cells": ["CD19", "MS4A1", "CD79A"],
    "NK_cells": ["NKG7", "GNLY", "NCAM1"],
    "Monocytes": ["CD14", "FCGR3A", "LYZ"],
    "Dendritic_cells": ["FCER1A", "CD1C"],
    "Plasma_cells": ["IGHG1", "MZB1"],
    "Epithelial_cells": ["EPCAM", "KRT8", "KRT18"]
  }
}
```

**常用标记基因库**:

#### 免疫细胞

- **T cells**: CD3D, CD3E, CD4, CD8A
- **B cells**: CD19, MS4A1 (CD20), CD79A
- **NK cells**: NKG7, GNLY, NCAM1 (CD56)
- **Monocytes**: CD14, FCGR3A (CD16)
- **Macrophages**: CD68, MSR1, MRC1 (CD206)

#### 脑细胞

- **Neurons**: RBFOX3 (NeuN), SYP, SNAP25
- **Astrocytes**: GFAP, AQP4, SLC1A2
- **Oligodendrocytes**: MOG, MBP, PLP1
- **Microglia**: TMEM119, P2RY12, CX3CR1

---

## 表观组学分析 (ChIP-seq/ATAC-seq)

**基础路径**: `/api/v1/epigenomics`

### 1. 序列比对

**端点**: `POST /epigenomics/align`

**描述**: 将 ChIP-seq 或 ATAC-seq 的 FASTQ 文件比对到参考基因组。

**请求体**:

```json
{
  "sample_id": 1,
  "input_files": [
    "/path/to/chip_R1.fastq.gz",
    "/path/to/chip_R2.fastq.gz"
  ],
  "reference_genome": "/path/to/bowtie2_index/hg38",
  "output_bam": "/path/to/output/aligned.bam",
  "aligner": "bowtie2",
  "threads": 8
}
```

**支持的比对工具**:

- `bowtie2`: 推荐用于 ChIP-seq/ATAC-seq (默认)
- `bwa`: 可选的比对工具

**响应**:

```json
{
  "workflow_id": 123,
  "status": "queued",
  "message": "bowtie2 alignment queued for execution"
}
```

### 2. Peak Calling

**端点**: `POST /epigenomics/peak-calling`

**描述**: 使用 MACS2/MACS3 从比对的 BAM 文件中检测 peaks。

**请求体 (带 control)**:

```json
{
  "sample_id": 1,
  "treatment_bam": "/path/to/chip_treatment.bam",
  "control_bam": "/path/to/chip_control.bam",
  "output_dir": "/path/to/peaks",
  "peak_caller": "macs2",
  "genome_size": "hs",
  "broad": false,
  "q_value": 0.05
}
```

**请求体 (不带 control, 用于 ATAC-seq)**:

```json
{
  "sample_id": 1,
  "treatment_bam": "/path/to/atac.bam",
  "output_dir": "/path/to/peaks",
  "peak_caller": "macs2",
  "genome_size": "hs",
  "nomodel": true,
  "shift": -100,
  "extsize": 200
}
```

**参数说明**:

- `genome_size`: 基因组大小
  - `"hs"`: Human - effective genome size 2.7e9 (MACS2 default)
  - `"mm"`: Mouse - effective genome size 1.87e9 (MACS2 default)
  - Note: These are effective genome sizes used by MACS2, not actual genome sizes
  - 或直接提供整数值
- `broad`: 是否检测 broad peaks (用于组蛋白标记如 H3K27me3)
- `q_value`: Q-value (FDR) 阈值
- `p_value`: P-value 阈值
- `nomodel`: 绕过模型构建 (ATAC-seq 推荐)
- `shift`: 移位大小
- `extsize`: 扩展大小

**响应**:

```json
{
  "workflow_id": 124,
  "status": "queued",
  "message": "macs2 peak calling queued for execution"
}
```

**输出文件**:

- `peaks_peaks.narrowPeak`: Narrow peaks (默认)
- `peaks_peaks.broadPeak`: Broad peaks (如果使用 --broad)
- `peaks_summits.bed`: Peak 顶点位置
- `peaks_peaks.xls`: 详细结果表格

### 3. 基序分析 (Motif Analysis)

**端点**: `POST /epigenomics/motif-analysis`

**描述**: 在 peak 区域中查找富集的 DNA 基序。

**请求体**:

```json
{
  "sample_id": 1,
  "peak_file": "/path/to/peaks_peaks.narrowPeak",
  "genome_fasta": "/path/to/reference/hg38.fa",
  "output_dir": "/path/to/motif_results",
  "tool": "homer",
  "motif_length": [8, 10, 12]
}
```

**支持的工具**:

- `homer`: HOMER findMotifsGenome.pl (推荐)
- `meme`: MEME Suite (规划中)

**响应**:

```json
{
  "workflow_id": 125,
  "status": "queued",
  "message": "homer motif analysis queued for execution"
}
```

**输出文件**:

- `homerResults.html`: 新发现基序的 HTML 报告
- `knownResults.html`: 已知基序匹配的 HTML 报告
- `homerMotifs.all.motifs`: 所有发现的基序

### 4. BigWig 生成

**端点**: `POST /epigenomics/bigwig`

**描述**: 生成用于基因组浏览器可视化的 BigWig 文件。

**请求体**:

```json
{
  "sample_id": 1,
  "input_bam": "/path/to/aligned.bam",
  "output_bigwig": "/path/to/coverage.bw",
  "genome_sizes": "/path/to/hg38.chrom.sizes",
  "normalize": true
}
```

**参数说明**:

- `normalize`: 是否进行 RPKM 归一化 (默认 true)
- `genome_sizes`: 染色体大小文件 (可从 UCSC 获取)

**响应**:

```json
{
  "workflow_id": 126,
  "status": "queued",
  "message": "BigWig generation queued for execution"
}
```

### 5. 完整流水线

**端点**: `POST /epigenomics/complete-pipeline`

**描述**: 运行完整的表观组学分析流水线，包括比对、peak calling、基序分析和 BigWig 生成。

**请求体**:

```json
{
  "sample_id": 1,
  "treatment_fastq": [
    "/path/to/treatment_R1.fastq.gz",
    "/path/to/treatment_R2.fastq.gz"
  ],
  "control_fastq": [
    "/path/to/control_R1.fastq.gz",
    "/path/to/control_R2.fastq.gz"
  ],
  "reference_genome": "/path/to/bowtie2_index/hg38",
  "genome_fasta": "/path/to/reference/hg38.fa",
  "genome_sizes": "/path/to/hg38.chrom.sizes",
  "output_dir": "/path/to/epigenomics_results",
  "aligner": "bowtie2",
  "peak_caller": "macs2",
  "genome_size": "hs",
  "run_motif_analysis": true,
  "generate_bigwig": true,
  "threads": 8
}
```

**响应**:

```json
{
  "workflow_id": 127,
  "status": "created",
  "message": "Complete epigenomics pipeline created. Use workflow engine for execution."
}
```

**注意**: 完整流水线需要使用工作流引擎 (如 Nextflow) 来协调各个步骤的顺序执行。

---

## 可视化模块

**基础路径**: `/api/v1/visualizations`

### 1. 火山图 (Volcano Plot)

**端点**: `POST /visualizations/volcano`

**描述**: 生成差异表达火山图数据。

**请求体**:

```json
{
  "results_file": "/path/to/deseq2_results.csv",
  "log2fc_col": "log2FoldChange",
  "pvalue_col": "pvalue",
  "padj_col": "padj",
  "gene_col": "gene",
  "fc_threshold": 1.0,
  "pval_threshold": 0.05
}
```

**响应示例**:

```json
{
  "x": [0.5, 1.2, -0.8, 2.3, ...],
  "y": [1.5, 3.2, 2.1, 5.4, ...],
  "genes": ["GENE1", "GENE2", "GENE3", ...],
  "colors": ["not_sig", "up", "down", "up", ...],
  "layout": {
    "title": "Volcano Plot",
    "xaxis": {"title": "Log2 Fold Change"},
    "yaxis": {"title": "-Log10(P-value)"}
  },
  "counts": {
    "up": 234,
    "down": 189,
    "not_sig": 12345
  }
}
```

### 2. PCA 图

**端点**: `POST /visualizations/pca`

**描述**: 生成 PCA 降维图数据。

**查询参数**:

- `h5ad_file`: AnnData 文件路径
- `color_by`: 着色依据的列名 (leiden, cell_type 等)
- `n_components`: 成分数量 (2 或 3)

**示例**:

```bash
curl -X POST "/api/v1/visualizations/pca?h5ad_file=/path/to/data.h5ad&color_by=leiden&n_components=2"
```

### 3. UMAP 图

**端点**: `POST /visualizations/umap`

**描述**: 生成 UMAP 降维图数据。

**请求体**:

```json
{
  "h5ad_file": "/path/to/processed.h5ad",
  "color_by": "cell_type"
}
```

**响应**:

```json
{
  "umap": [[1.2, 3.4], [2.1, 4.5], ...],
  "metadata": {
    "cell_type": ["T_cells", "B_cells", ...]
  },
  "cell_ids": ["AAACCTGAGCATCATC-1", ...],
  "layout": {...},
  "cluster_counts": {
    "T_cells": 1234,
    "B_cells": 567
  }
}
```

### 4. 热图 (Heatmap)

**端点**: `POST /visualizations/heatmap`

**描述**: 生成基因表达热图数据。

**请求体**:

```json
{
  "h5ad_file": "/path/to/annotated.h5ad",
  "gene_list": ["CD3D", "CD19", "CD14", "NKG7", "EPCAM"],
  "groupby": "cell_type",
  "standard_scale": "var"
}
```

### 5. IGV 基因组浏览器

**端点**: `GET /visualizations/igv`

**描述**: 获取 IGV.js 配置。

**查询参数**:

- `bam_file`: BAM 文件路径
- `region`: 基因组区域 (例如: chr1:1000-2000)
- `reference_genome`: 参考基因组 (hg38, hg19, mm10)

**示例**:

```bash
curl -X GET "/api/v1/visualizations/igv?bam_file=/path/to/aligned.bam&region=chr1:1000000-2000000&reference_genome=hg38"
```

### 6. QC 指标

**端点**: `GET /visualizations/quality-metrics`

**描述**: 获取单细胞 QC 指标分布。

**查询参数**:

- `h5ad_file`: AnnData 文件路径
- `metrics`: 指标列表 (n_genes_by_counts, total_counts, pct_counts_mt)

### 7. 支持的格式

**端点**: `GET /visualizations/formats`

**描述**: 获取支持的文件格式和可视化类型列表。

**响应**:

```json
{
  "formats": {
    "h5ad": {
      "description": "AnnData format for single-cell data",
      "visualizations": ["umap", "pca", "heatmap", "quality-metrics"]
    },
    "csv": {
      "description": "CSV/TSV tables for differential expression",
      "visualizations": ["volcano", "scatter"]
    },
    "bam": {
      "description": "Binary alignment format",
      "visualizations": ["igv", "coverage"]
    }
  }
}
```

---

## 工作流状态查询

所有分析任务提交后都会返回 `workflow_id`，可以使用以下端点查询状态：

**端点**: `GET /api/v1/workflows/{workflow_id}`

**响应**:

```json
{
  "id": 123,
  "name": "Cell Ranger Count",
  "workflow_type": "singlecell_cellranger",
  "status": "completed",
  "sample_id": 1,
  "input_files": {...},
  "output_files": {
    "filtered_h5": "/path/to/filtered_feature_bc_matrix.h5",
    "web_summary": "/path/to/web_summary.html"
  },
  "parameters": {...},
  "logs": "...",
  "error_message": null,
  "started_at": "2025-11-07T10:00:00Z",
  "completed_at": "2025-11-07T12:30:00Z"
}
```

**工作流状态**:

- `pending`: 等待执行
- `running`: 正在运行
- `completed`: 成功完成
- `failed`: 执行失败
- `cancelled`: 已取消

---

## 最佳实践

### 1. 基因组学分析流程

```bash
# 1. QC
curl -X POST "/api/v1/genomics/qc" -d '{...}'

# 2. 修剪
curl -X POST "/api/v1/genomics/trim" -d '{...}'

# 3. 比对
curl -X POST "/api/v1/genomics/align" -d '{...}'

# 4. 变异检测
curl -X POST "/api/v1/genomics/variant-calling" -d '{...}'

# 5. 注释
curl -X POST "/api/v1/genomics/annotate-variants" -d '{...}'
```

### 2. RNA-seq 分析流程

```bash
# 1. 定量
curl -X POST "/api/v1/transcriptomics/quantify" -d '{...}'

# 2. 差异表达
curl -X POST "/api/v1/transcriptomics/differential-expression" -d '{...}'

# 3. 可视化
curl -X POST "/api/v1/visualizations/volcano" -d '{...}'
```

### 3. 单细胞分析流程

```bash
# 1. Cell Ranger
curl -X POST "/api/v1/singlecell/cellranger" -d '{...}'

# 2. Scanpy 预处理
curl -X POST "/api/v1/singlecell/preprocess" -d '{...}'

# 3. 注释
curl -X POST "/api/v1/singlecell/annotate" -d '{...}'

# 4. UMAP 可视化
curl -X POST "/api/v1/visualizations/umap" -d '{...}'
```

### 4. ChIP-seq 分析流程

```bash
# 1. 比对
curl -X POST "/api/v1/epigenomics/align" \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{
    "sample_id": 1,
    "input_files": ["/data/chip_R1.fastq.gz", "/data/chip_R2.fastq.gz"],
    "reference_genome": "/ref/bowtie2/hg38",
    "output_bam": "/results/chip_aligned.bam",
    "aligner": "bowtie2",
    "threads": 8
  }'

# 2. Peak calling
curl -X POST "/api/v1/epigenomics/peak-calling" \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{
    "sample_id": 1,
    "treatment_bam": "/results/chip_aligned.bam",
    "control_bam": "/results/input_aligned.bam",
    "output_dir": "/results/peaks",
    "peak_caller": "macs2",
    "genome_size": "hs",
    "q_value": 0.05
  }'

# 3. 基序分析
curl -X POST "/api/v1/epigenomics/motif-analysis" \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{
    "sample_id": 1,
    "peak_file": "/results/peaks/peaks_peaks.narrowPeak",
    "genome_fasta": "/ref/hg38.fa",
    "output_dir": "/results/motifs",
    "tool": "homer"
  }'

# 4. 生成可视化文件
curl -X POST "/api/v1/epigenomics/bigwig" \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{
    "sample_id": 1,
    "input_bam": "/results/chip_aligned.bam",
    "output_bigwig": "/results/chip_coverage.bw",
    "genome_sizes": "/ref/hg38.chrom.sizes",
    "normalize": true
  }'
```

### 5. ATAC-seq 分析流程

```bash
# 1. 比对
curl -X POST "/api/v1/epigenomics/align" -d '{...}'

# 2. Peak calling (无 control, 使用 ATAC-seq 特定参数)
curl -X POST "/api/v1/epigenomics/peak-calling" \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{
    "sample_id": 1,
    "treatment_bam": "/results/atac_aligned.bam",
    "output_dir": "/results/atac_peaks",
    "peak_caller": "macs2",
    "genome_size": "hs",
    "nomodel": true,
    "shift": -100,
    "extsize": 200,
    "q_value": 0.01
  }'

# 3. 基序分析 (开放染色质区域的转录因子足迹)
curl -X POST "/api/v1/epigenomics/motif-analysis" -d '{...}'
```

---

## 错误处理

所有端点都返回标准的 HTTP 状态码：

- `200`: 成功
- `400`: 请求错误 (参数无效)
- `401`: 未认证
- `403`: 无权限
- `404`: 资源不存在
- `500`: 服务器错误

错误响应格式：

```json
{
  "detail": "Error message describing what went wrong"
}
```

---

## 认证

所有 API 端点都需要 JWT token 认证。在请求头中包含：

```
Authorization: Bearer <your_jwt_token>
```

获取 token:

```bash
curl -X POST "/api/v1/login/access-token" \
  -d "username=user@example.com&password=yourpassword"
```
