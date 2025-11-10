# Bioinformatics Format Conversion - Phase 1 Week 1

## 概述

Phase 1 Week 1 实现了生物信息学序列和区间格式的自动转换功能，支持以下格式：

### 序列格式（Sequence Formats）

- **FASTA** (`.fa`, `.fasta`, `.fna`)
- **FASTQ** (`.fq`, `.fastq`)
- **FASTQ.gz** (`.fastq.gz`, `.fq.gz`)

### 区间格式（Interval Formats）

- **BED** (`.bed`) - 支持 BED3, BED6, BED12 变体
- **bedGraph** (`.bg`, `.bedgraph`)
- **BigWig** (`.bw`, `.bigwig`)

## 转换路径

### 序列格式转换

```
FASTQ.gz ──decompress──> FASTQ ──convert──> FASTA ──export──> TSV
              │                    │
              └──────extract───────┴──────────────────────> TSV (with quality)
```

**支持的转换：**

- `FASTQ.gz → FASTQ`: 解压 gzip 压缩的 FASTQ 文件
- `FASTQ → FASTA`: 移除质量分数，仅保留序列
- `FASTQ → TSV`: 导出为表格（可选包含质量分数）
- `FASTA → TSV`: 导出为表格（ID, 描述, 序列, 长度）

### 区间格式转换

```
BigWig ──extract──> bedGraph ──convert──> BED ──export──> CSV
                         │                  │
                         └─────compress─────┘
```

**支持的转换：**

- `BigWig → bedGraph`: 提取染色体区间信号
- `bedGraph → BigWig`: 压缩为二进制格式（需要染色体大小文件）
- `bedGraph → BED`: 转换为 BED4 格式
- `BED → bedGraph`: 使用 score 列作为信号值
- `BED → CSV`: 导出为 CSV 表格
- `CSV → BED`: 从 CSV 导入 BED 格式

## 使用方法

### 1. 通过 FormatConverter 类

```python
from app.converters.format_converter import FormatConverter

converter = FormatConverter()

# 序列格式转换
result = converter.convert(
    source_path='input.fastq.gz',
    target_path='output.fasta',
    from_format='fastq.gz',
    to_format='fasta'
)

# 区间格式转换
result = converter.convert(
    source_path='input.bed',
    target_path='output.csv',
    from_format='bed',
    to_format='csv'
)

# BigWig to bedGraph (可选指定染色体)
result = converter.convert(
    source_path='input.bw',
    target_path='output.bedgraph',
    from_format='bigwig',
    to_format='bedgraph',
    chromosome='chr1'  # 可选：仅提取指定染色体
)
```

### 2. 直接使用专门的转换器

#### 序列转换器

```python
from app.converters.sequence_converter import get_sequence_converter

converter = get_sequence_converter()

# FASTQ to FASTA
converter.convert_fastq_to_fasta('input.fastq', 'output.fasta')

# FASTQ to TSV (包含质量分数)
converter.convert_fastq_to_tsv('input.fastq', 'output.tsv', include_quality=True)

# 解压 FASTQ.gz
converter.decompress_fastq_gz('input.fastq.gz', 'output.fastq')

# 获取序列统计信息
stats = converter.get_sequence_stats('input.fasta', 'fasta')
# 返回: {'count': 100, 'total_length': 50000, 'avg_length': 500, ...}
```

#### 区间转换器

```python
from app.converters.interval_converter import get_interval_converter

converter = get_interval_converter()

# BED to CSV
converter.convert_bed_to_csv('input.bed', 'output.csv', bed_type='6')

# bedGraph to BigWig (需要染色体大小)
chrom_sizes = {'chr1': 248956422, 'chr2': 242193529}
converter.convert_bedgraph_to_bigwig('input.bedgraph', 'output.bw', chrom_sizes)

# 获取 BED 统计信息
stats = converter.get_bed_stats('input.bed')
# 返回: {'count': 1000, 'total_length': 500000, 'chromosomes': ['chr1', 'chr2'], ...}
```

## API 端点

### 格式检测

```bash
GET /api/formats/detect?file_path=/path/to/file.fastq

# 响应
{
  "format": "fastq",
  "extensions": [".fq", ".fastq"]
}
```

### 转换路径查询

```bash
GET /api/formats/conversion-path?from=fastq.gz&to=fasta

# 响应
{
  "path": ["fastq.gz", "fastq", "fasta"],
  "steps": 2
}
```

### 转换时间估算

```bash
POST /api/formats/estimate-time
{
  "file_path": "/path/to/file.fastq.gz",
  "from_format": "fastq.gz",
  "to_format": "fasta"
}

# 响应
{
  "estimated_seconds": 8.5,
  "file_size_gb": 1.7
}
```

### 执行转换

```bash
POST /api/formats/convert
{
  "source_path": "/path/to/input.fastq",
  "target_path": "/path/to/output.fasta",
  "from_format": "fastq",
  "to_format": "fasta"
}

# 响应
{
  "status": "completed",
  "source": "/path/to/input.fastq",
  "target": "/path/to/output.fasta",
  "steps": [...],
  "total_time": 3.2
}
```

## 转换时间估算

基于每 GB 文件的处理时间：

| 转换类型            | 时间/GB | 说明           |
| ------------------- | ------- | -------------- |
| `fastq.gz → fastq`  | 5 秒    | 解压缩         |
| `fastq → fasta`     | 3 秒    | 移除质量分数   |
| `fasta → tsv`       | 2 秒    | 导出表格       |
| `fastq → tsv`       | 3 秒    | 导出表格       |
| `bed → csv`         | 2 秒    | 导出表格       |
| `bedgraph → bed`    | 2 秒    | 格式转换       |
| `bed → bedgraph`    | 2 秒    | 格式转换       |
| `bigwig → bedgraph` | 15 秒   | 提取二进制数据 |
| `bedgraph → bigwig` | 20 秒   | 压缩为二进制   |

## 文件格式说明

### FASTA 格式

```
>sequence_id description
ATCGATCGATCGATCG
>another_sequence
GCTAGCTAGCTAGCTA
```

### FASTQ 格式

```
@read_id description
ATCGATCGATCGATCG
+
IIIIIIIIIIIIIIII
```

### BED3 格式

```
chr1    100    200
chr1    300    400
```

### BED6 格式

```
chr1    100    200    feature1    100    +
chr1    300    400    feature2    200    -
```

### BED12 格式

完整的 12 列 BED 格式，包含 blockCount, blockSizes, blockStarts 等信息。

### bedGraph 格式

```
chr1    100    200    1.5
chr1    200    300    2.3
```

### BigWig 格式

二进制压缩格式，用于高效存储大规模基因组信号数据（如 ChIP-seq, RNA-seq coverage）。

## 依赖项

新增的 Python 包（已添加到 `pyproject.toml`）：

```toml
"biopython>=1.81",      # 序列格式处理
"pybedtools>=0.9.0",    # BED 格式处理
"pyBigWig>=0.3.18"      # BigWig 格式处理
```

## 测试

运行测试套件：

```bash
# 运行所有生信格式转换测试
pytest tests/test_bioinformatics_converters.py -v

# 仅运行序列转换器测试
pytest tests/test_bioinformatics_converters.py::TestSequenceConverter -v

# 仅运行区间转换器测试
pytest tests/test_bioinformatics_converters.py::TestIntervalConverter -v

# 运行集成测试
pytest tests/test_bioinformatics_converters.py::TestFormatConverterIntegration -v
```

## 后续开发计划

### Phase 1 Week 2 (下一步)

- **比对格式**: SAM, BAM
- **变异格式**: VCF, BCF
- **注释格式**: GTF, GFF3

### Phase 1 Week 3

- **表达矩阵**: MTX, 10X HDF5
- **Python 数据**: h5, npy, npz
- **R 数据**: RData
- **压缩格式**: gzip, bgzip, zip, tar

参见 [`docs/FORMAT_CONVERSION_ROADMAP.md`](./FORMAT_CONVERSION_ROADMAP.md) 获取完整计划。

## 常见问题

### Q: BigWig 转换为什么需要染色体大小文件？

A: BigWig 是索引格式，需要预先知道每个染色体的长度。可以：

1. 提供 `chrom_sizes` 字典
2. 从参考基因组 `.fai` 文件读取
3. 使用 UCSC 标准染色体大小

### Q: 如何处理大文件？

A: 所有转换器使用流式处理或分块读取：

- FASTQ 解压：1MB 块
- BED 转换：Pandas 批处理
- BigWig 提取：按染色体处理

### Q: 支持哪些 BED 变体？

A: 支持 BED3 (3 列), BED6 (6 列), BED12 (12 列)，会自动检测列数。

### Q: 质量分数如何编码？

A: 使用标准 Phred+33 编码（Sanger/Illumina 1.8+）。平均质量分数通过 ASCII 值减 33 计算。

## 相关文档

- [生信格式分析](./BIOINFORMATICS_FORMATS_ANALYSIS.md) - 52 种格式详细分析
- [格式转换路线图](./FORMAT_CONVERSION_ROADMAP.md) - 3 阶段实现计划
- [格式转换系统](./FORMAT_CONVERSION_SYSTEM.md) - 原有系统文档
- [测试文档](./TESTING.md) - 测试指南
