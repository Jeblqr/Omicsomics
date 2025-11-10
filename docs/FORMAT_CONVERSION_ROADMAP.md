# 格式转换系统路线图 (Format Conversion System Roadmap)

## 📊 当前状态 (Current Status)

### ✅ 已完成 - 基础版 (TODO #17)

**支持格式**: 7 种通用格式

- CSV, TSV, Excel, JSON (表格)
- RDS (R), h5ad (Python), pickle (Python)

**功能**:

- ✅ 自动转换模式 (Pipeline Editor)
- ✅ 手动转换模式 (Data Browser)
- ✅ 8 个 API 端点
- ✅ 完整文档

**文件**:

- `backend/app/converters/format_converter.py`
- `backend/app/converters/pipeline_auto_converter.py`
- `backend/app/models/format_conversion.py`
- `backend/app/api/routers/formats.py`
- `frontend/src/components/FormatConverterModal.tsx`
- `frontend/src/components/AutoConversionIndicator.tsx`
- `docs/FORMAT_CONVERSION_SYSTEM.md`
- `docs/FORMAT_CONVERSION_IMPLEMENTATION.md`

---

## 🎯 Phase 1: 生信格式自动转换 (TODO #18)

### 目标

添加 20 种生物信息学常用格式的自动转换支持

### 优先级

🔴 **高优先级** - 立即开始实现

### 预期工作量

2-3 周

### 格式清单 (20 种)

#### 1. 序列格式 (3 种)

- **FASTA** (.fa, .fasta, .fna) → TSV (序列列表)
- **FASTQ** (.fq, .fastq) → FASTA, TSV
- **FASTQ.gz** (.fastq.gz, .fq.gz) → FASTQ (解压)

**转换路径**:

```
FASTQ.gz → FASTQ → FASTA → TSV
```

**工具**: BioPython

```python
from Bio import SeqIO
# FASTQ to FASTA
SeqIO.convert("input.fastq", "fastq", "output.fasta", "fasta")
```

#### 2. 基因组区间格式 (5 种)

- **BED3** (.bed, 3 列) → TSV, CSV
- **BED6** (.bed, 6 列) → BED3, BED12, TSV
- **BED12** (.bed, 12 列) → BED6, GTF
- **bedGraph** (.bg, .bedgraph) → BED, WIG
- **BigWig** (.bw, .bigwig) → bedGraph, WIG

**转换路径**:

```
BigWig → bedGraph → BED → TSV/CSV
             ↓        ↓
           WIG     BED6 → BED12 → GTF
```

**工具**: pybedtools, pyBigWig

```python
import pyBigWig
bw = pyBigWig.open("input.bw")
# Convert to bedGraph
```

#### 3. 基因注释格式 (2 种)

- **GTF** (.gtf) → GFF3, BED12, TSV
- **GFF3** (.gff3) → GTF, BED, TSV

**转换路径**:

```
GTF ↔ GFF3 → BED → TSV
```

**工具**: pybedtools, BCBio.GFF

```python
from BCBio import GFF
# Parse GFF/GTF files
```

#### 4. 比对格式 (2 种)

- **SAM** (.sam) → BAM, TSV (统计)
- **BAM** (.bam) → SAM, BED (覆盖), TSV

**转换路径**:

```
BAM → SAM → TSV (统计表)
  ↓
BED (覆盖度)
```

**工具**: pysam, samtools

```python
import pysam
# BAM to SAM
pysam.view("-h", "-o", "output.sam", "input.bam")
```

#### 5. 变异格式 (2 种)

- **VCF** (.vcf) → BCF, BED, TSV, CSV
- **BCF** (.bcf) → VCF

**转换路径**:

```
BCF → VCF → TSV/CSV
        ↓
       BED (变异位置)
```

**工具**: pysam, bcftools

```python
import pysam
vcf = pysam.VariantFile("input.vcf")
# Parse variants
```

#### 6. 表达矩阵格式 (2 种)

- **MTX** (.mtx, Market Matrix) → CSV (稠密化), h5ad
- **10X HDF5** (.h5) → h5ad, MTX, CSV

**转换路径**:

```
10X HDF5 → h5ad (AnnData) → CSV
    ↓           ↓              ↓
  MTX        RDS (Seurat)    TSV
```

**工具**: scipy, h5py, anndata

```python
from scipy.io import mmread, mmwrite
import anndata as ad
# Read 10X data
adata = ad.read_10x_h5("input.h5")
```

#### 7. Python 数据格式 (3 种)

- **HDF5** (.h5, .hdf5) → CSV (数据集提取)
- **NPY** (.npy, NumPy 数组) → CSV, TSV
- **NPZ** (.npz, 压缩 NumPy) → 多个 NPY

**工具**: h5py, numpy

```python
import h5py
import numpy as np
# Read HDF5
f = h5py.File("input.h5", "r")
```

#### 8. R 数据格式 (1 种)

- **RData** (.RData, .rda) → RDS (单对象提取), CSV

**工具**: rpy2

```python
import rpy2.robjects as ro
# Load RData
ro.r['load']("input.RData")
```

#### 9. 压缩格式 (5 种)

- **gzip** (.gz) → 原格式 (解压)
- **bgzip** (.gz, bgzip 压缩) → 原格式, gzip
- **zip** (.zip) → 原文件 (解压)
- **tar** (.tar) → 原文件 (解包)
- **tar.gz** (.tar.gz, .tgz) → tar, 原文件

**工具**: Python built-in (gzip, zipfile, tarfile)

```python
import gzip
import zipfile
import tarfile
```

### 实施步骤

#### Week 1: 序列和区间格式

- [ ] 安装依赖: BioPython, pybedtools, pyBigWig
- [ ] 实现 FASTA/FASTQ 转换器
- [ ] 实现 BED/bedGraph/BigWig 转换器
- [ ] 单元测试
- [ ] 集成到 FormatConverter 类

#### Week 2: 比对和变异格式

- [ ] 安装依赖: pysam, samtools, bcftools
- [ ] 实现 SAM/BAM 转换器
- [ ] 实现 VCF/BCF 转换器
- [ ] 实现 GTF/GFF3 转换器
- [ ] 单元测试

#### Week 3: 表达矩阵和压缩格式

- [ ] 安装依赖: scipy, h5py, numpy, rpy2
- [ ] 实现 MTX/10X HDF5 转换器
- [ ] 实现 HDF5/NPY/NPZ 转换器
- [ ] 实现 RData 转换器
- [ ] 实现压缩格式处理
- [ ] 完整测试和文档

### 新增文件

#### Backend (3 个扩展文件)

1. `backend/app/converters/bioinformatics_converter.py` (新建)

   - BioinformaticsConverter 类
   - 序列、区间、比对格式转换方法

2. `backend/app/converters/genomics_converter.py` (新建)

   - GenomicsConverter 类
   - 变异、注释格式转换方法

3. `backend/app/converters/compression_handler.py` (新建)
   - CompressionHandler 类
   - 压缩文件自动检测和解压

#### Backend (1 个更新文件)

1. `backend/app/converters/format_converter.py` (更新)
   - 集成新转换器
   - 更新 FORMATS 字典
   - 更新转换路径逻辑

#### 依赖配置

`backend/pyproject.toml` (更新):

```toml
[tool.poetry.dependencies]
biopython = "^1.81"
pybedtools = "^0.9.0"
pyBigWig = "^0.3.18"
pysam = "^0.21.0"
rpy2 = "^3.5.0"
h5py = "^3.9.0"
```

### API 扩展

#### 新增端点

```http
GET /api/formats/bioinformatics
# 获取生信格式列表

POST /api/formats/decompress
# 解压缩文件

GET /api/formats/inspect/{file_path}
# 检查文件内容类型 (FASTA/FASTQ/BED/等)
```

### 文档更新

- 更新 `docs/FORMAT_CONVERSION_SYSTEM.md`
  - 添加 20 种新格式文档
  - 转换示例
  - 工具依赖说明
- 创建 `docs/BIOINFORMATICS_FORMATS_GUIDE.md`
  - 每种格式详细说明
  - 转换最佳实践

---

## 🎨 Phase 2: 交互式转换框架 (TODO #19)

### 目标

构建交互式转换向导，处理需要用户参数的复杂场景

### 优先级

🔴 **高优先级** - 紧随 Phase 1

### 预期工作量

3-4 周

### 核心组件 (6 个)

#### 1. ColumnMappingWizard

**用途**: 列名映射
**场景**: GWAS 汇总统计、基因表达矩阵

```tsx
<ColumnMappingWizard
  sourceColumns={["chr", "pos", "beta", "p"]}
  targetSchema={{
    chromosome: "required",
    position: "required",
    effect_size: "required",
    p_value: "required",
  }}
  onMapping={(mapping) => console.log(mapping)}
/>
```

#### 2. GenomeVersionSelector

**用途**: 基因组版本选择
**场景**: GWAS、BED、VCF 坐标转换

```tsx
<GenomeVersionSelector
  sourceVersion="hg19"
  targetOptions={["hg19", "hg38", "GRCh37", "GRCh38"]}
  onSelect={(version) => console.log(version)}
/>
```

#### 3. QCThresholdPanel

**用途**: 质量控制阈值
**场景**: 单细胞数据导入、蛋白质组学

```tsx
<QCThresholdPanel
  metrics={["min_genes", "min_cells", "max_percent_mito"]}
  defaults={{ min_genes: 200, min_cells: 3, max_percent_mito: 5 }}
  onUpdate={(thresholds) => console.log(thresholds)}
/>
```

#### 4. DataPreviewTable

**用途**: 数据预览
**场景**: 所有交互式转换

```tsx
<DataPreviewTable
  data={firstNRows}
  totalRows={1000000}
  highlighting={{
    columns: ["problematic_column"],
    message: "Missing values detected",
  }}
/>
```

#### 5. ConversionParameterForm

**用途**: 转换参数表单
**场景**: 归一化方法、单位转换

```tsx
<ConversionParameterForm
  parameters={[
    { name: "normalization", type: "select", options: ["TPM", "CPM", "RPKM"] },
    { name: "log_transform", type: "boolean", default: true },
  ]}
  onSubmit={(params) => console.log(params)}
/>
```

#### 6. ProgressTracker

**用途**: 转换进度追踪
**场景**: 大文件、多步转换

```tsx
<ProgressTracker
  steps={[
    { name: "Reading file", status: "completed" },
    { name: "Validating data", status: "running", progress: 45 },
    { name: "Converting format", status: "pending" },
  ]}
/>
```

### 应用场景 (10 个)

#### 1. GWAS 汇总统计标准化

**输入**: 各种 GWAS 结果文件
**交互**:

1. 上传文件 → 自动检测列
2. 列映射向导 (SNP → variant_id, Chr → chromosome, etc.)
3. 基因组版本选择 (hg19/hg38)
4. Effect 类型选择 (beta/OR/log(OR))
5. 单位转换 (自动或手动)
6. 预览标准化结果
7. 执行转换

**输出**: 标准化 GWAS TSV

#### 2. 基因表达矩阵标准化

**输入**: 各种工具输出的表达矩阵
**交互**:

1. 上传文件 → 检测维度
2. 方向确认 (gene×sample or sample×gene)
3. 基因 ID 类型 (Ensembl/Symbol/Entrez)
4. 样本元数据上传 (可选)
5. 归一化方法选择
6. 预览标准化矩阵
7. 执行转换

**输出**: Seurat RDS / h5ad / 标准 CSV

#### 3. 单细胞数据导入向导

**输入**: 10X, Drop-seq, Smart-seq 等
**交互**:

1. 文件类型识别 (matrix.mtx + barcodes + features)
2. 物种选择 (human/mouse/other)
3. 基因版本 (GRCh38/GRCm38)
4. QC 阈值设定 (min_genes, min_cells, max_mito%)
5. 预览细胞/基因过滤统计
6. 选择输出格式 (Seurat/h5ad/Loom)
7. 执行导入

**输出**: Seurat RDS / h5ad

#### 4. 蛋白质组学结果处理

**输入**: MaxQuant, Proteome Discoverer 等
**交互**:

1. 工具类型识别
2. 蛋白分组策略 (razor/shared)
3. FDR 阈值 (1%/5%)
4. 定量方法 (LFQ/iBAQ/intensity)
5. 缺失值处理 (filter/impute)
6. 预览蛋白矩阵
7. 执行转换

**输出**: 标准化蛋白矩阵 CSV

#### 5. 代谢组学数据标准化

**输入**: 各平台质谱数据
**交互**:

1. 代谢物 ID 映射 (HMDB/KEGG/PubChem)
2. 峰强度类型 (area/height)
3. 归一化方法 (median/quantile/internal standard)
4. 样本分组信息上传
5. 预览归一化结果
6. 执行转换

**输出**: 标准化代谢矩阵 CSV

#### 6. 表观遗传学数据转换

**输入**: Peak calling, methylation array
**交互**:

1. 数据类型 (ChIP-seq/ATAC-seq/Methylation)
2. 峰定义参数 (width, gap)
3. 基因组版本
4. 背景模型 (genomic/local)
5. 阈值设定 (p-value/FDR)
6. 预览峰分布
7. 执行转换

**输出**: BED peaks / BigWig coverage / 甲基化矩阵

#### 7. 多组学数据整合

**输入**: 多层组学数据
**交互**:

1. 数据层上传 (RNA/Protein/Metabolite)
2. 样本 ID 映射表
3. 特征对齐策略 (gene symbol/ID)
4. 归一化方法 (per-layer)
5. 缺失数据处理
6. 预览整合矩阵
7. 执行整合

**输出**: MultiAssayExperiment (R) / MuData (Python)

#### 8. 网络数据格式化

**输入**: 各种网络表示
**交互**:

1. 网络类型 (PPI/GRN/metabolic)
2. 节点 ID 类型
3. 边权重解释 (correlation/probability)
4. 方向性 (directed/undirected)
5. 预览网络统计
6. 执行转换

**输出**: Cytoscape format / GraphML / 邻接矩阵

#### 9. 临床数据标准化

**输入**: 电子病历、临床试验数据
**交互**:

1. PHI 脱敏确认
2. 日期格式统一
3. 分类变量编码
4. 连续变量单位
5. 缺失值编码
6. 预览标准化数据
7. 执行转换

**输出**: CDISC 标准 / CSV / REDCap

#### 10. 测序质控报告解析

**输入**: FastQC, MultiQC HTML/JSON
**交互**:

1. 报告类型识别
2. 指标选择
3. 阈值设定 (pass/warn/fail)
4. 样本分组
5. 预览汇总表
6. 执行解析

**输出**: CSV 指标表 / JSON 结构化

### API 设计

```http
# 启动交互式转换
POST /api/formats/interactive/start
{
  "file_path": "/data/gwas_summary.txt",
  "scenario": "gwas_standardization"
}
Response: {
  "session_id": "uuid",
  "detected_columns": [...],
  "suggested_mapping": {...},
  "next_step": "column_mapping"
}

# 预览转换结果
POST /api/formats/interactive/preview
{
  "session_id": "uuid",
  "parameters": {...}
}
Response: {
  "preview_data": [...],
  "warnings": [...],
  "estimated_time": 15.5
}

# 执行转换
POST /api/formats/interactive/convert
{
  "session_id": "uuid",
  "parameters": {...},
  "target_path": "/data/gwas_standardized.tsv"
}
Response: {
  "conversion_id": 123,
  "status": "running"
}
```

### 实施步骤

#### Week 1-2: 核心组件开发

- [ ] ColumnMappingWizard 组件
- [ ] GenomeVersionSelector 组件
- [ ] QCThresholdPanel 组件
- [ ] DataPreviewTable 组件
- [ ] ConversionParameterForm 组件
- [ ] ProgressTracker 组件

#### Week 3: 场景 1-5 实现

- [ ] GWAS 标准化向导
- [ ] 表达矩阵标准化
- [ ] 单细胞导入向导
- [ ] 蛋白质组学处理
- [ ] 代谢组学标准化

#### Week 4: 场景 6-10 实现

- [ ] 表观遗传学转换
- [ ] 多组学整合
- [ ] 网络数据格式化
- [ ] 临床数据标准化
- [ ] QC 报告解析

### 新增文件

#### Frontend (7 个组件)

1. `frontend/src/components/conversion/ColumnMappingWizard.tsx`
2. `frontend/src/components/conversion/GenomeVersionSelector.tsx`
3. `frontend/src/components/conversion/QCThresholdPanel.tsx`
4. `frontend/src/components/conversion/DataPreviewTable.tsx`
5. `frontend/src/components/conversion/ConversionParameterForm.tsx`
6. `frontend/src/components/conversion/ProgressTracker.tsx`
7. `frontend/src/components/InteractiveConversionModal.tsx`

#### Backend (2 个新文件)

1. `backend/app/converters/interactive_converter.py`
2. `backend/app/api/routers/interactive_formats.py`

---

## 🔵 Phase 3: 补充格式 (TODO #20)

### 目标

添加 15 种中优先级补充格式

### 优先级

🟡 **中优先级** - Phase 1 和 2 完成后

### 预期工作量

2-3 周

### 格式清单

- Parquet, Feather (列式存储)
- Loom (单细胞层次化)
- CRAM (压缩比对)
- PAF (轻量比对)
- MAF (突变注释)
- WIG (覆盖度)
- BigBed (二进制区间)
- GenBank (注释)
- mzML, mzXML, MGF (质谱)
- 等

---

## 📈 成功指标

### 覆盖率

- ✅ 基础版: 7 种格式
- 🎯 Phase 1: 27 种格式 (7+20)
- 🎯 Phase 2: 37 种场景 (27+10)
- 🎯 Phase 3: 52 种格式 (37+15)

### 用户体验

- 支持 90%+常见生信格式
- 自动识别成功率 >95%
- 交互式转换<5 步完成
- 转换时间<30s (1GB 数据)

### 工作流影响

- 减少手动转换时间 80%
- 避免格式错误 95%
- 跨工具数据流无缝衔接

---

## 📅 时间表

| Phase      | 内容          | 工期   | 开始日期   | 结束日期   |
| ---------- | ------------- | ------ | ---------- | ---------- |
| ✅ 基础版  | 7 种通用格式  | -      | 2025-01-10 | 2025-01-10 |
| 🔴 Phase 1 | 20 种生信格式 | 2-3 周 | 2025-01-11 | 2025-02-01 |
| 🔴 Phase 2 | 10 种交互场景 | 3-4 周 | 2025-02-02 | 2025-03-02 |
| 🟡 Phase 3 | 15 种补充格式 | 2-3 周 | 2025-03-03 | 2025-03-24 |

**总计**: 约 2.5-3 个月完成全部格式转换系统

---

## 📚 相关文档

- [格式分析清单](./BIOINFORMATICS_FORMATS_ANALYSIS.md)
- [格式转换系统使用指南](./FORMAT_CONVERSION_SYSTEM.md)
- [格式转换实现总结](./FORMAT_CONVERSION_IMPLEMENTATION.md)
- [工具变体系统](./TOOL_VARIANT_SYSTEM_SUMMARY.md)
- [跨运行时集成](./CROSS_RUNTIME_INTEGRATION.md)

---

**版本**: 2.0.0  
**创建日期**: 2025-01-10  
**下一步**: 开始 Phase 1 实现
