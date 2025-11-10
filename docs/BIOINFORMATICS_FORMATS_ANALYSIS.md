# 生物信息学数据格式完整清单 (Bioinformatics Data Formats Comprehensive List)

## 分类体系

### 类别 A: 自动转换格式 (Automatic Conversion)

**特点**: 格式明确、转换规则固定、无需用户交互

### 类别 B: 交互式转换格式 (Interactive Conversion)

**特点**: 需要用户指定参数、列映射、或数据解释

---

## 📊 类别 A: 自动转换格式

### 1. 通用表格格式 (General Tabular Formats)

| 格式    | 扩展名      | 描述       | 转换目标                    | 优先级      |
| ------- | ----------- | ---------- | --------------------------- | ----------- |
| CSV     | .csv        | 逗号分隔   | TSV, Excel, JSON, RDS, h5ad | ✅ 已实现   |
| TSV     | .tsv, .txt  | Tab 分隔   | CSV, Excel, JSON            | ✅ 已实现   |
| Excel   | .xlsx, .xls | Excel 表格 | CSV, TSV, JSON              | ✅ 已实现   |
| JSON    | .json       | JSON 格式  | CSV, TSV, Excel             | ✅ 已实现   |
| Parquet | .parquet    | 列式存储   | CSV, TSV                    | 🔴 高优先级 |
| Feather | .feather    | 轻量二进制 | CSV, TSV                    | 🟡 中优先级 |

### 2. R 数据格式 (R Data Formats)

| 格式  | 扩展名       | 描述     | 转换目标         | 优先级      |
| ----- | ------------ | -------- | ---------------- | ----------- |
| RDS   | .rds, .RDS   | R 单对象 | CSV, h5ad        | ✅ 已实现   |
| RData | .RData, .rda | R 多对象 | RDS (单对象提取) | 🔴 高优先级 |
| Robj  | .Robj        | R 对象   | RDS              | 🟡 中优先级 |

### 3. Python 数据格式 (Python Data Formats)

| 格式   | 扩展名        | 描述           | 转换目标         | 优先级      |
| ------ | ------------- | -------------- | ---------------- | ----------- |
| pickle | .pkl, .pickle | Python 序列化  | CSV, JSON        | ✅ 已实现   |
| h5ad   | .h5ad         | AnnData 单细胞 | CSV, RDS         | ✅ 已实现   |
| h5     | .h5, .hdf5    | HDF5 通用      | CSV (数据集提取) | 🔴 高优先级 |
| npy    | .npy          | NumPy 数组     | CSV, TSV         | 🔴 高优先级 |
| npz    | .npz          | NumPy 压缩     | 多个 npy         | 🔴 高优先级 |

### 4. 基因组序列格式 (Genomic Sequence Formats)

| 格式        | 扩展名            | 描述         | 转换目标            | 优先级      |
| ----------- | ----------------- | ------------ | ------------------- | ----------- |
| FASTA       | .fa, .fasta, .fna | 序列格式     | TSV (序列列表)      | 🔴 高优先级 |
| FASTQ       | .fq, .fastq       | 带质量序列   | TSV, FASTA (去质量) | 🔴 高优先级 |
| Multi-FASTA | .mfa              | 多序列 FASTA | 单 FASTA (拆分)     | 🟡 中优先级 |
| FASTQ.gz    | .fastq.gz, .fq.gz | 压缩 FASTQ   | FASTQ (解压)        | 🔴 高优先级 |

### 5. 基因组区间格式 (Genomic Interval Formats)

| 格式     | 扩展名         | 描述          | 转换目标         | 优先级      |
| -------- | -------------- | ------------- | ---------------- | ----------- |
| BED      | .bed           | 基因组区间    | TSV, CSV, GFF    | 🔴 高优先级 |
| BED6     | .bed           | 6 列标准 BED  | BED3, BED12      | 🔴 高优先级 |
| BED12    | .bed           | 12 列扩展 BED | BED6, GTF        | 🔴 高优先级 |
| bedGraph | .bg, .bedgraph | 覆盖度数据    | WIG, BED         | 🔴 高优先级 |
| WIG      | .wig           | 可变步长覆盖  | bedGraph, BigWig | 🟡 中优先级 |
| BigWig   | .bw, .bigwig   | 二进制覆盖    | WIG, bedGraph    | 🔴 高优先级 |
| BigBed   | .bb, .bigbed   | 二进制 BED    | BED              | 🟡 中优先级 |

### 6. 基因注释格式 (Gene Annotation Formats)

| 格式    | 扩展名    | 描述         | 转换目标    | 优先级      |
| ------- | --------- | ------------ | ----------- | ----------- |
| GTF     | .gtf      | 基因转录本   | GFF3, BED12 | 🔴 高优先级 |
| GFF     | .gff      | 基因特征     | GTF, BED    | 🔴 高优先级 |
| GFF3    | .gff3     | GFF 版本 3   | GTF, BED    | 🔴 高优先级 |
| GenBank | .gb, .gbk | GenBank 格式 | FASTA, GFF  | 🟡 中优先级 |

### 7. 比对格式 (Alignment Formats)

| 格式 | 扩展名 | 描述       | 转换目标        | 优先级      |
| ---- | ------ | ---------- | --------------- | ----------- |
| SAM  | .sam   | 文本比对   | BAM, TSV (统计) | 🔴 高优先级 |
| BAM  | .bam   | 二进制比对 | SAM, CRAM       | 🔴 高优先级 |
| CRAM | .cram  | 压缩比对   | BAM, SAM        | 🟡 中优先级 |
| PAF  | .paf   | 轻量比对   | TSV, SAM        | 🟡 中优先级 |

### 8. 变异格式 (Variant Formats)

| 格式 | 扩展名 | 描述       | 转换目标      | 优先级      |
| ---- | ------ | ---------- | ------------- | ----------- |
| VCF  | .vcf   | 变异调用   | TSV, CSV, BED | 🔴 高优先级 |
| BCF  | .bcf   | 二进制 VCF | VCF           | 🔴 高优先级 |
| MAF  | .maf   | 突变注释   | VCF, TSV      | 🟡 中优先级 |

### 9. 表达矩阵格式 (Expression Matrix Formats)

| 格式     | 扩展名 | 描述         | 转换目标           | 优先级      |
| -------- | ------ | ------------ | ------------------ | ----------- |
| MTX      | .mtx   | 稀疏矩阵     | CSV (稠密化), h5ad | 🔴 高优先级 |
| 10X HDF5 | .h5    | 10X Genomics | h5ad, MTX          | 🔴 高优先级 |
| Loom     | .loom  | 层次化存储   | h5ad, CSV          | 🟡 中优先级 |

### 10. 质谱数据格式 (Mass Spectrometry Formats)

| 格式  | 扩展名 | 描述         | 转换目标     | 优先级      |
| ----- | ------ | ------------ | ------------ | ----------- |
| mzML  | .mzml  | 质谱 XML     | CSV (峰列表) | 🟡 中优先级 |
| mzXML | .mzxml | 旧版质谱 XML | mzML         | 🟡 中优先级 |
| MGF   | .mgf   | Mascot 格式  | mzML, CSV    | 🟡 中优先级 |

### 11. 压缩格式 (Compressed Formats)

| 格式   | 扩展名        | 描述       | 转换目标      | 优先级      |
| ------ | ------------- | ---------- | ------------- | ----------- |
| gzip   | .gz           | GZIP 压缩  | 原格式 (解压) | 🔴 高优先级 |
| bgzip  | .gz (bgzip)   | 可索引压缩 | 原格式, gzip  | 🔴 高优先级 |
| zip    | .zip          | ZIP 归档   | 原文件 (解压) | 🔴 高优先级 |
| tar    | .tar          | TAR 归档   | 原文件 (解包) | 🔴 高优先级 |
| tar.gz | .tar.gz, .tgz | 压缩归档   | tar, 原文件   | 🔴 高优先级 |

---

## 📝 类别 B: 交互式转换格式

### 1. GWAS 汇总统计 (GWAS Summary Statistics)

**格式**: 多种变体（无统一标准）

**常见列**:

- SNP/variant ID
- Chromosome, Position
- Effect/Other allele
- Effect size (beta/OR)
- P-value
- Sample size
- MAF

**需要交互**:

- ✅ 列名映射 (不同研究使用不同列名)
- ✅ 坐标系统选择 (hg19/hg38/hg37)
- ✅ Effect allele 识别
- ✅ 单位转换 (beta/OR/log(OR))

**转换目标**:

- 标准化 GWAS 格式
- VCF (用于下游分析)
- BED (位置信息)

### 2. 基因表达矩阵 (Gene Expression Matrix)

**格式**: 各种工具输出的矩阵

**需要交互**:

- ✅ 行列方向确定 (gene×sample or sample×gene)
- ✅ 基因 ID 类型识别 (Ensembl/Symbol/Entrez)
- ✅ 样本元数据关联
- ✅ 批次效应处理选择

**转换目标**:

- Seurat 对象 (RDS)
- AnnData (h5ad)
- CSV 标准化格式

### 3. 代谢组学数据 (Metabolomics Data)

**格式**: 多种质谱平台输出

**需要交互**:

- ✅ 代谢物 ID 映射 (HMDB/KEGG/PubChem)
- ✅ 峰强度类型选择
- ✅ 归一化方法选择
- ✅ 样本分组信息

**转换目标**:

- 标准化代谢矩阵
- mzTab 格式
- CSV

### 4. 蛋白质组学搜库结果 (Proteomics Search Results)

**格式**: MaxQuant, Proteome Discoverer, Mascot 等

**需要交互**:

- ✅ 蛋白质分组策略
- ✅ FDR 阈值设置
- ✅ 定量方法选择 (LFQ/iBAQ/intensity)
- ✅ 缺失值处理

**转换目标**:

- 标准化蛋白矩阵
- mzTab
- CSV

### 5. 单细胞计数矩阵 (Single-cell Count Matrix)

**格式**: 10X, Drop-seq, Smart-seq 等

**需要交互**:

- ✅ 文件类型识别 (matrix.mtx + barcodes + features)
- ✅ 物种选择
- ✅ 基因版本选择
- ✅ QC 阈值设定

**转换目标**:

- Seurat 对象
- h5ad (scanpy)
- Loom

### 6. 表观遗传学数据 (Epigenomics Data)

**格式**: Peak calling, methylation array 等

**需要交互**:

- ✅ 峰定义参数
- ✅ 基因组版本
- ✅ 背景模型选择
- ✅ 差异甲基化阈值

**转换目标**:

- BED (peaks)
- BigWig (coverage)
- 甲基化矩阵 (CSV)

### 7. 多组学整合数据 (Multi-omics Data)

**格式**: 多层次数据组合

**需要交互**:

- ✅ 样本 ID 映射
- ✅ 特征对齐策略
- ✅ 归一化方法
- ✅ 缺失数据处理

**转换目标**:

- MultiAssayExperiment (R)
- MuData (Python)
- 标准化多表

### 8. 网络/通路数据 (Network/Pathway Data)

**格式**: 各种网络表示

**需要交互**:

- ✅ 节点 ID 类型
- ✅ 边权重解释
- ✅ 方向性确定
- ✅ 网络类型 (PPI/GRN/metabolic)

**转换目标**:

- Cytoscape 格式
- GraphML
- 邻接矩阵 (CSV)

### 9. 临床数据 (Clinical Data)

**格式**: 电子病历、临床试验数据

**需要交互**:

- ✅ PHI 脱敏确认
- ✅ 日期格式统一
- ✅ 分类变量编码
- ✅ 连续变量单位

**转换目标**:

- CDISC 标准
- CSV 标准化
- REDCap 格式

### 10. 测序质控报告 (Sequencing QC Reports)

**格式**: FastQC, MultiQC 等 HTML/JSON

**需要交互**:

- ✅ 指标选择
- ✅ 阈值设定
- ✅ 样本分组

**转换目标**:

- CSV (指标表)
- JSON (结构化)
- HTML (可视化)

---

## 🎯 实施优先级

### Phase 1: 高优先级自动转换 (立即实现)

**包含格式**: 20 种

1. **序列格式**: FASTA, FASTQ, FASTQ.gz
2. **区间格式**: BED (3/6/12), bedGraph, BigWig
3. **注释格式**: GTF, GFF3
4. **比对格式**: SAM, BAM
5. **变异格式**: VCF, BCF
6. **表达格式**: MTX, 10X HDF5
7. **Python 格式**: h5, npy, npz
8. **R 格式**: RData
9. **压缩格式**: gzip, bgzip, zip, tar, tar.gz

**预计工作量**: 2-3 周

### Phase 2: 交互式转换框架 (随后实现)

**包含场景**: 10 种

1. GWAS 汇总统计标准化
2. 基因表达矩阵标准化
3. 单细胞数据导入向导
4. 蛋白质组学结果处理
5. 代谢组学数据标准化
6. 表观遗传学数据转换
7. 多组学数据整合
8. 网络数据格式化
9. 临床数据标准化
10. QC 报告解析

**预计工作量**: 3-4 周

### Phase 3: 中优先级自动转换 (后续补充)

**包含格式**: 15 种

- Parquet, Feather
- Loom, CRAM, PAF
- MAF, WIG, BigBed
- GenBank, mzML, mzXML, MGF
- 等

**预计工作量**: 2-3 周

---

## 📦 转换能力矩阵扩展

### 自动转换路径 (Auto Conversion Paths)

#### 序列格式转换

```
FASTQ.gz → FASTQ → FASTA → TSV (序列列表)
   ↓          ↓        ↓
 FASTQ     FASTA    CSV
```

#### 基因组区间转换

```
BigWig → bedGraph → BED → TSV/CSV
   ↑         ↓        ↓
  WIG      BED6    BED12 → GTF
```

#### 比对格式转换

```
CRAM → BAM → SAM → TSV (统计表)
              ↓
            BED (覆盖度)
```

#### 变异格式转换

```
BCF → VCF → TSV/CSV
       ↓
      BED (变异位置)
```

#### 表达矩阵转换

```
10X HDF5 → h5ad (AnnData) → CSV
    ↓           ↓              ↓
  MTX        RDS (Seurat)    TSV
```

### 交互式转换流程

#### GWAS 标准化流程

```
用户上传 → 列映射界面 → 坐标系统选择 → 单位转换 → 标准化输出
         ↓                ↓                ↓            ↓
      预览50行      hg19/hg38/hg37    beta/OR     标准GWAS TSV
```

#### 单细胞导入流程

```
文件识别 → 物种选择 → 基因版本 → QC阈值 → 对象创建
   ↓           ↓          ↓         ↓         ↓
10X/h5ad   human/mouse  GRCh38   min_cells  Seurat/h5ad
```

---

## 🛠️ 技术实现策略

### 自动转换工具依赖

| 格式转换        | 依赖工具        | 语言          |
| --------------- | --------------- | ------------- |
| FASTQ 处理      | BioPython       | Python        |
| BED/GTF/GFF     | pybedtools      | Python        |
| BAM/SAM         | pysam, samtools | Python/Binary |
| VCF/BCF         | pysam, bcftools | Python/Binary |
| BigWig/bedGraph | pyBigWig        | Python        |
| HDF5            | h5py            | Python        |
| MTX             | scipy           | Python        |
| RData           | rpy2            | Python/R      |

### 交互式转换 UI 组件

1. **ColumnMappingWizard** - 列映射向导
2. **GenomeVersionSelector** - 基因组版本选择器
3. **QCThresholdPanel** - QC 阈值面板
4. **DataPreviewTable** - 数据预览表格
5. **ConversionParameterForm** - 转换参数表单
6. **ProgressTracker** - 转换进度追踪器

---

## 📈 预期效果

### 覆盖率提升

- **当前**: 7 种格式 (通用表格 + Python/R)
- **Phase 1**: 27 种格式 (+20 种生信格式)
- **Phase 2**: 37 种场景 (+10 种交互式)
- **Phase 3**: 52 种格式 (+15 种补充格式)

### 用户体验改善

- ✅ 支持 90%+ 常见生物信息学格式
- ✅ 自动识别文件类型
- ✅ 智能推荐转换路径
- ✅ 交互式参数配置
- ✅ 实时转换预览

### 工作流加速

- ⚡ 减少手动格式转换时间 80%
- ⚡ 避免格式不兼容错误 95%
- ⚡ 支持跨工具无缝数据流
- ⚡ 自动化最佳实践应用

---

**文档版本**: 1.0.0  
**创建日期**: 2025-01-10  
**下一步**: 更新 TODO 列表，开始 Phase 1 实现
