# 组学分析工具变体库（Omics Analysis Tool Variants）

本文档定义了组学分析流程中所有关键步骤的工具变体选择。

## 工具变体系统架构

每个分析功能（Function）都有多个工具变体（Tool Variants），用户可以根据：

- 数据类型（RNA-seq、单细胞、蛋白质组等）
- 分析目标（速度、精度、可视化等）
- 个人偏好（熟悉的工具、编程语言）

自由选择最适合的工具。

---

## 1. 数据预处理（Data Preprocessing）

### 1.1 质量控制（Quality Control）

| Tool Variant  | Runtime | Language | Data Type           | Strengths            | Use Case                |
| ------------- | ------- | -------- | ------------------- | -------------------- | ----------------------- |
| **Seurat QC** | R       | R        | Single-cell RNA-seq | 集成化、可视化丰富   | 单细胞 RNA-seq 标准流程 |
| **scanpy QC** | Python  | Python   | Single-cell RNA-seq | 速度快、内存效率高   | 大规模单细胞数据        |
| **FastQC**    | Binary  | Java     | Raw sequencing      | 通用性强、报告详细   | 原始测序数据评估        |
| **MultiQC**   | Python  | Python   | Multiple data types | 整合多个 QC 报告     | 批量样本 QC 汇总        |
| **Picard QC** | Binary  | Java     | Genomic alignment   | 深度统计、标准化指标 | BAM/SAM 文件 QC         |
| **RSeQC**     | Python  | Python   | RNA-seq             | RNA-seq 专用指标     | 转录组特异性 QC         |

### 1.2 过滤和归一化（Filtering & Normalization）

| Tool Variant         | Runtime | Language | Method                               | Strengths         |
| -------------------- | ------- | -------- | ------------------------------------ | ----------------- |
| **Seurat Normalize** | R       | R        | Log-normalization, SCTransform       | 单细胞标准方法    |
| **scanpy Normalize** | Python  | Python   | Log-normalization, Pearson residuals | 速度快、灵活      |
| **DESeq2 Normalize** | R       | R        | Median-of-ratios                     | Bulk RNA-seq 标准 |
| **edgeR Normalize**  | R       | R        | TMM normalization                    | 差异表达分析      |
| **limma-voom**       | R       | R        | Voom transformation                  | RNA-seq 线性建模  |
| **scran Normalize**  | R       | R        | Deconvolution                        | 单细胞池化归一化  |

### 1.3 批次效应校正（Batch Effect Correction）

| Tool Variant           | Runtime | Language | Method                   | Strengths            |
| ---------------------- | ------- | -------- | ------------------------ | -------------------- |
| **Harmony**            | R       | R        | Iterative clustering     | 速度快、保留生物差异 |
| **BBKNN**              | Python  | Python   | Graph-based              | 单细胞集成           |
| **scVI**               | Python  | Python   | Deep learning            | 非线性校正           |
| **Combat**             | R       | R        | Empirical Bayes          | Bulk 数据经典方法    |
| **MNN**                | R       | R        | Mutual nearest neighbors | 保留罕见细胞类型     |
| **Seurat Integration** | R       | R        | Anchors-based            | 多数据集整合         |

---

## 2. 降维与可视化（Dimensionality Reduction & Visualization）

### 2.1 降维（Dimensionality Reduction）

| Tool Variant       | Runtime | Language | Method               | Strengths            |
| ------------------ | ------- | -------- | -------------------- | -------------------- |
| **Seurat PCA**     | R       | R        | PCA                  | 线性、可解释         |
| **scanpy PCA**     | Python  | Python   | PCA                  | 高效、可扩展         |
| **UMAP (R)**       | R       | R        | UMAP                 | 非线性、保留全局结构 |
| **UMAP (Python)**  | Python  | Python   | UMAP                 | 原始实现、速度快     |
| **t-SNE (R)**      | R       | R        | t-SNE                | 局部结构保留         |
| **t-SNE (Python)** | Python  | Python   | t-SNE                | Barnes-Hut 加速      |
| **PHATE**          | Python  | Python   | Manifold learning    | 保留轨迹结构         |
| **DiffusionMap**   | R       | R        | Diffusion pseudotime | 发育轨迹             |

### 2.2 聚类（Clustering）

| Tool Variant          | Runtime | Language | Method               | Strengths            |
| --------------------- | ------- | -------- | -------------------- | -------------------- |
| **Seurat Clustering** | R       | R        | Louvain/Leiden       | 社区检测、分辨率可调 |
| **scanpy Clustering** | Python  | Python   | Louvain/Leiden       | 高效、可扩展         |
| **SC3**               | R       | R        | Consensus clustering | 稳健性高             |
| **DBSCAN**            | Python  | Python   | Density-based        | 识别噪声点           |
| **Hierarchical**      | R       | R        | Hierarchical         | 可视化层次结构       |
| **k-means**           | Python  | Python   | Centroid-based       | 简单快速             |

### 2.3 细胞类型注释（Cell Type Annotation）

| Tool Variant   | Runtime | Language | Method            | Strengths        |
| -------------- | ------- | -------- | ----------------- | ---------------- |
| **SingleR**    | R       | R        | Reference-based   | 自动化、准确     |
| **CellTypist** | Python  | Python   | Machine learning  | 大规模参考数据   |
| **SCINA**      | R       | R        | Marker-based      | 用户定义标记基因 |
| **Azimuth**    | R       | R        | Reference mapping | 人类细胞图谱     |
| **scmap**      | R       | R        | Nearest neighbor  | 快速投影         |
| **Garnett**    | R       | R        | Marker hierarchy  | 层次化注释       |

---

## 3. 差异分析（Differential Analysis）

### 3.1 差异表达分析（Differential Expression）

| Tool Variant           | Runtime | Language | Method            | Strengths           |
| ---------------------- | ------- | -------- | ----------------- | ------------------- |
| **DESeq2**             | R       | R        | Negative binomial | Bulk RNA-seq 金标准 |
| **edgeR**              | R       | R        | Negative binomial | 灵活、多种模型      |
| **limma**              | R       | R        | Linear modeling   | 微阵列和 RNA-seq    |
| **Seurat FindMarkers** | R       | R        | Wilcoxon/t-test   | 单细胞快速标记      |
| **scanpy rank_genes**  | Python  | Python   | t-test/Wilcoxon   | 高效、并行化        |
| **MAST**               | R       | R        | Hurdle model      | 单细胞零膨胀模型    |
| **scDE**               | R       | R        | Bayesian          | 单细胞贝叶斯方法    |

### 3.2 富集分析（Enrichment Analysis）

| Tool Variant        | Runtime | Language | Method                | Strengths        |
| ------------------- | ------- | -------- | --------------------- | ---------------- |
| **clusterProfiler** | R       | R        | ORA/GSEA              | 全面、可视化好   |
| **GSEA**            | Binary  | Java     | Gene Set Enrichment   | 经典方法         |
| **enrichR**         | R       | R        | Multiple databases    | 数据库丰富       |
| **gseapy**          | Python  | Python   | GSEA Python           | Python 生态集成  |
| **DAVID**           | Web API | -        | Functional annotation | Web 服务、更新快 |
| **Metascape**       | Web API | -        | Meta-analysis         | 整合分析         |

### 3.3 通路分析（Pathway Analysis）

| Tool Variant       | Runtime | Language | Database          | Strengths       |
| ------------------ | ------- | -------- | ----------------- | --------------- |
| **pathview**       | R       | R        | KEGG              | KEGG 通路可视化 |
| **ReactomePA**     | R       | R        | Reactome          | Reactome 数据库 |
| **graphite**       | R       | R        | Multiple          | 多数据库整合    |
| **piano**          | R       | R        | Gene set analysis | 方向性通路分析  |
| **GSEApy Enrichr** | Python  | Python   | Multiple          | Python 通路分析 |

---

## 4. 高级分析（Advanced Analysis）

### 4.1 轨迹推断（Trajectory Inference）

| Tool Variant  | Runtime | Language | Method                 | Strengths    |
| ------------- | ------- | -------- | ---------------------- | ------------ |
| **Monocle 3** | R       | R        | UMAP + principal graph | 复杂分支轨迹 |
| **Slingshot** | R       | R        | Cluster-based          | 简单直观     |
| **PAGA**      | Python  | Python   | Graph abstraction      | 全局拓扑     |
| **Palantir**  | Python  | Python   | Diffusion maps         | 分化潜能     |
| **Velocyto**  | Python  | Python   | RNA velocity           | 方向性推断   |
| **scVelo**    | Python  | Python   | Dynamical modeling     | RNA 速率模型 |

### 4.2 基因调控网络（Gene Regulatory Networks）

| Tool Variant   | Runtime  | Language | Method            | Strengths       |
| -------------- | -------- | -------- | ----------------- | --------------- |
| **SCENIC**     | R/Python | Both     | Regulon analysis  | TF-target 网络  |
| **GENIE3**     | R        | R        | Random forest     | 无监督、通用    |
| **GRNBoost2**  | Python   | Python   | Gradient boosting | 速度快          |
| **CellOracle** | Python   | Python   | Perturbation      | 扰动预测        |
| **Pando**      | R        | R        | Multiome          | ATAC + RNA 整合 |

### 4.3 细胞通讯（Cell-Cell Communication）

| Tool Variant    | Runtime | Language | Method              | Strengths         |
| --------------- | ------- | -------- | ------------------- | ----------------- |
| **CellChat**    | R       | R        | Network analysis    | 网络可视化强      |
| **CellPhoneDB** | Python  | Python   | Statistical testing | 配体-受体数据库大 |
| **NicheNet**    | R       | R        | Ligand activity     | 配体活性预测      |
| **LIANA**       | Python  | Python   | Multiple methods    | 整合多种方法      |
| **Connectome**  | R       | R        | Weighted networks   | 加权网络          |

---

## 5. 蛋白质组学（Proteomics）

### 5.1 质谱数据处理（MS Data Processing）

| Tool Variant            | Runtime | Language | Method              | Strengths    |
| ----------------------- | ------- | -------- | ------------------- | ------------ |
| **MaxQuant**            | Binary  | C#       | Label-free/SILAC    | 工业标准     |
| **Proteome Discoverer** | Binary  | -        | Thermo data         | Thermo 专用  |
| **MSFragger**           | Binary  | Java     | Ultra-fast search   | 速度极快     |
| **OpenMS**              | Binary  | C++      | Open-source         | 开源、可定制 |
| **Skyline**             | Binary  | C#       | Targeted proteomics | 靶向定量     |

### 5.2 蛋白质组统计分析（Proteomics Statistical Analysis）

| Tool Variant | Runtime | Language | Method                | Strengths    |
| ------------ | ------- | -------- | --------------------- | ------------ |
| **MSstats**  | R       | R        | Linear mixed model    | 统计严谨     |
| **limma**    | R       | R        | Linear modeling       | 灵活、成熟   |
| **proDA**    | R       | R        | Probabilistic dropout | 处理缺失值   |
| **Perseus**  | Binary  | C#       | MaxQuant 配套         | 整合分析     |
| **DEqMS**    | R       | R        | Variance modeling     | 蛋白质组专用 |

---

## 6. 表观基因组学（Epigenomics）

### 6.1 ChIP-seq 分析（ChIP-seq Analysis）

| Tool Variant   | Runtime | Language | Method               | Strengths    |
| -------------- | ------- | -------- | -------------------- | ------------ |
| **MACS2**      | Python  | Python   | Peak calling         | 经典方法     |
| **HOMER**      | Binary  | Perl     | Motif analysis       | Motif 发现强 |
| **ChIPseeker** | R       | R        | Annotation           | 峰注释可视化 |
| **DiffBind**   | R       | R        | Differential binding | 差异结合分析 |
| **SICER**      | Python  | Python   | Broad peaks          | 宽峰检测     |

### 6.2 ATAC-seq 分析（ATAC-seq Analysis）

| Tool Variant | Runtime | Language | Method              | Strengths        |
| ------------ | ------- | -------- | ------------------- | ---------------- |
| **ArchR**    | R       | R        | Single-cell ATAC    | 单细胞 ATAC 专用 |
| **Signac**   | R       | R        | Seurat framework    | 多组学整合       |
| **SnapATAC** | R       | R        | Dimension reduction | 大规模单细胞     |
| **chromVAR** | R       | R        | TF activity         | 转录因子活性     |
| **MACS2**    | Python  | Python   | Peak calling        | 通用峰检测       |

### 6.3 DNA 甲基化（DNA Methylation）

| Tool Variant  | Runtime | Language | Method                   | Strengths   |
| ------------- | ------- | -------- | ------------------------ | ----------- |
| **minfi**     | R       | R        | 450K/EPIC array          | 微阵列标准  |
| **RnBeads**   | R       | R        | Comprehensive            | 全流程分析  |
| **methylKit** | R       | R        | Bisulfite-seq            | WGBS 分析   |
| **Bismark**   | Binary  | Perl     | Alignment                | BS-seq 比对 |
| **DSS**       | R       | R        | Differential methylation | 差异甲基化  |

---

## 7. 基因组学（Genomics）

### 7.1 变异检测（Variant Calling）

| Tool Variant          | Runtime | Language | Method          | Strengths    |
| --------------------- | ------- | -------- | --------------- | ------------ |
| **GATK**              | Binary  | Java     | HaplotypeCaller | 金标准、严谨 |
| **FreeBayes**         | Binary  | C++      | Bayesian        | 快速、灵活   |
| **SAMtools/BCFtools** | Binary  | C        | Pileup-based    | 轻量、快速   |
| **Platypus**          | Binary  | Python/C | Local assembly  | 处理复杂区域 |
| **DeepVariant**       | Binary  | Python   | Deep learning   | AI 驱动      |

### 7.2 变异注释（Variant Annotation）

| Tool Variant   | Runtime | Language | Database      | Strengths    |
| -------------- | ------- | -------- | ------------- | ------------ |
| **ANNOVAR**    | Binary  | Perl     | Multiple      | 数据库全     |
| **SnpEff**     | Binary  | Java     | Multiple      | 功能预测     |
| **VEP**        | Binary  | Perl     | Ensembl       | Ensembl 集成 |
| **Funcotator** | Binary  | Java     | GATK 配套     | GATK 生态    |
| **CADD**       | Web API | -        | Pathogenicity | 致病性评分   |

### 7.3 CNV 检测（Copy Number Variation）

| Tool Variant   | Runtime | Language | Method           | Strengths    |
| -------------- | ------- | -------- | ---------------- | ------------ |
| **inferCNV**   | R       | R        | Single-cell      | 单细胞 CNV   |
| **CopywriteR** | R       | R        | Off-target reads | DNA-seq CNV  |
| **CNVkit**     | Python  | Python   | Targeted/WGS     | 灵活、可视化 |
| **GISTIC2**    | Binary  | MATLAB   | Recurrent CNV    | 肿瘤显著 CNV |
| **FACETS**     | R       | R        | Tumor purity     | 肿瘤纯度校正 |

---

## 8. 代谢组学（Metabolomics）

### 8.1 代谢物鉴定（Metabolite Identification）

| Tool Variant      | Runtime | Language | Method               | Strengths      |
| ----------------- | ------- | -------- | -------------------- | -------------- |
| **XCMS**          | R       | R        | Peak detection       | LC-MS 标准     |
| **MZmine**        | Binary  | Java     | Feature detection    | GUI 友好       |
| **MS-DIAL**       | Binary  | C#       | Untargeted           | 非靶向代谢组   |
| **MetaboAnalyst** | Web API | R        | Web-based            | Web 平台、易用 |
| **SIRIUS**        | Binary  | Java     | Structure prediction | 结构预测       |

### 8.2 代谢通路分析（Metabolic Pathway Analysis）

| Tool Variant      | Runtime | Language | Database         | Strengths  |
| ----------------- | ------- | -------- | ---------------- | ---------- |
| **MetaboAnalyst** | Web API | R        | KEGG/HMDB        | 整合分析   |
| **Mummichog**     | Python  | Python   | Network analysis | 网络富集   |
| **IMPaLA**        | Web API | -        | Integrated       | 整合多组学 |
| **FELLA**         | R       | R        | Graph-based      | 图论方法   |

---

## 9. 多组学整合（Multi-omics Integration）

### 9.1 多组学数据整合（Multi-omics Data Integration）

| Tool Variant   | Runtime  | Language | Method                | Strengths      |
| -------------- | -------- | -------- | --------------------- | -------------- |
| **MOFA**       | R/Python | Both     | Factor analysis       | 无监督因子分析 |
| **mixOmics**   | R        | R        | PLS-based             | 多种集成方法   |
| **DIABLO**     | R        | R        | Multi-block PLS       | 监督学习       |
| **Seurat WNN** | R        | R        | Weighted NN           | 多模态单细胞   |
| **totalVI**    | Python   | Python   | Deep learning         | 蛋白+转录组    |
| **MultiVI**    | Python   | Python   | Variational inference | 多模态 VAE     |

---

## 10. 机器学习与预测（Machine Learning & Prediction）

### 10.1 特征选择（Feature Selection）

| Tool Variant     | Runtime | Language | Method             | Strengths  |
| ---------------- | ------- | -------- | ------------------ | ---------- |
| **Boruta**       | R       | R        | Random Forest      | 特征重要性 |
| **mRMR**         | R       | R        | Mutual information | 最小冗余   |
| **scikit-learn** | Python  | Python   | Multiple methods   | 灵活、全面 |
| **glmnet**       | R       | R        | Lasso/Ridge        | 正则化回归 |

### 10.2 分类与预测（Classification & Prediction）

| Tool Variant              | Runtime  | Language | Method            | Strengths      |
| ------------------------- | -------- | -------- | ----------------- | -------------- |
| **caret**                 | R        | R        | Multiple models   | R 机器学习框架 |
| **scikit-learn**          | Python   | Python   | Multiple models   | Python ML 标准 |
| **XGBoost**               | R/Python | Both     | Gradient boosting | 性能强         |
| **Random Forest**         | R/Python | Both     | Ensemble          | 稳健           |
| **SVM**                   | R/Python | Both     | Support vector    | 小样本         |
| **Deep Learning (Keras)** | Python   | Python   | Neural networks   | 复杂模式       |

---

## 11. 可视化（Visualization）

### 11.1 静态可视化（Static Visualization）

| Tool Variant        | Runtime | Language | Output  | Strengths         |
| ------------------- | ------- | -------- | ------- | ----------------- |
| **ggplot2**         | R       | R        | PDF/PNG | Publication-ready |
| **ComplexHeatmap**  | R       | R        | PDF/PNG | 复杂注释热图      |
| **EnrichedHeatmap** | R       | R        | PDF     | 基因组区域热图    |
| **pheatmap**        | R       | R        | PDF/PNG | 简单热图          |
| **seaborn**         | Python  | Python   | PNG     | 统计可视化        |
| **matplotlib**      | Python  | Python   | PNG     | Python 基础绘图   |

### 11.2 交互式可视化（Interactive Visualization）

| Tool Variant  | Runtime  | Language | Output  | Strengths      |
| ------------- | -------- | -------- | ------- | -------------- |
| **plotly**    | Python/R | Both     | HTML    | 交互式图表     |
| **Shiny**     | R        | R        | Web App | R 交互应用     |
| **Dash**      | Python   | Python   | Web App | Python 仪表板  |
| **Bokeh**     | Python   | Python   | HTML    | 交互式可视化   |
| **CellxGene** | Python   | Python   | Web     | 单细胞数据浏览 |

### 11.3 网络可视化（Network Visualization）

| Tool Variant   | Runtime  | Language | Method           | Strengths   |
| -------------- | -------- | -------- | ---------------- | ----------- |
| **igraph**     | R/Python | Both     | Network analysis | 图论分析    |
| **networkx**   | Python   | Python   | Network analysis | Python 网络 |
| **Cytoscape**  | Binary   | Java     | Interactive      | GUI 强大    |
| **visNetwork** | R        | R        | Interactive web  | R 交互网络  |

---

## 实现策略

### 数据库设计

```sql
CREATE TABLE tool_functions (
    id UUID PRIMARY KEY,
    function_name VARCHAR(255) UNIQUE NOT NULL,  -- 如 "quality_control"
    display_name VARCHAR(255) NOT NULL,          -- 如 "Quality Control"
    category VARCHAR(100) NOT NULL,              -- 如 "preprocessing"
    description TEXT,
    data_types VARCHAR[] NOT NULL,               -- 适用的数据类型
    created_at TIMESTAMP DEFAULT NOW()
);

CREATE TABLE tool_variants (
    id UUID PRIMARY KEY,
    function_id UUID REFERENCES tool_functions(id),
    tool_id VARCHAR(255) NOT NULL,
    tool_name VARCHAR(255) NOT NULL,
    runtime VARCHAR(50) NOT NULL,                -- r, python, binary
    language VARCHAR(50),
    method VARCHAR(255),
    strengths TEXT,
    use_case TEXT,
    popularity_score INTEGER DEFAULT 0,
    tool_definition JSONB NOT NULL,
    created_at TIMESTAMP DEFAULT NOW(),
    UNIQUE(function_id, tool_id)
);

CREATE TABLE tool_compatibility (
    id UUID PRIMARY KEY,
    tool_variant_id UUID REFERENCES tool_variants(id),
    input_format VARCHAR(50) NOT NULL,
    output_format VARCHAR(50) NOT NULL,
    conversion_required BOOLEAN DEFAULT false
);
```

### API Endpoints

```
# 获取所有分析功能
GET /api/functions
Response: [{"id": "...", "name": "quality_control", "display_name": "Quality Control", ...}]

# 获取某功能的所有工具变体
GET /api/functions/{function_name}/variants
Response: [{"tool_id": "seurat_qc", "runtime": "r", ...}, {"tool_id": "scanpy_qc", ...}]

# 比较多个工具变体
POST /api/tools/compare
Body: {"tool_ids": ["seurat_qc", "scanpy_qc", "fastqc"]}
Response: {comparison table with features, pros/cons, performance}

# 获取推荐工具
GET /api/tools/recommend?function=quality_control&data_type=single_cell&runtime_preference=python
Response: {"recommended": ["scanpy_qc", "squidpy_qc"], "reasons": [...]}

# 检查工具兼容性
GET /api/tools/compatibility?from_tool=seurat_qc&to_tool=scanpy_hvg
Response: {"compatible": true, "conversion_needed": {"format": "rds_to_h5ad", "time_estimate": "15s"}}
```

### 前端 UI 组件

1. **FunctionSelector** - 选择分析功能
2. **ToolVariantPicker** - 为选定功能选择工具变体
3. **ToolComparisonModal** - 对比多个工具
4. **ToolRecommendationPanel** - 显示推荐工具
5. **CompatibilityChecker** - 检查工具链兼容性

---

## 成功指标

1. **覆盖率**：支持 50+ 分析功能，每个功能 3-5 个工具选项
2. **灵活性**：用户可自由混搭不同运行时的工具
3. **智能性**：推荐系统准确率 > 85%
4. **兼容性**：自动格式转换成功率 > 95%
5. **性能**：工具切换不增加 > 10% 执行时间
