# 数据流程与功能需求文档

## 1. 现状分析

### 1.1 当前数据流问题

- ❌ **Run 完成后数据去向不明确**：用户不知道输出文件存储在哪里
- ❌ **缺少数据管理中心**：没有统一的地方查看和管理所有产生的数据
- ❌ **可视化与数据脱节**：Pipeline 中的可视化工具无法直接访问已完成的 Run 数据
- ❌ **数据合并困难**：没有工具来合并多个实验的结果
- ❌ **数据版本混乱**：同一个分析的多次运行结果无法追溯

### 1.2 用户痛点

1. **数据丢失感**："我的分析结果去哪了？"
2. **重复下载**：每次可视化都要重新下载数据
3. **手动合并**：需要手动下载、合并多个 CSV 文件
4. **无法复用**：之前的分析结果无法直接用于新的 Pipeline
5. **存储混乱**：不知道哪些数据可以删除，哪些需要保留

---

## 2. 完整用户流程设计

### 2.1 标准工作流程

```
┌─────────────────────────────────────────────────────────────────┐
│                     用户完整工作流程                              │
└─────────────────────────────────────────────────────────────────┘

1. 数据上传阶段
   ├─ 用户上传原始数据（FASTQ, VCF, CSV等）
   ├─ 数据存储在 Projects/{project_id}/raw_data/
   └─ 自动加密、元数据提取、预览生成

2. Pipeline 构建阶段
   ├─ 在 Custom Pipelines 页面设计分析流程
   ├─ 配置各个工具的参数
   └─ 验证 Pipeline 有效性

3. Pipeline 执行阶段
   ├─ 提交 Run（关联到特定 Project）
   ├─ 实时监控执行状态（RunExecutionVisualizer）
   └─ 每个工具生成中间输出和最终输出

4. 数据产出阶段 ⭐ 核心问题
   ├─ 输出数据自动保存到 Projects/{project_id}/runs/{run_id}/outputs/
   ├─ 每个工具的输出独立存储
   ├─ 生成数据索引和元数据
   └─ ❌ 当前：用户不知道如何访问这些数据

5. 数据管理阶段 ⭐ 缺失功能
   ├─ ❌ 缺少：数据浏览器（Data Browser）
   ├─ ❌ 缺少：数据集管理（Dataset Manager）
   └─ ❌ 缺少：数据版本控制

6. 数据分析阶段 ⭐ 缺失功能
   ├─ ❌ 缺少：直接在平台内可视化历史数据
   ├─ ❌ 缺少：数据比较工具
   └─ ❌ 缺少：交互式探索界面

7. 数据整合阶段 ⭐ 缺失功能
   ├─ ❌ 缺少：多文件合并工具
   ├─ ❌ 缺少：数据转换工具
   └─ ❌ 缺少：批量导出功能

8. 结果分享阶段
   ├─ 下载最终结果
   ├─ 生成分析报告
   └─ 导出到外部工具
```

---

## 3. 核心功能需求

### 3.1 数据浏览器（Data Browser）★★★★★

#### 目标

为用户提供一个集中的地方查看、管理和操作所有项目数据。

#### 功能需求

##### 3.1.1 数据分类视图

```
Data Browser 页面结构：
├─ 📁 Raw Data（原始数据）
│  ├─ 显示所有上传的原始文件
│  ├─ 文件预览（前100行）
│  └─ 元数据（上传时间、大小、格式、MD5）
│
├─ 🔬 Run Outputs（运行输出）
│  ├─ 按 Run 组织
│  ├─ 每个 Run 显示所有输出文件
│  ├─ 按工具分组（FastQC → 质控报告、Trim → 清洗后数据）
│  └─ 状态标记（新数据、已查看、已下载、已删除）
│
├─ 📊 Datasets（数据集）⭐ 新概念
│  ├─ 用户创建的数据集合
│  ├─ 可包含原始数据、Run 输出、合并数据
│  ├─ 支持标签、描述、版本
│  └─ 可直接用于新的 Pipeline
│
└─ 🗑️ Archived（归档数据）
   ├─ 软删除的数据
   ├─ 30天后永久删除
   └─ 可恢复
```

##### 3.1.2 数据操作

- **预览**:
  - CSV/TSV: 表格预览（分页，前 1000 行）
  - JSON: 树形结构展示
  - 图片: 缩略图和全图查看
  - 文本: 语法高亮
- **下载**:
  - 单文件下载
  - 批量打包下载（ZIP）
  - 选择性下载（勾选多个文件）
- **删除**:
  - 软删除（移到归档）
  - 批量删除
  - 确认对话框
- **重命名**:
  - 修改显示名称（保留原始文件名）
  - 添加标签
- **元数据编辑**:
  - 添加描述
  - 设置类别
  - 标记重要性

##### 3.1.3 数据搜索与过滤

- **搜索**:
  - 按文件名搜索
  - 按内容搜索（CSV 列名、JSON 键）
  - 按日期范围搜索
  - 按大小搜索
- **过滤**:
  - 按文件类型（CSV、VCF、FASTQ 等）
  - 按产生来源（Run ID、工具名）
  - 按状态（新建、已修改、已使用）
  - 按标签

##### 3.1.4 数据统计

```
数据概览卡片：
┌──────────────────────────────────────────┐
│ 📊 Project Data Summary                  │
├──────────────────────────────────────────┤
│ Total Files: 247                         │
│ Total Size: 15.3 GB                      │
│ Raw Data: 45 files (8.2 GB)            │
│ Run Outputs: 156 files (6.8 GB)        │
│ Datasets: 12 (342 MB)                   │
│ Archived: 34 files (254 MB)            │
└──────────────────────────────────────────┘
```

---

### 3.2 数据集管理（Dataset Manager）★★★★★

#### 目标

允许用户创建、管理和复用数据集，实现数据的组织和版本控制。

#### 核心概念

**Dataset** = 一组相关数据文件的集合，具有：

- 唯一名称和描述
- 包含的文件列表（引用，非复制）
- 元数据（创建时间、作者、标签）
- 版本历史
- 可直接用作 Pipeline 输入

#### 功能需求

##### 3.2.1 创建数据集

```typescript
// 创建数据集的方式
方式1: 从 Run 输出创建
  - 在 Run Details 页面点击"Create Dataset from Run"
  - 自动包含该 Run 的所有输出
  - 建议命名：{pipeline_name}_{run_date}

方式2: 手动选择文件创建
  - 在 Data Browser 勾选多个文件
  - 点击"Create Dataset"
  - 输入名称、描述、标签

方式3: 从现有数据集克隆
  - 复制数据集结构
  - 可修改包含的文件
  - 创建新版本
```

##### 3.2.2 数据集版本管理

```
Dataset: "Proteomics_Batch1_Processed"
├─ v1.0 (2025-01-01)
│  ├─ 初始版本
│  └─ 包含 50 个样本的定量结果
│
├─ v1.1 (2025-01-05)
│  ├─ 添加批次校正后的数据
│  └─ 新增 10 个文件
│
└─ v2.0 (2025-01-10)
   ├─ 移除低质量样本
   └─ 重新归一化
```

##### 3.2.3 数据集操作

- **查看**: 展示所有包含的文件
- **编辑**: 添加/移除文件、修改元数据
- **复制**: 创建副本用于不同分析
- **合并**: 合并多个数据集（见 3.3）
- **导出**: 打包下载所有文件
- **删除**: 软删除（不删除原始文件）

##### 3.2.4 数据集关联

- **关联 Run**: 显示该数据集被哪些 Run 使用
- **关联 Pipeline**: 显示兼容的 Pipeline
- **依赖追踪**: 追溯数据来源链

---

### 3.3 数据合并工具（Data Merger）★★★★★

#### 目标

提供强大的工具来合并多个数据文件，支持各种组学数据格式。

#### 功能需求

##### 3.3.1 文件合并模式

###### 模式 1: 垂直合并（Vertical/Row Concatenation）

```
适用场景: 多批次样本合并
┌─────────────┐
│ File 1      │  Sample1, Sample2 (50 rows)
├─────────────┤
│ File 2      │  Sample3, Sample4 (50 rows)
├─────────────┤
│ File 3      │  Sample5, Sample6 (50 rows)
└─────────────┘
        ↓ 合并
┌─────────────┐
│ Merged      │  Sample1-6 (150 rows)
└─────────────┘

配置:
- 要求: 列名必须一致
- 选项:
  - 处理重复行（去重/保留）
  - 添加来源列（标记来自哪个文件）
  - 重置索引
```

###### 模式 2: 水平合并（Horizontal/Column Concatenation）

```
适用场景: 添加新的测量维度
┌──────┬──────┐
│Gene  │Expr  │  (100 genes × 2 cols)
└──────┴──────┘
        +
┌──────┬──────┐
│Gene  │Prot  │  (100 genes × 2 cols)
└──────┴──────┘
        ↓ 合并
┌──────┬──────┬──────┐
│Gene  │Expr  │Prot  │  (100 genes × 3 cols)
└──────┴──────┴──────┘

配置:
- 连接键: 选择用于匹配的列（Gene ID）
- 连接方式:
  - Inner join（交集）
  - Left join（保留左表）
  - Outer join（并集）
- 处理缺失值: 填充策略（0, NA, mean）
```

###### 模式 3: 智能合并（Smart Merge）

```
适用场景: 复杂的多源数据整合
示例: 整合基因组、转录组、蛋白质组数据
┌─────────────────────────────────────────┐
│ Genomics (VCF) → 提取变异位点          │
│ Transcriptomics (CSV) → 基因表达量     │
│ Proteomics (CSV) → 蛋白丰度            │
└─────────────────────────────────────────┘
        ↓ 智能匹配（通过 Gene Symbol）
┌─────────────────────────────────────────┐
│ Gene | Variants | Expression | Protein  │
└─────────────────────────────────────────┘

配置:
- 自动检测文件格式
- 推荐连接键
- 数据类型转换
- 列重命名建议
```

##### 3.3.2 合并界面设计

```
Data Merger 页面:
┌─────────────────────────────────────────────────────┐
│ 📂 Select Files to Merge                            │
├─────────────────────────────────────────────────────┤
│ [+] Add Files                                        │
│                                                      │
│ Selected Files (3):                                  │
│ ☑ proteomics_batch1.csv (1.2 MB, 1000 rows)       │
│ ☑ proteomics_batch2.csv (1.5 MB, 1200 rows)       │
│ ☑ proteomics_batch3.csv (980 KB, 950 rows)        │
│                                                      │
│ Preview:                                             │
│ ┌─────────────────────────────────────────────┐    │
│ │ Sample  │ Protein │ Abundance │ Batch       │    │
│ │ S001    │ P12345  │ 1234.5    │ Batch1      │    │
│ │ ...                                          │    │
│ └─────────────────────────────────────────────┘    │
└─────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────┐
│ ⚙️ Merge Configuration                              │
├─────────────────────────────────────────────────────┤
│ Merge Mode:  ◉ Vertical  ○ Horizontal  ○ Smart     │
│                                                      │
│ [Vertical Mode Options]                              │
│ ☑ Add source column (column name: "Batch")         │
│ ☑ Check for duplicate rows                          │
│ ☐ Reset row index                                   │
│                                                      │
│ Column Mapping:                                      │
│ File1.col → Merged.col                              │
│ "Protein_ID" → "Protein"  [Auto-detected ✓]        │
│ "Abundance" → "Abundance" [Auto-detected ✓]        │
│                                                      │
│ Output Options:                                      │
│ Name: proteomics_all_batches                        │
│ Format: ○ CSV  ◉ TSV  ○ Excel  ○ Parquet          │
│ ☑ Save as Dataset                                   │
│ ☑ Compress (gzip)                                   │
└─────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────┐
│ 📊 Merge Preview                                     │
├─────────────────────────────────────────────────────┤
│ Input Summary:                                       │
│ - Total files: 3                                     │
│ - Total rows: 3,150                                  │
│ - Total columns: 4 (after deduplication)            │
│                                                      │
│ Output Preview:                                      │
│ - Expected rows: 3,150                              │
│ - Expected columns: 5 (original 4 + 1 batch col)   │
│ - Estimated size: 3.8 MB                            │
│ - Warnings: None                                     │
│                                                      │
│ [Preview Output (first 10 rows)] [Merge] [Cancel]  │
└─────────────────────────────────────────────────────┘
```

##### 3.3.3 高级功能

- **批量合并**: 一次合并 10+ 个文件
- **分组合并**: 先按条件分组，再分别合并
- **条件过滤**: 合并时同时过滤数据
- **数据转换**: 合并时进行单位转换、标准化
- **质量检查**: 自动检测数据质量问题

---

### 3.4 数据可视化集成（Visualization Integration）★★★★☆

#### 目标

让用户可以直接在平台内对历史数据进行可视化，无需重新运行 Pipeline。

#### 功能需求

##### 3.4.1 快速可视化

```
在 Data Browser 中:
任意数据文件 → 右键 → "Quick Visualize"
    ↓
自动检测数据类型和结构
    ↓
推荐合适的可视化类型:
- 表格数据 → Heatmap, PCA, Volcano, Box Plot
- 时间序列 → Line Chart
- 分类数据 → Bar Chart, Pie Chart
    ↓
打开交互式可视化面板
- 调整参数
- 切换图表类型
- 导出高清图片
```

##### 3.4.2 可视化工作区（Visualization Workspace）

```
新页面: "Visualizations"
功能:
- 多面板布局（Grid Layout）
- 拖拽添加图表
- 每个面板连接一个数据源
- 实时交互（点击、缩放、筛选）
- 保存为"Visualization Dashboard"
- 分享给团队成员

示例布局:
┌──────────────────┬──────────────────┐
│ PCA Plot         │ Heatmap          │
│ (Dataset: Expr)  │ (Dataset: Expr)  │
├──────────────────┼──────────────────┤
│ Volcano Plot     │ Box Plot         │
│ (Dataset: DE)    │ (Dataset: Expr)  │
└──────────────────┴──────────────────┘
```

##### 3.4.3 对比可视化

```
功能: 在一个视图中对比多个数据集
示例:
- 对比不同批次的 QC 指标
- 对比不同处理条件的表达量
- 对比不同时间点的变化

实现:
- 选择 2-4 个数据集
- 选择可视化类型（Overlay, Side-by-side）
- 自动对齐坐标轴
- 颜色编码区分数据集
```

---

### 3.5 Pipeline 与数据集集成 ★★★★☆

#### 目标

让数据集可以无缝用作 Pipeline 的输入，实现数据复用。

#### 功能需求

##### 3.5.1 数据集作为输入源

```
在 Pipeline Editor 中:
Input 节点 → 配置面板
    ↓
数据源选择:
◉ Upload new file
◉ Select from Project Data  ⭐ 新功能
◉ Select from Dataset  ⭐ 新功能
    ↓
[如果选择 "Select from Dataset"]
    ↓
显示当前项目的所有数据集列表
- 显示数据集名称、描述、文件数量
- 筛选兼容的数据集（根据 Pipeline 要求）
- 预览数据集内容
    ↓
选择数据集 → 自动关联到 Input 节点
```

##### 3.5.2 Run 输出自动保存

```
Pipeline 执行完成后:
1. 自动创建数据集
   名称: {pipeline_name}_run_{run_id}
   包含: 所有输出文件

2. 关联到 Project
   路径: Projects/{project_id}/datasets/{dataset_id}

3. 通知用户
   "Run completed! Dataset created: {dataset_name}"
   [View Dataset] [Visualize] [Use in New Pipeline]
```

##### 3.5.3 数据血缘追踪（Data Lineage）

```
显示数据的完整来源链:

Raw Data (upload_20250101.fastq)
    ↓ Run #123 (QC Pipeline)
QC Output (qc_passed.fastq) → Dataset "QC_Batch1"
    ↓ Run #124 (Alignment Pipeline)
BAM File (aligned.bam) → Dataset "Aligned_Batch1"
    ↓ Run #125 (Variant Calling)
VCF File (variants.vcf) → Dataset "Variants_Batch1"
    ↓ Manual Merge
Combined VCF (all_variants.vcf) → Dataset "Final_Variants"

功能:
- 点击任意数据查看来源
- 追溯到原始上传文件
- 查看中间处理步骤
- 重现分析流程
```

---

### 3.6 数据导出与分享 ★★★☆☆

#### 功能需求

##### 3.6.1 批量导出

- 选择多个文件/数据集
- 打包为 ZIP 或 TAR.GZ
- 包含元数据文件（JSON）
- 生成 README（说明数据来源）

##### 3.6.2 格式转换导出

- CSV ↔ TSV ↔ Excel
- VCF → BED
- FASTQ → FASTA
- 自定义分隔符

##### 3.6.3 云端导出

- 导出到 AWS S3
- 导出到 Google Drive
- 导出到本地文件系统（Docker volume）

---

## 4. 技术实现要点

### 4.1 数据存储结构

```
MinIO Bucket: omicsomics-data
├─ projects/
│  └─ {project_id}/
│     ├─ raw_data/
│     │  └─ {file_id}/{filename}
│     ├─ runs/
│     │  └─ {run_id}/
│     │     ├─ inputs/
│     │     ├─ outputs/
│     │     │  └─ {tool_id}/{output_file}
│     │     └─ metadata.json
│     └─ datasets/
│        └─ {dataset_id}/
│           ├─ manifest.json (文件列表)
│           └─ metadata.json
```

### 4.2 数据库 Schema

```sql
-- 数据集表
CREATE TABLE datasets (
    id UUID PRIMARY KEY,
    project_id UUID REFERENCES projects(id),
    name VARCHAR(255) NOT NULL,
    description TEXT,
    version VARCHAR(50),
    created_by UUID REFERENCES users(id),
    created_at TIMESTAMP,
    updated_at TIMESTAMP,
    tags TEXT[],
    file_count INTEGER,
    total_size BIGINT
);

-- 数据集文件关联表
CREATE TABLE dataset_files (
    id UUID PRIMARY KEY,
    dataset_id UUID REFERENCES datasets(id),
    file_id UUID,
    file_path TEXT,
    file_type VARCHAR(50),
    file_size BIGINT,
    added_at TIMESTAMP,
    source_run_id UUID REFERENCES runs(id) NULL
);

-- 数据血缘表
CREATE TABLE data_lineage (
    id UUID PRIMARY KEY,
    child_file_id UUID,
    parent_file_id UUID,
    transformation_type VARCHAR(100),
    run_id UUID REFERENCES runs(id),
    created_at TIMESTAMP
);
```

### 4.3 API 端点

```
# 数据浏览器
GET    /api/projects/{project_id}/data/browser
GET    /api/projects/{project_id}/data/files
GET    /api/projects/{project_id}/data/files/{file_id}/preview
DELETE /api/projects/{project_id}/data/files/{file_id}

# 数据集管理
GET    /api/projects/{project_id}/datasets
POST   /api/projects/{project_id}/datasets
GET    /api/datasets/{dataset_id}
PUT    /api/datasets/{dataset_id}
DELETE /api/datasets/{dataset_id}
POST   /api/datasets/{dataset_id}/files

# 数据合并
POST   /api/data/merge
POST   /api/data/merge/preview

# 数据可视化
POST   /api/data/visualize
GET    /api/visualizations/{viz_id}

# 数据导出
POST   /api/data/export
GET    /api/data/export/{export_id}/download
```

---

## 5. 优先级与开发阶段

### Phase 1: 数据可见性（必须）★★★★★

- [ ] 数据浏览器基础版
- [ ] Run 输出文件列表展示
- [ ] 文件预览功能
- [ ] 文件下载功能

### Phase 2: 数据管理（核心）★★★★★

- [ ] 数据集创建和管理
- [ ] 数据集版本控制
- [ ] 数据搜索和过滤
- [ ] 数据统计面板

### Phase 3: 数据整合（重要）★★★★☆

- [ ] 数据合并工具（垂直合并）
- [ ] 数据合并工具（水平合并）
- [ ] 批量操作支持
- [ ] 合并预览功能

### Phase 4: 高级功能（增强）★★★☆☆

- [ ] 快速可视化
- [ ] 可视化工作区
- [ ] 数据血缘追踪
- [ ] Pipeline 与数据集集成

### Phase 5: 企业功能（可选）★★☆☆☆

- [ ] 数据分享和权限
- [ ] 云端导出
- [ ] 数据归档策略
- [ ] 数据质量报告

---

## 6. 成功指标

### 用户体验指标

- 用户能在 **3 次点击内** 找到任何历史数据
- 数据合并操作 **10 分钟内** 完成（包括配置）
- **90%** 的数据操作无需下载到本地
- **零** 数据丢失感投诉

### 技术指标

- 数据浏览器加载时间 < 2 秒
- 文件预览响应时间 < 1 秒
- 支持 **10,000+** 文件的项目
- 合并操作支持 **100MB+** 文件

---

## 7. 用户故事

### 故事 1: 多批次实验合并

```
作为一名生物信息学家
我想合并 3 批次的蛋白质组学数据
以便进行统一的差异分析

步骤:
1. 进入 Data Browser
2. 筛选 "Proteomics" 标签的文件
3. 选择 3 个批次的定量结果文件
4. 点击 "Merge Files"
5. 选择 "Vertical Merge"，添加 "Batch" 列
6. 预览合并结果
7. 保存为新数据集 "Proteomics_All_Batches"
8. 在新 Pipeline 中使用该数据集
```

### 故事 2: 历史数据可视化

```
作为一名研究员
我想查看上个月完成的 RNA-seq 分析结果
无需重新运行整个 Pipeline

步骤:
1. 进入 Data Browser → Run Outputs
2. 找到 Run #123 (RNA-seq)
3. 点击 "View Outputs"
4. 选择 "gene_expression.csv"
5. 点击 "Quick Visualize"
6. 系统推荐: Heatmap, PCA
7. 选择 Heatmap
8. 调整聚类参数
9. 导出高清 PNG
```

### 故事 3: 数据复用

```
作为一名数据分析师
我想使用上次 QC 后的干净数据
进行新的下游分析

步骤:
1. 创建新 Pipeline
2. 添加 Input 节点
3. 选择 "Select from Dataset"
4. 选择 "QC_Passed_Batch1" 数据集
5. 继续构建 Pipeline
6. 提交 Run（自动使用数据集中的文件）
```

---

## 8. 风险与挑战

### 技术挑战

1. **大文件处理**: 如何高效合并 GB 级文件
2. **并发访问**: 多用户同时访问同一数据
3. **存储优化**: 避免数据冗余存储
4. **性能优化**: 快速索引和搜索大量文件

### 解决方案

- 使用流式处理（streaming）处理大文件
- 数据集使用引用而非复制
- 实现增量索引更新
- 使用 Redis 缓存文件元数据

---

## 9. 总结

### 核心价值

1. **数据可见性**: 用户清楚知道数据在哪
2. **数据管理**: 统一的数据管理界面
3. **数据复用**: 避免重复计算和存储
4. **工作效率**: 减少下载、上传、手动合并的时间

### 下一步

1. 创建详细的 UI 设计稿
2. 设计数据库 Schema
3. 实现 API 端点
4. 开发前端组件
5. 集成到现有系统

---

**文档版本**: 1.0  
**创建日期**: 2025-01-10  
**作者**: System  
**状态**: 待评审和实现
