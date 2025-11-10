# 跨运行时工具集成系统（Cross-Runtime Tool Integration）

## 1. 问题现状

### 1.1 当前系统的局限性

**❌ 可视化工具缺失**

- 只有内置的简单图表（Line, Bar, Scatter 等）
- 缺少专业可视化工具集成：
  - R: ggplot2, EnrichedHeatmap, pheatmap, ComplexHeatmap
  - Python: seaborn, plotly, matplotlib, bokeh
  - JavaScript: D3.js, ECharts, Plotly.js

**❌ 工具选择单一**

- 每个功能只能用一个固定的工具
- 用户无法根据需求选择不同的工具包
- 例如 QC 阶段：
  - 某些用户习惯 Seurat (R)
  - 某些用户习惯 scanpy (Python)
  - 某些用户习惯 FastQC (binary)

**❌ 跨语言障碍**

- R 工具和 Python 工具无法在同一个 pipeline 中混用
- 数据格式不兼容：
  - R: RDS, RData
  - Python: pickle, h5ad (AnnData)
  - 通用: CSV, TSV, JSON
- 环境隔离导致数据传递困难

### 1.2 用户痛点

**场景 1: 单细胞分析流程**

```
用户想要的流程：
1. QC: Seurat (R) - 因为熟悉 Seurat 的 QC 参数
2. 找可变基因: scanpy (Python) - 因为 scanpy 的算法更快
3. 聚类: Seurat (R) - 因为喜欢 Seurat 的聚类方法
4. 热图可视化: EnrichedHeatmap (R) - 因为需要特定的热图样式

当前系统：❌ 无法实现，因为无法混用 R 和 Python 工具
```

**场景 2: 多组学整合**

```
用户想要的流程：
1. 蛋白质组 QC: MaxQuant (binary) → R 脚本可视化
2. 转录组 QC: FastQC (binary) → MultiQC (Python) 整合报告
3. 代谢组 QC: Python 脚本 → ggplot2 (R) 可视化

当前系统：❌ 无法实现，因为缺少运行时切换和数据格式转换
```

**场景 3: 自定义可视化**

```
用户想要：
- 用 ggplot2 创建 publication-ready 的图表
- 用 ComplexHeatmap 创建复杂的注释热图
- 用 plotly 创建交互式图表

当前系统：❌ 只能用内置的简单图表，无法使用专业可视化库
```

## 2. 解决方案设计

### 2.1 核心架构：多运行时容器化执行

```
┌─────────────────────────────────────────────────────────────┐
│                    Pipeline Execution Engine                 │
└─────────────────────────────────────────────────────────────┘
                              ↓
        ┌─────────────────────┴─────────────────────┐
        ↓                     ↓                      ↓
┌───────────────┐    ┌───────────────┐    ┌───────────────┐
│  R Runtime    │    │ Python Runtime│    │ Binary Runtime│
│  Container    │    │  Container    │    │  Container    │
├───────────────┤    ├───────────────┤    ├───────────────┤
│ - Seurat      │    │ - scanpy      │    │ - FastQC      │
│ - ggplot2     │    │ - seaborn     │    │ - samtools    │
│ - EnrichedHM  │    │ - plotly      │    │ - bcftools    │
│ - pheatmap    │    │ - anndata     │    │ - bedtools    │
└───────────────┘    └───────────────┘    └───────────────┘
        ↓                     ↓                      ↓
┌─────────────────────────────────────────────────────────────┐
│              Data Format Conversion Layer                    │
│  CSV ↔ TSV ↔ RDS ↔ h5ad ↔ pickle ↔ JSON ↔ VCF ↔ BED       │
└─────────────────────────────────────────────────────────────┘
```

### 2.2 工具选择系统（Tool Variant Selection）

**概念：为每个功能提供多个工具选项**

```json
{
  "function": "quality_control",
  "display_name": "Quality Control",
  "variants": [
    {
      "tool_id": "seurat_qc",
      "name": "Seurat QC",
      "runtime": "r",
      "language": "R",
      "description": "Seurat-based quality control for single-cell data",
      "inputs": ["count_matrix"],
      "outputs": ["qc_report", "filtered_matrix"],
      "parameters": {
        "min_features": 200,
        "max_features": 5000,
        "mt_percent_threshold": 10
      }
    },
    {
      "tool_id": "scanpy_qc",
      "name": "scanpy QC",
      "runtime": "python",
      "language": "Python",
      "description": "scanpy-based quality control for single-cell data",
      "inputs": ["count_matrix"],
      "outputs": ["qc_report", "filtered_h5ad"],
      "parameters": {
        "min_genes": 200,
        "max_genes": 5000,
        "pct_counts_mt": 10
      }
    },
    {
      "tool_id": "fastqc",
      "name": "FastQC",
      "runtime": "binary",
      "language": "Java",
      "description": "Quality control for raw sequencing data",
      "inputs": ["fastq_file"],
      "outputs": ["html_report", "zip_results"],
      "parameters": {
        "threads": 4
      }
    }
  ]
}
```

**UI 交互：**

```
Pipeline Editor 中添加 QC 节点时：

┌─────────────────────────────────────┐
│  Select Tool for Quality Control    │
├─────────────────────────────────────┤
│  ○ Seurat QC (R)                   │
│    Single-cell RNA-seq QC using    │
│    Seurat package                   │
│                                     │
│  ○ scanpy QC (Python)              │
│    Single-cell RNA-seq QC using    │
│    scanpy package                   │
│                                     │
│  ○ FastQC (Binary)                 │
│    Raw sequencing data QC          │
└─────────────────────────────────────┘
       [Confirm] [Compare Features]
```

### 2.3 数据格式自动转换

**转换矩阵：**

| From \ To | CSV | RDS | h5ad | pickle | JSON |
| --------- | --- | --- | ---- | ------ | ---- |
| CSV       | -   | ✅  | ✅   | ✅     | ✅   |
| RDS       | ✅  | -   | ✅   | ❌     | ✅   |
| h5ad      | ✅  | ✅  | -    | ✅     | ✅   |
| pickle    | ✅  | ❌  | ✅   | -      | ✅   |
| JSON      | ✅  | ✅  | ✅   | ✅     | -    |

**自动转换引擎：**

```python
class DataFormatConverter:
    """自动数据格式转换"""

    def convert(self, input_file: str, from_format: str, to_format: str) -> str:
        """
        自动转换数据格式

        Args:
            input_file: 输入文件路径
            from_format: 源格式 (csv, rds, h5ad, pickle, json)
            to_format: 目标格式

        Returns:
            转换后的文件路径
        """
        if from_format == to_format:
            return input_file

        converter = self._get_converter(from_format, to_format)
        return converter.convert(input_file)

    def auto_detect_format(self, file_path: str) -> str:
        """自动检测文件格式"""
        # 基于文件扩展名和内容检测
        pass
```

**Pipeline 中的自动转换：**

```
用户构建的 Pipeline:
[Seurat QC (R)] → [scanpy Variable Genes (Python)] → [EnrichedHeatmap (R)]
      ↓ RDS            ↓ h5ad                           ↓ RDS

系统自动插入转换节点:
[Seurat QC (R)] → [RDS→h5ad] → [scanpy] → [h5ad→RDS] → [EnrichedHM (R)]
```

### 2.4 专业可视化工具集成

**2.4.1 R 可视化工具**

```r
# ggplot2 工具定义
{
  "tool_id": "ggplot2_scatter",
  "name": "ggplot2 Scatter Plot",
  "runtime": "r",
  "category": "visualization",
  "script": "
    library(ggplot2)
    data <- read.csv('{{input_file}}')
    p <- ggplot(data, aes(x={{x_column}}, y={{y_column}}, color={{color_column}})) +
      geom_point(size={{point_size}}, alpha={{alpha}}) +
      theme_{{theme}}() +
      labs(title='{{title}}', x='{{x_label}}', y='{{y_label}}')
    ggsave('{{output_file}}', p, width={{width}}, height={{height}}, dpi={{dpi}})
  ",
  "parameters": {
    "x_column": "gene1",
    "y_column": "gene2",
    "color_column": "cluster",
    "point_size": 2,
    "alpha": 0.7,
    "theme": "bw",
    "title": "Gene Expression",
    "width": 10,
    "height": 8,
    "dpi": 300
  }
}
```

```r
# EnrichedHeatmap 工具定义
{
  "tool_id": "enriched_heatmap",
  "name": "EnrichedHeatmap",
  "runtime": "r",
  "category": "visualization",
  "script": "
    library(EnrichedHeatmap)
    library(circlize)
    mat <- read.table('{{input_matrix}}', header=TRUE, row.names=1)
    col_fun <- colorRamp2(c({{min_value}}, 0, {{max_value}}), c('blue', 'white', 'red'))
    EnrichedHeatmap(mat, col=col_fun, name='{{legend_name}}',
                    column_title='{{title}}',
                    top_annotation={{top_annotation}},
                    show_row_names={{show_row_names}})
  ",
  "parameters": {
    "input_matrix": "expression_matrix.txt",
    "min_value": -2,
    "max_value": 2,
    "legend_name": "Expression",
    "title": "Gene Expression Heatmap",
    "show_row_names": true
  }
}
```

**2.4.2 Python 可视化工具**

```python
# seaborn 工具定义
{
  "tool_id": "seaborn_heatmap",
  "name": "seaborn Heatmap",
  "runtime": "python",
  "category": "visualization",
  "script": "
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt

    data = pd.read_csv('{{input_file}}', index_col=0)
    plt.figure(figsize=({{width}}, {{height}}))
    sns.heatmap(data, cmap='{{colormap}}', center={{center}},
                annot={{annot}}, fmt='{{fmt}}',
                cbar_kws={'label': '{{cbar_label}}'})
    plt.title('{{title}}')
    plt.savefig('{{output_file}}', dpi={{dpi}}, bbox_inches='tight')
  ",
  "parameters": {
    "input_file": "matrix.csv",
    "width": 10,
    "height": 8,
    "colormap": "RdBu_r",
    "center": 0,
    "annot": false,
    "fmt": ".2f",
    "cbar_label": "Expression",
    "title": "Heatmap",
    "dpi": 300
  }
}
```

```python
# plotly 交互式图表
{
  "tool_id": "plotly_3d_scatter",
  "name": "Plotly 3D Scatter",
  "runtime": "python",
  "category": "visualization",
  "script": "
    import pandas as pd
    import plotly.express as px

    df = pd.read_csv('{{input_file}}')
    fig = px.scatter_3d(df, x='{{x_col}}', y='{{y_col}}', z='{{z_col}}',
                        color='{{color_col}}', size='{{size_col}}',
                        title='{{title}}', hover_data={{hover_cols}})
    fig.write_html('{{output_file}}')
  ",
  "parameters": {
    "x_col": "PC1",
    "y_col": "PC2",
    "z_col": "PC3",
    "color_col": "cluster",
    "title": "3D PCA Plot"
  }
}
```

### 2.5 运行时容器管理

**Docker 容器配置：**

```yaml
# R Runtime Container
r-runtime:
  image: rocker/tidyverse:4.3.0
  packages:
    - Seurat
    - ggplot2
    - pheatmap
    - ComplexHeatmap
    - EnrichedHeatmap
    - dplyr
    - tidyr
  volumes:
    - ./data:/data
    - ./output:/output

# Python Runtime Container
python-runtime:
  image: python:3.11-slim
  packages:
    - scanpy
    - seaborn
    - plotly
    - anndata
    - pandas
    - numpy
    - scipy
  volumes:
    - ./data:/data
    - ./output:/output

# Binary Runtime Container
binary-runtime:
  image: ubuntu:22.04
  tools:
    - fastqc
    - samtools
    - bcftools
    - bedtools
  volumes:
    - ./data:/data
    - ./output:/output
```

**执行引擎：**

```python
class RuntimeExecutor:
    """多运行时执行引擎"""

    def execute_tool(self, tool_def: dict, inputs: dict) -> dict:
        """
        执行工具

        Args:
            tool_def: 工具定义（包含 runtime 信息）
            inputs: 输入文件和参数

        Returns:
            输出文件路径字典
        """
        runtime = tool_def['runtime']

        if runtime == 'r':
            return self._execute_r_tool(tool_def, inputs)
        elif runtime == 'python':
            return self._execute_python_tool(tool_def, inputs)
        elif runtime == 'binary':
            return self._execute_binary_tool(tool_def, inputs)
        else:
            raise ValueError(f"Unsupported runtime: {runtime}")

    def _execute_r_tool(self, tool_def: dict, inputs: dict) -> dict:
        """在 R 容器中执行工具"""
        script = self._render_script(tool_def['script'], inputs)

        # 运行 R 脚本
        result = subprocess.run(
            ['docker', 'run', '--rm',
             '-v', f"{self.data_dir}:/data",
             '-v', f"{self.output_dir}:/output",
             'r-runtime',
             'Rscript', '-e', script],
            capture_output=True, text=True
        )

        return self._parse_outputs(result, tool_def['outputs'])
```

## 3. 实现计划

### 3.1 Phase 1: 可视化工具库（★★★★★）

**目标：集成专业可视化工具**

1. **R 可视化工具集成**

   - ggplot2 (散点图、箱线图、小提琴图、密度图)
   - pheatmap (热图)
   - ComplexHeatmap (复杂注释热图)
   - EnrichedHeatmap (富集热图)

2. **Python 可视化工具集成**

   - seaborn (统计可视化)
   - plotly (交互式图表)
   - matplotlib (基础绘图)

3. **可视化工具 JSON 定义**
   - 标准化工具定义格式
   - 参数配置界面
   - 预览和导出功能

**交付物：**

- VisualizationToolsPage（可视化工具库页面）
- 15+ 专业可视化工具定义
- 可视化工具选择和配置界面

### 3.2 Phase 2: 工具选择系统（★★★★★）

**目标：为每个功能提供多个工具选项**

1. **工具变体（Variant）系统**

   - 为同一功能定义多个工具实现
   - 例如：QC 功能 → Seurat QC / scanpy QC / FastQC

2. **工具比较功能**

   - 对比不同工具的特性、参数、输出
   - 帮助用户选择最适合的工具

3. **工具推荐引擎**
   - 基于输入数据类型推荐工具
   - 基于分析目标推荐工具组合

**交付物：**

- Tool Variant Selection UI
- Tool Comparison Modal
- Tool Recommendation System

### 3.3 Phase 3: 数据格式转换（★★★★★）

**目标：自动转换不同工具间的数据格式**

1. **格式转换器库**

   - CSV ↔ RDS (R)
   - CSV ↔ h5ad (Python AnnData)
   - RDS ↔ h5ad
   - JSON ↔ all formats

2. **自动转换插入**

   - Pipeline 执行时自动检测格式不匹配
   - 自动插入转换节点
   - 显示转换路径给用户

3. **格式兼容性检查**
   - 验证工具链中的格式兼容性
   - 警告不支持的转换
   - 建议替代方案

**交付物：**

- DataFormatConverter class
- Auto-conversion pipeline middleware
- Format compatibility checker

### 3.4 Phase 4: 多运行时支持（★★★★☆）

**目标：支持 R、Python、Binary 工具混合执行**

1. **容器化运行时**

   - R runtime Docker 镜像（预装 Seurat、ggplot2 等）
   - Python runtime Docker 镜像（预装 scanpy、seaborn 等）
   - Binary runtime Docker 镜像（预装 FastQC、samtools 等）

2. **运行时切换**

   - 根据工具自动选择运行时
   - 运行时间隔离和资源管理
   - 跨运行时数据传递

3. **依赖管理**
   - R 包安装和版本管理
   - Python 包安装和版本管理
   - 二进制工具安装

**交付物：**

- Multi-runtime Docker images
- RuntimeExecutor class
- Dependency management system

### 3.5 Phase 5: 高级功能（★★★☆☆）

1. **自定义脚本支持**

   - 用户上传 R/Python 脚本作为工具
   - 脚本参数自动解析
   - 脚本验证和沙箱执行

2. **工具市场**

   - 社区贡献的工具库
   - 工具评分和评论
   - 一键安装第三方工具

3. **性能优化**
   - 容器缓存和复用
   - 并行执行优化
   - 大数据格式转换优化

## 4. 技术实现细节

### 4.1 数据库 Schema 扩展

```sql
-- 工具变体表
CREATE TABLE tool_variants (
    id UUID PRIMARY KEY,
    function_name VARCHAR(255) NOT NULL,  -- 功能名称，如 "quality_control"
    tool_id VARCHAR(255) NOT NULL,        -- 工具 ID，如 "seurat_qc"
    tool_name VARCHAR(255) NOT NULL,      -- 显示名称
    runtime VARCHAR(50) NOT NULL,         -- r, python, binary
    language VARCHAR(50),                 -- R, Python, Java, etc.
    description TEXT,
    tool_definition JSONB NOT NULL,       -- 完整的工具定义
    popularity_score INTEGER DEFAULT 0,
    created_at TIMESTAMP DEFAULT NOW(),
    UNIQUE(function_name, tool_id)
);

-- 格式转换规则表
CREATE TABLE format_converters (
    id UUID PRIMARY KEY,
    from_format VARCHAR(50) NOT NULL,
    to_format VARCHAR(50) NOT NULL,
    converter_script TEXT NOT NULL,       -- 转换脚本
    runtime VARCHAR(50) NOT NULL,         -- 使用哪个运行时执行转换
    created_at TIMESTAMP DEFAULT NOW(),
    UNIQUE(from_format, to_format)
);

-- 运行时镜像表
CREATE TABLE runtime_images (
    id UUID PRIMARY KEY,
    runtime_name VARCHAR(50) NOT NULL,    -- r, python, binary
    docker_image VARCHAR(255) NOT NULL,
    installed_packages JSONB,             -- 已安装的包列表
    version VARCHAR(50),
    created_at TIMESTAMP DEFAULT NOW()
);
```

### 4.2 API Endpoints

```python
# 工具变体 API
GET    /api/tools/variants?function=quality_control
  # 获取某功能的所有工具变体

POST   /api/tools/variants/compare
  # 比较多个工具变体
  Body: {"tool_ids": ["seurat_qc", "scanpy_qc"]}

GET    /api/tools/recommend?data_type=single_cell&analysis_goal=qc
  # 获取工具推荐

# 格式转换 API
POST   /api/data/convert
  # 转换数据格式
  Body: {
    "file_path": "data.csv",
    "from_format": "csv",
    "to_format": "h5ad"
  }

GET    /api/data/formats/compatibility?from=rds&to=h5ad
  # 检查格式兼容性

# 运行时 API
GET    /api/runtimes
  # 获取所有可用运行时

POST   /api/runtimes/{runtime_name}/execute
  # 在指定运行时执行脚本
  Body: {
    "script": "...",
    "inputs": {...}
  }
```

## 5. 用户体验设计

### 5.1 工具选择界面

```
Pipeline Editor - 添加 QC 节点:

┌──────────────────────────────────────────────────────────┐
│  Quality Control Tools                         [Compare] │
├──────────────────────────────────────────────────────────┤
│  ┌─────────────────────────────────────────────────────┐ │
│  │ ● Seurat QC (R)                     [Most Popular] │ │
│  │   Single-cell RNA-seq quality control              │ │
│  │   Runtime: R 4.3.0 | Package: Seurat 5.0          │ │
│  │   ★★★★★ (1,234 uses)                              │ │
│  └─────────────────────────────────────────────────────┘ │
│                                                           │
│  ┌─────────────────────────────────────────────────────┐ │
│  │ ○ scanpy QC (Python)                               │ │
│  │   Fast quality control for large datasets          │ │
│  │   Runtime: Python 3.11 | Package: scanpy 1.9      │ │
│  │   ★★★★☆ (856 uses)                                │ │
│  └─────────────────────────────────────────────────────┘ │
│                                                           │
│  ┌─────────────────────────────────────────────────────┐ │
│  │ ○ FastQC (Binary)                                  │ │
│  │   Raw sequencing data quality control              │ │
│  │   Runtime: Java | Binary: FastQC 0.12             │ │
│  │   ★★★★☆ (643 uses)                                │ │
│  └─────────────────────────────────────────────────────┘ │
│                                                           │
│  [Add Custom Tool]          [Cancel]  [Add to Pipeline]  │
└──────────────────────────────────────────────────────────┘
```

### 5.2 格式转换可视化

```
Pipeline 执行预览:

┌──────────────────────────────────────────────────────────┐
│  Pipeline: Single-cell Analysis                          │
├──────────────────────────────────────────────────────────┤
│                                                           │
│  [Seurat QC (R)]                                         │
│        ↓ output: qc_results.rds                          │
│        ↓                                                  │
│   🔄 [Auto-convert: RDS → h5ad]  ⚠️ Added automatically │
│        ↓                                                  │
│  [scanpy Variable Genes (Python)]                        │
│        ↓ output: variable_genes.h5ad                     │
│        ↓                                                  │
│   🔄 [Auto-convert: h5ad → RDS]  ⚠️ Added automatically │
│        ↓                                                  │
│  [EnrichedHeatmap (R)]                                   │
│        ↓ output: heatmap.pdf                             │
│                                                           │
│  ℹ️ Conversion nodes were added automatically to ensure  │
│     compatibility between different runtimes.            │
│                                                           │
│  [View Conversion Details]      [Cancel]  [Run Pipeline] │
└──────────────────────────────────────────────────────────┘
```

## 6. 示例用例

### 6.1 混合运行时 Pipeline

```json
{
  "pipeline_name": "Multi-runtime Single-cell Analysis",
  "nodes": [
    {
      "id": "1",
      "tool_id": "seurat_qc",
      "runtime": "r",
      "parameters": {
        "min_features": 200,
        "max_features": 5000
      }
    },
    {
      "id": "2",
      "tool_id": "auto_convert",
      "from_format": "rds",
      "to_format": "h5ad"
    },
    {
      "id": "3",
      "tool_id": "scanpy_highly_variable_genes",
      "runtime": "python",
      "parameters": {
        "n_top_genes": 2000
      }
    },
    {
      "id": "4",
      "tool_id": "auto_convert",
      "from_format": "h5ad",
      "to_format": "csv"
    },
    {
      "id": "5",
      "tool_id": "ggplot2_heatmap",
      "runtime": "r",
      "parameters": {
        "color_scheme": "RdBu"
      }
    }
  ]
}
```

### 6.2 可视化工具使用

```r
# 用户可以直接使用 EnrichedHeatmap
{
  "tool_id": "enriched_heatmap_custom",
  "name": "Gene Region Enrichment Heatmap",
  "runtime": "r",
  "script": "
    library(EnrichedHeatmap)
    library(GenomicRanges)

    # 加载数据
    mat <- read.table('{{input_matrix}}', header=T, row.names=1)

    # 创建热图
    EnrichedHeatmap(
      mat,
      col = colorRamp2(c(-2, 0, 2), c('blue', 'white', 'red')),
      name = 'Enrichment',
      column_title = '{{title}}',
      top_annotation = HeatmapAnnotation(
        enriched = anno_enriched()
      )
    )
  "
}
```

## 7. 成功指标

### 7.1 功能指标

- ✅ 支持 10+ R 可视化工具
- ✅ 支持 10+ Python 可视化工具
- ✅ 支持 5+ 数据格式自动转换
- ✅ 支持 3 种运行时（R, Python, Binary）
- ✅ 工具变体选择 < 3 次点击
- ✅ 格式转换时间 < 30 秒（1GB 数据）

### 7.2 用户体验指标

- ✅ 用户可以自由混用不同语言的工具
- ✅ 90% 的格式转换自动完成
- ✅ 用户满意度 > 4.5/5
- ✅ 工具切换成功率 > 95%

## 8. 风险与挑战

### 8.1 技术风险

**容器资源消耗**

- 问题：多个运行时容器同时运行消耗大量资源
- 解决：容器池管理，复用和缓存容器实例

**格式转换数据丢失**

- 问题：某些格式转换可能丢失元数据
- 解决：转换前验证，支持无损转换路径，提供警告

**运行时版本冲突**

- 问题：不同工具需要不同版本的 R/Python
- 解决：使用独立的 Docker 镜像，版本隔离

### 8.2 用户体验风险

**工具选择困难**

- 问题：太多工具选项让用户困惑
- 解决：智能推荐、流行度排序、工具对比功能

**执行时间增加**

- 问题：格式转换和容器启动增加执行时间
- 解决：并行执行、转换缓存、容器预热

## 9. 总结

通过实现跨运行时工具集成系统，我们可以：

1. **打破语言障碍**：用户可以在同一个 pipeline 中混用 R、Python、Binary 工具
2. **灵活选择工具**：为每个功能提供多个工具选项，用户选择最适合的
3. **专业可视化**：集成 ggplot2、EnrichedHeatmap、seaborn 等专业可视化库
4. **自动化转换**：系统自动处理不同工具间的数据格式转换
5. **容器化执行**：隔离运行环境，避免依赖冲突

这将使 Omicsomics 成为真正的**多语言、多工具、多组学整合平台**。
