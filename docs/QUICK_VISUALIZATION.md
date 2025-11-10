# Quick Data Visualization

## 概述

快速数据可视化功能允许用户直接从数据浏览器一键生成数据可视化，无需编写代码或配置复杂的管道。系统会自动分析数据结构并推荐最合适的可视化方式。

## 功能特性

### 1. 智能数据类型检测

系统自动识别以下数据类型：

- **Numeric Matrix**: 主要包含数值的矩阵数据
- **Expression Matrix**: 基因表达矩阵或类似的组学数据（大量数值列）
- **Time Series**: 包含时间维度的数据
- **Categorical**: 主要包含分类变量的数据
- **Genomic Intervals**: 基因组坐标数据（chr, start, end）
- **Mixed**: 混合数值和分类数据
- **Unknown**: 无法确定的数据类型

### 2. 支持的图表类型

#### 2.1 Histogram（直方图）

- **用途**: 显示单个数值变量的分布
- **配置**: 选择列、设置分箱数
- **最适合**: 探索数据分布、识别异常值

#### 2.2 Scatter Plot（散点图）

- **用途**: 显示两个数值变量之间的关系
- **配置**: 选择 X 轴和 Y 轴列
- **最适合**: 探索变量相关性、识别模式

#### 2.3 Line Chart（折线图）

- **用途**: 显示趋势或时间序列
- **配置**: 选择要绘制的列（最多 5 条线）
- **最适合**: 时间序列数据、趋势分析

#### 2.4 Bar Chart（条形图）

- **用途**: 比较分类变量的计数
- **配置**: 选择分类列
- **最适合**: 类别比较、频次统计

#### 2.5 Box Plot（箱线图）

- **用途**: 显示数值变量的分布和异常值
- **配置**: 选择多个数值列（最多 10 个）
- **最适合**: 比较多个变量的分布、检测异常值

#### 2.6 Heatmap（热图）

- **用途**: 用颜色可视化矩阵数据
- **配置**: 无需配置（自动选择数值列）
- **最适合**: 表达矩阵、相关性矩阵、高维数据

#### 2.7 Violin Plot（小提琴图）

- **用途**: 显示分布形状和密度估计
- **配置**: 选择数值列
- **最适合**: 比较分布形状、查看数据密度

#### 2.8 Correlation Matrix（相关性矩阵）

- **用途**: 显示变量之间的成对相关性
- **配置**: 无需配置（自动选择数值列）
- **最适合**: 探索变量关系、特征选择

## 使用方法

### 方法 1: 从数据浏览器（推荐）

1. 打开 **Data Browser** 页面
2. 找到想要可视化的文件（支持 CSV, TSV, XLSX）
3. 点击文件行中的 📊 按钮
4. 系统自动分析文件并显示建议的可视化方式
5. 选择一个推荐的图表类型或自定义配置
6. 点击 "Generate Visualization" 生成图表

### 方法 2: 通过 API

```python
import requests

# 1. 分析文件
response = requests.post(
    "http://localhost:8000/api/quick-visualizer/analyze",
    json={"file_path": "/path/to/data.csv"}
)

analysis = response.json()
print(f"Data type: {analysis['data_type']}")
print(f"Suggestions: {analysis['suggestions']}")

# 2. 生成可视化
response = requests.post(
    "http://localhost:8000/api/quick-visualizer/visualize",
    json={
        "file_path": "/path/to/data.csv",
        "chart_type": "scatter",
        "options": {
            "x_column": "age",
            "y_column": "expression"
        }
    }
)

result = response.json()
if result['success']:
    # 图像为 base64 编码的 PNG
    image_data = result['image']  # data:image/png;base64,iVBORw0KG...
```

## API 端点

### POST /api/quick-visualizer/analyze

分析文件并返回数据类型和推荐的可视化方式。

**请求体**:

```json
{
  "file_path": "/absolute/path/to/file.csv"
}
```

**响应**:

```json
{
  "success": true,
  "data_type": "expression_matrix",
  "n_rows": 1000,
  "n_columns": 50,
  "limited": false,
  "column_info": [
    {
      "name": "gene_id",
      "dtype": "object",
      "n_unique": 1000,
      "n_missing": 0,
      "is_numeric": false,
      "is_categorical": true
    },
    {
      "name": "sample1",
      "dtype": "float64",
      "n_unique": 950,
      "n_missing": 5,
      "is_numeric": true,
      "is_categorical": false,
      "min": 0.1,
      "max": 15.7,
      "mean": 5.2,
      "median": 4.8
    }
  ],
  "suggestions": [
    {
      "chart_type": "heatmap",
      "title": "Heatmap",
      "description": "Visualize numeric data as heatmap",
      "options": {},
      "priority": 1
    }
  ]
}
```

### POST /api/quick-visualizer/visualize

生成指定类型的可视化图表。

**请求体**:

```json
{
  "file_path": "/absolute/path/to/file.csv",
  "chart_type": "scatter",
  "options": {
    "x_column": "age",
    "y_column": "expression"
  }
}
```

**响应**:

```json
{
  "success": true,
  "image": "data:image/png;base64,iVBORw0KG...",
  "chart_type": "scatter"
}
```

### GET /api/quick-visualizer/chart-types

列出所有可用的图表类型及其配置要求。

**响应**:

```json
[
  {
    "type": "histogram",
    "name": "Histogram",
    "description": "Distribution of a numeric variable",
    "required_options": [],
    "optional_options": ["column", "bins"]
  }
]
```

### GET /api/quick-visualizer/data-types

列出所有可检测的数据类型。

## 技术实现

### 后端架构

**文件**: `backend/app/services/quick_visualizer.py`

- **QuickVisualizer**: 主服务类
- **DataType**: 数据类型枚举
- **ChartType**: 图表类型枚举
- 基于 pandas 进行数据处理
- 使用 matplotlib + seaborn 生成图表
- 图表返回为 base64 编码的 PNG 图像

### 智能检测算法

```python
def _detect_data_type(df):
    # 1. 检查表达矩阵模式（>10个数值列 且 >80%为数值）
    if len(numeric_cols) > 10 and ratio > 0.8:
        return EXPRESSION_MATRIX

    # 2. 检查时间序列（日期列或时间关键字）
    if has_date_columns or has_time_keywords:
        return TIME_SERIES

    # 3. 检查基因组区间（chr/start/end关键字）
    if genomic_keywords_count >= 2:
        return GENOMIC_INTERVALS

    # 4. 检查数值矩阵（>70%数值列）
    if numeric_ratio > 0.7:
        return NUMERIC_MATRIX

    # 5. 检查分类数据（>50%分类列）
    if categorical_ratio > 0.5:
        return CATEGORICAL

    # 6. 混合类型
    if has_both_numeric_and_categorical:
        return MIXED
```

### 前端组件

**文件**: `frontend/src/components/QuickVisualizerModal.tsx`

- React + Material-UI 模态对话框
- 自动分析和显示建议
- 交互式图表配置
- 实时图像预览
- 与 Data Browser 集成

## 性能考虑

1. **数据限制**: 自动限制为前 10,000 行以提高性能
2. **列限制**: 热图和相关性矩阵限制列数（20-50 列）
3. **图像格式**: PNG 格式，DPI=100，适合快速预览
4. **缓存**: 建议在生产环境中添加图像缓存

## 支持的文件格式

- CSV (`.csv`)
- TSV (`.tsv`)
- TXT (`.txt`, tab-separated)
- Excel (`.xlsx`)

## 限制和注意事项

1. **文件大小**: 建议 < 100MB，大文件会被截断到前 10,000 行
2. **列数**: 热图和相关性矩阵建议 < 50 列
3. **数据质量**: 缺失值会被自动处理（dropna）
4. **图表定制**: 当前仅支持基本配置，高级定制请使用专业可视化工具

## 与其他功能的集成

### Data Browser

- 在文件列表中显示 📊 按钮（仅限 CSV/TSV/XLSX）
- 一键打开可视化模态框

### Pipeline Builder

- 未来可将快速可视化结果保存为管道节点
- 可作为 QC 步骤使用

### Visualization Workspace

- 未来可将快速可视化添加到仪表板
- 支持多图表组合

## 示例用例

### 用例 1: 探索基因表达数据

```
文件: gene_expression.csv (1000 genes × 20 samples)

自动检测: Expression Matrix
推荐图表:
1. Heatmap - 查看整体表达模式
2. Correlation Matrix - 查看样本相关性
3. Box Plot - 查看样本分布
```

### 用例 2: 质量控制检查

```
文件: qc_metrics.csv (多个QC指标)

自动检测: Numeric Matrix
推荐图表:
1. Box Plot - 比较各指标分布
2. Histogram - 查看单个指标分布
3. Scatter - 检查指标间关系
```

### 用例 3: 时间序列分析

```
文件: timecourse_data.csv (时间点 × 测量值)

自动检测: Time Series
推荐图表:
1. Line Chart - 显示趋势
2. Scatter - 查看时间点关系
```

## 未来改进

- [ ] 添加更多图表类型（PCA, UMAP, 聚类树等）
- [ ] 支持图表交互（缩放、选择、过滤）
- [ ] 保存可视化配置为模板
- [ ] 导出高分辨率图像（PDF, SVG）
- [ ] 批量可视化生成
- [ ] 与 R/Python 专业可视化包集成
- [ ] 添加统计检验和注释
- [ ] 支持 3D 可视化

## 故障排查

### 问题 1: "Visualization libraries not available"

**原因**: matplotlib/seaborn 未安装

**解决**:

```bash
cd backend
pip install matplotlib seaborn
```

### 问题 2: 图表生成失败

**可能原因**:

- 数据格式不兼容
- 列名包含特殊字符
- 缺失值过多

**解决**: 检查错误消息，确保数据格式正确

### 问题 3: 生成速度慢

**优化建议**:

- 减少数据行数（使用采样）
- 限制可视化的列数
- 使用更简单的图表类型（如 histogram 代替 heatmap）

## 参考资料

- [Matplotlib 文档](https://matplotlib.org/)
- [Seaborn 文档](https://seaborn.pydata.org/)
- [Pandas 可视化](https://pandas.pydata.org/docs/user_guide/visualization.html)
