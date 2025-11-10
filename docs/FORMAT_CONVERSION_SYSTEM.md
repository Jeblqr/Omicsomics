# 数据格式转换系统 (Data Format Conversion System)

## 概述 (Overview)

数据格式转换系统提供**自动转换**和**手动转换**两种模式，支持跨运行时工具集成（R、Python、Binary）。

### 支持的格式 (Supported Formats)

| 格式   | 扩展名        | 运行时 | 用途               |
| ------ | ------------- | ------ | ------------------ |
| CSV    | .csv          | 通用   | 表格数据           |
| TSV    | .tsv, .txt    | 通用   | Tab 分隔表格       |
| Excel  | .xlsx, .xls   | 通用   | Excel 文件         |
| JSON   | .json         | 通用   | 结构化数据         |
| RDS    | .rds, .RDS    | R      | R 数据对象         |
| h5ad   | .h5ad         | Python | AnnData 单细胞数据 |
| pickle | .pkl, .pickle | Python | Python 对象序列化  |

### 转换能力矩阵 (Conversion Matrix)

```
CSV ↔ TSV ↔ Excel ↔ JSON
 ↕     ↕      ↕       ↕
RDS   -      -       -
 ↕     ↕      ↕       ↕
h5ad  -      -       -
 ↕     ↕      ↕       ↕
pickle -     -       -
```

所有格式都可以通过 CSV 作为中介进行转换。

---

## 1. 自动转换模式 (Auto Conversion Mode)

### 工作原理

在 Pipeline Editor 中构建管道时，系统会自动检测工具间的格式不匹配，并提供自动插入转换节点的选项。

### 使用流程

#### 步骤 1: 构建管道

在 Pipeline Editor 中连接不同运行时的工具：

```
[Seurat QC (R)] → [scanpy HVG (Python)] → [EnrichedHeatmap (R)]
    ↓ RDS           ↓ h5ad                ↓ RDS
```

#### 步骤 2: 格式检测

系统自动检测格式不匹配并显示警告面板：

![Auto Conversion Indicator]

- 显示需要转换的连接数量
- 展示每个转换的路径 (如: RDS → CSV → h5ad)
- 预估转换时间
- 标注多步转换警告 (⚠️)

#### 步骤 3: 应用转换

点击 **"自动插入转换节点"** 按钮，系统会：

1. 在格式不匹配的连接间插入转换节点
2. 转换节点以黄色 (amber) 高亮显示
3. 连接线改为动画效果，表示数据转换

#### 步骤 4: 执行管道

运行管道时，转换节点会自动执行格式转换。

### API 端点

#### 分析管道格式兼容性

```http
POST /api/pipelines/{pipeline_id}/analyze-formats
```

**响应示例:**

```json
{
  "conversions_needed": [
    {
      "source_node": "node_1",
      "target_node": "node_2",
      "source_format": "rds",
      "target_format": "h5ad",
      "conversion_path": ["rds", "csv", "h5ad"],
      "estimated_seconds": 25.0,
      "num_steps": 2
    }
  ],
  "num_conversions": 1,
  "warnings": [],
  "total_estimated_time": 25.0
}
```

#### 自动插入转换节点

```http
POST /api/pipelines/{pipeline_id}/insert-conversions
```

**响应:** 返回修改后的管道定义，包含新插入的转换节点。

---

## 2. 手动转换模式 (Manual Conversion Mode)

### 使用场景

- 在 Data Browser 中转换单个文件
- 批量转换多个文件
- 准备特定格式的数据用于工具

### 使用流程

#### 方式 A: Data Browser 右键菜单

1. 在 Data Browser 中找到要转换的文件
2. 右键点击文件，选择 **"转换格式 (Convert Format)"**
3. 在弹出的 Format Converter Modal 中：
   - 源文件路径自动填充
   - 源格式自动检测
   - 选择目标格式
   - 目标文件路径自动生成
   - 查看转换路径和预估时间
4. 点击 **"开始转换 (Convert)"**

#### 方式 B: 独立转换页面

1. 导航到 **工具 (Tools) → 格式转换器 (Format Converter)**
2. 手动输入源文件路径和目标文件路径
3. 选择源格式和目标格式
4. 查看转换预览
5. 执行转换

#### 方式 C: 批量转换

```bash
# 使用 API 进行批量转换
curl -X POST http://localhost:8000/api/formats/convert \
  -H "Content-Type: application/json" \
  -d '{
    "source_path": "/data/file1.rds",
    "target_path": "/data/file1.h5ad",
    "from_format": "rds",
    "to_format": "h5ad"
  }'
```

### API 端点

#### 获取支持的格式

```http
GET /api/formats/supported
```

**响应:**

```json
[
  {
    "format_id": "csv",
    "extensions": [".csv"],
    "mime_type": "text/csv",
    "can_convert_to": ["tsv", "excel", "json", "rds", "h5ad", "pickle"],
    "can_convert_from": ["tsv", "excel", "json", "rds", "h5ad", "pickle"]
  }
]
```

#### 检测文件格式

```http
POST /api/formats/detect
```

**请求:**

```json
{
  "file_path": "/data/sample.csv"
}
```

**响应:**

```json
{
  "file_path": "/data/sample.csv",
  "format": "csv",
  "extensions": [".csv"],
  "mime_type": "text/csv"
}
```

#### 获取转换路径

```http
POST /api/formats/conversion-path
```

**请求:**

```json
{
  "from_format": "rds",
  "to_format": "h5ad"
}
```

**响应:**

```json
{
  "from_format": "rds",
  "to_format": "h5ad",
  "conversion_path": ["rds", "csv", "h5ad"],
  "num_steps": 2
}
```

#### 预估转换时间

```http
POST /api/formats/estimate
```

**请求:**

```json
{
  "file_path": "/data/sample.rds",
  "from_format": "rds",
  "to_format": "h5ad"
}
```

**响应:**

```json
{
  "file_path": "/data/sample.rds",
  "file_size_bytes": 1073741824,
  "file_size_mb": 1024.0,
  "from_format": "rds",
  "to_format": "h5ad",
  "estimated_seconds": 45.0,
  "estimated_minutes": 0.75
}
```

#### 执行转换

```http
POST /api/formats/convert
```

**请求:**

```json
{
  "source_path": "/data/sample.rds",
  "target_path": "/data/sample.h5ad",
  "from_format": "rds",
  "to_format": "h5ad",
  "parameters": {
    "obs": { "cell_type": ["A", "B", "C"] },
    "var": { "gene_name": ["Gene1", "Gene2"] }
  }
}
```

**响应:**

```json
{
  "id": 123,
  "source_path": "/data/sample.rds",
  "target_path": "/data/sample.h5ad",
  "from_format": "rds",
  "to_format": "h5ad",
  "conversion_path": ["rds", "csv", "h5ad"],
  "status": "pending",
  "conversion_mode": "manual",
  "created_at": "2025-01-10T12:00:00Z"
}
```

#### 查询转换状态

```http
GET /api/formats/conversions/{conversion_id}
```

**响应:**

```json
{
  "id": 123,
  "status": "completed",
  "duration_seconds": 42.5,
  "target_size_bytes": 1048576000
}
```

#### 列出转换历史

```http
GET /api/formats/conversions?status=completed&limit=10
```

---

## 3. 转换性能 (Conversion Performance)

### 预估时间 (每 GB)

| 转换         | 时间 (s/GB) |
| ------------ | ----------- |
| CSV ↔ TSV    | 2           |
| CSV ↔ JSON   | 10          |
| CSV ↔ Excel  | 15          |
| CSV ↔ pickle | 8           |
| CSV ↔ RDS    | 20          |
| CSV ↔ h5ad   | 25          |

### 性能优化建议

1. **大文件转换**: 对于 >10GB 的文件，建议在后台任务中执行
2. **批量转换**: 使用并行处理多个文件
3. **中间格式缓存**: 系统会缓存中间转换结果，避免重复转换
4. **流式处理**: 对超大文件使用分块读写

---

## 4. 前端组件 (Frontend Components)

### AutoConversionIndicator

在 Pipeline Editor 中显示格式不匹配警告。

**使用:**

```tsx
import AutoConversionIndicator from "@/components/AutoConversionIndicator";

<AutoConversionIndicator
  conversions={conversions}
  warnings={warnings}
  totalEstimatedTime={totalTime}
  onApplyConversions={handleApply}
  onDismiss={handleDismiss}
/>;
```

### FormatConverterModal

手动格式转换弹窗。

**使用:**

```tsx
import FormatConverterModal from "@/components/FormatConverterModal";

<FormatConverterModal
  isOpen={isOpen}
  onClose={handleClose}
  initialFilePath="/data/sample.csv"
  onConversionComplete={(targetPath) => {
    console.log("Converted to:", targetPath);
  }}
/>;
```

---

## 5. 后端实现 (Backend Implementation)

### FormatConverter 类

核心转换引擎。

**主要方法:**

- `detect_format(file_path)`: 检测文件格式
- `get_conversion_path(from_format, to_format)`: 获取转换路径
- `estimate_conversion_time(file_path, from_format, to_format)`: 预估时间
- `convert(source_path, target_path, ...)`: 执行转换

### PipelineAutoConverter 类

管道自动转换分析器。

**主要方法:**

- `analyze_pipeline(pipeline)`: 分析管道格式兼容性
- `insert_conversion_nodes(pipeline)`: 插入转换节点

### 数据库模型

#### FormatConversion

记录所有转换操作。

**字段:**

- `source_path`, `source_format`
- `target_path`, `target_format`
- `conversion_path`: 转换路径 (JSON)
- `conversion_mode`: 'auto' 或 'manual'
- `status`: 'pending', 'running', 'completed', 'failed'
- `duration_seconds`, `error_message`

#### ConversionRule

存储转换规则和脚本。

**字段:**

- `from_format`, `to_format`
- `method`: 'pandas', 'r_script', 'python_script', 'binary'
- `script_path`, `command_template`
- `avg_time_per_gb`, `success_rate`

---

## 6. 转换脚本示例 (Conversion Scripts)

### CSV → RDS (R)

```r
data <- read.csv("/path/to/input.csv", stringsAsFactors = FALSE)
saveRDS(data, "/path/to/output.rds")
```

### RDS → CSV (R)

```r
data <- readRDS("/path/to/input.rds")
write.csv(data, "/path/to/output.csv", row.names = FALSE)
```

### CSV → h5ad (Python)

```python
import pandas as pd
import anndata as ad

df = pd.read_csv("/path/to/input.csv", index_col=0)
adata = ad.AnnData(X=df.values)
adata.obs_names = df.index.astype(str)
adata.var_names = df.columns.astype(str)
adata.write_h5ad("/path/to/output.h5ad")
```

### h5ad → CSV (Python)

```python
import anndata as ad
import pandas as pd

adata = ad.read_h5ad("/path/to/input.h5ad")
df = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)
df.to_csv("/path/to/output.csv")
```

---

## 7. 故障排除 (Troubleshooting)

### 常见问题

#### Q: 转换失败，显示 "R conversion failed"

**A:** 确保安装了 R 和所需的包。检查 Rscript 是否在系统 PATH 中。

#### Q: h5ad 转换错误

**A:** 确保安装了 `anndata` 包：

```bash
pip install anndata
```

#### Q: 转换时间过长

**A:** 检查文件大小。对于 >10GB 的文件，转换可能需要几分钟。考虑：

- 使用更高效的格式 (如 pickle 而非 Excel)
- 分批处理大文件
- 在服务器上执行转换而非本地

#### Q: 管道中自动转换节点显示错误

**A:** 检查：

1. 工具的输入/输出格式定义是否正确
2. 格式转换路径是否存在
3. 后端 API 是否正常运行

---

## 8. 最佳实践 (Best Practices)

### 性能优化

1. **选择合适的格式**

   - 小文件 (<100MB): 使用 CSV/JSON（易读）
   - 大文件 (>1GB): 使用 pickle/h5ad（高效）
   - R 分析: 使用 RDS
   - Python 分析: 使用 h5ad/pickle

2. **避免不必要的转换**

   - 尽量在同一运行时内完成连续操作
   - 规划好工具链，减少运行时切换

3. **缓存中间结果**
   - 系统会自动缓存转换结果
   - 重用已转换的文件，避免重复转换

### 数据完整性

1. **验证转换结果**

   - 检查文件大小是否合理
   - 对比源文件和目标文件的行列数
   - 抽查数据值是否一致

2. **备份重要数据**
   - 转换前备份原始文件
   - 使用版本控制跟踪数据变化

---

## 9. 开发者指南 (Developer Guide)

### 添加新格式支持

#### 步骤 1: 更新 FORMATS 字典

```python
# backend/app/converters/format_converter.py
FORMATS = {
    'new_format': {
        'extensions': ['.nfmt'],
        'mime_type': 'application/x-new-format'
    }
}
```

#### 步骤 2: 添加转换时间预估

```python
CONVERSION_TIME_ESTIMATES = {
    ('csv', 'new_format'): 30,
    ('new_format', 'csv'): 30,
}
```

#### 步骤 3: 实现转换方法

```python
def _convert_csv_to_new_format(self, source_path: str, target_path: str):
    # 实现转换逻辑
    pass

def _convert_new_format_to_csv(self, source_path: str, target_path: str):
    # 实现转换逻辑
    pass
```

#### 步骤 4: 更新 \_convert_direct 方法

```python
def _convert_direct(self, source_path: str, target_path: str, ...):
    # ...
    elif from_format == 'csv' and to_format == 'new_format':
        self._convert_csv_to_new_format(source_path, target_path)
    elif from_format == 'new_format' and to_format == 'csv':
        self._convert_new_format_to_csv(source_path, target_path)
```

### 测试

```python
# tests/test_format_converter.py
def test_new_format_conversion():
    converter = FormatConverter()

    # 测试格式检测
    format_id = converter.detect_format('/path/to/file.nfmt')
    assert format_id == 'new_format'

    # 测试转换路径
    path = converter.get_conversion_path('csv', 'new_format')
    assert path == ['csv', 'new_format']

    # 测试转换
    result = converter.convert(
        source_path='/path/to/input.csv',
        target_path='/path/to/output.nfmt'
    )
    assert result['success'] == True
```

---

## 10. 未来规划 (Future Plans)

### 短期 (Q1 2025)

- [ ] 支持更多生物信息学格式 (VCF, BAM, FASTQ)
- [ ] 流式转换支持大文件 (>100GB)
- [ ] 并行批量转换
- [ ] 转换质量检查和报告

### 中期 (Q2-Q3 2025)

- [ ] 云端转换 (AWS S3, Google Cloud Storage)
- [ ] 增量转换 (只转换变化部分)
- [ ] 智能格式推荐 (基于后续工具需求)
- [ ] 转换性能优化 (GPU 加速)

### 长期 (Q4 2025+)

- [ ] 自定义转换脚本上传
- [ ] 转换市场 (社区贡献转换器)
- [ ] 机器学习优化转换参数
- [ ] 实时转换预览

---

## 相关文档

- [Tool Variant Selection System](./TOOL_VARIANT_SYSTEM_SUMMARY.md)
- [Cross-Runtime Integration](./CROSS_RUNTIME_INTEGRATION.md)
- [Pipeline Builder Guide](./PIPELINE_BUILDER.md)
- [Data Management](./DATA_MANAGEMENT_GUIDE.md)

---

**最后更新:** 2025-01-10
**版本:** 1.0.0
