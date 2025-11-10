# Dataset Manager

## 概述

数据集管理器（Dataset Manager）是一个用于组织和管理相关数据文件的系统，提供版本控制、文件来源追踪（lineage）和元数据管理功能。它允许用户将多个文件组织成逻辑数据集，并跟踪数据的处理历史。

## 核心概念

### 1. Dataset（数据集）

数据集是相关文件的逻辑集合，包含：

- **基本信息**: 名称、描述、数据类型
- **文件**: 包含的所有文件及其元数据
- **版本**: 版本号和历史记录
- **标签**: 用于分类和搜索的标签
- **状态**: active（活动）、archived（归档）、deleted（删除）

### 2. Version Control（版本控制）

- 每个数据集有版本号（从 1 开始）
- 创建新版本时复制文件列表和标签
- 版本之间通过 parent-child 关系连接
- 可以查看所有历史版本

### 3. File Lineage（文件来源）

追踪数据处理历史：

- **Source File**: 输入文件
- **Output File**: 输出文件
- **Operation**: 执行的操作（merge、conversion、pipeline 等）
- **Tool**: 使用的工具名称和版本
- **Parameters**: 操作参数

### 4. Tags（标签）

- 用于组织和分类数据集
- 每个标签有名称和颜色
- 支持多标签搜索

### 5. File Integrity（文件完整性）

- **MD5 Hash**: 快速校验
- **SHA256 Hash**: 安全校验
- 自动计算并存储文件哈希值

## 数据库架构

### Dataset Table

```sql
CREATE TABLE datasets (
    id VARCHAR PRIMARY KEY,
    project_id VARCHAR NOT NULL,
    name VARCHAR NOT NULL,
    description TEXT,
    version INTEGER DEFAULT 1,
    data_type VARCHAR,
    file_format VARCHAR,
    metadata JSON,
    status VARCHAR DEFAULT 'active',
    is_public BOOLEAN DEFAULT FALSE,
    created_at TIMESTAMP,
    updated_at TIMESTAMP,
    created_by VARCHAR,
    parent_id VARCHAR  -- For versioning
);
```

### DatasetFileEntry Table

```sql
CREATE TABLE dataset_file_entries (
    id VARCHAR PRIMARY KEY,
    file_path VARCHAR NOT NULL,
    file_name VARCHAR NOT NULL,
    file_type VARCHAR,
    file_size INTEGER,
    role VARCHAR,  -- 'primary', 'metadata', 'supplementary'
    description TEXT,
    md5_hash VARCHAR,
    sha256_hash VARCHAR,
    added_at TIMESTAMP
);
```

### DatasetLineage Table

```sql
CREATE TABLE dataset_lineage (
    id VARCHAR PRIMARY KEY,
    dataset_id VARCHAR NOT NULL,
    operation_type VARCHAR,
    operation_id VARCHAR,
    source_file_id VARCHAR,
    output_file_id VARCHAR,
    operation_params JSON,
    tool_name VARCHAR,
    tool_version VARCHAR,
    executed_at TIMESTAMP
);
```

### Tag Table

```sql
CREATE TABLE tags (
    id VARCHAR PRIMARY KEY,
    name VARCHAR UNIQUE NOT NULL,
    color VARCHAR DEFAULT '#3b82f6',
    description TEXT,
    created_at TIMESTAMP
);
```

## API 端点

### 数据集管理

#### POST /api/datasets/

创建新数据集

**请求体**:

```json
{
  "project_id": "project-123",
  "name": "RNA-seq Batch 1",
  "description": "First batch of RNA-seq samples",
  "data_type": "transcriptomics",
  "file_format": "FASTQ",
  "metadata": {
    "organism": "Homo sapiens",
    "tissue": "liver"
  },
  "tags": ["rna-seq", "batch1"]
}
```

**响应**:

```json
{
  "id": "dataset-456",
  "name": "RNA-seq Batch 1",
  "version": 1,
  "status": "active",
  "file_count": 0,
  "tag_names": ["rna-seq", "batch1"],
  "created_at": "2025-01-10T10:00:00Z"
}
```

#### GET /api/datasets/{dataset_id}

获取数据集详情

#### GET /api/datasets/

列出数据集

**查询参数**:

- `project_id`: 过滤项目
- `data_type`: 过滤数据类型
- `status`: 过滤状态（active/archived/deleted）
- `tags`: 逗号分隔的标签列表
- `search`: 搜索关键词
- `skip`: 分页偏移
- `limit`: 每页数量

#### PUT /api/datasets/{dataset_id}

更新数据集

#### DELETE /api/datasets/{dataset_id}

删除数据集（软删除）

#### POST /api/datasets/{dataset_id}/archive

归档数据集

### 版本控制

#### POST /api/datasets/{dataset_id}/versions

创建新版本

**请求体**:

```json
{
  "description": "Updated with reprocessed samples"
}
```

#### GET /api/datasets/{dataset_id}/versions

获取所有版本

### 文件管理

#### POST /api/datasets/{dataset_id}/files

添加文件到数据集

**请求体**:

```json
{
  "file_path": "/path/to/sample1.fastq",
  "role": "primary",
  "description": "Sample 1 raw reads",
  "compute_hash": true
}
```

**响应**:

```json
{
  "id": "file-789",
  "file_name": "sample1.fastq",
  "file_path": "/path/to/sample1.fastq",
  "file_type": ".fastq",
  "file_size": 1234567890,
  "role": "primary",
  "md5_hash": "5d41402abc4b2a76b9719d911017c592",
  "sha256_hash": "2c26b46b68ffc68ff99b453c1d30413413422d706483bfa0f98a5e886266e7ae",
  "added_at": "2025-01-10T10:05:00Z"
}
```

#### DELETE /api/datasets/{dataset_id}/files/{file_id}

从数据集移除文件

#### GET /api/datasets/{dataset_id}/files

列出数据集中的所有文件

### 来源追踪

#### POST /api/datasets/{dataset_id}/lineage

添加来源记录

**请求体**:

```json
{
  "operation_type": "conversion",
  "source_file_id": "file-789",
  "output_file_id": "file-790",
  "operation_id": "job-123",
  "operation_params": {
    "from_format": "FASTQ",
    "to_format": "BAM",
    "quality_threshold": 20
  },
  "tool_name": "BWA-MEM",
  "tool_version": "0.7.17"
}
```

#### GET /api/datasets/{dataset_id}/lineage

获取数据集的来源历史

#### GET /api/datasets/files/{file_id}/lineage

追踪文件的完整来源树

**响应**:

```json
{
  "file": {
    "id": "file-790",
    "name": "sample1.bam",
    "path": "/path/to/sample1.bam"
  },
  "created_by": [
    {
      "operation": "conversion",
      "tool": "BWA-MEM",
      "executed_at": "2025-01-10T11:00:00Z",
      "source_files": [
        {
          "id": "file-789",
          "name": "sample1.fastq"
        }
      ]
    }
  ],
  "used_in": [
    {
      "operation": "variant_calling",
      "tool": "GATK",
      "executed_at": "2025-01-10T12:00:00Z",
      "output_files": [
        {
          "id": "file-791",
          "name": "sample1.vcf"
        }
      ]
    }
  ]
}
```

### 标签管理

#### POST /api/datasets/{dataset_id}/tags?tag_name=experiment1

添加标签到数据集

#### DELETE /api/datasets/{dataset_id}/tags/{tag_id}

从数据集移除标签

#### GET /api/datasets/tags/list

列出所有标签

#### POST /api/datasets/tags

创建新标签

**请求体**:

```json
{
  "name": "high-priority",
  "color": "#ef4444",
  "description": "High priority datasets"
}
```

### 统计

#### GET /api/datasets/stats/summary?project_id=project-123

获取数据集统计

**响应**:

```json
{
  "total_datasets": 50,
  "active": 45,
  "archived": 5,
  "by_data_type": {
    "genomics": 20,
    "proteomics": 15,
    "transcriptomics": 10,
    "metabolomics": 5
  }
}
```

## 使用场景

### 场景 1: 创建和管理 RNA-seq 数据集

```python
import requests

API_BASE = "http://localhost:8000/api/datasets"

# 1. 创建数据集
dataset = requests.post(f"{API_BASE}/", json={
    "project_id": "proj-001",
    "name": "RNA-seq Liver Samples",
    "description": "RNA-seq from mouse liver tissue",
    "data_type": "transcriptomics",
    "tags": ["rna-seq", "mouse", "liver"]
}).json()

dataset_id = dataset["id"]

# 2. 添加FASTQ文件
files = [
    "/data/sample1_R1.fastq.gz",
    "/data/sample1_R2.fastq.gz",
    "/data/sample2_R1.fastq.gz",
    "/data/sample2_R2.fastq.gz"
]

for file_path in files:
    requests.post(f"{API_BASE}/{dataset_id}/files", json={
        "file_path": file_path,
        "role": "primary",
        "compute_hash": True
    })

# 3. 记录处理步骤
requests.post(f"{API_BASE}/{dataset_id}/lineage", json={
    "operation_type": "quality_control",
    "tool_name": "FastQC",
    "tool_version": "0.11.9",
    "operation_params": {
        "threads": 4
    }
})

# 4. 查看数据集
result = requests.get(f"{API_BASE}/{dataset_id}").json()
print(f"Dataset: {result['name']}")
print(f"Files: {result['file_count']}")
print(f"Tags: {result['tag_names']}")
```

### 场景 2: 版本控制和更新

```python
# 处理完成后创建新版本
new_version = requests.post(
    f"{API_BASE}/{dataset_id}/versions",
    json={"description": "Added quality-filtered reads"}
).json()

new_version_id = new_version["id"]

# 添加过滤后的文件
requests.post(f"{API_BASE}/{new_version_id}/files", json={
    "file_path": "/data/sample1_filtered.fastq.gz",
    "role": "primary",
    "description": "Quality filtered reads"
})

# 记录过滤操作
requests.post(f"{API_BASE}/{new_version_id}/lineage", json={
    "operation_type": "filtering",
    "tool_name": "Trimmomatic",
    "tool_version": "0.39",
    "operation_params": {
        "quality_threshold": 20,
        "min_length": 36
    }
})

# 查看所有版本
versions = requests.get(f"{API_BASE}/{dataset_id}/versions").json()
for v in versions:
    print(f"Version {v['version']}: {v['description']}")
```

### 场景 3: 文件来源追踪

```python
# 获取某个文件的完整处理历史
file_id = "file-123"
lineage = requests.get(f"{API_BASE}/files/{file_id}/lineage").json()

print(f"File: {lineage['file']['name']}")
print("\nCreated by:")
for op in lineage['created_by']:
    print(f"  - {op['operation']} using {op['tool']}")

print("\nUsed in:")
for op in lineage['used_in']:
    print(f"  - {op['operation']} using {op['tool']}")
    print(f"    Outputs: {[f['name'] for f in op['output_files']]}")
```

### 场景 4: 搜索和过滤

```python
# 搜索包含"liver"的数据集
results = requests.get(f"{API_BASE}/", params={
    "search": "liver",
    "data_type": "transcriptomics",
    "status": "active",
    "limit": 10
}).json()

for ds in results:
    print(f"{ds['name']} (v{ds['version']}) - {ds['file_count']} files")

# 按标签搜索
tagged_datasets = requests.get(f"{API_BASE}/", params={
    "tags": "rna-seq,high-priority"
}).json()
```

## 前端界面

### Dataset Manager 页面

**URL**: `/dataset-manager`

**功能**:

1. **数据集卡片视图**

   - 显示数据集名称、版本、描述
   - 文件数量和标签
   - 操作按钮（查看、归档、删除）

2. **过滤和搜索**

   - 文本搜索
   - 状态过滤（Active/Archived）
   - 数据类型过滤

3. **创建数据集对话框**

   - 名称、描述、数据类型
   - 初始标签

4. **数据集详情对话框**
   - **Files 标签页**: 文件列表，显示名称、类型、大小、角色
   - **Lineage 标签页**: 处理历史，显示操作、工具、时间
   - **Versions 标签页**: 版本列表，显示版本号、文件数、状态

## 与其他功能的集成

### 1. Data Browser 集成

```javascript
// 从Data Browser创建数据集
const selectedFiles = getSelectedFiles();
const dataset = await createDataset({
  name: "New Dataset",
  project_id: currentProject.id,
});

// 添加选中的文件
for (const file of selectedFiles) {
  await addFileToDataset(dataset.id, file.path);
}
```

### 2. Data Merger 集成

```javascript
// 合并操作后自动记录来源
const mergeResult = await mergeFiles({
  files: [file1, file2],
  mode: "horizontal",
});

// 记录到数据集
await addLineage(dataset.id, {
  operation_type: "merge",
  source_file_id: [file1.id, file2.id],
  output_file_id: mergeResult.file_id,
  tool_name: "DataMerger",
  operation_params: {
    mode: "horizontal",
    join_type: "inner",
  },
});
```

### 3. Pipeline Builder 集成

```javascript
// 管道运行后创建输出数据集
const run = await executePipeline(pipeline);

const outputDataset = await createDataset({
  name: `${pipeline.name} - Run ${run.id}`,
  description: `Output from pipeline run ${run.id}`,
  data_type: pipeline.output_type,
});

// 添加所有输出文件
for (const output of run.outputs) {
  await addFileToDataset(outputDataset.id, output.file_path, {
    role: output.type,
    description: output.description,
  });

  // 记录来源
  await addLineage(outputDataset.id, {
    operation_type: "pipeline",
    operation_id: run.id,
    tool_name: pipeline.name,
    operation_params: run.parameters,
  });
}
```

### 4. Quick Visualizer 集成

```javascript
// 从数据集文件直接可视化
const files = await getDatasetFiles(dataset.id);
const primaryFile = files.find((f) => f.role === "primary");

// 打开可视化
openQuickVisualizer(primaryFile.file_path, primaryFile.file_name);
```

## 最佳实践

### 1. 数据集组织

- **按项目分组**: 每个研究项目一个或多个数据集
- **使用描述性名称**: "RNA-seq_LiverSamples_Batch1" 而不是 "Dataset1"
- **添加完整描述**: 包含实验条件、样本信息等
- **使用标签**: 便于搜索和分类

### 2. 版本管理

- **何时创建新版本**:
  - 数据重新处理
  - 添加新样本
  - 重要参数改变
- **版本描述**: 说明此版本的变化

### 3. 文件角色

- **primary**: 主要数据文件
- **metadata**: 元数据文件（样本信息、实验设计）
- **supplementary**: 辅助文件（QC 报告、日志）

### 4. 来源追踪

- **记录所有操作**: 即使是简单的格式转换
- **包含参数**: 确保可重现性
- **记录工具版本**: 便于追溯

### 5. 文件完整性

- **始终计算哈希**: 用于验证文件完整性
- **定期验证**: 检查文件是否损坏

## 故障排查

### 问题 1: 文件哈希计算失败

**原因**:

- 文件权限不足
- 文件不存在
- 磁盘 I/O 错误

**解决**:

```python
# 跳过哈希计算
add_file(dataset_id, file_path, compute_hash=False)

# 或手动设置哈希
import hashlib

def compute_md5(file_path):
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()
```

### 问题 2: 数据库关系错误

**症状**: 外键约束错误

**解决**: 确保正确的创建顺序

```python
# 正确顺序
1. 创建项目
2. 创建数据集
3. 添加文件
4. 添加来源记录
```

### 问题 3: 大量文件的性能问题

**优化**:

```python
# 批量添加文件
from sqlalchemy.orm import Session

def batch_add_files(db: Session, dataset_id: str, files: List[str]):
    file_entries = []
    for file_path in files:
        entry = DatasetFileEntry(
            file_path=file_path,
            file_name=Path(file_path).name,
            # ... other fields
        )
        file_entries.append(entry)

    # 批量插入
    db.bulk_save_objects(file_entries)
    db.commit()
```

## 未来改进

- [ ] 自动数据集快照
- [ ] 数据集导出（ZIP with metadata）
- [ ] 数据集模板系统
- [ ] 与云存储集成（S3, Azure Blob）
- [ ] 数据集访问权限管理
- [ ] 数据集统计和摘要视图
- [ ] 自动来源追踪（从 Pipeline 自动记录）
- [ ] 数据集比较工具
- [ ] 数据集合并和拆分
- [ ] 高级搜索（全文搜索、元数据查询）

## 参考资料

- [SQLAlchemy Documentation](https://docs.sqlalchemy.org/)
- [Data Provenance](https://en.wikipedia.org/wiki/Data_lineage)
- [Semantic Versioning](https://semver.org/)
