# Omicsomics - 统一组学分析平台

一个面向研究与临床的 Web 平台，支持常见组学数据（基因组学、转录组学、单细胞、表观、蛋白质组、代谢组、宏基因组等）的统一接收、处理和分析。

## 特性

### 已实现的核心功能 ✅

1. **用户认证与授权**

   - JWT token 认证
   - 用户注册和登录
   - 基于角色的访问控制

2. **项目管理**

   - 创建、查看、更新、删除项目
   - 项目级别的权限控制
   - 项目元数据管理

3. **样本管理**

   - 样本的 CRUD 操作
   - 灵活的 JSON 元数据支持
   - 样本与项目的关联

4. **文件存储**

   - 基于 MinIO 的对象存储
   - S3 兼容 API
   - 文件上传和下载
   - 预签名 URL 支持

5. **工作流执行**

   - Nextflow 流水线集成
   - FastQC 质量控制
   - 异步任务执行
   - 工作流状态跟踪
   - 日志记录

6. **质量控制（QC）**
   - FastQC 支持
   - 批量 QC 分析
   - QC 结果存储和查询

### 🧬 组学分析模块 ✅

7. **基因组学分析 (WGS/WES)**

   - 预处理: FastQC, fastp/Trimmomatic
   - 对齐: BWA-MEM, Bowtie2, Minimap2
   - 变异检测: GATK4 HaplotypeCaller, FreeBayes, DeepVariant
   - 变异注释: VEP, SnpEff, ANNOVAR
   - API 端点: 6 个 (`/qc`, `/trim`, `/align`, `/variant-calling`, `/annotate-variants`, `/complete-pipeline`)

8. **转录组学分析 (bulk RNA-seq)**

   - 对齐/定量: STAR, HISAT2, Salmon, Kallisto
   - count 矩阵生成: featureCounts
   - 差异表达: DESeq2, edgeR, limma-voom
   - 富集分析: GSEA (规划中)
   - API 端点: 3 个 (`/quantify`, `/count-matrix`, `/differential-expression`)

9. **单细胞分析 (scRNA-seq)**

   - 预处理: Cell Ranger count
   - 质控与标准化: Scanpy pipeline (EmptyDrops, SCTransform)
   - 降维与聚类: PCA, UMAP, Leiden clustering
   - 批次校正: Seurat integration, Harmony
   - 细胞注释: Marker-based annotation
   - API 端点: 4 个 (`/cellranger`, `/preprocess`, `/integrate`, `/annotate`)

10. **表观组学分析 (ChIP-seq/ATAC-seq)** ✨ NEW

    - 对齐: Bowtie2, BWA
    - Peak calling: MACS2/MACS3 (narrow/broad peaks)
    - 基序分析: HOMER
    - 信号可视化: BigWig 生成
    - API 端点: 5 个 (`/align`, `/peak-calling`, `/motif-analysis`, `/bigwig`, `/complete-pipeline`)

11. **蛋白质组学分析 (LC-MS/MS)** ✨ NEW

    - 原始文件转换: ThermoRawFileParser
    - 蛋白鉴定/定量: MaxQuant
    - 快速肽段搜索: MSFragger
    - Label-free 定量: LFQ
    - API 端点: 5 个 (`/convert-raw`, `/maxquant`, `/msfragger`, `/lfq-quantification`, `/complete-pipeline`)

12. **代谢组学分析 (LC-MS/GC-MS)** ✨ NEW

    - 特征检测: XCMS, MZmine
    - 谱图注释: GNPS, MS-DIAL
    - 定量归一化: median, quantile, PQN
    - API 端点: 4 个 (`/feature-detection`, `/spectral-annotation`, `/quantification`, `/complete-pipeline`)

13. **多组学整合 (Multi-omics Integration)** ✨ NEW

    - 无监督整合: MOFA2 (Multi-Omics Factor Analysis)
    - 有监督整合: DIABLO (生物标志物发现)
    - 样本匹配: 跨组学数据集匹配
    - 通路富集: 多层次功能分析
    - API 端点: 5 个 (`/mofa2`, `/diablo`, `/pathway-enrichment`, `/match-samples`, `/complete-pipeline`)

14. **可视化与交互**
    - 火山图 (Volcano plot) - 差异表达可视化
    - UMAP/PCA - 降维可视化
    - 热图 (Heatmap) - 基因表达矩阵
    - IGV.js - 基因组浏览器集成
    - QC 指标分布图
    - API 端点: 7 个 (Plotly 格式数据导出)

**总计**: 8 大组学模块，39 个 API 端点组 🎉

## 技术栈

- **后端**: FastAPI 0.121.0 + Python 3.11
- **数据库**: PostgreSQL 18.0 (AsyncIO 支持)
- **对象存储**: MinIO
- **ORM**: SQLAlchemy 2.0 (async)
- **迁移**: Alembic
- **认证**: JWT (python-jose)
- **密码哈希**: bcrypt

## 快速开始

详细部署指南请参见 [DEPLOYMENT.md](DEPLOYMENT.md)

### 1. 环境准备

```bash
micromamba create -n omicsomics-dev python=3.11
micromamba activate omicsomics-dev
micromamba install -n omicsomics-dev postgresql
cd backend && pip install -e .
```

### 2. 初始化

```bash
# 数据库
initdb -D local_db_data
pg_ctl -D local_db_data -l postgresql.log start
createdb omicsomics
cd backend && alembic upgrade head
```

### 3. 启动服务（3 个终端）

```bash
# 终端1: PostgreSQL
pg_ctl -D local_db_data start

# 终端2: MinIO
./scripts/start_minio.sh

# 终端3: FastAPI
cd backend
export SECRET_KEY="your-key" DATABASE_URL="postgresql+asyncpg://jeblqr@localhost/omicsomics"
uvicorn app.main:app --host 127.0.0.1 --port 8001 --reload
```

### 4. 访问

- **API 文档**: http://127.0.0.1:8001/docs
- **MinIO 控制台**: http://127.0.0.1:9003 (minioadmin/minioadmin123)

## API 示例

```bash
# 注册
curl -X POST "http://localhost:8001/api/v1/register" \
  -H "Content-Type: application/json" \
  -d '{"email": "user@example.com", "password": "pass123", "full_name": "User"}'

# 登录
TOKEN=$(curl -s -X POST "http://localhost:8001/api/v1/login/access-token" \
  -H "Content-Type: application/x-www-form-urlencoded" \
  -d "username=user@example.com&password=pass123" | jq -r '.access_token')

# 创建项目
curl -X POST "http://localhost:8001/api/v1/projects/" \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{"name": "RNA-Seq", "description": "Transcriptomics"}'

# 创建样本
curl -X POST "http://localhost:8001/api/v1/samples/" \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{"name": "Sample1", "project_id": 1, "metadata_": {"tissue": "liver"}}'

# 上传文件
curl -X POST "http://localhost:8001/api/v1/files/upload" \
  -H "Authorization: Bearer $TOKEN" \
  -F "file=@data.fastq" -F "sample_id=1" -F "file_type=fastq"

# 运行QC
curl -X POST "http://localhost:8001/api/v1/qc/fastqc" \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{"sample_id": 1, "file_ids": [1]}'
```

## 项目结构

```
Omicsomics/
├── backend/              # FastAPI 后端
│   ├── app/
│   │   ├── api/         # API 路由
│   │   ├── models/      # 数据库模型
│   │   ├── schemas/     # Pydantic schemas
│   │   ├── services/    # 业务逻辑
│   │   ├── storage/     # S3 客户端
│   │   └── workflows/   # 工作流执行器
│   └── alembic/         # 数据库迁移
├── scripts/             # 启动脚本
├── bin/                 # 二进制文件(MinIO)
├── local_db_data/       # PostgreSQL 数据
└── local_minio_data/    # MinIO 存储
```

## 数据模型

```
User → Projects → Samples → Files
                         └→ Workflows
```

## 开发

```bash
# 数据库迁移
alembic revision --autogenerate -m "description"
alembic upgrade head

# 运行测试
cd backend && pytest
```

## 安全

**生产环境务必**:

1. 更改 `SECRET_KEY`
2. 更改 MinIO 凭据
3. 启用 HTTPS
4. 使用强密码
5. 定期备份

## 路线图

### MVP ✅

- ✅ 认证、项目、样本、文件、工作流、QC

### v1.0 (计划)

- 单细胞、蛋白质组、可视化、MultiQC

### v2.0 (未来)

- 多组学整合、ML、实时协作、插件系统

## 许可

Apache v2.0 License

---

**详细文档**: [DEPLOYMENT.md](DEPLOYMENT.md) | **API**: http://127.0.0.1:8001/docs
