# Omicsomics - 统一组学分析平台

> 一个面向研究与临床的 Web 平台，支持常见组学数据（基因组学、转录组学、单细胞、表观、蛋白质组、代谢组、宏基因组等）的统一接收、处理和分析。

[![CI/CD Pipeline](https://github.com/Jeblqr/Omicsomics/actions/workflows/ci.yml/badge.svg)](https://github.com/Jeblqr/Omicsomics/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/Python-3.11+-blue.svg)](https://www.python.org/)
[![FastAPI](https://img.shields.io/badge/FastAPI-0.115+-green.svg)](https://fastapi.tiangolo.com/)
[![React](https://img.shields.io/badge/React-18-61dafb.svg)](https://reactjs.org/)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![codecov](https://codecov.io/gh/Jeblqr/Omicsomics/branch/main/graph/badge.svg)](https://codecov.io/gh/Jeblqr/Omicsomics)

## 特性

### 已实现的核心功能 ✅

1. **用户认证与授权**

---

- JWT token 认证

## 📋 目录 - 用户注册和登录

- 基于角色的访问控制

- [特性](#特性)

- [技术栈](#技术栈)2. **项目管理**

- [快速开始](#快速开始)

- [API 文档](#api-文档) - 创建、查看、更新、删除项目

- [前端使用](#前端使用) - 项目级别的权限控制

- [部署指南](#部署指南) - 项目元数据管理

- [测试](#测试)

- [项目结构](#项目结构)3. **样本管理**

- [常见问题](#常见问题)

  - 样本的 CRUD 操作

--- - 灵活的 JSON 元数据支持

- 样本与项目的关联

## ✨ 特性

4. **文件存储**

### 后端模块 (9 个核心模块，44 个 API 端点)

- 基于 MinIO 的对象存储

| 模块 | 功能 | API 端点 | 状态 | - S3 兼容 API

|------|------|---------|------| - 文件上传和下载

| **Genomics** | WGS/WES 分析：QC、比对、变异检测 | 6 | ✅ | - 预签名 URL 支持

| **Transcriptomics** | RNA-seq：定量、差异表达、富集 | 3 | ✅ |

| **Single-cell** | scRNA-seq：QC、聚类、轨迹分析 | 4 | ✅ |5. **工作流执行**

| **Epigenomics** | ChIP-seq、ATAC-seq、DNA 甲基化 | 5 | ✅ |

| **Proteomics** | 肽段鉴定、蛋白定量、差异分析 | 5 | ✅ | - Nextflow 流水线集成

| **Metabolomics** | 特征检测、代谢物注释、定量 | 4 | ✅ | - FastQC 质量控制

| **Multi-omics** | MOFA2、DIABLO 多组学整合 | 5 | ✅ | - 异步任务执行

| **GWAS** | PLINK QC、关联分析、MTAG | 5 | ✅ | - 工作流状态跟踪

| **Visualizations** | 火山图、热图、PCA、UMAP | 7 | ✅ | - 日志记录

### 前端界面 (8 个分析页面)6. **质量控制（QC）**

- FastQC 支持

- ✅ **GWAS 分析** - 质控、关联分析、跨性状 MTAG - 批量 QC 分析

- ✅ **多组学整合** - MOFA2 无监督、DIABLO 监督学习 - QC 结果存储和查询

- ✅ **代谢组学** - XCMS 检测、GNPS 注释、定量归一化

- ✅ **蛋白质组学** - Mascot 搜索、iTRAQ 定量、差异分析### 🧬 组学分析模块 ✅

- ✅ **表观基因组学** - ChIP-seq、ATAC-seq、甲基化分析

- ✅ **单细胞分析** - QC 过滤、Louvain 聚类、拟时序 7. **基因组学分析 (WGS/WES)**

- ✅ **基因组学** - FastQC、BWA 比对、GATK 变异检测

- ✅ **转录组学** - Salmon 定量、DESeq2 差异表达、GSEA - 预处理: FastQC, fastp/Trimmomatic

  - 对齐: BWA-MEM, Bowtie2, Minimap2

--- - 变异检测: GATK4 HaplotypeCaller, FreeBayes, DeepVariant

- 变异注释: VEP, SnpEff, ANNOVAR

## 🛠 技术栈 - API 端点: 6 个 (`/qc`, `/trim`, `/align`, `/variant-calling`, `/annotate-variants`, `/complete-pipeline`)

### 后端 8. **转录组学分析 (bulk RNA-seq)**

- **FastAPI 0.121.0** - 高性能异步 Web 框架

- **SQLAlchemy 2.0** - 异步 ORM - 对齐/定量: STAR, HISAT2, Salmon, Kallisto

- **PostgreSQL 15** - 关系型数据库 - count 矩阵生成: featureCounts

- **asyncpg** - 异步数据库驱动 - 差异表达: DESeq2, edgeR, limma-voom

- **MinIO** - S3 兼容对象存储 - 富集分析: GSEA (规划中)

- **JWT + Passlib** - 身份认证 - API 端点: 3 个 (`/quantify`, `/count-matrix`, `/differential-expression`)

- **Alembic** - 数据库迁移

- **pytest + httpx** - 测试框架 9. **单细胞分析 (scRNA-seq)**

### 前端 - 预处理: Cell Ranger count

- **React 18** + **TypeScript 5** - UI 框架 - 质控与标准化: Scanpy pipeline (EmptyDrops, SCTransform)

- **Vite** - 构建工具 - 降维与聚类: PCA, UMAP, Leiden clustering

- **shadcn/ui** - UI 组件库 - 批次校正: Seurat integration, Harmony

- **Tailwind CSS** - 样式框架 - 细胞注释: Marker-based annotation

- **Lucide React** - 图标库 - API 端点: 4 个 (`/cellranger`, `/preprocess`, `/integrate`, `/annotate`)

- **React Router** - 路由管理

10. **表观组学分析 (ChIP-seq/ATAC-seq)** ✨ NEW

### 生物信息学工具

- PLINK 1.9+, BWA/STAR, GATK 4, Salmon/Kallisto - 对齐: Bowtie2, BWA

- Scanpy, MACS2, XCMS, MOFA2 - Peak calling: MACS2/MACS3 (narrow/broad peaks)

  - 基序分析: HOMER

--- - 信号可视化: BigWig 生成

    - API 端点: 5 个 (`/align`, `/peak-calling`, `/motif-analysis`, `/bigwig`, `/complete-pipeline`)

## 🚀 快速开始

11. **蛋白质组学分析 (LC-MS/MS)** ✨ NEW

### 使用 Docker (推荐)

    - 原始文件转换: ThermoRawFileParser

````bash - 蛋白鉴定/定量: MaxQuant

# 1. 克隆仓库    - 快速肽段搜索: MSFragger

git clone https://github.com/Jeblqr/Omicsomics.git    - Label-free 定量: LFQ

cd Omicsomics    - API 端点: 5 个 (`/convert-raw`, `/maxquant`, `/msfragger`, `/lfq-quantification`, `/complete-pipeline`)



# 2. 一键启动（使用便捷脚本）12. **代谢组学分析 (LC-MS/GC-MS)** ✨ NEW

./scripts/dev-start.sh

    - 特征检测: XCMS, MZmine

# 3. 访问应用    - 谱图注释: GNPS, MS-DIAL

# Backend API: http://localhost:8000/docs    - 定量归一化: median, quantile, PQN

# Frontend: http://localhost:5173    - API 端点: 4 个 (`/feature-detection`, `/spectral-annotation`, `/quantification`, `/complete-pipeline`)

# MinIO Console: http://localhost:9001

```13. **多组学整合 (Multi-omics Integration)** ✨ NEW



### 本地开发    - 无监督整合: MOFA2 (Multi-Omics Factor Analysis)

    - 有监督整合: DIABLO (生物标志物发现)

#### 后端设置    - 样本匹配: 跨组学数据集匹配

    - 通路富集: 多层次功能分析

```bash    - API 端点: 5 个 (`/mofa2`, `/diablo`, `/pathway-enrichment`, `/match-samples`, `/complete-pipeline`)

# 1. 创建环境

conda create -n omicsomics python=3.1114. **可视化与交互**

conda activate omicsomics    - 火山图 (Volcano plot) - 差异表达可视化

    - UMAP/PCA - 降维可视化

# 2. 安装依赖    - 热图 (Heatmap) - 基因表达矩阵

cd backend && pip install -e .    - IGV.js - 基因组浏览器集成

    - QC 指标分布图

# 3. 启动数据库和存储    - API 端点: 7 个 (Plotly 格式数据导出)

docker-compose up -d db minio

**总计**: 8 大组学模块，39 个 API 端点组 🎉

# 4. 配置环境变量

export DATABASE_URL="postgresql+asyncpg://omics_user:omics_pass@localhost:5432/omicsomics"## 技术栈

export SECRET_KEY="your-secret-key-here"

export MINIO_ENDPOINT="localhost:9000"- **后端**: FastAPI 0.121.0 + Python 3.11

export MINIO_ACCESS_KEY="minioadmin"- **数据库**: PostgreSQL 18.0 (AsyncIO 支持)

export MINIO_SECRET_KEY="minioadmin"- **对象存储**: MinIO

- **ORM**: SQLAlchemy 2.0 (async)

# 5. 初始化数据库- **迁移**: Alembic

alembic upgrade head- **认证**: JWT (python-jose)

- **密码哈希**: bcrypt

# 6. 启动服务

uvicorn app.main:app --host 0.0.0.0 --port 8000 --reload## 快速开始

````

详细部署指南请参见 [DEPLOYMENT.md](DEPLOYMENT.md)

#### 前端设置

### 1. 环境准备

````bash

cd frontend```bash

npm installmicromamba create -n omicsomics-dev python=3.11

npm run devmicromamba activate omicsomics-dev

# 访问 http://localhost:5173micromamba install -n omicsomics-dev postgresql

```cd backend && pip install -e .

````

---

### 2. 初始化

## 📚 API 文档

````bash

访问 **http://localhost:8000/docs** 查看完整的 Swagger UI 文档。# 数据库

initdb -D local_db_data

### 快速示例pg_ctl -D local_db_data -l postgresql.log start

createdb omicsomics

```bashcd backend && alembic upgrade head

# 1. 注册用户```

curl -X POST "http://localhost:8000/api/v1/register" \

  -H "Content-Type: application/json" \### 3. 启动服务（3 个终端）

  -d '{"email": "user@example.com", "password": "Pass123", "full_name": "John"}'

```bash

# 2. 登录获取 Token# 终端1: PostgreSQL

TOKEN=$(curl -X POST "http://localhost:8000/api/v1/login/access-token" \pg_ctl -D local_db_data start

  -H "Content-Type: application/x-www-form-urlencoded" \

  -d "username=user@example.com&password=Pass123" | jq -r '.access_token')# 终端2: MinIO

./scripts/start_minio.sh

# 3. 创建项目

curl -X POST "http://localhost:8000/api/v1/projects/" \# 终端3: FastAPI

  -H "Authorization: Bearer $TOKEN" \cd backend

  -H "Content-Type: application/json" \export SECRET_KEY="your-key" DATABASE_URL="postgresql+asyncpg://jeblqr@localhost/omicsomics"

  -d '{"name": "My Project", "description": "GWAS Study"}'uvicorn app.main:app --host 127.0.0.1 --port 8001 --reload

````

# 4. 运行 GWAS QC

curl -X POST "http://localhost:8000/api/v1/gwas/qc" \### 4. 访问

-H "Authorization: Bearer $TOKEN" \

-H "Content-Type: application/json" \- **API 文档**: http://127.0.0.1:8001/docs

-d '{- **MinIO 控制台**: http://127.0.0.1:9003 (minioadmin/minioadmin123)

    "sample_id": "sample_001",

    "bfile": "path/to/data",## API 示例

    "maf": 0.01,

    "geno": 0.02,```bash

    "mind": 0.02,# 注册

    "hwe": 0.000001curl -X POST "http://localhost:8001/api/v1/register" \

}' -H "Content-Type: application/json" \

````-d '{"email": "user@example.com", "password": "pass123", "full_name": "User"}'



### API 端点总览# 登录

TOKEN=$(curl -s -X POST "http://localhost:8001/api/v1/login/access-token" \

| 端点 | 功能 | 认证 |  -H "Content-Type: application/x-www-form-urlencoded" \

|------|------|------|  -d "username=user@example.com&password=pass123" | jq -r '.access_token')

| `/api/v1/register` | 用户注册 | ❌ |

| `/api/v1/login` | 用户登录 | ❌ |# 创建项目

| `/api/v1/projects` | 项目管理 | ✅ |curl -X POST "http://localhost:8001/api/v1/projects/" \

| `/api/v1/samples` | 样本管理 | ✅ |  -H "Authorization: Bearer $TOKEN" \

| `/api/v1/genomics/*` | 基因组学分析 (6个端点) | ✅ |  -H "Content-Type: application/json" \

| `/api/v1/transcriptomics/*` | 转录组学分析 (3个端点) | ✅ |  -d '{"name": "RNA-Seq", "description": "Transcriptomics"}'

| `/api/v1/single-cell/*` | 单细胞分析 (4个端点) | ✅ |

| `/api/v1/epigenomics/*` | 表观基因组学 (5个端点) | ✅ |# 创建样本

| `/api/v1/proteomics/*` | 蛋白质组学 (5个端点) | ✅ |curl -X POST "http://localhost:8001/api/v1/samples/" \

| `/api/v1/metabolomics/*` | 代谢组学 (4个端点) | ✅ |  -H "Authorization: Bearer $TOKEN" \

| `/api/v1/multiomics/*` | 多组学整合 (5个端点) | ✅ |  -H "Content-Type: application/json" \

| `/api/v1/gwas/*` | GWAS 分析 (5个端点) | ✅ |  -d '{"name": "Sample1", "project_id": 1, "metadata_": {"tissue": "liver"}}'

| `/api/v1/visualizations/*` | 可视化 (7个端点) | ✅ |

# 上传文件

---curl -X POST "http://localhost:8001/api/v1/files/upload" \

  -H "Authorization: Bearer $TOKEN" \

## 💻 前端使用  -F "file=@data.fastq" -F "sample_id=1" -F "file_type=fastq"



### 功能特点# 运行QC

curl -X POST "http://localhost:8001/api/v1/qc/fastqc" \

- ✅ 现代化 UI 设计（shadcn/ui + Tailwind CSS）  -H "Authorization: Bearer $TOKEN" \

- ✅ 8 个专业分析页面，24+ 标签页  -H "Content-Type: application/json" \

- ✅ 实时表单验证和错误提示  -d '{"sample_id": 1, "file_ids": [1]}'

- ✅ 加载状态管理```

- ✅ JWT 认证集成

- ✅ 响应式布局## 项目结构



### 使用流程```

Omicsomics/

1. **登录** → 输入邮箱和密码├── backend/              # FastAPI 后端

2. **创建项目** → 填写项目信息│   ├── app/

3. **上传样本** → 添加样本数据│   │   ├── api/         # API 路由

4. **选择分析** → 从 8 个模块中选择│   │   ├── models/      # 数据库模型

5. **配置参数** → 设置分析参数│   │   ├── schemas/     # Pydantic schemas

6. **提交任务** → 后台处理│   │   ├── services/    # 业务逻辑

7. **查看结果** → JSON 格式展示│   │   ├── storage/     # S3 客户端

│   │   └── workflows/   # 工作流执行器

---│   └── alembic/         # 数据库迁移

├── scripts/             # 启动脚本

## 🐳 部署指南├── bin/                 # 二进制文件(MinIO)

├── local_db_data/       # PostgreSQL 数据

### Docker Compose 部署└── local_minio_data/    # MinIO 存储

````

````bash

# 1. 克隆仓库## 数据模型

git clone https://github.com/Jeblqr/Omicsomics.git

cd Omicsomics```

User → Projects → Samples → Files

# 2. 配置环境变量                         └→ Workflows

cp .env.example .env```

# 编辑 .env，设置生产密钥

## 开发

# 3. 启动服务

docker-compose up -d```bash

# 数据库迁移

# 4. 初始化数据库alembic revision --autogenerate -m "description"

docker-compose exec backend alembic upgrade headalembic upgrade head

````

# 运行测试

### 环境变量配置 cd backend && pytest

````

```bash

# .env 文件示例## 安全

DATABASE_URL=postgresql+asyncpg://user:pass@db:5432/omicsomics

SECRET_KEY=your-super-secret-key-change-this**生产环境务必**:

MINIO_ENDPOINT=minio:9000

MINIO_ACCESS_KEY=minioadmin1. 更改 `SECRET_KEY`

MINIO_SECRET_KEY=minioadmin2. 更改 MinIO 凭据

MINIO_BUCKET=omicsomics3. 启用 HTTPS

```4. 使用强密码

5. 定期备份

### 生产环境建议

## 路线图

- ✅ 使用强密码和密钥

- ✅ 启用 HTTPS (Nginx + Let's Encrypt)### MVP ✅

- ✅ 配置防火墙规则

- ✅ 定期备份数据库- ✅ 认证、项目、样本、文件、工作流、QC

- ✅ 监控系统资源

### v1.0 (计划)

---

- 单细胞、蛋白质组、可视化、MultiQC

## 🧪 测试

### v2.0 (未来)

### 运行后端测试

- 多组学整合、ML、实时协作、插件系统

```bash

cd backend## 许可



# 方式1: 使用 pytestApache v2.0 License

pytest

---

# 方式2: 使用测试脚本（推荐）

cd .. && ./scripts/run-tests.sh**详细文档**: [DEPLOYMENT.md](DEPLOYMENT.md) | **API**: http://127.0.0.1:8001/docs


# 查看测试覆盖率
pytest --cov=app --cov-report=html
````

### 测试脚本功能

`./scripts/run-tests.sh` 提供交互式菜单：

```
1) 运行所有测试
2) 快速测试 (auth + projects)
3) GWAS 测试
4) 认证测试
5) 项目管理测试
6) 自定义测试路径
```

### 当前测试状态

```
✅ 11 个测试通过 (Auth, Projects, Health)
⏳ 54 个测试待数据 (分析模块需要示例数据)
```

---

## 📁 项目结构

```
Omicsomics/
├── backend/                    # FastAPI 后端
│   ├── app/
│   │   ├── api/               # API 路由 (44个端点)
│   │   ├── models/            # 数据库模型
│   │   ├── schemas/           # Pydantic 模式
│   │   ├── pipelines/         # 分析管道 (9个模块)
│   │   ├── dependencies.py    # 依赖注入
│   │   ├── database.py        # 数据库连接
│   │   └── main.py           # 应用入口
│   ├── alembic/              # 数据库迁移
│   ├── tests/                # 测试文件
│   └── requirements/         # 依赖列表
│
├── frontend/                  # React 前端
│   ├── src/
│   │   ├── components/ui/    # shadcn/ui 组件
│   │   ├── pages/
│   │   │   └── analysis/     # 8个分析页面
│   │   ├── hooks/            # 自定义 Hooks
│   │   └── main.tsx         # 入口文件
│   ├── package.json
│   └── vite.config.ts
│
├── infrastructure/           # Docker 配置
│   ├── docker-compose.yml
│   └── Dockerfile
│
├── scripts/                 # 实用脚本
│   ├── dev-start.sh        # 开发环境一键启动
│   └── run-tests.sh        # 测试运行器
│
└── README.md               # 项目文档
```

---

## ❓ 常见问题

### 1. 数据库连接失败

```bash
# 检查 PostgreSQL 状态
docker ps | grep postgres

# 查看日志
docker logs omicsomics-db

# 重启
docker-compose restart db
```

### 2. MinIO 无法访问

```bash
# 访问控制台
# http://localhost:9001
# 用户名: minioadmin
# 密码: minioadmin

# 初始化 bucket
python scripts/init_minio.py
```

### 3. 前端 CORS 错误

在 `backend/app/main.py` 中检查 CORS 配置：

```python
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173"],  # 确保包含前端 URL
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
```

### 4. Token 过期

Token 默认 30 分钟过期，可在 `.env` 中调整：

```bash
ACCESS_TOKEN_EXPIRE_MINUTES=60
```

---

## 📊 项目统计

- **代码量**: 18,500+ 行
  - 后端 Python: 15,000 行
  - 前端 TypeScript: 3,500 行
- **API 端点**: 44 个
- **前端页面**: 8 个核心分析页面
- **测试覆盖**: 17% (基础设施完成)

---

## 🤝 贡献

欢迎贡献！请遵循以下步骤：

1. Fork 本仓库
2. 创建特性分支 (`git checkout -b feature/AmazingFeature`)
3. 提交更改 (`git commit -m 'Add AmazingFeature'`)
4. 推送分支 (`git push origin feature/AmazingFeature`)
5. 开启 Pull Request

---

## 📄 许可证

Apache 2.0 License

---

## 🙏 致谢

- [FastAPI](https://fastapi.tiangolo.com/)
- [React](https://reactjs.org/)
- [shadcn/ui](https://ui.shadcn.com/)
- [SQLAlchemy](https://www.sqlalchemy.org/)
- [MinIO](https://min.io/)

---

## 📞 联系

- **GitHub**: https://github.com/Jeblqr/Omicsomics
- **Issues**: https://github.com/Jeblqr/Omicsomics/issues

---

**版本**: 0.85.0 | **状态**: 🚀 活跃开发中 | **更新**: 2025-01-07
