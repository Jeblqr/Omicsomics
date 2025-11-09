# Omicsomics 项目结构说明

## 📁 目录结构总览

```
Omicsomics/
├── backend/              # FastAPI 后端服务
├── frontend/             # Vue 3 + TypeScript 前端
├── infrastructure/       # Docker 和基础设施配置
├── scripts/             # 开发和测试脚本
├── docs/                # 项目文档
├── test_data/           # 测试样本数据
├── workflows/           # CWL 和 Nextflow 工作流定义
├── bin/                 # 二进制可执行文件
├── data/                # 运行时数据目录
├── local_db_data/       # PostgreSQL 本地数据
└── [配置文件]           # 根目录配置文件
```

## 🔧 核心目录详解

### 1. `backend/` - 后端服务

```
backend/
├── app/                  # 应用代码
│   ├── api/             # API 路由和端点
│   │   ├── v1/          # API v1 版本
│   │   │   ├── auth.py           # 认证和授权
│   │   │   ├── projects.py       # 项目管理
│   │   │   ├── data_files.py     # 数据文件上传
│   │   │   ├── pipelines.py      # 流程管理
│   │   │   ├── runs.py           # 运行管理
│   │   │   ├── visualizations.py # 可视化
│   │   │   └── ...
│   │   └── deps.py      # 依赖注入
│   ├── converters/      # 数据格式转换器
│   │   ├── genomics/    # 基因组学 (VCF, BED, GFF)
│   │   ├── transcriptomics/  # 转录组学 (CSV, TSV, counts)
│   │   ├── proteomics/  # 蛋白质组学 (mzML, MGF, CSV)
│   │   └── base.py      # 转换器基类
│   ├── core/            # 核心功能
│   │   ├── async_processor.py  # 异步文件处理
│   │   ├── pipeline_executor.py # 流程执行引擎
│   │   └── security.py         # 安全工具
│   ├── crud/            # 数据库 CRUD 操作
│   ├── models/          # SQLAlchemy ORM 模型
│   ├── schemas/         # Pydantic 数据模型
│   ├── main.py          # FastAPI 应用入口
│   ├── database.py      # 数据库配置
│   ├── settings.py      # 应用设置
│   └── celery_app.py    # Celery 配置
├── alembic/             # 数据库迁移
│   ├── versions/        # 迁移版本
│   └── env.py           # Alembic 配置
├── tests/               # 测试代码
├── Dockerfile           # Docker 镜像定义
├── pyproject.toml       # Python 项目配置
└── alembic.ini          # Alembic 配置文件
```

**关键功能**:

- ✅ RESTful API (FastAPI)
- ✅ 异步文件处理 (Celery)
- ✅ 多格式数据转换
- ✅ 流程执行引擎
- ✅ 用户认证 (JWT)
- ✅ 对象存储 (MinIO)

### 2. `frontend/` - 前端应用

```
frontend/
├── src/
│   ├── assets/          # 静态资源
│   ├── components/      # Vue 组件
│   │   ├── common/      # 通用组件
│   │   ├── dashboard/   # 仪表板组件
│   │   ├── projects/    # 项目管理组件
│   │   ├── pipelines/   # 流程管理组件
│   │   ├── runs/        # 运行管理组件
│   │   └── visualizations/  # 数据可视化组件
│   ├── views/           # 页面视图
│   │   ├── Dashboard.vue
│   │   ├── Projects.vue
│   │   ├── Pipelines.vue
│   │   ├── Runs.vue
│   │   └── ...
│   ├── stores/          # Pinia 状态管理
│   │   ├── auth.ts
│   │   ├── projects.ts
│   │   ├── pipelines.ts
│   │   └── ...
│   ├── router/          # Vue Router 配置
│   ├── services/        # API 服务层
│   │   └── api.ts       # Axios 实例和 API 调用
│   ├── types/           # TypeScript 类型定义
│   ├── App.vue          # 根组件
│   └── main.ts          # 应用入口
├── public/              # 公共静态文件
├── Dockerfile           # Docker 镜像定义
├── vite.config.ts       # Vite 配置
├── tsconfig.json        # TypeScript 配置
└── package.json         # Node.js 依赖
```

**关键功能**:

- ✅ Vue 3 Composition API
- ✅ TypeScript 类型安全
- ✅ Pinia 状态管理
- ✅ Element Plus UI 框架
- ✅ 响应式设计
- ✅ 实时数据更新

### 3. `infrastructure/` - 基础设施

```
infrastructure/
├── docker-compose.example.yml  # Docker Compose 模板
└── [docker-compose.yml]        # 本地配置 (gitignore)
```

**服务组件**:

- `db` - PostgreSQL 15 数据库
- `minio` - 对象存储服务
- `redis` - 缓存和消息队列
- `backend` - FastAPI 应用
- `celery-worker` - 异步任务处理
- `frontend` - Vue 前端
- `cloudflared` - Cloudflare Tunnel (可选)

### 4. `scripts/` - 脚本工具

```
scripts/
├── dev-start.sh                  # 开发环境启动
├── start_all.sh                  # 启动所有服务
├── start_frontend.sh             # 启动前端
├── start_minio.sh                # 启动 MinIO
├── init_minio.py                 # 初始化 MinIO
├── quick_test.py                 # 快速测试
├── test_async_processing.py     # 测试异步处理
├── test_pipeline_execution.py   # 测试流程执行
├── verify_fixes.py               # 验证修复
└── ...                           # 其他测试脚本
```

### 5. `docs/` - 文档

```
docs/
├── README.md                      # 文档索引
├── ARCHITECTURE.md                # 架构设计
├── DEPLOYMENT.md                  # 部署指南
├── PROJECT_CLEANUP_REPORT.md     # 清理报告
├── TESTING.md                     # 测试指南
├── api/                           # API 文档
├── architecture/                  # 架构图表
└── archive/                       # 归档的旧文档
```

## 🗂️ 配置文件

### 根目录配置

| 文件                         | 用途         | Git 追踪 |
| ---------------------------- | ------------ | -------- |
| `.gitignore`                 | Git 忽略规则 | ✅       |
| `.env.example`               | 环境变量模板 | ✅       |
| `.env`                       | 实际环境变量 | ❌       |
| `README.md`                  | 项目说明     | ✅       |
| `IMPLEMENTATION_COMPLETE.md` | 实现报告     | ✅       |
| `SECURITY.md`                | 安全策略     | ✅       |

### 便捷脚本（根目录）

| 脚本              | 功能                 | 用途     |
| ----------------- | -------------------- | -------- |
| `docker-start.sh` | 启动所有 Docker 服务 | 生产环境 |
| `docker-stop.sh`  | 停止所有服务         | 管理     |
| `docker-logs.sh`  | 查看服务日志         | 调试     |
| `docker-shell.sh` | 进入容器 shell       | 调试     |

## 📊 数据目录

### 运行时数据

```
data/               # 应用数据（MinIO 或本地）
local_db_data/      # PostgreSQL 数据（Docker volume）
test_data/          # 测试样本
  ├── sample_genomics.vcf
  ├── sample_proteomics.csv
  └── sample_transcriptomics.csv
```

### 工作流定义

```
workflows/
├── cwl/            # Common Workflow Language 定义
│   ├── genomics/
│   ├── transcriptomics/
│   └── proteomics/
└── nextflow/       # Nextflow 流程定义
    └── ...
```

## 🔐 安全和敏感文件

### ❌ 不应提交到 Git

```
.env                              # 环境变量（包含密钥）
infrastructure/docker-compose.yml # 本地配置（包含密钥）
*.pem, *.key, *.crt              # 证书文件
__pycache__/                     # Python 缓存
node_modules/                    # Node.js 依赖
*.log                            # 日志文件
local_db_data/                   # 数据库数据
data/                            # 应用数据
```

### ✅ 应提交到 Git

```
.env.example                      # 环境变量模板
infrastructure/docker-compose.example.yml  # Docker 配置模板
backend/app/                      # 应用代码
frontend/src/                     # 前端代码
docs/                             # 文档
scripts/                          # 脚本（不含密钥）
test_data/sample_*.{vcf,csv}     # 样本数据
```

## 🎯 重要目录说明

### 需要创建的目录（首次运行）

这些目录在运行时自动创建，不需要提交到 git：

```bash
data/              # MinIO 数据
local_db_data/     # PostgreSQL 数据
backend/app/__pycache__/  # Python 缓存
frontend/node_modules/    # Node.js 依赖
```

### 持久化数据

```
local_db_data/     # PostgreSQL - 用户、项目、文件元数据
data/              # MinIO - 上传的数据文件、结果文件
```

## 📝 配置文件管理

### 首次设置

1. **复制环境变量模板**:

   ```bash
   cp .env.example .env
   # 编辑 .env，填入真实值
   ```

2. **复制 Docker Compose 配置**:

   ```bash
   cp infrastructure/docker-compose.example.yml infrastructure/docker-compose.yml
   # 编辑 docker-compose.yml，配置密钥
   ```

3. **更新密钥**:
   - `POSTGRES_PASSWORD` - PostgreSQL 密码
   - `MINIO_ROOT_PASSWORD` - MinIO 密码
   - `CLOUDFLARED_TOKEN` - Cloudflare Tunnel token (可选)

### 版本控制策略

- ✅ **提交**: 所有 `.example` 配置文件
- ❌ **不提交**: 包含真实密钥的配置文件
- 📝 **文档化**: 在 README 中说明配置方法

## 🔄 开发工作流

### 日常开发

```bash
# 1. 启动服务
./docker-start.sh

# 2. 查看日志
./docker-logs.sh backend
./docker-logs.sh frontend

# 3. 进入容器调试
./docker-shell.sh backend

# 4. 停止服务
./docker-stop.sh
```

### 代码目录

- **后端开发**: `backend/app/`
- **前端开发**: `frontend/src/`
- **数据库迁移**: `backend/alembic/versions/`
- **测试代码**: `backend/tests/`, `scripts/test_*.py`

## 📚 相关文档

- [README.md](../README.md) - 项目介绍和快速开始
- [ARCHITECTURE.md](ARCHITECTURE.md) - 架构设计文档
- [DEPLOYMENT.md](DEPLOYMENT.md) - 部署指南
- [SECURITY.md](../SECURITY.md) - 安全策略
- [TESTING.md](TESTING.md) - 测试指南

---

**注意**: 此文档描述的是清理后的项目结构。确保遵循安全最佳实践，不要提交敏感信息到 git。
