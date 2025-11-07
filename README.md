# Omicsomics - 统一组学分析平台

一个面向研究与临床的Web平台，支持常见组学数据（基因组学、转录组学、单细胞、表观、蛋白质组、代谢组、宏基因组等）的统一接收、处理和分析。

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
   - 样本的CRUD操作
   - 灵活的JSON元数据支持
   - 样本与项目的关联

4. **文件存储**
   - 基于MinIO的对象存储
   - S3兼容API
   - 文件上传和下载
   - 预签名URL支持

5. **工作流执行**
   - Nextflow流水线集成
   - FastQC质量控制
   - 异步任务执行
   - 工作流状态跟踪
   - 日志记录

6. **质量控制（QC）**
   - FastQC支持
   - 批量QC分析
   - QC结果存储和查询

## 技术栈

- **后端**: FastAPI 0.121.0 + Python 3.11
- **数据库**: PostgreSQL 18.0 (AsyncIO支持)
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

### 3. 启动服务（3个终端）

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

MIT License

---

**详细文档**: [DEPLOYMENT.md](DEPLOYMENT.md) | **API**: http://127.0.0.1:8001/docs
