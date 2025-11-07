# Omicsomics 开发环境启动指南

本项目包含完整的多组学分析平台，可从源码直接部署。

## 前置要求

1. **Micromamba 环境**

   ```bash
   micromamba create -n omicsomics-dev python=3.11
   micromamba activate omicsomics-dev
   ```

2. **安装 Python 依赖**

   ```bash
   cd backend
   pip install -e .
   ```

3. **PostgreSQL**  
   本地已通过 micromamba 安装：
   ```bash
   micromamba install -n omicsomics-dev postgresql
   ```

## 启动服务

### 方式一：使用启动脚本（推荐）

在**3 个独立终端**中分别运行：

#### 终端 1: PostgreSQL

```bash
cd /home/jeblqr/data1/projects/Omicsomics
pg_ctl -D local_db_data -l postgresql.log start
# 查看状态
pg_isready
```

#### 终端 2: MinIO

```bash
cd /home/jeblqr/data1/projects/Omicsomics
./scripts/start_minio.sh
# 应该看到 "API: http://127.0.0.1:9002"
```

#### 终端 3: FastAPI Backend

```bash
cd /home/jeblqr/data1/projects/Omicsomics/backend
micromamba activate omicsomics-dev
export SECRET_KEY="your-secret-key-here-change-in-production"
export DATABASE_URL="postgresql+asyncpg://jeblqr@localhost/omicsomics"
uvicorn app.main:app --host 127.0.0.1 --port 8001 --reload
```

### 方式二：一键启动（后台运行）

⚠️ 注意：这将在后台启动所有服务

```bash
cd /home/jeblqr/data1/projects/Omicsomics
./scripts/start_all.sh
```

## 首次初始化

1. **创建数据库**（仅首次）

   ```bash
   cd /home/jeblqr/data1/projects/Omicsomics
   micromamba run -n omicsomics-dev createdb omicsomics
   ```

2. **运行数据库迁移**

   ```bash
   cd backend
   micromamba run -n omicsomics-dev alembic upgrade head
   ```

3. **初始化 MinIO Bucket**  
   这会在首次上传文件时自动完成，也可以手动运行：
   ```bash
   micromamba run -n omicsomics-dev python scripts/init_minio.py
   ```

## 验证服务

```bash
# 检查 PostgreSQL
pg_isready

# 检查 FastAPI
curl http://127.0.0.1:8001/health

# 检查 MinIO
curl http://127.0.0.1:9002/minio/health/live
```

## 访问接口

- **API 文档**: http://127.0.0.1:8001/docs
- **MinIO 控制台**: http://127.0.0.1:9003
  - 用户名: `minioadmin`
  - 密码: `minioadmin123`

## 停止服务

```bash
# 停止 FastAPI (Ctrl+C 或)
lsof -ti:8001 | xargs kill

# 停止 MinIO (Ctrl+C 或)
lsof -ti:9002 | xargs kill

# 停止 PostgreSQL
pg_ctl -D local_db_data stop
```

## 目录结构

```
Omicsomics/
├── backend/               # FastAPI后端
│   ├── app/              # 应用代码
│   ├── alembic/          # 数据库迁移
│   └── pyproject.toml    # Python依赖
├── bin/                  # 二进制文件（MinIO）
├── local_db_data/        # PostgreSQL数据
├── local_minio_data/     # MinIO对象存储
└── scripts/              # 启动脚本
    ├── start_all.sh     # 一键启动
    ├── start_minio.sh   # MinIO启动
    └── init_minio.py    # MinIO初始化
```

## 常见问题

### 端口被占用

```bash
# 检查占用
lsof -i :8001  # FastAPI
lsof -i :9002  # MinIO
lsof -i :5432  # PostgreSQL

# 终止进程
lsof -ti:PORT | xargs kill -9
```

### 数据库连接失败

```bash
# 确保PostgreSQL正在运行
pg_isready

# 重启PostgreSQL
pg_ctl -D local_db_data restart
```

### MinIO 连接失败

确保 MinIO 在独立终端运行并保持前台输出，或者查看日志：

```bash
tail -f minio.log
```
