# Omicsomics 配置和部署指南

## 🎯 项目结构

```
Omicsomics/
├── manage.sh                 # 统一管理脚本 (启动/停止/日志等)
├── .env                      # 环境变量配置 (本地，不提交)
├── .env.example              # 环境变量模板
├── infrastructure/
│   └── docker-compose.yml    # Docker Compose 配置 (唯一)
├── backend/                  # FastAPI 后端
├── frontend/                 # Vue.js 前端
├── scripts/                  # 测试脚本
└── docs/                     # 文档

⚠️  已删除冗余文件:
  - 多余的 docker-compose 文件 (example, local, bak)
  - 所有 .sh 脚本 (统一使用 manage.sh)
  - bin/minio 二进制 (使用 Docker 镜像)
```

## 🚀 快速开始

### 1. 配置环境变量

```bash
# 复制模板
cp .env.example .env

# 编辑配置
nano .env
```

必填配置项：

```bash
# 数据库
POSTGRES_PASSWORD=your_secure_password

# MinIO
MINIO_ROOT_PASSWORD=your_minio_password

# Cloudflared (可选)
CLOUDFLARED_TOKEN=your_token
```

### 2. 启动服务

```bash
# 一键启动
./manage.sh start

# 查看状态
./manage.sh status

# 查看日志
./manage.sh logs
```

### 3. 访问服务

- **前端**: http://localhost:5173
- **后端 API**: http://localhost:8001
- **API 文档**: http://localhost:8001/docs
- **MinIO 控制台**: http://localhost:9001

## 📦 容器架构

所有服务运行在独立容器中：

| 容器名                   | 服务          | 端口       | 说明           |
| ------------------------ | ------------- | ---------- | -------------- |
| omicsomics-db            | PostgreSQL 15 | 5432       | 数据库         |
| omicsomics-redis         | Redis 7       | 6379       | 缓存和消息队列 |
| omicsomics-minio         | MinIO         | 9000, 9001 | 对象存储       |
| omicsomics-backend       | FastAPI       | 8001       | 后端 API       |
| omicsomics-celery-worker | Celery        | -          | 异步任务       |
| omicsomics-frontend      | Vue.js        | 5173       | 前端界面       |
| omicsomics-cloudflared   | Cloudflared   | -          | Tunnel (可选)  |

## 🛠️ 管理命令

### 基本操作

```bash
./manage.sh start      # 启动所有服务
./manage.sh stop       # 停止所有服务
./manage.sh restart    # 重启所有服务
./manage.sh status     # 查看服务状态
```

### 日志和调试

```bash
./manage.sh logs              # 查看所有日志
./manage.sh logs backend      # 查看后端日志
./manage.sh logs frontend     # 查看前端日志
./manage.sh shell backend     # 进入后端容器
```

### 数据管理

```bash
./manage.sh init       # 初始化数据库和 MinIO
./manage.sh clean      # 清理所有容器和数据 (⚠️ 危险)
```

### 测试

```bash
./manage.sh test       # 运行测试
```

## 🔧 手动操作 (高级)

如果不使用 `manage.sh`，可以直接使用 docker compose：

```bash
cd infrastructure

# 启动
docker compose up -d

# 停止
docker compose down

# 查看日志
docker compose logs -f backend

# 进入容器
docker compose exec backend bash

# 重启单个服务
docker compose restart backend

# 查看资源使用
docker compose stats
```

## 🔐 环境变量

所有敏感配置都通过 `.env` 文件管理：

```bash
# 数据库
POSTGRES_USER=postgres
POSTGRES_PASSWORD=changeme
POSTGRES_DB=omicsomics
POSTGRES_PORT=5432

# MinIO
MINIO_ROOT_USER=minioadmin
MINIO_ROOT_PASSWORD=changeme123
MINIO_BUCKET=omicsomics
MINIO_PORT=9000
MINIO_CONSOLE_PORT=9001

# Redis
REDIS_PORT=6379

# 后端
BACKEND_PORT=8001
LOG_LEVEL=INFO

# 前端
FRONTEND_PORT=5173
VITE_API_BASE_URL=http://localhost:8001

# Cloudflared (可选)
CLOUDFLARED_TOKEN=your_token_here
```

## 🔄 开发工作流

### 后端开发

```bash
# 启动服务
./manage.sh start

# 后端代码在 backend/app/ 中修改
# 支持热重载，修改后自动生效

# 查看日志
./manage.sh logs backend
```

### 前端开发

```bash
# 启动服务
./manage.sh start

# 前端代码在 frontend/src/ 中修改
# 支持热重载，修改后自动生效

# 查看日志
./manage.sh logs frontend
```

### 数据库迁移

```bash
cd infrastructure

# 生成迁移
docker compose exec backend alembic revision --autogenerate -m "description"

# 运行迁移
docker compose exec backend alembic upgrade head

# 回滚迁移
docker compose exec backend alembic downgrade -1
```

## 🐛 故障排查

### 服务无法启动

```bash
# 查看详细日志
./manage.sh logs

# 检查容器状态
./manage.sh status

# 重启服务
./manage.sh restart
```

### 端口冲突

编辑 `.env` 文件修改端口：

```bash
POSTGRES_PORT=5433
BACKEND_PORT=8002
FRONTEND_PORT=5174
```

### 清理重启

```bash
# 停止服务
./manage.sh stop

# 清理所有数据 (⚠️ 会删除数据)
./manage.sh clean

# 重新初始化
./manage.sh init
./manage.sh start
```

### 数据库连接失败

```bash
# 检查数据库容器
docker ps | grep omicsomics-db

# 查看数据库日志
./manage.sh logs db

# 重启数据库
cd infrastructure
docker compose restart db
```

## 📚 相关文档

- [README.md](README.md) - 项目介绍
- [DEPLOYMENT.md](docs/DEPLOYMENT.md) - 详细部署文档
- [ARCHITECTURE.md](docs/ARCHITECTURE.md) - 系统架构
- [PROJECT_STRUCTURE.md](docs/PROJECT_STRUCTURE.md) - 项目结构

## ⚠️ 重要说明

1. **环境文件**: `.env` 文件包含敏感信息，已在 `.gitignore` 中排除
2. **数据持久化**: 数据存储在 Docker volumes 中
3. **网络隔离**: 所有容器在 `omicsomics-network` 网络中通信
4. **健康检查**: 关键服务配置了健康检查，确保依赖顺序
5. **热重载**: 开发模式下，代码修改会自动重载

## 🎓 最佳实践

1. **使用管理脚本**: 优先使用 `./manage.sh` 而不是直接 docker compose
2. **查看日志**: 出问题先看日志 `./manage.sh logs`
3. **环境变量**: 不同环境使用不同的 `.env` 文件
4. **定期备份**: 备份 PostgreSQL 数据和 MinIO 数据
5. **版本控制**: 不要提交 `.env` 文件到 git
