# Omicsomics 部署指南

本文档详细说明如何部署 Omicsomics 平台的各种场景。

---

## 目录

- [系统要求](#系统要求)
- [Docker 部署（推荐）](#docker-部署推荐)
- [本地开发环境](#本地开发环境)
- [生产环境部署](#生产环境部署)
- [环境变量配置](#环境变量配置)
- [数据库迁移](#数据库迁移)
- [故障排除](#故障排除)

---

## 系统要求

### 最低配置

- **CPU**: 4 核
- **内存**: 8GB RAM
- **存储**: 50GB 可用空间
- **操作系统**: Linux (Ubuntu 20.04+), macOS, Windows (WSL2)

### 推荐配置

- **CPU**: 8 核或更多
- **内存**: 16GB RAM 或更多
- **存储**: 100GB+ SSD
- **操作系统**: Ubuntu 22.04 LTS

### 软件依赖

- Docker 20.10+
- Docker Compose 2.0+
- Git 2.30+
- (可选) Python 3.11+ for 本地开发
- (可选) Node.js 18+ for 前端开发

---

## Docker 部署（推荐）

### 快速开始

```bash
# 1. 克隆仓库
git clone https://github.com/Jeblqr/Omicsomics.git
cd Omicsomics

# 2. 启动所有服务
cd infrastructure
docker compose up -d

# 3. 检查服务状态
docker compose ps

# 4. 查看日志
docker compose logs -f
```

### 服务访问

启动成功后，可以通过以下地址访问：

- **前端应用**: http://localhost:5173
- **后端 API 文档**: http://localhost:8001/docs
- **MinIO 控制台**: http://localhost:9001 (minioadmin/minioadmin)
- **PostgreSQL**: localhost:5432

### 停止和重启服务

```bash
# 停止所有服务
docker compose down

# 停止并删除数据卷（警告：会删除所有数据）
docker compose down -v

# 重启特定服务
docker compose restart backend
docker compose restart frontend

# 查看特定服务日志
docker compose logs -f backend
docker compose logs -f frontend
```

### 更新服务

```bash
# 拉取最新代码
git pull

# 重新构建并启动
docker compose up -d --build
```

---

## 本地开发环境

### 后端设置

```bash
# 1. 创建 Python 虚拟环境
cd backend
python3.11 -m venv venv
source venv/bin/activate  # Linux/macOS
# 或 Windows: venv\Scripts\activate

# 2. 安装依赖
pip install -e .
pip install -r requirements/dev.txt

# 3. 配置环境变量
cp .env.example .env
# 编辑 .env 文件设置必要的变量

# 4. 启动数据库和 MinIO（使用 Docker）
cd ../infrastructure
docker compose up -d db minio redis

# 5. 运行数据库迁移
cd ../backend
alembic upgrade head

# 6. 启动开发服务器
uvicorn app.main:app --host 0.0.0.0 --port 8001 --reload
```

### 前端设置

```bash
# 1. 安装依赖
cd frontend
npm install

# 2. 配置环境变量
cp .env.example .env.local
# 编辑 .env.local 设置 API URL

# 3. 启动开发服务器
npm run dev

# 4. 访问 http://localhost:5173
```

### 使用便捷脚本

```bash
# 一键启动开发环境
./scripts/dev-start.sh

# 运行测试
./scripts/run-tests.sh

# 初始化 MinIO bucket
python scripts/init_minio.py
```

---

## 生产环境部署

### 准备工作

1. **服务器配置**

   - 确保服务器满足推荐配置
   - 开放必要的端口（80, 443）
   - 配置防火墙规则

2. **域名和 SSL**

   - 注册域名并配置 DNS
   - 获取 SSL 证书（Let's Encrypt 推荐）

3. **数据备份**
   - 配置自动备份策略
   - 定期测试恢复流程

### Docker Compose 生产配置

```bash
# 1. 创建生产环境配置
cp infrastructure/docker-compose.yml infrastructure/docker-compose.prod.yml

# 2. 编辑生产配置
vim infrastructure/docker-compose.prod.yml
```

关键修改项：

```yaml
version: "3.8"

services:
  backend:
    restart: always
    environment:
      - DEBUG=False
      - SECRET_KEY=${SECRET_KEY} # 使用强密钥
      - DATABASE_URL=${DATABASE_URL}
    deploy:
      resources:
        limits:
          cpus: "4"
          memory: 4G

  frontend:
    restart: always
    build:
      context: ../frontend
      dockerfile: Dockerfile.prod # 生产构建
    deploy:
      resources:
        limits:
          cpus: "2"
          memory: 2G

  db:
    restart: always
    volumes:
      - postgres_data:/var/lib/postgresql/data
      - ./backups:/backups # 备份目录
    deploy:
      resources:
        limits:
          cpus: "2"
          memory: 4G

  # 添加 Nginx 反向代理
  nginx:
    image: nginx:alpine
    restart: always
    ports:
      - "80:80"
      - "443:443"
    volumes:
      - ./nginx.conf:/etc/nginx/nginx.conf
      - ./ssl:/etc/nginx/ssl
    depends_on:
      - frontend
      - backend
```

### Nginx 配置示例

创建 `infrastructure/nginx.conf`:

```nginx
events {
    worker_connections 1024;
}

http {
    upstream backend {
        server backend:8001;
    }

    upstream frontend {
        server frontend:5173;
    }

    server {
        listen 80;
        server_name your-domain.com;
        return 301 https://$server_name$request_uri;
    }

    server {
        listen 443 ssl http2;
        server_name your-domain.com;

        ssl_certificate /etc/nginx/ssl/cert.pem;
        ssl_certificate_key /etc/nginx/ssl/key.pem;

        # 前端
        location / {
            proxy_pass http://frontend;
            proxy_set_header Host $host;
            proxy_set_header X-Real-IP $remote_addr;
        }

        # 后端 API
        location /api/ {
            proxy_pass http://backend;
            proxy_set_header Host $host;
            proxy_set_header X-Real-IP $remote_addr;
            proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
            proxy_set_header X-Forwarded-Proto $scheme;
        }

        # API 文档
        location /docs {
            proxy_pass http://backend;
        }

        location /openapi.json {
            proxy_pass http://backend;
        }
    }
}
```

### 启动生产环境

```bash
# 使用生产配置启动
docker compose -f infrastructure/docker-compose.prod.yml up -d

# 查看状态
docker compose -f infrastructure/docker-compose.prod.yml ps

# 查看日志
docker compose -f infrastructure/docker-compose.prod.yml logs -f
```

---

## 环境变量配置

### 后端环境变量 (`.env`)

```bash
# === 应用配置 ===
DEBUG=False
SECRET_KEY=your-super-secret-key-change-this-in-production
ALGORITHM=HS256
ACCESS_TOKEN_EXPIRE_MINUTES=30

# === 数据库配置 ===
DATABASE_URL=postgresql+asyncpg://omics_user:strong_password@db:5432/omicsomics
DATABASE_POOL_SIZE=20
DATABASE_MAX_OVERFLOW=10

# === MinIO 对象存储 ===
MINIO_ENDPOINT=minio:9000
MINIO_ACCESS_KEY=your_minio_access_key
MINIO_SECRET_KEY=your_minio_secret_key
MINIO_BUCKET=omicsomics-data
MINIO_SECURE=false  # 生产环境设为 true

# === Redis 配置 ===
REDIS_URL=redis://redis:6379/0
CELERY_BROKER_URL=redis://redis:6379/0
CELERY_RESULT_BACKEND=redis://redis:6379/0

# === 日志配置 ===
LOG_LEVEL=INFO
LOG_FORMAT=json

# === CORS 配置 ===
CORS_ORIGINS=http://localhost:5173,https://your-domain.com

# === 文件上传限制 ===
MAX_UPLOAD_SIZE=10737418240  # 10GB in bytes
```

### 前端环境变量 (`.env.local`)

```bash
VITE_API_URL=http://localhost:8001
VITE_API_TIMEOUT=30000
VITE_MAX_FILE_SIZE=10737418240
```

### 安全建议

**生产环境必须修改:**

1. `SECRET_KEY` - 使用强随机字符串
2. 数据库密码
3. MinIO 凭据
4. 关闭 DEBUG 模式

**生成安全密钥:**

```bash
# Python
python -c "import secrets; print(secrets.token_urlsafe(32))"

# OpenSSL
openssl rand -base64 32
```

---

## 数据库迁移

### 创建新迁移

```bash
cd backend
alembic revision --autogenerate -m "描述变更内容"
```

### 应用迁移

```bash
# 升级到最新版本
alembic upgrade head

# 升级到特定版本
alembic upgrade <revision_id>

# 降级一个版本
alembic downgrade -1

# 查看迁移历史
alembic history

# 查看当前版本
alembic current
```

### Docker 环境中的迁移

```bash
# 在运行的容器中执行迁移
docker compose exec backend alembic upgrade head

# 或在启动时自动执行（在 docker-compose.yml 中配置）
services:
  backend:
    command: >
      sh -c "alembic upgrade head &&
             uvicorn app.main:app --host 0.0.0.0 --port 8001"
```

---

## 故障排除

### 常见问题

#### 1. 数据库连接失败

**症状**: `psycopg2.OperationalError: could not connect to server`

**解决方案**:

```bash
# 检查数据库容器状态
docker compose ps db

# 查看数据库日志
docker compose logs db

# 重启数据库
docker compose restart db

# 确认数据库 URL 正确
echo $DATABASE_URL
```

#### 2. MinIO 无法访问

**症状**: `botocore.exceptions.ConnectionError`

**解决方案**:

```bash
# 检查 MinIO 状态
docker compose ps minio

# 访问 MinIO 控制台
# http://localhost:9001

# 初始化 bucket
python scripts/init_minio.py

# 检查环境变量
echo $MINIO_ENDPOINT
echo $MINIO_ACCESS_KEY
```

#### 3. 前端 CORS 错误

**症状**: `Access to XMLHttpRequest has been blocked by CORS policy`

**解决方案**:

在 `backend/app/main.py` 中检查 CORS 配置:

```python
from fastapi.middleware.cors import CORSMiddleware

app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173"],  # 添加前端 URL
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
```

#### 4. 端口占用

**症状**: `bind: address already in use`

**解决方案**:

```bash
# 查找占用端口的进程
lsof -i :5173  # 前端
lsof -i :8001  # 后端
lsof -i :5432  # PostgreSQL

# 杀死进程
kill -9 <PID>

# 或修改 docker-compose.yml 中的端口映射
```

#### 5. Docker 磁盘空间不足

**症状**: `no space left on device`

**解决方案**:

```bash
# 清理未使用的镜像
docker system prune -a

# 清理未使用的卷
docker volume prune

# 清理构建缓存
docker builder prune

# 查看磁盘使用
docker system df
```

### 性能优化

#### 数据库优化

```sql
-- 创建索引加速查询
CREATE INDEX idx_projects_user_id ON projects(user_id);
CREATE INDEX idx_samples_project_id ON samples(project_id);
CREATE INDEX idx_files_sample_id ON files(sample_id);

-- 分析表统计信息
ANALYZE projects;
ANALYZE samples;
ANALYZE files;
```

#### 后端优化

```python
# backend/app/database.py

# 增加连接池大小
engine = create_async_engine(
    DATABASE_URL,
    pool_size=20,  # 增加连接池
    max_overflow=40,
    pool_pre_ping=True,
    echo=False,  # 生产环境关闭 SQL 日志
)
```

#### 前端优化

```javascript
// 启用生产构建
npm run build

// 使用 nginx 提供静态文件
// 启用 gzip 压缩
// 配置浏览器缓存
```

---

## 监控和日志

### 日志查看

```bash
# 查看所有服务日志
docker compose logs -f

# 查看特定服务日志
docker compose logs -f backend

# 查看最近 100 行日志
docker compose logs --tail=100 backend

# 导出日志到文件
docker compose logs backend > backend.log
```

### 健康检查

```bash
# 后端健康检查
curl http://localhost:8001/healthz

# 预期响应: {"status":"ok"}
```

### 性能监控

可以集成以下工具进行监控：

- **Prometheus** - 指标收集
- **Grafana** - 可视化仪表板
- **Sentry** - 错误追踪
- **ELK Stack** - 日志聚合

---

## 备份和恢复

### 数据库备份

```bash
# 创建备份
docker compose exec db pg_dump -U omics_user omicsomics > backup_$(date +%Y%m%d).sql

# 恢复备份
cat backup_20251109.sql | docker compose exec -T db psql -U omics_user omicsomics
```

### MinIO 数据备份

```bash
# 使用 mc (MinIO Client)
mc mirror local/minio-data s3/backup-bucket
```

### 自动化备份脚本

```bash
#!/bin/bash
# backup.sh

DATE=$(date +%Y%m%d_%H%M%S)
BACKUP_DIR="/backups"

# 数据库备份
docker compose exec -T db pg_dump -U omics_user omicsomics | gzip > $BACKUP_DIR/db_$DATE.sql.gz

# 保留最近 30 天的备份
find $BACKUP_DIR -name "db_*.sql.gz" -mtime +30 -delete

echo "Backup completed: db_$DATE.sql.gz"
```

添加到 crontab:

```bash
# 每天凌晨 2 点执行备份
0 2 * * * /path/to/backup.sh
```

---

## 更多资源

- [Docker 文档](https://docs.docker.com/)
- [PostgreSQL 文档](https://www.postgresql.org/docs/)
- [FastAPI 文档](https://fastapi.tiangolo.com/)
- [Nginx 文档](https://nginx.org/en/docs/)

---

**需要帮助？** 请查看 [FAQ](../docs/FAQ.md) 或提交 [Issue](https://github.com/Jeblqr/Omicsomics/issues)
