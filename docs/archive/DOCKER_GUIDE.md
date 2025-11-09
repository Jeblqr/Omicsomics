# Docker 部署指南

## 🐳 使用 Docker 运行 Omicsomics

所有服务都已容器化,方便部署和管理。

## 快速开始

### 1. 启动所有服务

```bash
./docker-start.sh
```

这会启动:

- ✅ PostgreSQL 数据库 (端口 5432)
- ✅ MinIO 对象存储 (端口 9000, 9001)
- ✅ Backend API (端口 8001)
- ✅ Frontend (端口 5173)
- ✅ Cloudflare Tunnel (如果配置)

### 2. 停止所有服务

```bash
./docker-stop.sh
```

### 3. 查看日志

```bash
# 查看后端日志
./docker-logs.sh backend

# 查看前端日志
./docker-logs.sh frontend

# 查看数据库日志
./docker-logs.sh db
```

### 4. 进入容器

```bash
# 进入后端容器
./docker-shell.sh

# 在容器内可以运行命令
# 例如: alembic upgrade head
```

## 📋 服务访问

启动成功后:

| 服务             | URL                        | 说明         |
| ---------------- | -------------------------- | ------------ |
| **前端**         | http://localhost:5173      | React 应用   |
| **后端 API**     | http://localhost:8001      | FastAPI 服务 |
| **API 文档**     | http://localhost:8001/docs | Swagger UI   |
| **MinIO 控制台** | http://localhost:9001      | 对象存储管理 |
| **数据库**       | localhost:5432             | PostgreSQL   |

**MinIO 登录:**

- 用户名: `minio`
- 密码: `minio123`

**数据库登录:**

- 用户名: `postgres`
- 密码: `postgres`
- 数据库: `omicsomics`

## 🔧 开发模式

### 代码热重载

Docker Compose 已配置代码卷挂载:

```yaml
volumes:
  - ../backend/app:/app/app:delegated
  - ../frontend:/usr/src/app:delegated
```

修改代码后:

- **Backend**: 自动重载 (uvicorn --reload)
- **Frontend**: 自动重载 (Vite HMR)

### 数据库迁移

```bash
# 进入后端容器
./docker-shell.sh

# 创建新迁移
alembic revision --autogenerate -m "migration message"

# 应用迁移
alembic upgrade head

# 回滚
alembic downgrade -1
```

## 🏗️ 架构说明

### 服务依赖关系

```
frontend → backend → db
                  → minio
```

- Frontend 依赖 Backend
- Backend 依赖 Database 和 MinIO
- 所有服务在同一 Docker 网络中通信

### 数据持久化

使用 Docker volumes:

- `pgdata`: PostgreSQL 数据
- `minio-data`: MinIO 对象存储数据

数据不会因为容器重启而丢失。

### 环境变量

Backend 环境变量 (在 `docker-compose.yml` 中配置):

```yaml
DATABASE_URL: postgresql+asyncpg://postgres:postgres@db:5432/omicsomics
OBJECT_STORAGE_ENDPOINT: http://minio:9000
OBJECT_STORAGE_ACCESS_KEY: minio
OBJECT_STORAGE_SECRET_KEY: minio123
OBJECT_STORAGE_BUCKET: omicsomics
```

## 🐛 故障排除

### 检查容器状态

```bash
cd infrastructure
docker compose ps
```

### 查看所有日志

```bash
cd infrastructure
docker compose logs
```

### 重启单个服务

```bash
cd infrastructure
docker compose restart backend
```

### 重新构建镜像

```bash
cd infrastructure
docker compose build --no-cache backend
docker compose up -d backend
```

### 清理所有数据

```bash
cd infrastructure
docker compose down -v  # ⚠️ 会删除所有数据!
```

### 端口冲突

如果端口被占用,修改 `infrastructure/docker-compose.yml`:

```yaml
ports:
  - "8002:8001" # 使用 8002 代替 8001
```

## 📦 生产部署

### 1. 修改环境变量

创建 `.env` 文件:

```bash
# Database
DATABASE_URL=postgresql+asyncpg://prod_user:secure_password@db:5432/omicsomics

# Object Storage
OBJECT_STORAGE_ACCESS_KEY=production_key
OBJECT_STORAGE_SECRET_KEY=production_secret

# API
SECRET_KEY=your-super-secret-key-here
```

### 2. 禁用开发功能

修改 `docker-compose.yml`:

```yaml
backend:
  command: ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8001"]
  # 移除 --reload 参数
```

### 3. 配置反向代理

使用 Nginx 或 Caddy 作为反向代理:

```nginx
server {
    listen 80;
    server_name your-domain.com;

    location / {
        proxy_pass http://localhost:5173;
    }

    location /api {
        proxy_pass http://localhost:8001;
    }
}
```

### 4. 启用 HTTPS

使用 Let's Encrypt + Certbot 或 Cloudflare Tunnel (已配置)。

## 🔄 更新部署

```bash
# 拉取最新代码
git pull

# 重新构建并启动
./docker-stop.sh
./docker-start.sh
```

## 📊 监控

查看资源使用:

```bash
docker stats
```

查看容器日志大小:

```bash
docker compose logs --tail=100 backend
```

## 🆘 常见问题

### Q: 构建失败?

```bash
# 清理 Docker 缓存
docker system prune -a
# 重新构建
./docker-start.sh
```

### Q: 数据库连接失败?

检查数据库是否就绪:

```bash
docker compose exec db pg_isready -U postgres
```

### Q: MinIO 初始化失败?

手动创建 bucket:

```bash
# 访问 http://localhost:9001
# 登录后创建 bucket: omicsomics
```

### Q: 前端无法连接后端?

检查 CORS 配置和环境变量:

```bash
./docker-logs.sh backend | grep CORS
```

## 📝 开发工作流

1. **启动服务**: `./docker-start.sh`
2. **修改代码**: 编辑器中修改,自动重载
3. **查看日志**: `./docker-logs.sh backend`
4. **测试 API**: http://localhost:8001/docs
5. **停止服务**: `./docker-stop.sh`

---

**优势:**

- ✅ 环境一致性
- ✅ 快速部署
- ✅ 易于迁移
- ✅ 隔离依赖
- ✅ 简化配置
