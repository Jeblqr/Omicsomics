# 项目启动状态

## ✅ 所有服务已成功启动并运行

### Docker 容器状态
```bash
cd infrastructure && docker compose ps
```

所有 5 个容器正在运行：
- ✅ **backend** - FastAPI 后端 (端口 8001)
- ✅ **frontend** - React + Vite 前端 (端口 5173)
- ✅ **db** - PostgreSQL 15 数据库 (端口 5432)
- ✅ **minio** - MinIO 对象存储 (端口 9000-9001)
- ✅ **cloudflared** - Cloudflare 隧道

### 数据库迁移已完成
- 创建了初始迁移文件: `5e8f44ef46d2_initial_migration.py`
- 已成功运行迁移，创建了所有数据库表：
  - users
  - projects
  - samples
  - files
  - alembic_version

### 测试的 API 端点

#### ✅ 健康检查
```bash
curl http://localhost:8001/healthz
# {"status":"ok"}
```

#### ✅ 用户注册
```bash
curl -X POST http://localhost:8001/api/v1/register \
  -H "Content-Type: application/json" \
  -d '{"email":"test@example.com","password":"testpass123","full_name":"Test User"}'
# {"email":"test@example.com","full_name":"Test User","id":1}
```

#### ✅ 用户登录
```bash
curl -X POST http://localhost:8001/api/v1/login/access-token \
  -H "Content-Type: application/x-www-form-urlencoded" \
  -d "username=test@example.com&password=testpass123"
# {"access_token":"eyJ...","token_type":"bearer"}
```

### 前端配置
- **访问地址**: http://localhost:5173/
- **API 代理**: 配置为代理 `/api` 到后端
- **容器间通信**: 使用 Docker 服务名 `backend:8001`

### 配置更新
1. ✅ `vite.config.ts` - 支持环境变量 `VITE_API_TARGET`
2. ✅ `docker-compose.yml` - 添加 `VITE_API_TARGET=http://backend:8001`

## 如何使用

### 启动所有服务
```bash
cd infrastructure
docker compose up -d
```

### 查看日志
```bash
# 后端日志
docker logs infrastructure-backend-1 -f

# 前端日志
docker logs infrastructure-frontend-1 -f

# 数据库日志
docker logs infrastructure-db-1 -f
```

### 停止所有服务
```bash
cd infrastructure
docker compose down
```

### 访问应用
- **前端应用**: http://localhost:5173/
- **后端 API**: http://localhost:8001/
- **API 文档**: http://localhost:8001/docs
- **MinIO 控制台**: http://localhost:9001/ (用户名: minio, 密码: minio123)

## 测试建议

### 1. 在浏览器中测试
1. 访问 http://localhost:5173/
2. 注册新用户
3. 登录
4. 创建项目
5. 测试各个分析页面

### 2. 使用 pytest 测试后端
```bash
export TEST_DATABASE_URL="postgresql+asyncpg://postgres:postgres@localhost:5432/omicsomics"
cd backend
pytest tests/ -v
```

## 故障排查

### 前端无法访问后端
1. 检查容器是否运行: `docker compose ps`
2. 检查后端日志: `docker logs infrastructure-backend-1`
3. 从前端容器测试连接: `docker exec infrastructure-frontend-1 wget -O - http://backend:8001/healthz`

### 数据库连接错误
1. 检查数据库是否健康: `docker compose ps` (应该显示 "healthy")
2. 运行迁移: `docker exec infrastructure-backend-1 alembic upgrade head`

### 端口已被占用
- 修改 `docker-compose.yml` 中的端口映射
- 或停止占用端口的其他服务

## 开发模式（本地开发）

如果想在本地开发而不使用 Docker：

### 后端
```bash
cd backend
# 使用 micromamba 环境
export DATABASE_URL="postgresql+asyncpg://postgres:postgres@localhost:5432/omicsomics"
micromamba run -n omicsomics-dev uvicorn app.main:app --reload --port 8001
```

### 前端
```bash
cd frontend
npm run dev
# 访问 http://localhost:5174/
```

## 下一步
- [ ] 准备样本数据文件
- [ ] 测试文件上传功能
- [ ] 测试各个分析模块
- [ ] 配置生产环境部署
