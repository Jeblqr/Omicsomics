# 本地配置说明

## 📁 本地配置文件

以下文件包含真实的密钥和配置，已保存在本地，**不会提交到 git**：

### 1. `.env`
环境变量配置文件，包含：
- 数据库密码
- MinIO 凭证
- Cloudflared token

**使用方法**：
```bash
# 已自动创建，包含从 docker-compose.local.yml 提取的配置
# 如需修改，直接编辑此文件
nano .env
```

### 2. `infrastructure/docker-compose.local.yml`
完整的 Docker Compose 配置，包含所有真实密钥。

**使用方法**：
```bash
# 启动服务时使用本地配置
cd infrastructure
docker compose -f docker-compose.local.yml up -d

# 或使用根目录的便捷脚本（需修改脚本使其使用 local 配置）
```

## 🔒 Git 保护

`.gitignore` 已配置排除以下文件：
```gitignore
.env
.env.local
infrastructure/docker-compose.local.yml
```

## 🔄 与模板文件的关系

| 模板文件 (提交到git) | 本地文件 (不提交) | 说明 |
|---------------------|------------------|------|
| `.env.example` | `.env` | 环境变量配置 |
| `docker-compose.yml` | `docker-compose.local.yml` | Docker 配置 |

### 配置文件对比

**docker-compose.yml** (模板，使用环境变量):
```yaml
environment:
  POSTGRES_PASSWORD: ${POSTGRES_PASSWORD:-changeme}
  CLOUDFLARED_TOKEN: ${CLOUDFLARED_TOKEN}  # 从环境读取
```

**docker-compose.local.yml** (实际配置，包含真实值):
```yaml
environment:
  POSTGRES_PASSWORD: postgres  # 真实密码
  CLOUDFLARED_TOKEN: eyJhIjoi...  # 真实 token
```

## 🚀 使用说明

### 首次设置（已完成）
✅ 已自动从原 docker-compose.yml 提取配置
✅ 已创建 .env 和 docker-compose.local.yml
✅ 已更新 .gitignore

### 日常使用

**方式 1：使用本地配置文件**
```bash
cd infrastructure
docker compose -f docker-compose.local.yml up -d
```

**方式 2：使用环境变量 + 模板**
```bash
# .env 文件会自动被 docker compose 读取
cd infrastructure
docker compose up -d
```

### 更新配置

如需修改密钥：
1. 编辑 `.env` 或 `infrastructure/docker-compose.local.yml`
2. 重启服务：`docker compose restart`

## ⚠️ 重要提醒

1. **永远不要提交** `.env` 或 `docker-compose.local.yml` 到 git
2. 如果需要在其他机器部署，从 `.env.example` 和 `docker-compose.yml` 复制并填入真实值
3. 定期备份 `.env` 和 `docker-compose.local.yml`（本地备份，不要上传到公共位置）

## 🔐 安全性

- ✅ 真实密钥仅存在于本地文件
- ✅ Git 仓库中只有模板文件
- ✅ `.gitignore` 保护敏感文件
- ✅ 已创建 `SECURITY.md` 说明安全策略

查看更多安全信息：[SECURITY.md](../SECURITY.md)
