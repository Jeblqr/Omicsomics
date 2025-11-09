# Troubleshooting Guide

## 已修复的问题

### 1. ✅ 后端启动失败 (AttributeError: SINGLE_CELL)

**症状**: 后端容器启动失败，日志显示 `AttributeError: SINGLE_CELL`

**原因**: `backend/app/schemas/tools.py` 中的 `ToolCategory` enum 缺少 `SINGLE_CELL` 定义

**解决方案**:
```python
# backend/app/schemas/tools.py
class ToolCategory(str, Enum):
    # ... 其他类别 ...
    SINGLE_CELL = "single_cell"  # 已添加
    CUSTOM = "custom"
```

**状态**: ✅ 已修复并推送到 git (commit 4460313)

---

### 2. ✅ 项目名称唯一性约束问题

**症状**: 创建同名项目时返回 500 错误：`duplicate key value violates unique constraint "ix_projects_name"`

**原因**: `Project` 模型中 `name` 字段设置了 `unique=True`，不允许不同用户创建同名项目

**解决方案**:
1. 移除模型中的 unique 约束：
```python
# backend/app/models/project.py
name: Mapped[str] = mapped_column(String(255), index=True)  # 移除 unique=True
```

2. 创建并应用数据库迁移：
```bash
alembic revision --autogenerate -m "Remove unique constraint from project name"
alembic upgrade head
```

**状态**: ✅ 已修复，数据库迁移已应用

---

### 3. ✅ 前端代理配置问题

**症状**: 前端无法连接到后端，出现 `ENOTFOUND backend` 错误

**原因**: Vite proxy 的 `changeOrigin: false` 导致某些请求无法正确转发

**解决方案**:
```typescript
// frontend/vite.config.ts
proxy: {
  '/api': {
    target: process.env.VITE_API_TARGET || 'http://localhost:8001',
    changeOrigin: true,  // 改为 true
    secure: false,
    rewrite: (path) => path,
  }
}
```

**状态**: ✅ 已修复并推送到 git (commit 4460313)

---

## 当前状态

### ✅ 正常工作的功能

1. **项目管理**
   - 创建项目 ✅
   - 列出项目 ✅
   - 更新项目 ✅
   - 删除项目 ✅
   - 支持同名项目（不同用户）✅

2. **文件上传**
   - 文件上传到 MinIO ✅
   - 加密存储 ✅
   - 下载文件 ✅
   - 解密下载 ✅

3. **Pipeline Templates**
   - 列出模板 ✅
   - 按类别过滤 ✅
   - 获取模板详情 ✅

4. **后端服务**
   - 正常启动 ✅
   - API 文档可访问 ✅
   - 数据库连接正常 ✅
   - MinIO 存储正常 ✅

### ⚠️ 已知小问题

1. **OPTIONS 请求返回 400** (custom-pipelines)
   - 影响: 部分 CORS 预检请求失败
   - 优先级: 低
   - 原因: 可能是客户端发送了不正确的预检请求
   - 不影响实际功能使用

2. **历史错误日志**
   - 某些旧的 500 错误日志仍在日志中
   - 这些是修复前的错误，已不再发生
   - 可通过 `docker restart infrastructure-backend-1` 清除日志

---

## 故障排查步骤

### 问题: 创建项目失败

1. **检查后端日志**:
```bash
docker logs infrastructure-backend-1 --tail 50 | grep -i error
```

2. **检查数据库约束**:
```bash
docker exec infrastructure-db-1 psql -U postgres -d omicsomics -c "\d projects"
```

3. **验证迁移状态**:
```bash
docker exec infrastructure-backend-1 alembic current
docker exec infrastructure-backend-1 alembic history
```

### 问题: 文件上传失败

1. **检查 MinIO 状态**:
```bash
docker logs infrastructure-minio-1 --tail 20
```

2. **验证 bucket 存在**:
```bash
docker exec infrastructure-minio-1 mc ls local/
```

3. **检查上传日志**:
```bash
docker logs infrastructure-backend-1 | grep -i upload
```

### 问题: 前端无法连接后端

1. **检查容器网络**:
```bash
docker network inspect infrastructure_default
```

2. **测试容器间连接**:
```bash
docker exec infrastructure-frontend-1 ping -c 2 backend
docker exec infrastructure-frontend-1 wget -O- http://backend:8001/docs
```

3. **检查代理配置**:
```bash
cat frontend/vite.config.ts | grep -A 10 proxy
```

---

## 维护建议

### 定期检查

1. **查看后端健康状态**:
```bash
curl http://localhost:8001/healthz
```

2. **监控日志错误**:
```bash
docker logs infrastructure-backend-1 --since 1h | grep -i error
```

3. **检查数据库大小**:
```bash
docker exec infrastructure-db-1 psql -U postgres -d omicsomics -c "SELECT pg_size_pretty(pg_database_size('omicsomics'));"
```

### 性能优化

1. **清理旧日志**:
```bash
docker logs infrastructure-backend-1 > /dev/null 2>&1
```

2. **清理 Python 缓存**:
```bash
docker exec infrastructure-backend-1 find /app -type d -name __pycache__ -exec rm -rf {} +
```

3. **重启所有服务**:
```bash
cd infrastructure && docker-compose restart
```

---

## 联系与支持

如果遇到新问题：

1. 查看本文档的已知问题
2. 检查 GitHub Issues
3. 查看后端日志：`docker logs infrastructure-backend-1`
4. 查看前端日志：`docker logs infrastructure-frontend-1`
5. 提供详细的错误信息和复现步骤

---

**最后更新**: 2025-11-09
**维护者**: Jeblqr
**状态**: 所有核心功能正常运行 ✅
