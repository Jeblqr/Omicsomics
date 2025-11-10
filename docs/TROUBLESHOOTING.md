# Troubleshooting Guide

## 常见问题与解决方案

### 🔴 Pipelines 页面显示空白

**症状：**
- 点击 Pipelines 模块后显示 "Loading..."
- 然后整个页面变成纯白色
- 没有任何错误提示

**根本原因：**
TypeScript 接口定义与 API 返回的数据结构不匹配，导致 React 组件渲染时崩溃。

**原始问题：**
```typescript
// 错误的接口定义
interface PipelineTemplate {
  steps: string[];  // ❌ 实际是对象数组，不是字符串数组
}
```

**API 实际返回：**
```json
{
  "steps": [
    {
      "name": "Quality Control",
      "tool": "fastqc",
      "version": "0.11.9",
      "parameters": {...}
    }
  ]
}
```

**解决方案：**
更新 TypeScript 接口以匹配 API 响应（已在 commit 54984d9 修复）：
```typescript
interface PipelineStep {
  name: string;
  tool: string;
  version: string;
  parameters: Record<string, any>;
}

interface PipelineTemplate {
  steps: PipelineStep[];  // ✅ 正确的类型
}
```

---

### 🔴 文字显示模糊或不可见

**症状：**
- 警告消息显示为白色或浅色文字
- 在浅色背景上难以阅读

**解决方案：**
为所有文本元素添加明确的颜色样式：
- 标题：`color: '#212529'` (深色)
- 描述：`color: '#6c757d'` (中灰色)
- 警告：`color: '#856404'` (棕黄色)

---

### 🔴 Pipelines 页面显示 "Not authenticated"

**症状：**
- 未登录时访问 Pipelines 页面
- 显示红色错误消息

**解决方案：**
1. 登录账号：
   - Email: `demo@omicsomics.com`
   - Password: `demo123456`
2. 刷新页面

---

### 🔴 文件上传失败

**症状：**
- 点击 "Upload File" 后无响应
- 或显示错误消息

**检查清单：**
1. ✅ 确保已选择项目（Project Switcher）
2. ✅ 确保已登录
3. ✅ 检查文件大小是否合理
4. ✅ 检查后端日志：`docker compose logs backend --tail=50`
5. ✅ 检查 MinIO 服务是否运行：`docker compose ps minio`

---

## 调试工具

### 查看前端日志
```bash
cd infrastructure
docker compose logs frontend --tail=100
```

### 查看后端日志
```bash
cd infrastructure
docker compose logs backend --tail=100
```

### 检查所有服务状态
```bash
cd infrastructure
docker compose ps
```

### 测试 API 端点
```bash
# 获取认证 token
TOKEN=$(curl -s -X POST http://localhost:8001/api/v1/auth/login/access-token \
  -H "Content-Type: application/x-www-form-urlencoded" \
  -d "username=demo@omicsomics.com&password=demo123456" | jq -r '.access_token')

# 测试 pipelines 端点
curl -s -H "Authorization: Bearer $TOKEN" http://localhost:8001/api/v1/pipelines/ | jq '.[0]'

# 测试 projects 端点
curl -s -H "Authorization: Bearer $TOKEN" http://localhost:8001/api/v1/projects/ | jq '.'
```

### 运行手动测试脚本
```bash
./scripts/test_frontend_manual.sh
```

---

## 浏览器调试

### 打开开发者工具
按 `F12` 或右键 → "检查"

### 查看控制台错误
1. 打开 Console 标签
2. 查找红色错误消息
3. 检查是否有网络请求失败（Network 标签）

### 清除缓存
强制刷新：`Ctrl + Shift + R` (Linux/Windows) 或 `Cmd + Shift + R` (Mac)

---

## 数据库问题

### 检查数据库连接
```bash
docker compose exec db psql -U omicsomics -c "SELECT version();"
```

### 查看数据库表
```bash
docker compose exec db psql -U omicsomics -c "\dt"
```

### 运行数据库迁移
```bash
./manage.sh db:migrate
```

### 检查迁移状态
```bash
./manage.sh db:status
```

---

## 重启服务

### 重启前端
```bash
cd infrastructure
docker compose restart frontend
```

### 重启后端
```bash
cd infrastructure
docker compose restart backend
```

### 重启所有服务
```bash
cd infrastructure
docker compose restart
```

### 完全重新部署
```bash
./manage.sh stop
./manage.sh clean
./manage.sh start
```

---

## 已知问题记录

### ✅ 已修复 - Pipelines 页面空白 (2025-11-10)
- **Commit:** 54984d9
- **问题：** TypeScript 接口不匹配导致组件崩溃
- **修复：** 更新 PipelineTemplate 接口

### ✅ 已修复 - 文字颜色问题 (2025-11-10)
- **Commit:** 1a1fa98
- **问题：** 多个页面文字显示模糊
- **修复：** 为所有文本元素添加明确颜色

### ✅ 已修复 - 认证错误提示 (2025-11-10)
- **Commit:** 1a1fa98
- **问题：** 401 错误显示不友好
- **修复：** 添加 "Please log in" 提示消息

---

## 获取帮助

如果问题仍然存在：
1. 检查 GitHub Issues：https://github.com/Jeblqr/Omicsomics/issues
2. 查看项目文档：`docs/` 目录
3. 查看架构文档：`docs/ARCHITECTURE.md`
4. 查看 API 文档：http://localhost:8001/docs
