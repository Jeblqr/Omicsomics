# 异步文件处理文档

## 概述

Omicsomics 平台支持异步文件处理，适用于大文件或长时间运行的转换任务。异步处理使用 Celery + Redis 实现，避免阻塞 API 请求。

## 架构组件

### 1. Celery Worker
- 独立的后台进程，处理文件转换任务
- 配置在 `docker-compose.yml` 中的 `celery-worker` 服务
- 使用 Redis 作为消息代理和结果后端

### 2. 任务队列
- `file_processing`: 默认队列，处理常规文件
- `large_files`: 专用队列，处理大文件（>100MB）

### 3. 任务类型
- `process_file_async`: 标准文件处理任务
- `process_large_file_async`: 大文件优化处理任务

## 使用方法

### 通过 API 上传文件并异步处理

```bash
curl -X POST "http://localhost:8001/data/upload" \
  -H "Authorization: Bearer $TOKEN" \
  -F "project_id=1" \
  -F "file=@large_file.csv" \
  -F "process_file=true" \
  -F "async_processing=true" \
  -F "sample_id=sample_001"
```

响应示例：
```json
{
  "id": 123,
  "filename": "large_file.csv",
  "processing": {
    "processed": false,
    "async": true,
    "task_id": "a1b2c3d4-e5f6-7890-abcd-ef1234567890",
    "status": "submitted",
    "message": "File processing submitted to background queue"
  }
}
```

### 查询任务状态

```bash
curl -X GET "http://localhost:8001/data/task/{task_id}/status" \
  -H "Authorization: Bearer $TOKEN"
```

响应示例（进行中）：
```json
{
  "state": "PROGRESS",
  "status": "Converting to unified format",
  "progress": 60
}
```

响应示例（完成）：
```json
{
  "state": "SUCCESS",
  "status": "Processing complete",
  "progress": 100,
  "result": {
    "success": true,
    "processed_file_id": 124,
    "omics_type": "transcriptomics",
    "record_count": 15000
  }
}
```

## 本地开发设置

### 1. 启动 Redis
```bash
# 使用 Docker Compose
docker-compose up -d redis

# 或手动启动
redis-server
```

### 2. 启动 Celery Worker
```bash
cd backend
celery -A app.celery_app worker --loglevel=info --concurrency=2
```

### 3. 启动后端服务
```bash
cd backend
uvicorn app.main:app --reload --port 8001
```

### 4. 测试异步处理
```bash
python scripts/test_async_file_processing.py
```

## Docker Compose 部署

完整的 Docker Compose 配置已包含所有必要服务：

```bash
# 启动所有服务（包括 Celery worker）
./docker-start.sh

# 查看 Celery worker 日志
./docker-logs.sh celery-worker

# 查看所有服务状态
docker-compose ps
```

## 监控和调试

### 查看 Celery Worker 日志
```bash
# Docker Compose
docker-compose logs -f celery-worker

# 本地
# Worker 日志会输出到控制台
```

### Flower（可选监控工具）
如需图形化监控界面，可添加 Flower：

```bash
# 安装
pip install flower

# 启动
celery -A app.celery_app flower --port=5555
```

访问 http://localhost:5555 查看任务状态、worker 状态等。

## 配置选项

### 环境变量

在 `.env` 文件或 Docker Compose 环境中配置：

```env
# Redis/Celery 配置
CELERY_BROKER_URL=redis://redis:6379/0
CELERY_RESULT_BACKEND=redis://redis:6379/0

# 任务超时配置
CELERY_TASK_TIME_LIMIT=3600        # 1小时硬超时
CELERY_TASK_SOFT_TIME_LIMIT=3300   # 55分钟软超时
```

### Celery 配置

在 `backend/app/celery_app.py` 中调整：

```python
celery_app.conf.update(
    worker_prefetch_multiplier=1,     # 每次预取1个任务
    worker_max_tasks_per_child=50,    # worker重启前最多处理50个任务
    task_track_started=True,          # 跟踪任务开始状态
)
```

## 性能优化建议

1. **小文件（<10MB）**: 使用同步处理（`async_processing=false`）
2. **中等文件（10-100MB）**: 使用异步处理，默认队列
3. **大文件（>100MB）**: 使用异步处理，large_files 队列

### 自动选择策略

可在上传端点添加自动检测逻辑：

```python
# 根据文件大小自动决定是否异步处理
file_size_mb = len(content) / (1024 * 1024)
auto_async = file_size_mb > 10

if process_file and (async_processing or auto_async):
    # 使用 Celery 异步处理
    ...
```

## 故障排查

### 问题：任务一直处于 PENDING 状态
- 检查 Redis 是否运行：`redis-cli ping`
- 检查 Celery worker 是否运行：`docker-compose ps celery-worker`
- 查看 worker 日志：`docker-compose logs celery-worker`

### 问题：任务失败但没有错误信息
- 查看完整的 worker 日志
- 检查文件格式是否支持
- 验证数据库连接

### 问题：Worker 内存占用过高
- 减少 `worker_prefetch_multiplier`
- 降低 `worker_max_tasks_per_child`
- 增加 worker 实例数量

## 相关文件

- `backend/app/celery_app.py` - Celery 应用配置
- `backend/app/tasks/file_tasks.py` - 文件处理任务
- `backend/app/api/routers/data.py` - 上传和任务状态 API
- `infrastructure/docker-compose.yml` - Docker 服务配置
- `scripts/test_async_file_processing.py` - 测试脚本

## 后续计划

- [ ] 添加任务取消功能
- [ ] 实现 WebSocket 实时进度推送
- [ ] 批量文件处理支持
- [ ] 任务优先级队列
- [ ] 任务失败重试策略
