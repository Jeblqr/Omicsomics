# Omicsomics - 快速参考

## 🚀 快速启动 (Docker 方式 - 推荐)

```bash
# 启动所有服务
./docker-start.sh

# 访问
# - 前端: http://localhost:5173
# - 后端: http://localhost:8001/docs
# - MinIO: http://localhost:9001 (minio/minio123)

# 停止
./docker-stop.sh

# 查看日志
./docker-logs.sh backend
```

## 📊 已实现功能

### ✅ 数据管理

- 项目创建/删除 (支持级联删除)
- 文件上传/下载 (支持加密存储)
- 样本管理
- 自动文件处理和格式转换

### ✅ 文件格式支持

| 组学类型       | 输入格式                         | 输出格式         |
| -------------- | -------------------------------- | ---------------- |
| **基因组学**   | VCF, BED, GTF, GFF, FASTA, FASTQ | VCF, BED, PLINK  |
| **转录组学**   | CSV, TSV, Excel                  | CSV, TSV, DESeq2 |
| **蛋白质组学** | mzML, mzXML, MGF, CSV, TSV       | CSV, TSV, MGF    |
| **代谢组学**   | mzData, mzML, CDF, CSV, TSV      | CSV, TSV         |

### ✅ 统一数据格式

- 自动检测文件类型
- 转换为统一 JSON 格式
- 保留原始文件
- 支持双向转换

### ✅ Pipeline 功能

- Pipeline 创建/编辑
- Run 执行 (异步后台任务)
- 进度跟踪
- 结果存储

### ⚠️ 已知问题

- Run 执行可能需要进一步测试
- 前端 Run 页面可能需要优化

## 🔧 开发

### 本地开发 (不使用 Docker)

```bash
# 后端
cd backend
micromamba activate omicsomics-dev
uvicorn app.main:app --reload

# 前端
cd frontend
npm install
npm run dev
```

### Docker 开发

```bash
# 代码会自动热重载
# 只需修改文件,容器会自动更新

# 进入后端容器
./docker-shell.sh
```

## 🧪 测试

### 测试文件处理

```bash
cd scripts
python test_new_converters.py
```

### 验证转换器

```bash
./scripts/verify_new_converters.sh
```

### API 测试

```bash
# 启动后访问
http://localhost:8001/docs
```

## 📝 API 示例

### 上传文件 (自动处理)

```bash
curl -X POST "http://localhost:8001/data/upload" \
  -H "Authorization: Bearer $TOKEN" \
  -F "project_id=1" \
  -F "file=@test.vcf" \
  -F "organism=Homo sapiens" \
  -F "reference_genome=hg38"
```

### 获取处理后的数据

```bash
curl "http://localhost:8001/data/{datafile_id}/processed" \
  -H "Authorization: Bearer $TOKEN"
```

## 🐛 故障排除

### Docker 问题

```bash
# 查看容器状态
cd infrastructure
docker compose ps

# 查看日志
docker compose logs backend

# 重启服务
docker compose restart backend
```

### 数据库问题

```bash
# 检查数据库连接
docker compose exec db pg_isready -U postgres

# 运行迁移
./docker-shell.sh
alembic upgrade head
```

### 端口冲突

编辑 `infrastructure/docker-compose.yml` 修改端口映射

## 📚 文档

- `DOCKER_GUIDE.md` - Docker 部署完整指南
- `docs/FILE_PROCESSING.md` - 文件处理系统文档
- `NEW_CONVERTERS_SUMMARY.md` - 新转换器总结
- `API_FORMAT_STANDARDS.md` - API 格式标准

## 🎯 下一步

1. ✅ Docker 容器化完成
2. ⏳ 测试 Run 执行功能
3. ⏳ 优化前端 UI
4. ⏳ 添加更多测试用例
5. ⏳ 性能优化 (大文件处理)

## 📞 获取帮助

- 查看日志: `./docker-logs.sh backend`
- API 文档: http://localhost:8001/docs
- 进入容器: `./docker-shell.sh`
