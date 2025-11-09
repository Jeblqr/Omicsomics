# Omicsomics - 统一组学分析平台

> 一个现代化的 Web 平台，支持多种组学数据（基因组学、转录组学、蛋白质组学、代谢组学等）的统一接收、存储、处理和分析。

[![CI/CD](https://github.com/Jeblqr/Omicsomics/actions/workflows/ci.yml/badge.svg)](https://github.com/Jeblqr/Omicsomics/actions)
[![Python](https://img.shields.io/badge/Python-3.11+-blue.svg)](https://www.python.org/)
[![FastAPI](https://img.shields.io/badge/FastAPI-0.115+-green.svg)](https://fastapi.tiangolo.com/)
[![React](https://img.shields.io/badge/React-18-61dafb.svg)](https://reactjs.org/)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](LICENSE)

---

## ✨ 核心特性

### 🔬 多组学分析支持

- **基因组学**: WGS/WES 变异检测、注释
- **转录组学**: RNA-seq 差异表达分析
- **单细胞**: scRNA-seq 聚类与注释
- **蛋白质组学**: LC-MS/MS 定量分析
- **代谢组学**: 特征检测与代谢物注释
- **表观组学**: ChIP-seq/ATAC-seq Peak calling
- **多组学整合**: MOFA2、DIABLO 整合分析
- **GWAS**: 全基因组关联分析

### 🎨 功能亮点

- ✅ **数据可视化** - 热图、火山图、散点图、MA 图
- ✅ **Pipeline 模板库** - 10+ 预配置分析流程
- ✅ **Pipeline Runs** - 实时监控与日志查看
- ✅ **质量控制** - 自动化 QC 报告与评分
- ✅ **批量操作** - 支持批量上传、删除、导出
- ✅ **高级搜索** - 多作用域搜索与过滤
- ✅ **数据导出** - CSV、TSV、Excel、JSON 格式

### 🏗️ 技术架构

- **后端**: FastAPI + SQLAlchemy 2.0 (异步)
- **前端**: React 18 + TypeScript + TailwindCSS
- **数据库**: PostgreSQL 15
- **存储**: MinIO (S3 兼容)
- **任务队列**: Redis + Celery
- **容器化**: Docker + Docker Compose

---

## 🚀 快速开始

### 前置要求

- Docker & Docker Compose
- Git

### 一键部署（推荐）

````bash
```bash
# 1. 克隆项目
git clone https://github.com/Jeblqr/Omicsomics.git
cd Omicsomics

# 2. 配置环境变量
cp .env.example .env
# 编辑 .env 文件，填入你的配置

# 3. 启动服务
./manage.sh start
````

访问服务：

- **前端**: http://localhost:5173
- **后端 API**: http://localhost:8001/docs
- **MinIO**: http://localhost:9001

详细说明请查看 [QUICKSTART.md](QUICKSTART.md)

````

**默认凭据:**
- MinIO: `minioadmin` / `minioadmin`
- 首次使用请通过前端注册新用户

### 本地开发

详细的本地开发环境配置请参见 **[DEPLOYMENT.md](docs/DEPLOYMENT.md)**

---

## 📚 文档导航

| 文档 | 描述 |
|------|------|
| [DEPLOYMENT.md](docs/DEPLOYMENT.md) | 部署指南 - Docker、本地开发、生产环境配置 |
| [USER_GUIDE.md](docs/USER_GUIDE.md) | 用户手册 - 功能使用与操作指南 |
| [ARCHITECTURE.md](docs/ARCHITECTURE.md) | 架构文档 - 系统设计与技术细节 |
| [API.md](docs/api/README.md) | API 文档 - 端点说明与示例 |
| [IMPLEMENTATION_COMPLETE.md](IMPLEMENTATION_COMPLETE.md) | 实施报告 - 功能实现完成情况 |

---

## 📊 项目统计

- **后端代码**: ~15,000 行 Python
- **前端代码**: ~3,500 行 TypeScript
- **API 端点**: 47+ 个
- **Pipeline 模板**: 10 个
- **支持的组学类型**: 8+ 种

---

## 🎯 核心工作流

**典型使用场景:**

1. **RNA-seq 分析**
   - 上传 FASTQ 文件 → 选择 "RNA-seq Basic Analysis" 模板 → 配置参考基因组 → 执行分析 → 下载差异表达基因列表

2. **蛋白质组学 QC**
   - 上传蛋白组数据 → 选择 "Proteomics QC" 模板 → 查看质量报告 → 导出清洗后的数据

3. **单细胞分析**
   - 上传 CellRanger 输出 → 选择 "Single Cell RNA-seq" 模板 → 聚类与注释 → 可视化 UMAP

---

## 🧪 测试

```bash
# 运行后端测试
cd backend
pytest

# 运行集成测试
./scripts/run_integration_tests.sh

# 查看测试覆盖率
pytest --cov=app --cov-report=html
````

**当前测试状态:**

- ✅ 认证模块
- ✅ 项目管理
- ✅ 样本管理
- ✅ 文件上传
- ✅ Pipeline 模板
- ✅ QC 报告

---

## 📁 项目结构

```
Omicsomics/
├── backend/                 # FastAPI 后端
│   ├── app/
│   │   ├── api/            # API 路由
│   │   ├── models/         # 数据库模型
│   │   ├── services/       # 业务逻辑
│   │   └── main.py         # 应用入口
│   ├── alembic/            # 数据库迁移
│   └── tests/              # 测试文件
│
├── frontend/                # React 前端
│   ├── src/
│   │   ├── components/     # UI 组件
│   │   ├── pages/          # 页面
│   │   └── main.tsx        # 入口
│   └── package.json
│
├── infrastructure/          # Docker 配置
│   ├── docker-compose.yml
│   └── Dockerfile
│
├── docs/                    # 文档目录
│   ├── DEPLOYMENT.md
│   ├── USER_GUIDE.md
│   ├── ARCHITECTURE.md
│   └── api/
│
├── scripts/                 # 实用脚本
│   ├── dev-start.sh        # 开发环境启动
│   └── run-tests.sh        # 测试运行器
│
└── README.md               # 本文件
```

---

## 🤝 贡献指南

我们欢迎任何形式的贡献！

1. Fork 本仓库
2. 创建特性分支 (`git checkout -b feature/AmazingFeature`)
3. 提交更改 (`git commit -m 'Add AmazingFeature'`)
4. 推送分支 (`git push origin feature/AmazingFeature`)
5. 开启 Pull Request

详细贡献指南请参见 [CONTRIBUTING.md](CONTRIBUTING.md)

---

## 📄 许可证

本项目采用 Apache 2.0 许可证 - 详见 [LICENSE](LICENSE) 文件

---

## 🙏 致谢

本项目使用了以下优秀的开源项目:

- [FastAPI](https://fastapi.tiangolo.com/) - 现代化的 Python Web 框架
- [React](https://reactjs.org/) - 用户界面库
- [SQLAlchemy](https://www.sqlalchemy.org/) - Python SQL 工具包
- [MinIO](https://min.io/) - 高性能对象存储
- [PostgreSQL](https://www.postgresql.org/) - 关系型数据库
- [TailwindCSS](https://tailwindcss.com/) - CSS 框架
- [shadcn/ui](https://ui.shadcn.com/) - React 组件库

---

## 📞 支持与联系

- **文档**: [docs/](docs/)
- **Issues**: [GitHub Issues](https://github.com/Jeblqr/Omicsomics/issues)
- **Discussions**: [GitHub Discussions](https://github.com/Jeblqr/Omicsomics/discussions)

---

**版本**: 1.0.0 | **状态**: 🚀 生产就绪 | **更新**: 2025-11-09

**快速链接**: [部署指南](docs/DEPLOYMENT.md) | [用户手册](docs/USER_GUIDE.md) | [API 文档](http://localhost:8001/docs) | [实施报告](IMPLEMENTATION_COMPLETE.md)
