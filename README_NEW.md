# Omicsomics 🧬

一个全栈多组学数据分析平台，支持基因组学、转录组学、表观基因组学、蛋白质组学、代谢组学、GWAS 和多组学整合分析。

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.11-blue.svg)](https://python.org)
[![FastAPI](https://img.shields.io/badge/FastAPI-0.121-green.svg)](https://fastapi.tiangolo.com)
[![React](https://img.shields.io/badge/React-18-blue.svg)](https://reactjs.org)

## ✨ 特性

### 🔬 支持的分析类型

- **基因组学 (Genomics)**: BWA 比对、GATK 变异检测、SnpEff 注释、CNVkit 拷贝数分析
- **转录组学 (Transcriptomics)**: STAR 比对、RSEM 定量、DESeq2 差异表达
- **单细胞 (Single-cell)**: Scanpy 聚类、标记基因、轨迹推断
- **表观基因组学 (Epigenomics)**: MACS2 峰检测、HOMER motif 分析
- **蛋白质组学 (Proteomics)**: MaxQuant/MSFragger 搜库、LFQ 定量
- **代谢组学 (Metabolomics)**: XCMS 特征检测、GNPS 注释
- **GWAS**: PLINK 关联分析、LD 计算、PRS 评分、MTAG 跨性状分析
- **多组学整合 (Multi-omics)**: MOFA2、DIABLO
- **可视化**: 火山图、热图、PCA、UMAP、网络图

### 🎯 核心功能

- ✅ **44 个 API 端点** 覆盖完整的多组学分析流程
- ✅ **异步架构** 高性能并发处理
- ✅ **JWT 认证** 安全的用户管理
- ✅ **MinIO 存储** S3 兼容的对象存储
- ✅ **Docker 部署** 容器化部署方案
- ✅ **现代 UI** React + TypeScript + Tailwind CSS

## 🚀 快速开始

### 前置要求

- Docker & Docker Compose
- Python 3.11+
- Node.js 18+
- Micromamba (推荐) 或 Conda

### 1. 克隆仓库

```bash
git clone https://github.com/Jeblqr/Omicsomics.git
cd Omicsomics
```

### 2. 启动开发环境

```bash
# 使用快速启动脚本
./scripts/dev-start.sh
```

这将自动：

- 启动 PostgreSQL 数据库
- 启动 MinIO 对象存储
- 创建测试数据库
- 显示服务访问信息

### 3. 启动后端服务

```bash
cd backend
micromamba create -n omicsomics-dev python=3.11 -y
micromamba activate omicsomics-dev
pip install -e .
uvicorn app.main:app --reload
```

后端 API 将在 http://localhost:8001 运行

### 4. 启动前端服务

```bash
cd frontend
npm install
npm run dev
```

前端界面将在 http://localhost:5173 运行

## 🧪 运行测试

```bash
# 使用测试脚本（推荐）
./scripts/run-tests.sh

# 或手动运行
cd backend
export TEST_DATABASE_URL="postgresql+asyncpg://postgres:postgres@localhost:5432/omicsomics_test"
micromamba run -n omicsomics-dev pytest tests/ -v
```

**当前测试状态**: 11 通过 / 54 失败 (需要示例数据)

## 📁 项目结构

```
Omicsomics/
├── backend/                 # FastAPI 后端
│   ├── app/
│   │   ├── api/            # API 路由
│   │   ├── models/         # 数据库模型
│   │   ├── pipelines/      # 分析流程
│   │   ├── schemas/        # Pydantic 模式
│   │   └── services/       # 业务逻辑
│   └── tests/              # 测试套件
├── frontend/               # React 前端
│   └── src/
│       ├── components/     # UI 组件
│       ├── pages/          # 页面
│       └── lib/            # 工具函数
├── infrastructure/         # Docker 配置
│   └── docker-compose.yml
├── docs/                   # 文档
└── scripts/                # 实用脚本
```

## 🔧 技术栈

### 后端

- **框架**: FastAPI 0.121
- **数据库**: PostgreSQL 15 + SQLAlchemy 2.0
- **存储**: MinIO (S3 兼容)
- **认证**: JWT + Passlib
- **任务队列**: BackgroundTasks

### 前端

- **框架**: React 18 + TypeScript
- **构建**: Vite
- **UI**: Tailwind CSS + shadcn/ui
- **状态管理**: React Context
- **路由**: React Router

### 部署

- **容器化**: Docker + Docker Compose
- **反向代理**: Cloudflare Tunnel (可选)

## 📖 API 文档

启动后端后访问:

- Swagger UI: http://localhost:8001/docs
- ReDoc: http://localhost:8001/redoc

## 🌟 主要 API 端点

### 认证

- `POST /api/auth/register` - 用户注册
- `POST /api/auth/login` - 用户登录

### 基因组学

- `POST /api/genomics/alignment` - 序列比对
- `POST /api/genomics/variant-calling` - 变异检测
- `POST /api/genomics/annotation` - 变异注释

### GWAS

- `POST /api/gwas/qc` - PLINK 质控
- `POST /api/gwas/association` - 关联分析
- `POST /api/gwas/mtag` - MTAG 跨性状分析

### 多组学

- `POST /api/multiomics/mofa2` - MOFA2 整合
- `POST /api/multiomics/diablo` - DIABLO 整合

_更多端点请参见 API 文档_

## 🎨 前端页面

- **认证页面**: 登录/注册
- **仪表板**: 项目概览和工作流状态
- **数据管理**: 上传和管理数据文件
- **分析页面**:
  - GWAS 分析 (QC, 关联分析, MTAG)
  - 多组学整合 (MOFA2, DIABLO)
  - 代谢组学分析 (特征检测, 注释, 定量)
- **结果可视化**: 交互式图表和表格

## 🗺️ 路线图

### 近期目标

- [x] 完成核心分析模块
- [x] 实现 GWAS 分析
- [x] 添加多组学整合
- [ ] 完善前端 UI
- [ ] 添加示例数据集
- [ ] 完整的 Docker 部署

### 长期目标

- [ ] 批量任务处理
- [ ] 自定义工作流构建
- [ ] 实时结果预览
- [ ] 云端部署支持
- [ ] 插件系统

## 🤝 贡献

欢迎贡献！请查看 [CONTRIBUTING.md](CONTRIBUTING.md) 了解详情。

## 📄 许可证

本项目采用 MIT 许可证 - 详见 [LICENSE](LICENSE) 文件

## 🙏 致谢

- FastAPI 框架
- React 社区
- shadcn/ui 组件库
- 所有开源生物信息学工具的开发者

## 📧 联系方式

- GitHub Issues: https://github.com/Jeblqr/Omicsomics/issues
- Email: your.email@example.com

## 📊 项目状态

**当前版本**: 0.7.0 (开发中)

**完成度**: 70%

- Backend: 95% ✅
- Frontend: 40% 🔄
- Tests: 20% ⚠️
- Docs: 50% 📝

---

⭐ 如果这个项目对你有帮助，请给个 Star！
