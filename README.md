# Omicsomics — 统一组学分析平台（草案）

这是一个面向研究与临床的组学数据统一管理与分析平台草案，目标是以统一的数据模型支持常见组学（基因组、转录组、单细胞、表观、蛋白组、代谢组、宏基因组/微生物组及影像/空间组学），并提供可重复的分析流水线、交互式可视化、API/SDK 与合规控制。

## 目录

- 项目概述
- 主要功能
- 仓库结构（当前）
- 快速开始
- 如何使用 `outline.md`
- 开发与贡献指南
- 下一步建议
- 许可与联系方式

## 项目概述

详见 `outline.md`，该文件包含平台能力、支持的文件格式、数据模型、工作流建议、推荐工具清单、可视化组件、部署与实施路线。

## 主要功能（摘要）

- 文件上传与自动校验（支持分片上传、S3 直传）
- 元数据管理（项目/实验/样本/aliquot 模型）
- 标准化 QC（FastQC / MultiQC 等）
- 可重复的工作流执行（Nextflow/WDL/CWL 等）
- 支持多类组学分析模块（基因组、转录组、单细胞、蛋白组、代谢组等）
- 交互式可视化（UMAP、火山图、IGV 嵌入、网络视图等）
- API/SDK（Python / R / JS）和 CLI 支持
- 权限管理、审计与合规（GDPR/HIPAA 等）

## 仓库结构（当前）

- `backend/` — FastAPI 服务、SQLAlchemy 模型、Alembic 迁移、测试与 Dockerfile
- `frontend/` — React + Vite 前端（React Query、路由、可视化骨架）
- `workflows/` — Nextflow 与 CWL 原型流水线（FastQC → MultiQC）
- `infrastructure/` — `docker-compose.yml` 与本地依赖服务（Postgres、MinIO、API、前端）
- `docs/` — 架构与 API 文档占位
- `.github/workflows/ci.yml` — GitHub Actions：Python lint/test、Node lint/build
- `outline.md` — 详细平台蓝图（组学覆盖面、模块、下一步路线）

## 快速开始（开发者）

1. 克隆仓库：

   git clone <仓库-url>
   cd Omicsomics

2. 启动开发环境（Postgres + MinIO + API + 前端）：

   docker compose -f infrastructure/docker-compose.yml up --build

   - 后端默认监听 `http://localhost:8001`（可在 `.env` 或 compose 中覆盖 `API_PORT`）
   - 前端开发服务器位于 `http://localhost:5173`

3. 单独运行后端（可选）：

   ```bash
   cd backend
   python -m venv .venv
   source .venv/bin/activate
   pip install -e .[dev]
   uvicorn app.main:app --host 0.0.0.0 --port 8001 --reload
   ```

4. 运行测试与静态检查：

   ```bash
   cd backend
   pytest

   cd ../frontend
   npm install
   npm run lint
   npm run build
   ```

5. 阅读 `outline.md` 掌握整体设计、模块优先级与下一阶段任务。

## 如何使用 `outline.md`

`outline.md` 包含：支持的组学、文件格式、数据模型、推荐工具、流水线设计、可视化组件、分阶段实施路线与工程注意点。建议阅读顺序：

1. 总体能力与契约 — 明确输入/输出与成功标准
2. 数据模型与元数据规范 — 确定核心对象与字段
3. 工作流与执行层 — 选择工作流语言与容器策略
4. 各组学分析模块 — 选择优先实现的分析流水线

## 开发与贡献指南

- 提交前请先在本地创建分支并确保变更关联到 issue 或讨论。
- 文档或设计变更通过 PR，代码变更需要包含测试。
- 大体实现建议使用容器化（Docker）与工作流（Nextflow/CWL/WDL）。

## 下一步建议（可选）

- 扩展元数据/样本 JSON Schema，并与 Postgres/Alembic 模型保持同步
- 增加项目 CRUD 的更多端点（更新、删除）、权限模型与审计日志
- 对接工作流编排（Nextflow Tower / Argo Workflows）并在前端显示运行状态
- 丰富前端数据可视化（UMAP、火山图、网络图）与表格筛选能力

## 许可与联系方式

当前为私人/项目草案；请在仓库中添加 LICENSE 文件以明确许可。

---

要我接下来为你：

- 生成 JSON Schema？
- 还是先生成 MVP 工程骨架（包含简单后端 API 与上传示例）？
  请选择其中一个或告诉我你更具体的需求。
