# Omicsomics 项目完整任务清单（TODO_FULL）

此文档为项目从当前状态到可上线（release）所需的完整、连贯、可执行的任务清单。每项都包含目的、实现要点、关键文件路径、验收标准与测试步骤。请在每项完成后更新状态与测试结果。

> 说明：本清单分为“已完成”“进行中”“待办”三部分，并针对 Run pipeline 与文件处理系统给出详细接口契约和测试步骤。

---

## 一、总体说明

目标：将 Omicsomics 打造成可在容器化生产环境运行的完整平台，支持原始文件上传、压缩包解压预览、自动转换为统一格式（UnifiedData）、运行 pipeline（Run）并查询进度、日志与结果，具备完整的测试、CI/CD、部署与监控方案。

成功标准（上线条件）：

- 后端 API 功能齐全且通过端到端集成测试（上传 → 处理 → 检索 → 删除）。
- 前端用户流程顺畅：上传/预览/处理/下载／启动/监控 Run。
- 所有关键服务均在 Docker Compose 下可启动并健康运行。
- 至少一套 CI 自动化（构建、单元测试、集成测试）可在推送时触发。
- 基本的监控与告警能力到位（服务健康、任务失败告警）。

---

## 二、当前状态（快照）

- 已完成（重要）：API 格式标准化、Dashboard 修复、项目级联删除、Run 后台会话修复、FileProcessor 服务、转录组/蛋白质组/代谢组转换器、上传接口集成、processed 数据检索 endpoint、验证脚本。
- 进行中：Docker 化与管理脚本（镜像构建正在进行或刚完成）。
- 待办（关键）：端到端集成测试、前端上传/处理/Run UI、异步后台任务队列（针对大文件）、CI/CD、生产部署、监控、密钥与安全管理等。

---

## 三、详细任务列表（包含测试）

下面按优先级列出任务。每个任务包含：目的、文件/路径、实现要点、验收条件、测试步骤。

### 1. 已完成（记录与验证说明）

1.1 API 格式标准化

- 目的：统一后端 API 返回结构与 HTTP 状态码，避免前端解析差异。
- 主要文件：`backend/app/api/routers/*.py`、`API_FORMAT_STANDARDS.md`
- 验收：创建/删除/列表等常用接口返回符合规范（201/204/200），前端无需特殊分支处理。
- 验证：已通过手工与自动脚本验证（见 `scripts/test_api_formats.py`）。

  1.2 Dashboard hook 修复

- 目的：使用统一的 context hook 显示项目数量与数据。
- 主要文件：`frontend/src/pages/dashboard/*`, `frontend/src/contexts/ProjectsContext.tsx`
- 验收：增删项目后 Dashboard 刷新，展示正确的 count。

  1.3 项目级联删除

- 目的：删除项目时级联移除 Runs、DataFiles、Samples，避免 FK 约束阻塞。
- 主要文件：`backend/app/services/projects.py`
- 验收：删除成功，数据库不留 orphan 记录。验证脚本：`scripts/verify_fixes.py`。

  1.4 Run 后台 DB 会话修复

- 目的：使后台执行的 run 使用独立的 DB session，避免请求生命周期问题。
- 主要文件：`backend/app/services/pipeline_executor.py`、`backend/app/api/routers/runs.py`
- 验收：后台任务能更新进度并写入 DB（手工触发 run 后查看进度）。

  1.5 FileProcessor 与转换器

- 目的：实现自动格式检测与转换为统一格式（UnifiedData）。
- 主要文件：`backend/app/services/file_processor.py`, `backend/app/converters/*.py`（transcriptomics/proteomics/metabolomics）
- 验收：通过 `scripts/test_new_converters.py` 与 `scripts/verify_new_converters.sh` 验证。注意：部分二进制、netCDF 功能为占位或需要额外库（在 TODO 中列出）。

  1.6 上传 API 集成

---

## 附：来自 MERGED_SUMMARIES.md 的合并摘要

以下为自动合并的会话摘要（保留原文以便历史比对）：

### 集成测试摘要

- 状态：所有 7 个集成测试通过
- 覆盖：health、transcriptomics、genomics、proteomics、upload-no-process、cascade-delete、metadata-linking
- 关键修复：SQLAlchemy JSONB 字段需要 `flag_modified()`；CSV 智能检测（proteomics vs transcriptomics）

### 前端数据处理 UI 摘要

- 状态：已完成（Data Catalog 页面）
- 功能：上传时 `process_file` 选项、文件列表显示处理状态与 `omics_type` 标签、查看处理后数据模态、已处理文件多选并准备创建 Run（占位实现）
- 测试：自动化脚本 `scripts/test_frontend_processing.py` 验证上传 → 处理 → 检索流程

### 项目进度摘录

- 完成度：13/16（约 81%）
- 待办主要项：#14 异步后台处理（Celery/RQ）、#15 高级蛋白质组二进制解析、#16 CI/CD 流水线

### 运行与部署要点

- 后端/数据库/前端/MinIO 使用 Docker Compose 管理。
- 数据库（`db`）与后端/前端各自运行于独立容器，持久化使用 volumes（如 `pgdata`），支持单独升级而不影响数据。

### 已清理与建议

- 已将本次会话的概要整合进 `MERGED_SUMMARIES.md`（并在 `TODO_FULL.md` 中引用）。
- 建议：保留 `MERGED_SUMMARIES.md` 作为历史快照，`TODO_FULL.md` 继续作为主任务与运行决策中心。

- 目的：上传文件、保存 raw、可选处理并保存 processed JSON，metadata 链接 raw↔processed。
- 主要文件：`backend/app/api/routers/data.py`
- 验收：上传返回 processing info 并在 metadata 中包含 processed_file_id。

---

### 2. 进行中

2.1 Docker 化全栈（infrastructure）

- 目的：将后端/前端/数据库/MinIO 容器化，简化部署与迁移。
- 主要文件：`backend/Dockerfile`、`infrastructure/docker-compose.yml`、根脚本 `docker-start.sh`,`docker-stop.sh`,`docker-logs.sh`,`docker-shell.sh`。
- 验收：`./docker-start.sh` 能构建并启动服务，访问：`http://localhost:8001/docs` 与 `http://localhost:5173`。
- 测试：在宿主上运行 `./docker-start.sh`；若构建失败，查看 `./docker-logs.sh backend`。

---

### 3. 待办（按优先级）

> 说明：下列任务含详细验收与测试步骤，建议逐项完成并在完成后添加测试记录。

3.1 端到端集成测试（upload → process → retrieve）

- 目的：验证上传（多种格式）、后台处理、processed retrieval、metadata 链接、删除级联在真实或容器环境下工作。
- 关键文件：`scripts/test_api.sh` 或 `tests/integration/test_end_to_end.py`（建议使用 pytest + requests 或 httpx）。
- 实现要点：
  - 使用 Docker 启动全部服务（db/minio/backend/frontend）。
  - 脚本或 pytest：上传样例文件（VCF、CSV、MGF、ZIP），检查 raw 存储与返回的 raw_file_id；若 archive 则先检查 preview；触发处理并轮询直到处理完成；GET /data/{id}/processed 获取统一数据并断言字段。
- 验收条件：所有样例在默认配置下通过；脚本可在 CI 中作为集成测试运行。
- 测试示例（shell）：

```bash
# 启动服务
./docker-start.sh

# 运行集成测试（容器内或本机安装依赖）
pytest -q tests/integration/test_end_to_end.py
```

3.2 前端：文件上传 UI（upload, unzip, preview）

- 目的：提供用户上传界面，支持压缩包解压预览、选择内部文件并提交处理参数。
- 关键路径：`frontend/src/pages/upload/*`, `frontend/src/components/FilePreview.tsx`, `frontend/src/api/*`。
- 功能点：
  - 表单：project、file、process_file、sample_id、organism、reference_genome。
  - 若上传的是压缩包，则在后端返回 preview 列表，前端显示内部清单并允许选择文件。
  - 支持单文件/多文件选择（按项目需求），并支持简单的文件大小/类型校验。
- 验收：用户能上传压缩包并在 UI 中预览内部文件，能选择并提交处理请求，看到处理队列条目。
- 测试：E2E（Playwright / Cypress）脚本：上传 zip → 选择内文件 → 提交 → 轮询处理到完成 → 点击下载 processed JSON。

  3.3 后端：文件处理端到端（接收 archive，安全解压，处理选中文件）

- 目的：实现 archive 的安全解压与基于 `file_path_in_archive` 的选择性处理。
- 关键文件：`backend/app/api/routers/data.py`, `backend/app/services/file_processor.py`。
- 实现要点：
  - 安全解压（禁止路径穿越）：当解压时仅允许相对路径且限制解压目录。
  - 解压后返回文件清单（name、path、size、mime-type）给前端作为 preview。
  - 当前端选择文件并请求处理时，FileProcessor 从临时目录读取那个文件并处理；处理结果写入对象存储并创建 processed datafile。
  - 清理临时文件（成功或失败后）。
- 验收与测试：同 3.1 集成测试；另外增加单元测试覆盖解压边界条件（含恶意文件名）。

  3.4 前端与后端的 processing 集成（轮询/推送）

- 目的：当处理是同步或异步时，前端能获得处理状态并显示进度。若使用异步队列，还需实现 progress 上报与 websocket 或 long-polling。
- 实现方式（短期优先）：使用轮询 GET `/data/{raw_id}/processing_status`；长期：用 websocket 或 SSE 实时推送。
- 验收：前端在处理开始后显示 spinner/进度（若后端返回 progress），处理完成后提供下载链接或展示 unified data preview。

  3.5 Run pipeline 后端与验证

- 目的：提供稳定、可监控的 Run API，支持启动、停止、查询进度与日志，并存储结果文件。
- 关键文件：`backend/app/api/routers/runs.py`, `backend/app/services/pipeline_executor.py`, `backend/app/models/*.py`。
- API 契约参考（简要）：
  - POST `/runs` 或 `/runs/{id}/start` — 启动 Run，返回 202 + run_id（异步）或 200（同步）。
  - POST `/runs/{id}/stop` — 请求停止 Run。
  - GET `/runs/{id}/status` — 返回 {status, progress, started_at, finished_at}。
  - GET `/runs/{id}/logs` — 返回日志（可分页）。
- 实现要点：
  - 后台任务使用独立 DB 会话（已实现）。
  - 将 logs 以增量方式写入 DB 或对象存储（便于前端流式拉取）。
- 验收：能用 API 启动 Run（使用一个小的 pipeline），看到 progress 从 0→100，能获取日志并下载结果。
- 测试：编写 `tests/integration/test_run_pipeline.py`，模拟 pipeline 执行并断言状态更新与日志写入。

  3.6 Run pipeline 前端

- 目的：在 UI 中提供启动/停止/监控 Run 的能力，并可查看历史 Runs / 下载结果。
- 关键路径：`frontend/src/pages/runs/*`, `frontend/src/components/RunControls.tsx`。
- 功能：显示 Runs 列表、Run Details 页面（progress、logs、results）、Start/Stop 控件、连接 processed data files。
- 验收：在前端能启动 Run，实时看到进度条与滚动日志，处理完成能下载结果。
- 测试：E2E 测试，启动 Run 并等待结束，断言 UI 显示正确的状态与可下载文件。

  3.7 异步后台队列（大文件/长任务）

- 目的：将大文件/耗时转换迁移到任务队列，避免阻塞 Web 进程并支持更可靠的重试与进度上报。
- 可选方案：Celery (Redis/RabbitMQ)、RQ (Redis)、FastAPI BackgroundTasks with lightweight worker、Prefect/Argo workflows。
- 推荐：Celery + Redis（成熟且有监控/flower），或在轻量部署上使用 RQ。
- 实现要点：
  - 新增 `worker` 服务到 docker-compose（Redis/RabbitMQ 和 worker image）。
  - FileProcessor 或 pipeline_executor 提交任务到队列并返回 processing_id。
  - Worker 执行任务并定期写进度到 DB（或专门的 progress store），前端轮询或 websocket 订阅读取。
- 验收：提交 >100MB 文件时后端返回 202，并在 10–30s 内异步完成并可查询进度。
- 测试：使用受控环境模拟慢速转换并断言 progress 更新与任务重试逻辑。

  3.8 CI/CD、生产部署、监控、密钥管理

- CI/CD：GitHub Actions（或 GitLab CI），步骤：lint → unit tests → build docker images → push to registry → deploy to staging。
- 生产：编写 k8s manifests / helm chart 或 Docker Compose production 版本，包含 secrets 管理（Vault/Secrets Manager）、备份策略、bucket lifecycle。
- 监控：Prometheus + Grafana（采集 metrics）、ELK 或 Loki/Fluentd（日志）
- 密钥管理与审计：确保 secrets 不入代码库，使用环境变量或 secrets provider。对外接口启用 HTTPS，JWT secret 安全化，数据库账号权限最小化。

---

## 四、测试计划（详尽）

本节详细说明如何运行单元测试、集成测试与 E2E 测试。建议把测试集成进 CI。

### 单元测试

- 工具：pytest
- 位置：`backend/tests/unit/`、`frontend/tests/`（若有）
- 运行：在后端虚拟环境 / 容器中运行

```bash
# 在 backend 环境
cd backend
pip install -r requirements-dev.txt
pytest -q backend/tests/unit
```

### 集成测试（推荐放到 CI 的单独阶段）

- 工具：pytest + httpx / requests
- 目标：在真实依赖（Postgres + MinIO）或用 testcontainers 模拟的环境中运行。
- 实例：`tests/integration/test_end_to_end.py` 包含：
  - 启动 docker compose（或 CI 中的服务）
  - 上传示例文件（vcf/csv/mgf/zip）
  - 触发处理并轮询直到结束
  - 请求 GET /data/{id}/processed 并断言 unified schema

运行示例：

```bash
# 先用 docker-compose 启动服务（本地或 CI）
./docker-start.sh

# 运行集成测试
pytest -q tests/integration/test_end_to_end.py
```

### E2E 测试（前端）

- 工具：Playwright 或 Cypress
- 覆盖：上传 UI、文件选择、提交处理、处理状态、Run 启动/停止/日志

示例（Playwright）流程：

1. 登录到前端（若有认证）
2. 打开 Upload 页面 → 上传 ZIP → 看到 preview → 选择内部文件并提交
3. 等待处理完成 → 断言下载链接存在并可下载
4. 打开 Runs 页面 → 启动 Run → 观察进度条与日志 → 断言完成后可下载结果

### 性能测试（可选）

- 大文件上传（>100MB）性能与队列处理测试
- 工具：locust 或自制脚本

---

## 五、安全与合规（实施建议）

1. 不要将 secrets 写入 repo，使用 `.env` + Docker secrets 或云提供的 secret manager。
2. 在容器中运行时限制权限，避免 root 运行应用进程。
3. 对上传文件做严格检查（类型、大小、解压路径）并对临时文件夹定期清理。
4. 对 API 加入速率限制和权限校验（按 project / user）。
5. 对敏感 API 操作（删除项目、触发生产 pipeline）加入审计日志。

---

## 六、时间估计（粗略）

- 完成 Docker 化并稳定运行：0.5–1 天（取决于环境） — 已在进行中。
- 集成测试（编写 + 调试）：1–2 天。
- 前端 Upload UI 与预览：2–3 天。
- Run UI 与日志：1–2 天。
- 异步队列 PoC（Celery+Redis）：2–3 天。
- CI/CD 初版（Actions）：1–2 天。
- 完整生产化（监控、备份、k8s）：1–2 周（取决范围）。

注意：以上为单人估计。若并行开发或已有模板可复用，时间可缩短。

---

## 七、优先级建议（短期 roadmap）

1. 完成 Docker 构建与验证（保证任何人可一键启动开发环境）。
2. 实现后端解压/preview 与基于 `file_path_in_archive` 的处理 API（关键依赖）。
3. 编写集成测试覆盖 upload→process→retrieve。把测试加入 CI（跑在 staging）。
4. 实现前端 Upload/Preview/Process Integration 与 Run UI。与后端并行进行（后端先实现契约）。
5. 引入任务队列（针对大文件）并改造处理流程。
6. CI/CD、监控、生产化与文档收尾。

---

## 八、下一步（我可以代劳的动作）

请选择其中一个让我立即开始：

- A. 立即实现并测试后端解压/preview & 处理 API（优先推荐）。
- B. 立即开始前端 Upload/Preview UI（需后端 mock 或并行开发）。
- C. 先写集成测试与 CI 配置（便于后续每次变更自动回归）。

我会在开始后把代码修改、测试结果与调用示例（curl 或 Playwright 脚本）发给您，并把变更提交到仓库（或给出补丁）。

---

## 九、附录：关键 API 示例（快速参考）

- Upload with processing (multipart):

```bash
curl -X POST "http://localhost:8001/data/upload" \
  -H "Authorization: Bearer $TOKEN" \
  -F "project_id=1" \
  -F "file=@test_archive.zip" \
  -F "process_file=true" \
  -F "sample_id=S1" \
  -F "organism=Homo sapiens"
```

- Preview internal files (archive upload returns preview) — 前端可先请求 preview 然后再选择：

```json
{"preview": [{"name":"sample1.csv","path":"/tmp/abc/sample1.csv","size":12345,"mime":"text/csv"}, ...]}
```

- Trigger processing on selected internal file:

```bash
curl -X POST "http://localhost:8001/data/<raw_file_id>/process" \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{"file_path_in_archive":"sample1.csv","process_file":true,"sample_id":"S1"}'
```

- Get processed unified data:

```bash
curl -X GET "http://localhost:8001/data/<processed_file_id>/processed" \
  -H "Authorization: Bearer $TOKEN"
```

- Start a run:

```bash
curl -X POST "http://localhost:8001/runs" \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{"pipeline_id": 1, "inputs": [<datafile_ids>], "config": {...}}'
```

---

## 十、联系方式与协作约定

- 我会把每次实现/修改通过补丁或直接编辑（apply_patch）提交，并在任务管理器中同步状态。
- 每完成一个子任务，我会附上：修改文件列表、测试命令、输出摘录、下一步建议。

---

文件位置：`/home/jeblqr/data1/projects/Omicsomics/TODO_FULL.md`

请回复 `A` / `B` / `C` 或给出更具体的优先级指示，我会立即开始实现您选择的工作。
