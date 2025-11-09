# Omicsomics 系统架构

本文档描述 Omicsomics 平台的整体架构、技术选型和设计理念。

---

## 目录

- [系统概览](#系统概览)
- [技术栈](#技术栈)
- [架构设计](#架构设计)
- [数据流](#数据流)
- [数据库设计](#数据库设计)
- [API 设计](#api-设计)
- [安全机制](#安全机制)
- [性能优化](#性能优化)

---

## 系统概览

Omicsomics 是一个现代化的多组学数据分析平台，采用**前后端分离**架构，支持容器化部署。

### 核心特点

- **多组学支持**: 基因组学、转录组学、蛋白质组学、代谢组学等
- **微服务架构**: 模块化设计，易于扩展
- **异步处理**: 高并发，响应快速
- **容器化部署**: Docker Compose 一键部署
- **RESTful API**: 标准化接口设计

---

## 技术栈

### 后端

| 技术 | 版本 | 用途 |
|------|------|------|
| **Python** | 3.11+ | 主要编程语言 |
| **FastAPI** | 0.115+ | Web 框架 |
| **SQLAlchemy** | 2.0 | ORM (异步) |
| **PostgreSQL** | 15 | 关系型数据库 |
| **MinIO** | Latest | 对象存储 (S3 兼容) |
| **Redis** | 7 | 缓存 + 消息队列 |
| **Celery** | 5.3+ | 异步任务队列 |
| **Alembic** | Latest | 数据库迁移 |
| **Pydantic** | 2.0+ | 数据验证 |

### 前端

| 技术 | 版本 | 用途 |
|------|------|------|
| **React** | 18 | UI 框架 |
| **TypeScript** | 5 | 类型安全 |
| **Vite** | 5 | 构建工具 |
| **TailwindCSS** | 3 | CSS 框架 |
| **shadcn/ui** | Latest | UI 组件库 |
| **TanStack Query** | 5 | 数据获取与缓存 |
| **React Router** | 6 | 路由管理 |
| **Recharts** | 2 | 数据可视化 |
| **D3.js** | 7 | 高级可视化 |

### DevOps

| 技术 | 用途 |
|------|------|
| **Docker** | 容器化 |
| **Docker Compose** | 服务编排 |
| **GitHub Actions** | CI/CD |
| **Nginx** | 反向代理 (生产环境) |

---

## 架构设计

### 整体架构

```
┌────────────────────────────────────────────────────────────┐
│                        用户                                 │
└──────────────────┬─────────────────────────────────────────┘
                   │ HTTPS
┌──────────────────▼─────────────────────────────────────────┐
│                    Nginx (反向代理)                         │
└──────────┬────────────────────────────┬────────────────────┘
           │                            │
    ┌──────▼──────┐              ┌─────▼──────┐
    │  Frontend   │              │  Backend   │
    │  (React)    │              │ (FastAPI)  │
    │  Port 5173  │              │ Port 8001  │
    └─────────────┘              └─────┬──────┘
                                       │
                ┌──────────────────────┼──────────────────────┐
                │                      │                      │
         ┌──────▼──────┐      ┌───────▼──────┐      ┌───────▼──────┐
         │ PostgreSQL  │      │    MinIO     │      │    Redis     │
         │   (数据库)   │      │  (对象存储)   │      │    (缓存)    │
         │  Port 5432  │      │  Port 9000   │      │  Port 6379   │
         └─────────────┘      └──────────────┘      └───────┬──────┘
                                                             │
                                                      ┌──────▼──────┐
                                                      │   Celery    │
                                                      │  (任务队列)  │
                                                      └─────────────┘
```

### 分层架构

```
┌────────────────────────────────────────────┐
│          Presentation Layer                │
│     (React Components + UI)                │
└────────────────┬───────────────────────────┘
                 │ API Calls (HTTP/JSON)
┌────────────────▼───────────────────────────┐
│         Application Layer (Backend)        │
│  ┌──────────────────────────────────────┐  │
│  │   API Routers (FastAPI Endpoints)   │  │
│  └──────────────┬───────────────────────┘  │
│                 │                           │
│  ┌──────────────▼───────────────────────┐  │
│  │   Business Logic (Services)          │  │
│  │  - Pipeline Templates                │  │
│  │  - Quality Control                   │  │
│  │  - Search & Export                   │  │
│  └──────────────┬───────────────────────┘  │
│                 │                           │
│  ┌──────────────▼───────────────────────┐  │
│  │   Data Access Layer (Models)         │  │
│  │  - SQLAlchemy Models                 │  │
│  │  - Pydantic Schemas                  │  │
│  └──────────────┬───────────────────────┘  │
└─────────────────┼───────────────────────────┘
                  │
┌─────────────────▼───────────────────────────┐
│           Infrastructure Layer              │
│  ┌──────────┐  ┌────────┐  ┌────────────┐  │
│  │PostgreSQL│  │ MinIO  │  │   Redis    │  │
│  │  (数据)  │  │ (文件) │  │  (缓存)    │  │
│  └──────────┘  └────────┘  └────────────┘  │
└─────────────────────────────────────────────┘
```

---

## 数据流

### 文件上传流程

```
User → Frontend → Backend API → MinIO
                     ↓
                PostgreSQL (元数据)
                     ↓
                Celery Task (异步处理)
                     ↓
                Pipeline Execution
```

### Pipeline 执行流程

```
1. 用户选择模板
2. 配置参数
3. 提交 Run 请求
   ↓
4. Backend 创建 Run 记录
   ↓
5. Celery Worker 异步执行
   ↓
6. 实时更新状态 (Redis)
   ↓
7. 存储结果 (MinIO + PostgreSQL)
   ↓
8. 前端轮询获取状态
   ↓
9. 展示结果和日志
```

### 数据查询流程

```
User Request
    ↓
Frontend (React Query)
    ↓
API Endpoint
    ↓
Service Layer (Business Logic)
    ↓
SQLAlchemy Query (Async)
    ↓
PostgreSQL
    ↓
Response (JSON)
    ↓
Frontend Render
```

---

## 数据库设计

### ER 图

```
┌──────────┐      1      ∞  ┌──────────┐      1      ∞  ┌──────────┐
│  User    │◄──────────────►│ Project  │◄──────────────►│  Sample  │
└──────────┘                └──────────┘                └──────────┘
     1                            1                           1
     │                            │                           │
     │                            │                           │
     │ ∞                          │ ∞                         │ ∞
     ▼                            ▼                           ▼
┌──────────┐                ┌──────────┐                ┌──────────┐
│Pipeline  │                │   Run    │                │   File   │
│ Template │                └──────────┘                └──────────┘
└──────────┘                      1                           1
                                  │                           │
                                  │ ∞                         │ ∞
                                  ▼                           ▼
                            ┌──────────┐                ┌──────────┐
                            │Run Task  │                │QC Report │
                            └──────────┘                └──────────┘
```

### 核心表设计

#### users 表
```sql
CREATE TABLE users (
    id SERIAL PRIMARY KEY,
    email VARCHAR(255) UNIQUE NOT NULL,
    hashed_password VARCHAR(255) NOT NULL,
    full_name VARCHAR(255),
    is_active BOOLEAN DEFAULT TRUE,
    is_superuser BOOLEAN DEFAULT FALSE,
    created_at TIMESTAMP DEFAULT NOW(),
    updated_at TIMESTAMP DEFAULT NOW()
);
```

#### projects 表
```sql
CREATE TABLE projects (
    id SERIAL PRIMARY KEY,
    name VARCHAR(255) NOT NULL,
    description TEXT,
    user_id INTEGER REFERENCES users(id),
    created_at TIMESTAMP DEFAULT NOW(),
    updated_at TIMESTAMP DEFAULT NOW(),
    INDEX idx_user_id (user_id)
);
```

#### data_files 表
```sql
CREATE TABLE data_files (
    id SERIAL PRIMARY KEY,
    filename VARCHAR(255) NOT NULL,
    file_type VARCHAR(50),
    data_type VARCHAR(50),  -- genomics, transcriptomics, etc.
    file_path VARCHAR(500),  -- MinIO path
    file_size BIGINT,
    md5_hash VARCHAR(32),
    project_id INTEGER REFERENCES projects(id),
    status VARCHAR(50) DEFAULT 'uploaded',
    metadata JSONB,
    created_at TIMESTAMP DEFAULT NOW(),
    INDEX idx_project_id (project_id),
    INDEX idx_data_type (data_type),
    INDEX idx_status (status)
);
```

#### pipeline_runs 表
```sql
CREATE TABLE pipeline_runs (
    id SERIAL PRIMARY KEY,
    name VARCHAR(255) NOT NULL,
    template_id VARCHAR(100),
    status VARCHAR(50) DEFAULT 'pending',
    progress INTEGER DEFAULT 0,
    start_time TIMESTAMP,
    end_time TIMESTAMP,
    parameters JSONB,
    input_files JSONB,
    output_files JSONB,
    logs TEXT,
    error_message TEXT,
    user_id INTEGER REFERENCES users(id),
    project_id INTEGER REFERENCES projects(id),
    created_at TIMESTAMP DEFAULT NOW(),
    INDEX idx_status (status),
    INDEX idx_template_id (template_id),
    INDEX idx_user_id (user_id)
);
```

---

## API 设计

### RESTful API 规范

```
GET    /api/v1/resource/       # 列表
POST   /api/v1/resource/       # 创建
GET    /api/v1/resource/{id}   # 详情
PUT    /api/v1/resource/{id}   # 更新
DELETE /api/v1/resource/{id}   # 删除
```

### API 模块

| 模块 | 端点前缀 | 功能 |
|------|---------|------|
| 认证 | `/api/v1/auth` | 注册、登录、JWT |
| 项目 | `/api/v1/projects` | CRUD 项目 |
| 数据 | `/api/v1/data` | 文件上传、导出、QC |
| Pipeline | `/api/v1/pipelines` | 模板管理 |
| Runs | `/api/v1/runs` | 执行管理 |
| 工具 | `/api/v1/tools` | 工具配置 |

### 响应格式

**成功响应**:
```json
{
  "id": 1,
  "name": "Project A",
  "created_at": "2025-11-09T10:00:00Z"
}
```

**错误响应**:
```json
{
  "detail": "Resource not found",
  "status_code": 404
}
```

**列表响应**:
```json
{
  "items": [...],
  "total": 100,
  "page": 1,
  "size": 20
}
```

---

## 安全机制

### 认证与授权

- **JWT Token**: 基于 JSON Web Token 的无状态认证
- **密码加密**: Bcrypt 哈希算法
- **HTTPS**: 生产环境强制使用 HTTPS
- **CORS**: 配置跨域资源共享策略

### 数据安全

- **SQL 注入防护**: SQLAlchemy ORM 参数化查询
- **XSS 防护**: React 自动转义
- **CSRF 防护**: SameSite Cookie 策略
- **文件验证**: 上传文件类型和大小限制

### 访问控制

```python
# 基于用户的资源访问控制
async def check_project_access(
    project_id: int,
    user: User,
    db: AsyncSession
) -> bool:
    project = await db.get(Project, project_id)
    return project.user_id == user.id or user.is_superuser
```

---

## 性能优化

### 数据库优化

1. **索引优化**: 关键字段添加索引
2. **连接池**: SQLAlchemy 异步连接池
3. **查询优化**: N+1 问题避免，使用 joinedload
4. **分页**: 大数据集分页查询

### 缓存策略

1. **Redis 缓存**: 热数据缓存
2. **Frontend 缓存**: React Query 自动缓存
3. **浏览器缓存**: 静态资源缓存

### 异步处理

```python
# FastAPI 异步端点
@router.get("/data/")
async def list_data(
    db: AsyncSession = Depends(get_db)
) -> List[DataFile]:
    result = await db.execute(select(DataFile))
    return result.scalars().all()

# Celery 异步任务
@celery_app.task
def process_pipeline(run_id: int):
    # 耗时的分析流程
    ...
```

### 前端优化

1. **代码分割**: React.lazy + Suspense
2. **虚拟滚动**: 大列表优化
3. **防抖节流**: 搜索和滚动事件优化
4. **图片懒加载**: 延迟加载非关键资源

---

## 可扩展性

### 水平扩展

- **无状态后端**: 可部署多个实例
- **负载均衡**: Nginx upstream 配置
- **数据库主从**: PostgreSQL 读写分离
- **对象存储**: MinIO 分布式集群

### 垂直扩展

- **增加资源**: CPU、内存、磁盘
- **优化配置**: 数据库连接池、Worker 数量

---

## 监控与日志

### 日志级别

- **ERROR**: 错误信息
- **WARNING**: 警告信息  
- **INFO**: 常规信息
- **DEBUG**: 调试信息

### 监控指标

- **API 响应时间**: 平均响应时间
- **数据库查询**: 慢查询监控
- **任务队列**: Celery 任务状态
- **系统资源**: CPU、内存、磁盘使用率

---

## 更多资源

- [FastAPI 最佳实践](https://fastapi.tiangolo.com/tutorial/best-practices/)
- [React 性能优化](https://react.dev/learn/render-and-commit)
- [PostgreSQL 性能调优](https://www.postgresql.org/docs/current/performance-tips.html)
- [Docker 最佳实践](https://docs.docker.com/develop/dev-best-practices/)

---

**维护者**: Omicsomics Development Team  
**最后更新**: 2025-11-09
