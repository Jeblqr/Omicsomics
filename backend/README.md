# Omicsomics Backend

This directory contains the backend services for the Omicsomics unified omics analysis platform. The backend is implemented with FastAPI and designed to provide modular services that handle data ingestion, metadata management, workflow orchestration, and analytics.

## Structure

```
backend/
├── app/
│   ├── api/          # FastAPI routers, dependencies, and request handlers
│   ├── analytics/    # Analytics services (machine learning, statistics)
│   ├── auth/         # Authentication & authorization logic
│   ├── config/       # Configuration utilities
│   ├── datasets/     # Dataset registration, indexing, access control
│   ├── integrations/ # External service integrations (storage, workflow engines, etc.)
│   ├── logging/      # Logging configuration and helpers
│   ├── metadata/     # Metadata schemas and CRUD operations
│   ├── models/       # SQLAlchemy models or domain models
│   ├── pipelines/    # Pipeline definitions and orchestration helpers
│   ├── projects/     # Project management logic
│   ├── samples/      # Sample handling, sample groups
│   ├── schemas/      # Pydantic schemas
│   ├── search/       # Search/index services
│   ├── services/     # Business services (ingest, QC, reporting)
│   ├── storage/      # Object storage and filesystem abstractions
│   ├── tasks/        # Background tasks and message queue workers
│   ├── utils/        # Shared utilities (validators, ingestion helpers, etc.)
│   ├── visualizations/ # Visualization data preparation
│   ├── workflows/    # Workflow engine integrations
│   ├── __init__.py
│   └── main.py       # FastAPI app entry point
├── requirements/     # Optional extra requirements groups (e.g. prod.txt, dev.txt)
├── scripts/          # Helper scripts for migrations, bootstrapping, etc.
├── tests/            # Pytest-based tests
└── pyproject.toml    # Package metadata & dependencies
```

## Getting Started

1. Create and activate a virtual environment.

   ```bash
   python3 -m venv .venv
   source .venv/bin/activate
   ```

2. Install dependencies (with optional dev extras).

   ```bash
   pip install -e .[dev]
   ```

3. Apply database migrations (requires Postgres running, see `infrastructure/docker-compose.yml`).

   ```bash
   alembic upgrade head
   ```

4. Run development server (defaults to port 8001, configurable via `API_PORT`).

   ```bash
   uvicorn app.main:app --reload --port 8001
   ```

## Next Steps

- Extend CRUD endpoints (update/delete) and add authentication/authorization dependencies.
- Flesh out metadata, sample, and dataset models plus migrations.
- Integrate workflow orchestration clients (Nextflow Tower, Argo, Cromwell) under `app/workflows`.
- Build background task workers for ingestion/QC (Celery, Dramatiq, or FastAPI background tasks).
- Expand tests with fixtures covering multi-tenant scenarios and validation edge cases.

Refer to `outline.md` in the project root for the detailed architecture plan.
