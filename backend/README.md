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

3. Run development server.

   ```bash
   uvicorn app.main:app --reload --port 8000
   ```

## Next Steps

- Implement `app/main.py` to bootstrap FastAPI with routers and configuration.
- Define Pydantic schemas in `app/schemas` and SQLAlchemy models in `app/models`.
- Set up database migrations using Alembic under `backend/alembic` (to be added).
- Implement service modules: ingestion, QC, metadata, pipeline orchestration, analytics.
- Add integration tests under `backend/tests`.

Refer to `outline.md` in the project root for the detailed architecture plan.
