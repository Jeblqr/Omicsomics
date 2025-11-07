#!/bin/bash
cd /workspaces/codespaces-react/backend
export SECRET_KEY="test-secret-key"
export DATABASE_URL="postgresql+asyncpg://postgres:postgres@localhost:5432/omicsomics"
export OBJECT_STORAGE_ENDPOINT="http://localhost:9000"
export OBJECT_STORAGE_ACCESS_KEY="minioadmin"
export OBJECT_STORAGE_SECRET_KEY="minioadmin123"
exec /workspaces/codespaces-react/.venv/bin/uvicorn app.main:app --host 0.0.0.0 --port 8001 --reload
