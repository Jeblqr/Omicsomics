#!/bin/bash

cd /home/jeblqr/data1/projects/Omicsomics/backend

export SECRET_KEY="your-secret-key-change-in-prod"
export DATABASE_URL="postgresql+asyncpg://postgres:postgres@localhost:5432/omicsomics"
export OBJECT_STORAGE_ENDPOINT="http://localhost:9000"
export OBJECT_STORAGE_ACCESS_KEY="minioadmin"
export OBJECT_STORAGE_SECRET_KEY="minioadmin123"

# Activate micromamba environment and run
eval "$(micromamba shell hook --shell bash)"
micromamba activate omicsomics-dev

# Use poetry if available, otherwise use uvicorn directly
if command -v poetry &> /dev/null; then
    poetry run uvicorn app.main:app --host 0.0.0.0 --port 8001 --reload
else
    python3 -m uvicorn app.main:app --host 0.0.0.0 --port 8001 --reload
fi
