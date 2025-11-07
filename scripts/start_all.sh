#!/bin/bash
# Omicsomics - Complete startup script

set -e

PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$PROJECT_ROOT"

echo "=== Starting Omicsomics Services ==="

# 1. Check PostgreSQL
echo "[1/3] Checking PostgreSQL..."
if ! pg_isready -q; then
    echo "  Starting PostgreSQL..."
    pg_ctl -D local_db_data -l postgresql.log start
    sleep 2
fi
echo "  ✓ PostgreSQL is running"

# 2. Start MinIO
echo "[2/3] Starting MinIO..."
if lsof -ti:9002 > /dev/null 2>&1; then
    echo "  ✓ MinIO is already running on port 9002"
else
    nohup env MINIO_ROOT_USER=minioadmin MINIO_ROOT_PASSWORD=minioadmin123 \
        ./bin/minio server local_minio_data \
        --address "127.0.0.1:9002" \
        --console-address "127.0.0.1:9003" \
        > minio.log 2>&1 &
    sleep 3
    echo "  ✓ MinIO started on port 9002"
fi

# 3. Start FastAPI
echo "[3/3] Starting FastAPI backend..."
cd backend
if lsof -ti:8001 > /dev/null 2>&1; then
    echo "  Stopping existing FastAPI process..."
    lsof -ti:8001 | xargs kill -9 2>/dev/null || true
    sleep 1
fi

micromamba run -n omicsomics-dev sh -c \
    'SECRET_KEY="your-secret-key-here-change-in-production" DATABASE_URL="postgresql+asyncpg://jeblqr@localhost/omicsomics" uvicorn app.main:app --host 127.0.0.1 --port 8001 --reload' &

sleep 3
echo "  ✓ FastAPI started on port 8001"

echo ""
echo "=== All services are running ==="
echo "  - PostgreSQL: local_db_data/"
echo "  - MinIO API: http://127.0.0.1:9002"
echo "  - MinIO Console: http://127.0.0.1:9003"
echo "  - FastAPI: http://127.0.0.1:8001"
echo "  - API Docs: http://127.0.0.1:8001/docs"
echo ""
echo "Logs:"
echo "  - MinIO: $PROJECT_ROOT/minio.log"
echo "  - PostgreSQL: $PROJECT_ROOT/local_db_data/postgresql.log"
echo ""
