#!/bin/bash
set -e

echo "🔄 Starting services..."

# Ensure database is ready
until pg_isready -h db -U postgres > /dev/null 2>&1; do
  echo "⏳ Waiting for PostgreSQL..."
  sleep 1
done

# Ensure MinIO is ready
until curl -f http://minio:9000/minio/health/live > /dev/null 2>&1; do
  echo "⏳ Waiting for MinIO..."
  sleep 1
done

echo "✅ All services are ready!"
echo ""
echo "💡 Available services:"
echo "  - PostgreSQL: db:5432 (postgres/postgres)"
echo "  - MinIO API: minio:9000"
echo "  - MinIO Console: http://localhost:9001 (minioadmin/minioadmin)"
echo ""
