#!/bin/bash
set -e

echo "🚀 Setting up Omicsomics development environment..."

# Wait for services to be ready
echo "⏳ Waiting for PostgreSQL..."
until pg_isready -h db -U postgres; do
  sleep 1
done

echo "⏳ Waiting for MinIO..."
until curl -f http://minio:9000/minio/health/live > /dev/null 2>&1; do
  sleep 1
done

# Install backend dependencies
echo "📦 Installing backend dependencies..."
cd /workspace/backend
pip install -e ".[dev]"

# Install frontend dependencies
echo "📦 Installing frontend dependencies..."
cd /workspace/frontend
npm install

# Setup database
echo "🗄️ Setting up database..."
cd /workspace/backend
export DATABASE_URL="postgresql+asyncpg://postgres:postgres@db:5432/omicsomics"
alembic upgrade head

# Setup MinIO buckets
echo "🪣 Setting up MinIO buckets..."
pip install minio
python3 << 'EOF'
from minio import Minio
from minio.error import S3Error

client = Minio(
    "minio:9000",
    access_key="minioadmin",
    secret_key="minioadmin",
    secure=False
)

buckets = ["omicsomics-data", "omicsomics-results", "omicsomics-temp"]
for bucket in buckets:
    try:
        if not client.bucket_exists(bucket):
            client.make_bucket(bucket)
            print(f"✅ Created bucket: {bucket}")
        else:
            print(f"✅ Bucket already exists: {bucket}")
    except S3Error as e:
        print(f"❌ Error creating bucket {bucket}: {e}")
EOF

echo "✅ Development environment setup complete!"
echo ""
echo "🎯 Quick commands:"
echo "  - Run backend: cd backend && uvicorn app.main:app --reload --host 0.0.0.0 --port 8001"
echo "  - Run frontend: cd frontend && npm run dev"
echo "  - Run tests: cd backend && pytest"
echo "  - Run tests with coverage: cd backend && pytest --cov=app --cov-report=html"
echo ""
