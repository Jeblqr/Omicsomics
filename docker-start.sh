#!/bin/bash
# Docker-based startup script for Omicsomics

set -e

echo "🐳 Starting Omicsomics with Docker..."
echo "======================================"

cd infrastructure

# Check if Docker is running
if ! docker info > /dev/null 2>&1; then
    echo "❌ Docker is not running. Please start Docker first."
    exit 1
fi

# Stop existing containers
echo "🛑 Stopping existing containers..."
docker compose down

# Build and start services
echo "🔨 Building images..."
docker compose build

echo "🚀 Starting services..."
docker compose up -d db minio

echo "⏳ Waiting for database to be ready..."
sleep 5

# Wait for database
until docker compose exec -T db pg_isready -U postgres > /dev/null 2>&1; do
    echo "   Waiting for PostgreSQL..."
    sleep 2
done
echo "✅ Database is ready"

# Initialize MinIO
echo "📦 Initializing MinIO..."
sleep 3

# Run database migrations
echo "🔄 Running database migrations..."
docker compose run --rm backend alembic upgrade head || echo "⚠️  Migrations may have already been applied"

# Start backend
echo "🚀 Starting backend..."
docker compose up -d backend

# Wait for backend
echo "⏳ Waiting for backend to be ready..."
sleep 5
until curl -f http://localhost:8001/healthz > /dev/null 2>&1; do
    echo "   Waiting for backend... (checking /healthz)"
    sleep 2
done
echo "✅ Backend is ready"

# Start frontend
echo "🚀 Starting frontend..."
docker compose up -d frontend

echo ""
echo "======================================"
echo "✅ Omicsomics is now running!"
echo "======================================"
echo ""
echo "📊 Services:"
echo "   Frontend:  http://localhost:5173"
echo "   Backend:   http://localhost:8001"
echo "   API Docs:  http://localhost:8001/docs"
echo "   MinIO:     http://localhost:9001 (admin: minio / minio123)"
echo "   Database:  localhost:5432 (user: postgres / postgres)"
echo ""
echo "📝 Logs:"
echo "   docker-compose logs -f backend"
echo "   docker-compose logs -f frontend"
echo ""
echo "🛑 Stop all:"
echo "   docker-compose down"
echo ""
