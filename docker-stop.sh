#!/bin/bash
# Stop all Omicsomics Docker containers

echo "🛑 Stopping Omicsomics..."

cd infrastructure
docker compose down

echo "✅ All services stopped"
