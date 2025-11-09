#!/bin/bash
# View logs from Docker containers

SERVICE=${1:-backend}

cd infrastructure

echo "📝 Viewing logs for: $SERVICE"
echo "   (Press Ctrl+C to exit)"
echo ""

docker compose logs -f $SERVICE
