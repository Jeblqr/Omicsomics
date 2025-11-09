#!/bin/bash
# Open a shell in the backend container

cd infrastructure

echo "🐚 Opening shell in backend container..."
docker compose exec backend /bin/bash
