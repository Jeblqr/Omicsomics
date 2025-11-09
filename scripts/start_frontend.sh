#!/usr/bin/env bash
# Start Frontend Development Server

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
FRONTEND_DIR="$PROJECT_ROOT/frontend"

cd "$FRONTEND_DIR"

echo "🚀 Starting Omicsomics Frontend Development Server..."
echo "📁 Frontend directory: $FRONTEND_DIR"
echo ""

# Check if node_modules exists
if [ ! -d "node_modules" ]; then
  echo "📦 Installing dependencies..."
  npm install
  echo ""
fi

# Start development server
echo "🌐 Starting Vite dev server..."
echo "   Frontend will be available at: http://localhost:5173"
echo "   API backend should be running at: http://localhost:8001"
echo ""

npm run dev -- --host
